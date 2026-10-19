'use strict';
let evidence;
function cells(id,values){const row=document.createElement('tr');for(const value of values){const cell=document.createElement('td');cell.textContent=String(value);row.append(cell);}document.getElementById(id).append(row);}
function number(x){return typeof x==='number'?x.toPrecision(4):'not recorded';}
fetch('/api/research-evidence').then(r=>{if(!r.ok)throw Error('Evidence summary is unavailable');return r.json();}).then(data=>{
  evidence=data;
  document.getElementById('status').textContent=`Research snapshot: ${data.date}. ${data.scope}`;
  for(const c of data.axisymmetric)cells('axis',[c.name,number(c.final_time),c.core_detected_all_samples?'Retained':'Not retained',number(c.max_carrier_charge_drift),'Fixed boundaries; finite domain/time']);
  for(const c of data.linear)cells('linear',[c.m,number(c.T),number(c.maximum_gain),c.scope]);
  for(const c of data.nonlinear)cells('nonlinear',[c.name,c.ntheta,number(c.final_time),`${c.min_angles_with_core} / ${c.ntheta}`,number(c.max_relative_energy_drift),c.reason]);
  if(!data.nonlinear.length)document.getElementById('empty').textContent='No completed experimental nonlinear runs are included in this snapshot.';
  for(const text of data.checks){const item=document.createElement('li');item.textContent=text;document.getElementById('checks').append(item);}
}).catch(e=>{document.getElementById('status').textContent=e.message;});
document.getElementById('export').onclick=()=>{if(!evidence)return;const url=URL.createObjectURL(new Blob([JSON.stringify(evidence,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='vortex-research-evidence.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
