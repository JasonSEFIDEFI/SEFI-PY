'use strict';
let dataset;
const colors=['#ee9c68','#71d9ed'];
function plot(id,speed){
  const canvas=document.getElementById(id),rect=canvas.getBoundingClientRect(),ratio=window.devicePixelRatio||1;
  canvas.width=Math.round(rect.width*ratio);canvas.height=Math.round(rect.height*ratio);
  const ctx=canvas.getContext('2d');ctx.scale(ratio,ratio);
  const w=rect.width,h=rect.height,left=48,top=35,right=w-20,bottom=h-45,max=speed?1.1:6;
  const X=k=>left+(right-left)*k/3,Y=v=>bottom-(bottom-top)*v/max;
  ctx.font='13px system-ui';ctx.fillStyle='#e1eaf2';ctx.fillText(speed?'Gapless phase speed / light speed':'Mode angular frequency',left,20);
  for(let i=0;i<=5;i++){const v=max*i/5,y=Y(v);ctx.strokeStyle='#34455a';ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(right,y);ctx.stroke();ctx.fillText(v.toFixed(2),2,y+4);}
  for(let k=0;k<=3;k++)ctx.fillText(String(k),X(k)-3,bottom+20);
  ctx.fillText('Dimensionless wave number k',left,bottom+40);
  dataset.models.forEach((model,j)=>{for(const field of speed?['phase_speed']:['low','high']){
    ctx.strokeStyle=colors[j];ctx.lineWidth=2;ctx.setLineDash(field==='high'?[7,5]:[]);ctx.beginPath();
    dataset.k.forEach((k,i)=>{const x=X(k),y=Y(model[field][i]);if(i===0)ctx.moveTo(x,y);else ctx.lineTo(x,y);});ctx.stroke();
  }});ctx.setLineDash([]);
  const k=Number(document.getElementById('wave').value);ctx.strokeStyle='#bdc9d3';ctx.beginPath();ctx.moveTo(X(k),top);ctx.lineTo(X(k),bottom);ctx.stroke();
}
function render(){
  if(!dataset)return;
  const k=Number(document.getElementById('wave').value),i=Math.round(k/.0125);
  document.getElementById('waveValue').textContent=k.toFixed(3);
  document.getElementById('readout').textContent=dataset.models.map(m=>`${m.name}: frequency ${m.low[i].toFixed(4)}, phase speed ${m.phase_speed[i].toFixed(4)}`).join(' | ')+(k===0?' (speed is the k -> 0 limit)':'');
  plot('frequencies',false);plot('speeds',true);
}
fetch('/api/dispersion').then(r=>{if(!r.ok)throw Error('Comparison data could not be loaded');return r.json();}).then(data=>{
  dataset=data;render();document.getElementById('status').textContent=data.scope;
}).catch(e=>{document.getElementById('status').textContent=e.message;});
document.getElementById('wave').addEventListener('input',render);window.addEventListener('resize',render);
document.getElementById('download').onclick=()=>{if(!dataset)return;const url=URL.createObjectURL(new Blob([JSON.stringify(dataset,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='vortex-background-comparison.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
