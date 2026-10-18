// Dependency-free UI wiring smoke test using a minimal DOM/canvas harness.
// Does not substitute for visual verification in a real browser.
import test from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
test('UI initializes, drives simulation, logs protocol, validates settings and latches stop',async()=>{
const elements=new Map(),storage=new Map(),context=new Proxy({}, {get:(o,k)=>o[k]??(()=>{}),set:(o,k,v)=>(o[k]=v,true)});
class Element{constructor(tag='div'){this.tag=tag;this.children=[];this.value='';this.textContent='';this.style={};this.dataset={};this.checked=false;this.classList={toggle(){}};}set id(v){this._id=v;elements.set(v,this);}get id(){return this._id;}append(...nodes){this.children.push(...nodes);}replaceChildren(...nodes){this.children=nodes;}setAttribute(k,v){this[k]=v;}addEventListener(){}getContext(){return context;}getBoundingClientRect(){return {width:800,height:440};}click(){this.onclick?.();}setPointerCapture(){}}
const html=await readFile(new URL('./index.html',import.meta.url),'utf8');
for(const m of html.matchAll(/<([a-z]+)\b([^>]*\bid="([^"]+)"[^>]*)>/g)){const e=new Element(m[1]);e.id=m[3];e.value=m[2].match(/\bvalue="([^"]*)"/)?.[1]||'';e.checked=m[2].includes('checked');}
const tabs=['experiments','inventory','control','safety'].map(tab=>{const e=new Element('button');e.dataset.tab=tab;return e;});
globalThis.document={getElementById:id=>elements.get(id),createElement:tag=>new Element(tag),createTextNode:text=>text,querySelectorAll:selector=>selector==='[data-tab]'?tabs:selector==='.tab'?tabs.map(t=>elements.get(t.dataset.tab)):[]};
globalThis.localStorage={getItem:key=>storage.get(key)??null,setItem:(key,v)=>storage.set(key,v)};globalThis.devicePixelRatio=1;let frame;globalThis.requestAnimationFrame=fn=>{frame=fn;};
await import('./app.js');const el=id=>elements.get(id);assert.ok(el('config').children.length>=5);assert.ok(el('phases').textContent.includes('CH 000'));
el('c-trajectory').value='circle';el('apply').click();el('record').click();el('run').click();let now=performance.now();for(let i=0;i<100;i++)frame(now+=100);assert.ok(Number(el('time').textContent)>5);assert.ok(Number(el('max').textContent)>0);el('finish').click();assert.ok(storage.get('sefi-logs').includes('SIMULATED'));
el('stop').click();const frozen=el('time').textContent;el('run').click();frame(now+=100);assert.equal(el('time').textContent,frozen);assert.match(el('status').textContent,/latched/);
el('reset').click();el('holdDuration').value='1';el('settling').value='0';el('staticRun').click();for(let i=0;i<30;i++)frame(now+=100);assert.match(el('criteria').textContent,/Simulation: PASS/);assert.match(el('criteria').textContent,/Physical: NOT TESTED/);
el('c-amplitude').value='1';el('apply').click();assert.match(el('issues').textContent,/Amplitude exceeds/);
await el('probe').onclick();assert.match(el('connection').textContent,/"position": null/);assert.match(el('connection').textContent,/DISABLED/);
});
