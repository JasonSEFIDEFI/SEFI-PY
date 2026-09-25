import {useEffect,useRef,useState} from 'react';
import type {Mode} from '../content';
import {BINS,distribution,distinguishability,random,sample} from './physics';
import type {Experiment} from './physics';
interface Props {experiment:Experiment;mode:Mode;running:boolean;paths:number;spread:number;overlap:number;reset:number;burst:number;onCount:(n:number)=>void}
export default function ExperimentCanvas(props:Props){
 const canvas=useRef<HTMLCanvasElement>(null),latest=useRef(props),points=useRef<{bin:number;jitter:number}[]>([]),rng=useRef(random(421)),dist=useRef(distribution(props.experiment)),clock=useRef(0);latest.current=props;
 const [unsupported,setUnsupported]=useState(false);
 function add(n:number){const list=points.current;for(let i=0;i<n&&list.length<6000;i++)list.push({bin:sample(dist.current.cdf,rng.current()),jitter:rng.current()});latest.current.onCount(list.length);}
 useEffect(()=>{dist.current=distribution(props.experiment);points.current=[];rng.current=random(421);props.onCount(0);},[props.experiment.width,props.experiment.spacing,props.experiment.information,props.experiment.interaction,props.experiment.detector,props.reset,props.onCount]);
 useEffect(()=>{if(props.burst>0)add(1000);},[props.burst]); // Samples use the current, validated reference distribution.
 useEffect(()=>{
  const el=canvas.current;if(!el)return;const context=el.getContext('2d');if(!context){setUnsupported(true);return;}const ctx:CanvasRenderingContext2D=context;let raf=0,last=0,accumulator=0;let w=900,h=400;
  function resize(){const rect=el!.getBoundingClientRect();w=rect.width;h=rect.height;const scale=Math.min(devicePixelRatio||1,2);el!.width=Math.round(w*scale);el!.height=Math.round(h*scale);ctx!.setTransform(scale,0,0,scale,0,0);}
  const observer=new ResizeObserver(resize);observer.observe(el);resize();
  const line=(x:number,y:number,x2:number,y2:number,color:string,width=1)=>{ctx.strokeStyle=color;ctx.lineWidth=width;ctx.beginPath();ctx.moveTo(x,y);ctx.lineTo(x2,y2);ctx.stroke();};
  function render(now:number){
   const p=latest.current,dt=Math.min((now-last)/1000||0,.05);last=now;if(p.running){clock.current+=dt;accumulator+=dt*85;if(accumulator>=10){add(Math.floor(accumulator));accumulator%=1;}}
   const t=clock.current,D=distinguishability(p.experiment),sx=w*.08,bx=w*.31,ex=w*.78,cy=h*.5,span=h*.72,sep=(.16+p.experiment.spacing/400*.17)*h,slits=[cy-sep/2,cy+sep/2],slot=7+p.experiment.width/80*15;
   ctx.clearRect(0,0,w,h);const bg=ctx.createRadialGradient(w*.48,cy,10,w*.48,cy,w*.65);bg.addColorStop(0,'#13263e');bg.addColorStop(1,'#08111e');ctx.fillStyle=bg;ctx.fillRect(0,0,w,h);
   ctx.strokeStyle='#8dbada09';ctx.lineWidth=1;for(let x=0;x<w;x+=36)line(x,0,x,h,'#8dbada09');for(let y=0;y<h;y+=36)line(0,y,w,y,'#8dbada09');
   ctx.save();ctx.beginPath();ctx.rect(sx,cy-span/2,bx-sx,span);ctx.clip();for(let k=0;k<11;k++){const r=((k/11+t*.12)%1)*(bx-sx+80);ctx.strokeStyle='#72e5df55';ctx.lineWidth=1.1;ctx.beginPath();ctx.arc(sx,cy,r,-.95,.95);ctx.stroke();}ctx.restore();
   if(p.mode==='qm'||p.mode==='qft'){
    ctx.save();ctx.beginPath();ctx.rect(bx+4,30,ex-bx-5,h-60);ctx.clip();
    if(p.mode==='qft'){for(let x=bx+10;x<ex;x+=9)for(let y=35;y<h-35;y+=9){const a=Math.hypot(x-bx,y-slits[0]),b=Math.hypot(x-bx,y-slits[1]),z=(Math.cos(a*.09-t*3)+Math.cos(b*.09-t*3))/2;ctx.fillStyle=`rgba(74,150,255,${.02+Math.max(0,z)*.12})`;ctx.fillRect(x,y,5,5);}}
    slits.forEach((y,j)=>{for(let k=0;k<15;k++){const r=((k/15+t*.055)%1)*(ex-bx+140);ctx.strokeStyle=j===0?`rgba(98,226,219,${.24*(1-D*.5)})`:`rgba(156,132,255,${.26*(1-D*.5)})`;ctx.lineWidth=1;ctx.beginPath();ctx.arc(bx,y,r,-Math.PI/2,Math.PI/2);ctx.stroke();}});ctx.restore();
   }else{
    for(let k=0;k<p.paths;k++){const side=k%2,start=slits[side],f=k/Math.max(1,p.paths-1),spread=p.spread/100*(1-D*.65),end=cy+(f-.5)*span*.95,cp1=start+(f-.5)*span*spread,cp2=cy+(f-.5)*span*(.4+spread);
     ctx.beginPath();ctx.moveTo(sx,cy);ctx.quadraticCurveTo(bx*.8,cy,bx,start);ctx.bezierCurveTo(bx+(ex-bx)*.3,cp1,bx+(ex-bx)*.6,cp2,ex,end);ctx.strokeStyle=side===0?`rgba(115,222,232,${.09+p.overlap/100*.2})`:`rgba(177,144,250,${.1+p.overlap/100*.2})`;ctx.lineWidth=.8;ctx.stroke();
     const q=(t*.13+k/p.paths)%1,iq=1-q,x=iq**3*bx+3*iq*iq*q*(bx+(ex-bx)*.3)+3*iq*q*q*(bx+(ex-bx)*.6)+q**3*ex,y=iq**3*start+3*iq*iq*q*cp1+3*iq*q*q*cp2+q**3*end;ctx.fillStyle=side?'#c3adff':'#8ee9e0';ctx.globalAlpha=.6;ctx.beginPath();ctx.arc(x,y,1.3,0,Math.PI*2);ctx.fill();ctx.globalAlpha=1;
    }
    const glow=ctx.createRadialGradient((bx+ex)/2,cy,0,(bx+ex)/2,cy,h*.26);glow.addColorStop(0,`rgba(132,113,245,${p.overlap/100*.12})`);glow.addColorStop(1,'transparent');ctx.fillStyle=glow;ctx.fillRect(bx,0,ex-bx,h);
   }
   ctx.fillStyle='#1e354c';ctx.fillRect(bx-4,34,8,h-68);for(const sy of slits){ctx.clearRect(bx-5,sy-slot/2,10,slot);ctx.fillStyle='#83e6df';ctx.fillRect(bx-4,sy-slot/2,8,2);ctx.fillRect(bx-4,sy+slot/2-2,8,2);}
   if(p.experiment.detector&&D>0){for(const sy of slits){ctx.strokeStyle=`rgba(252,190,127,${.3+.6*D})`;ctx.lineWidth=1.2;ctx.strokeRect(bx+10,sy-11,16,22);if(p.mode==='qft'){ctx.setLineDash([3,5]);line(bx+26,sy,bx+65,sy-35,'#f5c18c88');ctx.setLineDash([]);}}}
   ctx.shadowBlur=20;ctx.shadowColor='#68e2df';ctx.fillStyle='#b5fff6';ctx.beginPath();ctx.arc(sx,cy,4,0,Math.PI*2);ctx.fill();ctx.shadowBlur=0;
   line(ex,cy-span/2,ex,cy+span/2,'#94acc244');ctx.fillStyle='#8fece399';for(const hit of points.current){const y=cy-span/2+hit.bin/(BINS-1)*span;ctx.fillRect(ex+3+hit.jitter*w*.045,y,1.5,1.5);}
   const max=Math.max(...dist.current.probabilities),curveX=w*.88;ctx.beginPath();dist.current.probabilities.forEach((v,i)=>{const x=curveX+v/max*w*.07,y=cy-span/2+i/(BINS-1)*span;i?ctx.lineTo(x,y):ctx.moveTo(x,y);});ctx.strokeStyle='#9eb0ed';ctx.lineWidth=1.5;ctx.stroke();
   ctx.font='10px system-ui';ctx.fillStyle='#9badc4';ctx.textAlign='center';ctx.fillText('SOURCE',sx,23);ctx.fillText('TWO SLITS',bx,23);ctx.fillText('DETECTIONS',ex+w*.025,23);ctx.fillText('EXPECTED',curveX+w*.025,h-13);ctx.fillText('−30 mm',ex,cy-span/2-9);ctx.fillText('+30 mm',ex,cy+span/2+15);ctx.textAlign='left';
   raf=requestAnimationFrame(render);
  }
  raf=requestAnimationFrame(render);return ()=>{cancelAnimationFrame(raf);observer.disconnect();};
 },[]);
 return <div className="canvas-shell"><canvas ref={canvas} aria-label="Illustrative two-slit apparatus, sampled detections and quantum reference distribution" role="img"/><span className="canvas-caption">{props.mode==='qm'?'Amplitude illustration':props.mode==='qft'?'Field analogy':'Proposed geometry illustration'} · apparatus not to scale</span>{unsupported&&<p className="canvas-fallback">Canvas is unavailable. The explanations and fringe-contrast values below remain accessible.</p>}</div>;
}
