/** Equal-illumination scalar Fraunhofer reference, normalized over the visible screen. */
export interface Experiment {width:number;spacing:number;information:number;interaction:number;detector:boolean}
export const DEFAULT:Experiment={width:30,spacing:180,information:0,interaction:1,detector:true};
export const BINS=401, HALF_SCREEN=.03, WAVELENGTH=550e-9, DISTANCE=1;
export function distinguishability(p:Experiment){return p.detector?p.information/100*p.interaction:0;}
export function visibility(p:Experiment){return Math.sqrt(Math.max(0,1-distinguishability(p)**2));}
export function sinc(x:number){return Math.abs(x)<1e-8?1-x*x/6:Math.sin(x)/x;}
export function intensity(y:number,p:Experiment){const s=y/Math.hypot(DISTANCE,y);return sinc(Math.PI*p.width*1e-6*s/WAVELENGTH)**2*(1+visibility(p)*Math.cos(2*Math.PI*p.spacing*1e-6*s/WAVELENGTH));}
export function distribution(p:Experiment){
 if(![p.width,p.spacing,p.information,p.interaction].every(Number.isFinite)||p.width<10||p.width>80||p.spacing<120||p.spacing>400||p.information<0||p.information>100||p.interaction<0||p.interaction>1)throw Error('Unsupported experiment parameters');
 const values=Array.from({length:BINS},(_,i)=>intensity((2*i/(BINS-1)-1)*HALF_SCREEN,p)),sum=values.reduce((a,b)=>a+b,0),probabilities=values.map(x=>x/sum);let n=0;const cdf=probabilities.map(x=>n+=x);cdf[BINS-1]=1;return {probabilities,cdf};
}
export function sample(cdf:number[],u:number){if(!Number.isFinite(u)||u<0||u>=1)throw Error('Sample must be in [0,1)');let a=0,b=cdf.length-1;while(a<b){const m=(a+b)>>>1;if(u<cdf[m])b=m;else a=m+1;}return a;}
export function random(seed:number){let x=seed>>>0;return ()=>{x=(Math.imul(1664525,x)+1013904223)>>>0;return x/4294967296;};}
