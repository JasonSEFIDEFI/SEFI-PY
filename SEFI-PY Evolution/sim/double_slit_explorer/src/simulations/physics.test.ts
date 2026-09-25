import { describe, expect, it } from 'vitest';
import { BINS, DEFAULT, DISTANCE, WAVELENGTH, distribution, distinguishability, intensity, random, sample, sinc, visibility } from './physics';

describe('quantum reference', () => {
 it('normalizes a symmetric nonnegative distribution across the allowed range', () => {
  for (const width of [10, 30, 80]) for (const spacing of [120, 180, 400]) for (const information of [0, 37, 100]) {
   const { probabilities, cdf } = distribution({ ...DEFAULT, width, spacing, information });
   expect(probabilities).toHaveLength(BINS);
   expect(probabilities.reduce((a,b)=>a+b,0)).toBeCloseTo(1,12);
   expect(cdf.at(-1)).toBe(1);
   expect(Math.min(...probabilities)).toBeGreaterThanOrEqual(0); expect(Math.max(...probabilities.map((p,i)=>Math.abs(p-probabilities[BINS-1-i])))).toBeLessThan(1e-12);
  }
 });
 it('saturates the ideal visibility bound and restores coherence when the marker is disabled', () => {
  for (const information of [0,25,50,75,100]) {const p={...DEFAULT,information}; expect(distinguishability(p)**2+visibility(p)**2).toBeCloseTo(1,14);}
  expect(visibility({...DEFAULT,information:100})).toBe(0);
  expect(visibility({...DEFAULT,information:100,detector:false})).toBe(1);
  expect(visibility({...DEFAULT,information:100,interaction:0})).toBe(1);
 });
 it('removes the interference term but preserves the diffraction envelope at full marking', () => {
  const p={...DEFAULT,information:100};
  for(const y of [-.022,-.006,0,.008,.024]) {const s=y/Math.hypot(DISTANCE,y); expect(intensity(y,p)).toBeCloseTo(sinc(Math.PI*p.width*1e-6*s/WAVELENGTH)**2,14);}
 });
 it('places coherent dark fringes at the predicted angular separation', () => {
  for(const spacing of [120,180,400]) {const s=WAVELENGTH/(2*spacing*1e-6),y=DISTANCE*s/Math.sqrt(1-s*s); expect(intensity(y,{...DEFAULT,spacing})).toBeLessThan(1e-12);}
 });
 it('narrows the envelope as slit width increases', () => {
  const s=WAVELENGTH/(80e-6), y=DISTANCE*s/Math.sqrt(1-s*s);
  expect(intensity(y,{...DEFAULT,width:80,information:100})).toBeLessThan(1e-12);
  expect(intensity(y,{...DEFAULT,width:10,information:100})).toBeGreaterThan(.9);
 });
 it('reproduces sampling and converges to broad reference probabilities', () => {
  const {cdf,probabilities}=distribution(DEFAULT), a=random(421),b=random(421);
  for(let i=0;i<500;i++)expect(sample(cdf,a())).toBe(sample(cdf,b()));
  const rng=random(123),counts=[0,0,0,0,0],expected=[0,0,0,0,0],n=100000;
  probabilities.forEach((p,i)=>expected[Math.min(4,Math.floor(i*5/BINS))]+=p);
  for(let i=0;i<n;i++)counts[Math.min(4,Math.floor(sample(cdf,rng())*5/BINS))]++;
  counts.forEach((x,i)=>expect(Math.abs(x/n-expected[i])).toBeLessThan(.007));
 });
 it('rejects unsupported parameters and invalid random draws', () => {
  expect(()=>distribution({...DEFAULT,width:0})).toThrow();
  expect(()=>distribution({...DEFAULT,information:NaN})).toThrow();
  expect(()=>sample([.5,1],1)).toThrow();
  expect(()=>sample([.5,1],-.1)).toThrow();
  expect(sinc(0)).toBe(1);
 });
});
