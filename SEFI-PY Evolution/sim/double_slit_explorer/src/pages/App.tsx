import { useCallback, useState } from 'react';
import { motion, useReducedMotion } from 'framer-motion';
import { chapters, links, sources } from '../content';
import type { Mode } from '../content';
import ExperimentCanvas from '../simulations/ExperimentCanvas';
import { DEFAULT, distinguishability, visibility } from '../simulations/physics';
import type { Experiment } from '../simulations/physics';
import Slider from '../components/Slider';
import { Identity, Coupling } from '../components/Concepts';
import AdvancedGeometry from '../components/AdvancedGeometry';

export default function App() {
  const reduced = useReducedMotion();
  const [mode, setMode] = useState<Mode>('qm');
  const [experiment, setExperiment] = useState<Experiment>({ ...DEFAULT });
  const [running, setRunning] = useState(() => !window.matchMedia('(prefers-reduced-motion: reduce)').matches);
  const [paths, setPaths] = useState(36);
  const [spread, setSpread] = useState(65);
  const [overlap, setOverlap] = useState(50);
  const [count, setCount] = useState(0);
  const [reset, setReset] = useState(0);
  const [burst, setBurst] = useState(0);
  const [notice, setNotice] = useState('');
  const onCount = useCallback((n: number) => setCount(n), []);
  const chapter = chapters.find(c => c.id === mode)!;
  const index = chapters.indexOf(chapter);
  const d = distinguishability(experiment);
  const v = visibility(experiment);
  const change = (key: keyof Experiment, value: number | boolean) => setExperiment(e => ({ ...e, [key]: value }));
  const save = () => {
    const canvas = document.querySelector('canvas');
    if (!canvas) return;
    canvas.toBlob(blob => {
      if (!blob) { setNotice('Your browser could not export this image.'); return; }
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a'); a.href = url; a.download = `double-slit-${mode}.png`; a.click();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
      setNotice('Image saved. Keep the model notes with the image when sharing.');
    });
  };
  return <>
    <a className="skip" href="#explore">Skip to the experiment</a>
    <header className="nav shell"><a className="brand" href="#top"><span className="brand-mark">◎</span><span>SEFI <b>Double Slit Explorer</b></span></a><nav aria-label="Main navigation"><a href="#explore">Explore</a><a href="#question">The bigger question</a><a href={links.github} target="_blank" rel="noreferrer">GitHub ↗</a></nav></header>
    <main>
      <section id="top" className="hero shell">
        <div className="hero-copy"><p className="eyebrow">ONE EXPERIMENT. FIVE PERSPECTIVES.</p><h1>What if the Double Slit Experiment is really a <em>geometry problem?</em></h1><p className="lead">Explore the most famous experiment in physics through multiple interpretations.</p><a className="primary" href="#explore">Begin Exploration <span>↗</span></a><p className="small hero-note">A five-minute journey · No physics background needed</p></div>
        <div className="hero-art" aria-hidden="true"><div className="orbital-glow"/><svg viewBox="0 0 500 500"><defs><linearGradient id="orbit"><stop stopColor="#59e3ee"/><stop offset="1" stopColor="#aa7ef7"/></linearGradient></defs><circle cx="250" cy="250" r="165" fill="none" stroke="#253956" strokeDasharray="2 12"/>{[0,40,80,120].map((angle,i) => <motion.ellipse key={angle} cx="250" cy="250" rx="175" ry={45+i*13} fill="none" stroke="url(#orbit)" strokeWidth="1.2" opacity={.7} initial={{rotate:angle}} animate={{rotate:running&&!reduced?[angle,angle+360]:angle}} transition={{duration:40+i*10,repeat:Infinity,ease:'linear'}} style={{transformOrigin:'250px 250px'}}/>)}<circle cx="250" cy="250" r="7" fill="#d9ffff"/></svg><span className="art-label">Different descriptions.<br/>One observable experiment.</span></div>
      </section>
      <div className="trust-strip shell"><span><i className="dot cyan"/>Established quantum physics</span><span><i className="dot violet"/>Clearly marked proposed interpretations</span><span>Interactive, not a proof of new physics</span></div>

      <section id="explore" className="exploration shell">
        <div className="section-heading"><div><p className="eyebrow">THE GUIDED EXPLORER</p><h2>Same experiment.<br/>A new lens at every step.</h2></div><p>Change a setting. Watch the pattern. Ask what changed—and what is only a change in interpretation.</p></div>
        <div className="chapter-nav" aria-label="Choose a chapter">{chapters.map((c,i) => <button key={c.id} className={mode===c.id?'active':''} aria-pressed={mode===c.id} onClick={()=>setMode(c.id)}><span>0{i+1}</span>{c.name}</button>)}</div>
        <div className="chapter-heading"><div><span className={`badge ${index>1?'proposal':''}`}>{index<2?'Established framework':'Proposed interpretation'}</span><p className="eyebrow chapter-number">CHAPTER 0{index+1} / 05</p><h3>{chapter.title}</h3><p>{chapter.body}</p></div><div className="step-controls"><button disabled={index===0} aria-label="Previous chapter" onClick={()=>setMode(chapters[index-1].id)}>←</button><button disabled={index===4} aria-label="Next chapter" onClick={()=>setMode(chapters[index+1].id)}>→</button></div></div>
        <div className="lab-grid"><div className="experiment-panel"><div className="panel-top"><span><i className="dot cyan"/>{index<2?'QUANTUM REFERENCE':'QUANTUM REFERENCE + INTERPRETIVE OVERLAY'}</span><span>{count.toLocaleString()} detections</span></div><ExperimentCanvas experiment={experiment} mode={mode} running={running} paths={paths} spread={spread} overlap={overlap} reset={reset} burst={burst} onCount={onCount}/><div className="transport"><button className="primary compact" onClick={()=>setRunning(!running)}>{running?'Ⅱ Pause':'▶ Play'}</button><button onClick={()=>setBurst(b=>b+1)}>+1,000 detections</button><button onClick={()=>setReset(r=>r+1)}>Clear screen</button><button onClick={save}>Save image ↓</button></div><p className="caption">The dots are sampled detections; the curve is the expected distribution. Every chapter uses the same quantum reference calculation. Accumulation stops at 6,000 dots.</p><p role="status" className="small">{notice}</p></div>
          <aside className="controls"><p className="eyebrow">YOUR EXPERIMENT</p><h4>Shape the pattern</h4><Slider label="Slit width · a" value={experiment.width} min={10} max={80} unit=" μm" onChange={x=>change('width',x)} help="Wider openings make the diffraction envelope narrower."/><Slider label="Slit separation · d" value={experiment.spacing} min={120} max={400} unit=" μm" onChange={x=>change('spacing',x)} help="Larger separation brings the interference fringes closer together."/><label className="toggle"><input type="checkbox" checked={experiment.detector} onChange={e=>change('detector',e.target.checked)}/><span>Which-path marker enabled</span></label><Slider label="Marker interaction · g" value={experiment.interaction} min={0} max={1} step={.05} onChange={x=>change('interaction',x)} help="An illustrative coupling setting, not a calibrated detector efficiency."/><div className="mini-readout"><span>Fringe visibility <b>{Math.round(v*100)}%</b></span><div><i style={{width:`${v*100}%`}}/></div><small>Predicted contrast in the ideal balanced model</small></div></aside></div>

        {mode==='qft'&&<div className="interpretation-box"><span className="badge">FIELD VIEW</span><h4>A detection is an interaction.</h4><p>The field texture is an analogy for a quantum field description, not a numerical QFT solution. The marker can correlate a path with another system. Losing interference does not require destroying the particle.</p><div className="transfer"><span>Path alternatives</span><span style={{opacity:.2+.8*d}}>···· information transfer ···· →</span><span>Marker system</span></div></div>}
        {mode==='gwfm'&&<div className="interpretation-box"><span className="badge proposal">ILLUSTRATIVE GEOMETRY ONLY</span><h4>Explore a proposed bundle of paths.</h4><p>These curves are a proposed way to picture possibilities. They are not observed trajectories. “Boundary selection” is a research interpretation; no new detector probabilities are derived here.</p><div className="three-controls"><Slider label="Visible path count" value={paths} min={12} max={80} onChange={setPaths}/><Slider label="Geometric spread" value={spread} min={20} max={100} unit="%" onChange={setSpread}/><Slider label="Visual overlap" value={overlap} min={0} max={100} unit="%" onChange={setOverlap}/></div><p className="small">These three controls change the drawing only. They do not alter the detection statistics.</p></div>}
        {mode==='sefi'&&<Identity observation={d} running={running}/>}
        {mode==='defi'&&<Coupling observation={d} running={running}/>}
        <div className="explanation-grid"><article><span className="card-number">01 / NOTICE</span><h4>What am I seeing?</h4><p>{chapter.seeing}</p></article><article><span className="card-number">02 / ESTABLISHED</span><h4>What physics says</h4><p>{chapter.conventional}</p></article><article className={index>1?'purple-card':''}><span className="card-number">03 / BOUNDARY</span><h4>{index>1?'What is being proposed':'What this does not show'}</h4><p>{chapter.proposal}</p></article></div>
        <details className="plain-language"><summary>✦ Explain Like I’m 12</summary><p>{chapter.plain}</p></details>

        <section className="observation" aria-labelledby="observation-title"><div className="observation-heading"><div><p className="eyebrow">ONE CONTROL. FIVE EXPLANATIONS.</p><h3 id="observation-title">How much information<br/>is collected?</h3><p>Move from indistinguishable paths toward a readable path record.</p></div><div className="observation-control"><Slider label="Path information setting" value={experiment.information} min={0} max={100} unit="%" onChange={x=>change('information',x)}/><div className="range-labels"><span>No path record</span><span>Fully distinguishable*</span></div><p className="small">*At full interaction with the marker enabled. This is an idealized distinguishability control, not a literal percentage of measured information.</p></div></div><div className="five-readings">
          <article><span>01 / QUANTUM MECHANICS</span><h4>{d<.1?'Alternatives interfere':d>.95?'Fringes disappear':'Contrast fades'}</h4><p>{d<.1?'The paths cannot be told apart, so their amplitudes interfere.':'A path record reduces the interference term. The diffraction envelope remains.'}</p><div className="fringe-swatch" style={{opacity:.25+.75*v}}/></article>
          <article><span>02 / QUANTUM FIELD THEORY</span><h4>{d<.1?'No readable marker':'Systems become correlated'}</h4><p>Interactions can leave a path record in another system. The detector pattern follows the same quantum reference.</p><div className="information-dots" style={{letterSpacing:`${2+10*d}px`}}>● · ● · ●</div></article>
          <article className="purple-card"><span>03 / GWFM · PROPOSED</span><h4>{d<.1?'Overlapping possibilities':'Boundary selection picture'}</h4><p>The narrowing path bundle is an interpretive picture of measurement, not a derived trajectory or a new prediction.</p><div className="bundle-swatch" style={{transform:`scaleY(${1-.65*d})`}}>≋ ≋ ≋</div></article>
          <article className="purple-card"><span>04 / SEFI · PROPOSED</span><h4>Identity in context</h4><p>The identity layers highlight how an interaction may be described. No measured identity observable is defined by these graphics.</p><div className="identity-swatch" style={{opacity:.4+.6*d}}>◎</div></article>
          <article className="purple-card"><span>05 / DEFI · PROPOSED</span><h4>Connected descriptions</h4><p>Geometry and identity respond together in an illustrative coupling. This is an analogy awaiting a physical dynamical law.</p><div className="coupling-swatch" style={{gap:`${28-20*d}px`}}><i/>↔<i/></div></article>
        </div></section>

        <details className="advanced"><summary><span>Advanced Research View<small>Equations, assumptions, and the limits of the model</small></span><span>＋</span></summary><div className="advanced-body"><h4>The calculation behind the screen</h4><p>For equal coherent slit illumination in the far field, with identical rectangular slits and an ideal path marker:</p><div className="equation">I(y) ∝ sinc²(π a sin θ / λ) [1 + V cos(2π d sin θ / λ)]</div><div className="equation">θ = arctan(y/L) · D = g p · V = √(1 − D²)</div><p>Here sinc(x) = sin(x)/x; p is the information slider divided by 100; g is the interaction control. D is set to zero when the marker is disabled. The D = gp mapping is an educational choice. The visibility relation uses an ideal pure marker and saturates the more general bound D² + V² ≤ 1.</p><div className="math-facts"><span>λ = 550 nm</span><span>L = 1 m</span><span>Screen: −30 to +30 mm</span><span>401 sampling bins</span></div><p>The plotted distribution is normalized conditional on arrival within this finite screen window. A seeded pseudo-random sampler creates individual dots. Changing a physical setting clears the accumulated sample; changing the chapter does not. The geometry graphic is schematic and is not drawn to the screen’s physical scale.</p><h4>Proposed geometry is not yet a predictive model</h4><p>Path density, visual overlap, identity layers, and geometry–identity coupling are display concepts here. A physical theory would need defined fields and observables, an action or evolution law, boundary conditions, a derivation of these probabilities, and tests against experiment. The app does not derive matter, spacetime, or a replacement for quantum mechanics.</p><AdvancedGeometry paths={paths} spread={spread} experiment={experiment}/><h4>Assumptions and omissions</h4><p>Monochromatic illumination, equal slit weights, far-field scalar diffraction, no extra phase shift, and an idealized marker. No finite detector resolution, background noise, multi-particle effects, QFT dynamics, or microscopic detector simulation is included.</p></div></details>
      </section>

      <section id="question" className="bigger shell"><p className="eyebrow">THE BIGGER QUESTION</p><h2>Could geometry help us<br/>understand <em>what matter is?</em></h2><p className="lead">A useful question is a beginning—not a conclusion. Can a proposed geometry reproduce quantum predictions? Can it support stable excitations? Can those excitations acquire the measured properties of particles?</p><div className="research-path"><span>Geometry</span><b>?</b><span>Particles</span><b>?</b><span>Matter</span></div><p>These are open research goals. This explorer makes the questions visible; it does not establish their answers.</p><div className="link-row"><a href={links.gwfm} target="_blank" rel="noreferrer">GWFM manuscripts ↗</a><a href={links.sefi} target="_blank" rel="noreferrer">SEFI manuscripts ↗</a><a href={links.defi} target="_blank" rel="noreferrer">DEFI manuscripts ↗</a><a href={links.engine} target="_blank" rel="noreferrer">SEFI-PY research ↗</a></div></section>
      <section className="sources shell"><details><summary>Sources & scientific context</summary><ul>{sources.map(s=><li key={s.url}><a href={s.url} target="_blank" rel="noreferrer">{s.label} ↗</a></li>)}</ul><p>Established sources support the quantum reference. Project manuscripts describe the proposed interpretations; they are not treated as established physics.</p></details></section>
    </main><footer className="shell"><span>SEFI Double Slit Explorer</span><p>An educational research companion. Curiosity, with clear boundaries.</p><a href="#top">Back to top ↑</a></footer>
  </>;
}
