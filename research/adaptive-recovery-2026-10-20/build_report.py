from pathlib import Path
import csv,json,statistics,hashlib,subprocess,sys
from xml.sax.saxutils import escape
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,KeepTogether
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.graphics.shapes import Drawing,Line,Rect,String,Circle
from reportlab.graphics import renderSVG

P=Path(__file__).resolve().parent;out=P/'results/evaluation'
read=lambda n:list(csv.DictReader((out/n).open()))
summary=read('fidelity_summary.csv');diff=read('paired_differences.csv');timing=read('timing.csv')
lock=json.loads((P/'selection.json').read_text());decision=json.loads((out/'decision.json').read_text());protocol=json.loads((P/'protocol.json').read_text())
methods=['fixed',lock['ew'],lock['rolling'],'mixture'];labels={'fixed':'Fixed prior',lock['ew']:'EW (rate 0.03)',lock['rolling']:'Rolling (128)','mixture':'Adaptive mixture','ml':'Likelihood only','no_recovery':'No recovery'}
scenarios=protocol['scenarios'];pretty={'stationary':'Stationary','mean_step':'Mean steps','rotation':'Rotation','joint_drift':'Joint drift','bursts':'Bursts'}
palette=[colors.HexColor(x) for x in ['#667b8c','#17856a','#bc7822','#4257a4']]
score=lambda s,m:float(next(r['mean'] for r in summary if r['scenario']==s and r['method']==m))
cost={m:statistics.median(float(r['microseconds_per_sample']) for r in timing if r['method']==m) for m in methods}
primary=decision['primary'];drift_ew=score('joint_drift',lock['ew']);drift_mix=score('joint_drift','mixture');drift_roll=score('joint_drift',lock['rolling'])
claim=f'In a frozen synthetic pilot-assisted displacement-recovery benchmark, a conventional adaptive mixture improved mean entanglement fidelity over a development-selected 128-sample rolling estimator by {primary["mean"]*100:.3f} percentage points in the prespecified drifting-noise scenario (paired seed-bootstrap 95% interval {primary["lo"]*100:.3f} to {primary["hi"]*100:.3f} points; 20 independent streams). It did not outperform the tested EW estimator in that scenario and was slower in the measured Python implementation. This does not establish a SEFI-specific advantage.'
claims='# Supported claim and exclusions\n\n'+claim+'\n\n## Do not claim\n\n- A new SEFI algorithm, theory validation, quantum advantage, or superiority to all conventional methods.\n- Hardware QEC, a protecting-code threshold, measured logical error reduction, or real-time feedback readiness.\n- Faster or better at matched computation cost: the candidate uses more measured computation.\n- Generalization beyond the five scripted synthetic noise families.\n- Independent preregistration or independent reproduction; the protocol freeze and validation were local.\n\n## Useful engineering observation\n\n'+f'The EW baseline achieved drift fidelity {drift_ew:.6f}, versus {drift_roll:.6f} for the rolling baseline. Its median implementation time was {cost[lock["ew"]]:.2f} versus {cost[lock["rolling"]]:.2f} microseconds per observation ({cost[lock["rolling"]]/cost[lock["ew"]]:.2f} times lower). This is a secondary descriptive comparison of implementations, not a new estimator or SEFI-specific speedup.\n'
(P/'CLAIMS.md').write_text(claims,encoding='utf-8')

figs=P/'figures';figs.mkdir(exist_ok=True)
def fidelity_chart():
 d=Drawing(490,225);left=45;bottom=50;width=425;height=145;low=.78;high=.98
 for v in [.80,.85,.90,.95]:
  y=bottom+(v-low)/(high-low)*height;d.add(Line(left,y,left+width,y,strokeColor=colors.HexColor('#dce3e8')));d.add(String(12,y-3,f'{v:.2f}',fontSize=8))
 for i,s in enumerate(scenarios):
  for j,m in enumerate(methods):
   x=left+i*85+j*15+6;y=bottom+(score(s,m)-low)/(high-low)*height
   d.add(Rect(x,bottom,12,y-bottom,fillColor=palette[j],strokeColor=None))
  d.add(String(left+i*85+32,34,pretty[s],textAnchor='middle',fontSize=8))
 for j,m in enumerate(methods):
  x=35+j*118;d.add(Rect(x,10,7,7,fillColor=palette[j],strokeColor=None));d.add(String(x+10,10,labels[m],fontSize=8))
 d.add(String(45,211,'Unconditional entanglement fidelity (truncated axis: 0.78-0.98)',fontSize=9))
 return d
def speed_chart():
 d=Drawing(490,215);left=55;bottom=45;width=380;height=140
 for x in [0,25,50,75,100,125]:
  xx=left+x/130*width;d.add(Line(xx,bottom,xx,bottom+height,strokeColor=colors.HexColor('#e1e7eb')));d.add(String(xx,bottom-14,str(x),textAnchor='middle',fontSize=8))
 for y in [.945,.95,.955,.96,.965]:
  yy=bottom+(y-.945)/.022*height;d.add(Line(left,yy,left+width,yy,strokeColor=colors.HexColor('#e1e7eb')));d.add(String(15,yy-3,f'{y:.3f}',fontSize=8))
 # Fixed prior is deliberately outside this detailed drift-quality window; annotate it explicitly.
 for j,m in enumerate(methods):
  value=score('joint_drift',m)
  if value<.945:continue
  xx=left+cost[m]/130*width;yy=bottom+(value-.945)/.022*height
  d.add(Circle(xx,yy,4,fillColor=palette[j],strokeColor=None));d.add(String(xx+7,yy+5,labels[m],fontSize=8,fillColor=palette[j]))
 d.add(String(245,11,'Median microseconds per observation; lower is faster',textAnchor='middle',fontSize=9))
 d.add(String(45,202,'Drift fidelity vs processing cost; fixed prior fidelity 0.798 is below plot',fontSize=9))
 return d
charts=[('fidelity_comparison',fidelity_chart()),('quality_vs_runtime',speed_chart())]
for name,d in charts:renderSVG.drawToFile(d,str(figs/(name+'.svg')))
styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='SmallBody',parent=styles['BodyText'],fontSize=9,leading=13,spaceAfter=7));styles.add(ParagraphStyle(name='TinyCell',fontSize=8,leading=10))
story=[]
def para(t,style='SmallBody'):return Paragraph(t,styles[style])
def add(t,style='SmallBody'):story.append(para(t,style))
def table(rows,widths):
 t=Table([[para(str(x),'TinyCell') for x in row] for row in rows],colWidths=widths,repeatRows=1)
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8eef4')),('VALIGN',(0,0),(-1,-1),'TOP'),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#9eafbe')),('LINEBELOW',(0,1),(-1,-1),.3,colors.HexColor('#dfe5eb'))]));story.append(t)
add('Adaptive displacement recovery','Title');add('A bounded benchmark of conventional estimators','Heading2');add('Jason Duran Dutton | Technical research report | October 20, 2026')
add('<b>Outcome:</b> The mixture passes the preselected rolling-baseline comparison, but does not beat the simpler EW estimator on drifting noise. No distinct SEFI advantage or speedup is established.')
add('The evaluation contains 100 independent scenario/seed cases, 240,000 noisy observations, and six methods. All cases completed. The five noise families are synthetic; no quantum hardware was used.')
story.append(charts[0][1]);add('Figure 1. Mean held-out fidelity. The axis is truncated to show small differences. Complete values, including controls, are retained in the CSV results. Fixed calibration is best for stationary noise and bursts among these tested methods.')
table([['Condition','Fixed','EW','Rolling','Mixture']]+[[pretty[s]]+[f'{score(s,m):.6f}' for m in methods] for s in scenarios],[100,90,90,90,90])
add('Scope: sixteen auxiliary probe modes sharing a displacement with an unprotected two-Fock-level payload encoding. Ideal displacement correction; leakage is counted as failure. This is not a stabilizer-QEC demonstration.')
story.append(PageBreak())
add('What was specified and tested','Title')
add('1. Physical and information boundary','Heading2')
add('Each noisy observation is y = beta + epsilon, where each real component of epsilon has variance 1/32. Every method receives the same 512 noisy calibration observations and current/past probe observations. Hidden displacements, future observations and change schedules are unavailable to the decoder.')
add('For residual displacement r, the score is exp(-|r| squared) times (1 - |r| squared / 2) squared. This is the exact entanglement-fidelity overlap for the specified oscillator encoding; it does not renormalize away leakage.')
add('2. Estimators and selection','Heading2')
add('Fixed, EW and rolling estimators use the conventional Gaussian posterior formula with estimated moments. Observation covariance is corrected by subtracting known measurement covariance and projecting onto positive matrices. Adaptive moment updates occur after the current estimate.')
add('The candidate mixes three EW experts (rates 0.01, 0.05, 0.15) using current-observation predictive likelihoods, with fixed-share weight smoothing 0.005. This is conventional empirical Bayes and expert weighting. No uniquely SEFI mathematical step has been derived.')
add('Development: four seeds across five equally weighted conditions selected EW rate 0.03 and rolling window 128. The rolling family had the higher aggregate development fidelity and became the locked primary comparator. The mixture was not tuned.')
add('Evaluation: twenty separate seeds, 2400 observations per condition and seed. Noise families cover stationary noise, mean steps, covariance rotation, joint mean/covariance drift and intermittent bursts. Development and evaluation use the same specified families, so unknown-family generalization is not tested.')
add('3. Statistical interpretation','Heading2')
add('The primary endpoint is mixture minus selected rolling-baseline mean fidelity on joint drift. A paired seed-cluster bootstrap with 10,000 resamples gives a gain of <b>0.210 percentage points</b>, 95% interval <b>0.177 to 0.244 points</b>. Remaining comparisons are descriptive, with no multiple-comparison correction.')
add('The EW baseline nevertheless has higher drift fidelity than the mixture: 0.949550 versus 0.948473. A positive primary result against a selected comparator does not imply dominance over every tested baseline.')
add('The protocol was frozen locally before inspecting development/evaluation outputs. It was not independently preregistered. Source/protocol hashes connect the saved results to that freeze. Neither estimator tuning nor selection used evaluation scores.')
story.append(PageBreak())
add('Speed, limitations and supported wording','Title');story.append(charts[1][1])
table([['Implementation','Median us / observation','Drift fidelity']]+[[labels[m],f'{cost[m]:.2f}',f'{score("joint_drift",m):.6f}'] for m in methods],[180,140,140])
add('Timing uses seven alternating-order complete-stream repetitions on this Windows/Python machine after warm-up. Initialization, prediction and updates are included; simulator and disk I/O are excluded. Results are implementation timings, not worst-case latency or controller deadlines. Microseconds per observation are whole-stream averages, not tail-latency measurements.')
add(f'The mixture is {cost["mixture"]/cost[lock["ew"]]:.2f} times slower than EW in this measurement. EW is {cost[lock["rolling"]]/cost[lock["ew"]]:.2f} times faster than rolling and also has higher drift fidelity. This secondary engineering observation is not a novel algorithmic result.')
add('Supported statement','Heading2');add(escape(claim))
add('Limits','Heading2');add('No measured hardware, actuator effects, photon loss recovery, protecting-code threshold, universal decoder superiority, matched-compute advantage, or new-theory validation is established. The tested set does not cover all strong adaptive methods. The mixture is worse than fixed calibration under stationary noise and bursts. Calibration resources and the shared-displacement assumption require experimental justification.')
story.append(PageBreak())
add('Reproducibility and route to use','Title')
add('Validation completed','Heading2');add('Seven scientific tests passed. The covariance projection and posterior formula match independent matrix calculations. Causal-prefix checks pass. An independent legacy physical-scoring implementation exactly reproduces all 600 saved seed/method scores. All 100 prescribed evaluation cases completed; no interrupted cases enter the summaries.')
add('Run python -m unittest -v test_adaptive, then benchmark.py development, benchmark.py evaluation and benchmark.py timing. NumPy 2.3.5 is the numerical dependency. Preserve the supplied results before rerunning. Full streams, calibration, truth and every method prediction are stored in the compressed NPZ archive.')
add('Immediate use','Heading2');add('Release this as a reproducible robustness benchmark and baseline-selection result. Request an independent reproduction focused on causal information access, the shared-displacement model, comparisons and physical scoring. The evidence supports choosing among these implementations, not promoting an established SEFI breakthrough.')
add('Publication route','Heading2');add('A dated GitHub package is the immediate deliverable. A later Zenodo archive can give a versioned DOI after author metadata/licensing review. No Zenodo deposit has been made. A reusable packaged subsystem could be assessed for pyOpenSci scope; this single analysis is not automatically eligible.')
add('JOSS requires an open-source license, sustained open development and demonstrated research use. The current SEFI-PY license prohibits commercial use without permission, so it does not meet the OSI open-source requirement. Licensing must be resolved before an open-source review route. No license was changed. This new benchmark alone is not a ready JOSS submission. No submission or outreach was performed.')
add('Official route references','Heading2')
for label,url in [('Zenodo GitHub integration','https://help.zenodo.org/docs/github/'),('pyOpenSci scope','https://www.pyopensci.org/software-peer-review/about/package-scope.html'),('JOSS submission requirements','https://joss.readthedocs.io/en/latest/submitting.html')]:add(f'<link href="{url}" color="#28558b">{label}</link>')
add('Prior work and integrity','Heading2');add('The physical_reference.py module is copied unchanged from the earlier SEFI-PY displacement benchmark, which demonstrated equivalence to a conventional Gaussian estimator. This report is a companion computational study and does not replace, withdraw or resubmit any manuscript. Original console equations and prior results are unchanged.')
add('Artifact status: computational technical report; author review and independent reproduction remain pending. No formal priority or novelty claim is made.')
def footer(c,doc):
 c.setFont('Helvetica',8);c.setFillColor(colors.HexColor('#627487'));c.drawString(48,28,'SEFI-PY | Adaptive recovery benchmark | 2026-10-20');c.drawRightString(564,28,str(doc.page))
SimpleDocTemplate(str(P/'Adaptive_Recovery_Report.pdf'),pagesize=(612,792),leftMargin=48,rightMargin=48,topMargin=40,bottomMargin=45,title='Adaptive displacement recovery: bounded benchmark',author='Jason Duran Dutton').build(story,onFirstPage=footer,onLaterPages=footer)
md='# Results: adaptive displacement recovery\n\n'+claim+'\n\n## Fidelity\n\n| Condition | Fixed | EW | Rolling | Mixture |\n|---|---:|---:|---:|---:|\n'
for s in scenarios:md+='| '+pretty[s]+' | '+' | '.join(f'{score(s,m):.6f}' for m in methods)+' |\n'
md+='\n## Timing\n\n| Method | Median microseconds per observation |\n|---|---:|\n'
for m in methods:md+=f'| {labels[m]} | {cost[m]:.2f} |\n'
md+='\nAll 100 cases completed; 240,000 observations, 20 independent seeds per condition. Seven scientific checks passed and 600 saved seed/method scores independently reproduced. Primary positive result does not establish superiority over the EW baseline or a SEFI-specific mechanism. Full details: Adaptive_Recovery_Report.pdf.\n'
(P/'RESULTS.md').write_text(md,encoding='utf-8')
print(json.dumps({'report':str(P/'Adaptive_Recovery_Report.pdf'),'claim':claim,'timing':cost}),flush=True)
