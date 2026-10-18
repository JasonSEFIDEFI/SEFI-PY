"""Populate the editable report and render one readable technical PDF."""
import csv
import html
import json
import re
from pathlib import Path
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle,Image,KeepTogether

ROOT=Path(__file__).resolve().parent


def populate():
    rows=list(csv.DictReader((ROOT/'results/table_1_all_results.csv').open()))
    selected=[x for x in rows if float(x['sigma'])==.2]
    g=next(x for x in selected if x['method']=='Geometry');ml=next(x for x in selected if x['method']=='ML')
    shift=next(x for x in csv.DictReader((ROOT/'results/distribution_shift.csv').open()) if x['method']=='Geometry')
    runtime=list(csv.DictReader((ROOT/'results/runtime.csv').open()))
    record=json.loads((ROOT/'results/run_record.json').read_text())
    def f(x):return f'{float(x):.5f}'
    primary=['| Method | Fe | 95% CI for Fe | State fidelity | Accuracy | Delta Fe vs geometry |',
             '|---|---|---|---|---|---|']
    for r in selected:
        primary.append(f"| {r['method']} | {f(r['entanglement_fidelity'])} | [{f(r['entanglement_fidelity_ci_low'])}, {f(r['entanglement_fidelity_ci_high'])}] | {f(r['state_fidelity'])} | {f(r['accuracy'])} | {f(r['paired_delta_fe_vs_geometry'])} |")
    noise=['| Sigma | No recovery Fe | ML Fe | Geometry / posterior Fe | Paired gain vs ML |',
           '|---|---|---|---|---|']
    for sigma in sorted({float(x['sigma']) for x in rows}):
        same=[x for x in rows if float(x['sigma'])==sigma]
        item=lambda n:next(x for x in same if x['method']==n)
        noise.append(f"| {sigma:g} | {f(item('No recovery')['entanglement_fidelity'])} | {f(item('ML')['entanglement_fidelity'])} | {f(item('Geometry')['entanglement_fidelity'])} | {f(item('Geometry')['paired_delta_fe_vs_ml'])} |")
    replacements={
        'PRIMARY_G':f(g['entanglement_fidelity']),'PRIMARY_ML':f(ml['entanglement_fidelity']),
        'DELTA':f(g['paired_delta_fe_vs_ml']),
        'DELTA_CI':f"[{f(g['paired_delta_ci_low'])}, {f(g['paired_delta_ci_high'])}]",
        'RELATIVE_REDUCTION':f"{100*float(g['relative_infidelity_reduction_vs_ml']):.2f}%",
        'PRIMARY_TABLE':'\n'.join(primary),'NOISE_TABLE':'\n'.join(noise),
        'SHIFT_G':f(shift['fe']),'SHIFT_DELTA':f(shift['paired_delta_vs_ml']),
        'SHIFT_CI':f"[{f(shift['delta_ci_low'])}, {f(shift['delta_ci_high'])}]",
        'EQUIVALENCE_ERROR':f"{record['controls']['conventional_estimate_error']:.2e}",
        'CUTOFF_ERROR':f"{record['physical_validation']['cutoff_48_max_overlap_error']:.2e}",
        'KRAUS_ERROR':f"{record['physical_validation']['kraus_completeness_max_error']:.2e}",
        'CHOI_MIN':f"{record['physical_validation']['minimum_choi_eigenvalue']:.2e}",
        'ML_RUNTIME':f"{float(next(r for r in runtime if r['method']=='ML')['median_microseconds_per_trial']):.5f}",
        'G_RUNTIME':f"{float(next(r for r in runtime if r['method']=='Geometry')['median_microseconds_per_trial']):.5f}"}
    source=(ROOT/'report_template.md').read_text(encoding='utf-8')
    for key,value in replacements.items():source=source.replace('{{'+key+'}}',value)
    if '{{' in source:raise ValueError('unfilled report field')
    (ROOT/'BENCHMARK_REPORT.md').write_text(source,encoding='utf-8')
    preliminary=source.split('Recommended wording: "',1)[1].split('"',1)[0]
    (ROOT/'NSF_PRELIMINARY_WORK.md').write_text('# Preliminary work: suggested NSF wording\n\n'+preliminary+'\n\nUse only with the benchmark report and its explicit limitations. This package is local preliminary computational work, not a published result.\n',encoding='utf-8')
    return source


def rich(text):
    s=html.escape(text)
    s=re.sub(r'`([^`]+)`',r'<font name="Courier">\1</font>',s)
    s=re.sub(r'\*\*([^*]+)\*\*',r'<b>\1</b>',s)
    return s


def pdf(source):
    w=504
    body=ParagraphStyle('body',fontName='Helvetica',fontSize=9.8,leading=13.2,spaceAfter=8,textColor=colors.HexColor('#233140'))
    title=ParagraphStyle('title',parent=body,fontName='Helvetica-Bold',fontSize=23,leading=27,spaceAfter=18,textColor=colors.HexColor('#005574'))
    heading=ParagraphStyle('heading',parent=body,fontName='Helvetica-Bold',fontSize=12,leading=15,spaceBefore=5,spaceAfter=10,keepWithNext=True,textColor=colors.HexColor('#005574'))
    cell=ParagraphStyle('cell',parent=body,fontSize=7.5,leading=10,spaceAfter=0)
    story=[];buffer=[]
    def flush():
        if buffer:story.append(Paragraph(rich(' '.join(buffer)),body));buffer.clear()
    lines=source.splitlines();i=0
    while i<len(lines):
        line=lines[i].strip()
        if not line:flush();i+=1;continue
        if line=='---PAGE---':flush();story.append(PageBreak());i+=1;continue
        if line.startswith('# '):flush();story.append(Paragraph(rich(line[2:]),title));i+=1;continue
        if line.startswith('## '):flush();story.append(Paragraph(rich(line[3:]),heading));i+=1;continue
        if line.startswith('!['):
            flush();path=line.split('](',1)[1].rstrip(')');image_path=ROOT/path
            iw,ih=ImageReader(str(image_path)).getSize();story.append(Image(str(image_path),width=w,height=w*ih/iw));story.append(Spacer(1,8));i+=1;continue
        if line.startswith('|'):
            flush();table_lines=[]
            while i<len(lines) and lines[i].strip().startswith('|'):
                current=lines[i].strip()
                if not re.fullmatch(r'[|\- :]+',current):table_lines.append([Paragraph(rich(x.strip()),cell) for x in current.strip('|').split('|')])
                i+=1
            cols=len(table_lines[0]);widths=[w/cols]*cols
            if cols==6:widths=[130,48,100,66,55,105]
            if cols==2:widths=[190,314]
            t=Table(table_lines,colWidths=widths,repeatRows=1,hAlign='LEFT')
            t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e3eff4')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8fa')]),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('LINEBELOW',(0,0),(-1,0),.5,colors.HexColor('#a4c4d2'))]))
            story.append(t);story.append(Spacer(1,10));continue
        buffer.append(line);i+=1
    flush()
    def footer(canvas,doc):
        canvas.saveState();canvas.setFont('Helvetica',8);canvas.setFillColor(colors.HexColor('#617184'))
        canvas.drawString(54,31,'Jason D. Dutton | Preliminary QIS benchmark | October 6, 2026')
        canvas.drawRightString(558,31,str(doc.page));canvas.restoreState()
    doc=SimpleDocTemplate(str(ROOT/'Benchmark_Report.pdf'),pagesize=(612,792),leftMargin=54,rightMargin=54,topMargin=45,bottomMargin=48,title='Pilot-assisted bosonic displacement recovery benchmark',author='Jason D. Dutton')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    print('Generated editable benchmark report, NSF paragraph, and Benchmark_Report.pdf.')


if __name__=='__main__':pdf(populate())
