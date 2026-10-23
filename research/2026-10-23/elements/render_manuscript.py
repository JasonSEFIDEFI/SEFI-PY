"""Render the editable Markdown manuscript with equation images and paged tables."""
from pathlib import Path
import hashlib
import html
import re
import shutil
import matplotlib
matplotlib.use('Agg')
from matplotlib import mathtext
from matplotlib.font_manager import FontProperties
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image,Table,TableStyle,KeepTogether,PageBreak

ROOT=Path(__file__).resolve().parent
TMP=ROOT.parents[1]/'tmp/pdfs/sefi-elements';TMP.mkdir(parents=True,exist_ok=True)
OUT=ROOT.parent/'pdf/SEFI_ELEMENTS_Working_Manuscript_2026-10-23.pdf';OUT.parent.mkdir(exist_ok=True)
for name,file in [('TimesCustom','times.ttf'),('TimesCustom-Bold','timesbd.ttf'),('TimesCustom-Italic','timesi.ttf'),('TimesCustom-BoldItalic','timesbi.ttf')]:
    pdfmetrics.registerFont(TTFont(name,'C:/Windows/Fonts/'+file))
pdfmetrics.registerFontFamily('TimesCustom',normal='TimesCustom',bold='TimesCustom-Bold',italic='TimesCustom-Italic',boldItalic='TimesCustom-BoldItalic')
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='MainText',fontName='TimesCustom',fontSize=10.3,leading=14.2,spaceAfter=7))
styles.add(ParagraphStyle(name='MainTitle',fontName='TimesCustom-Bold',fontSize=22,leading=25,spaceAfter=10,textColor=colors.HexColor('#16334b')))
styles.add(ParagraphStyle(name='SubTitle',fontName='TimesCustom-Bold',fontSize=15,leading=19,spaceAfter=12))
styles.add(ParagraphStyle(name='SectionTitle',fontName='TimesCustom-Bold',fontSize=12,leading=15,spaceBefore=13,spaceAfter=6,keepWithNext=True))
styles.add(ParagraphStyle(name='SubSectionTitle',fontName='TimesCustom-Bold',fontSize=10.8,leading=14,spaceBefore=9,spaceAfter=5,keepWithNext=True))
styles.add(ParagraphStyle(name='TableCell',fontName='TimesCustom',fontSize=8.3,leading=10.2))
styles.add(ParagraphStyle(name='CodeText',fontName='Courier',fontSize=8,leading=10.5,spaceAfter=5))
WIDTH=468

def normalize(s):return s.replace('\u2011','-').replace('\u2013','-').replace('\u2014','-')

def math_image(eq,size=10.3):
    eq=eq.strip().replace(r'\sqrt A',r'\sqrt{A}');key=hashlib.sha256((eq+str(size)).encode()).hexdigest()[:20]
    path=TMP/(key+'.png')
    mathtext.math_to_image('$'+eq+'$',str(path),prop=FontProperties(size=size),dpi=240,format='png')
    with PILImage.open(path) as img:w,h=img.size
    return path,w*72/240,h*72/240

def inline(s):
    s=normalize(s);parts=[]
    for i,part in enumerate(re.split(r'(\$[^$\n]+\$)',s)):
        if part.startswith('$') and part.endswith('$'):
            p,w,h=math_image(part[1:-1]);parts.append(f'<img src="{p.as_posix()}" width="{w:.2f}" height="{h:.2f}" valign="-3"/>')
        else:
            part=html.escape(part)
            part=re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',r'<link href="\2" color="#184d78">\1</link>',part)
            part=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',part)
            part=re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)',r'<i>\1</i>',part)
            part=re.sub(r'`([^`]+)`',r'<font name="Courier" size="9">\1</font>',part)
            parts.append(part)
    return ''.join(parts)

def cell(s):
    labels={'SEMF_known_shell_features':'SEMF + known shell features','SEFI_geometric_prototype':'Geometric prototype',
        'conventional_matched_kernel':'Matched conventional kernel','distance_stability_hypothesis':'Stability-distance hypothesis',
        'conventional_polynomial':'Conventional polynomial','conventional_kernel':'Conventional kernel',
        'training_median':'Training median','proton_extrapolation':'Proton extrapolation','neutron_frontier':'Neutron frontier',
        'discovery_year_proxy':'Discovery-year proxy','random_4200':'Random split 4200'}
    return Paragraph(inline(labels.get(s,s)),styles['TableCell'])

def footer(canvas,doc):
    canvas.saveState();canvas.setStrokeColor(colors.HexColor('#ccd5db'));canvas.line(54,39,558,39)
    canvas.setFont('TimesCustom',8);canvas.setFillColor(colors.HexColor('#506170'))
    canvas.drawString(54,27,'SEFI:ELEMENTS | Author-review working manuscript | 23 October 2026')
    canvas.drawRightString(558,27,str(doc.page));canvas.restoreState()

def main():
    lines=(ROOT/'MANUSCRIPT.md').read_text(encoding='utf-8').splitlines();story=[];i=0
    while i<len(lines):
        line=lines[i].strip()
        if not line:i+=1;continue
        if line.startswith('|'):
            table=[]
            while i<len(lines) and lines[i].strip().startswith('|'):
                fields=[v.strip() for v in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r'[:\-]+',v) for v in fields):table.append([cell(v) for v in fields])
                i+=1
            cols=len(table[0]);widths={5:[91,143,48,86,100],4:[168,65,118,117]}.get(cols,[WIDTH/cols]*cols)
            t=Table(table,colWidths=widths,repeatRows=1,hAlign='LEFT')
            t.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8eff4')),
                ('LINEBELOW',(0,0),(-1,0),.7,colors.HexColor('#8095a5')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8fa')]),
                ('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]))
            if len(table)>5:
                t.setStyle(TableStyle([('NOSPLIT',(0,j),(-1,min(j+3,len(table)-1))) for j in range(1,len(table),4)]))
            story.extend([t,Spacer(1,9)]);continue
        if line.startswith('$$'):
            eq=line[2:-2];p,w,h=math_image(eq,11)
            scale=min(1.,WIDTH/w);story.extend([Spacer(1,5),Image(str(p),width=w*scale,height=h*scale),Spacer(1,10)]);i+=1;continue
        if line.startswith('!['):
            match=re.match(r'!\[([^]]*)\]\(([^)]+)\)',line)
            p=ROOT/match[2]
            with PILImage.open(p) as img:w,h=img.size
            figure_width=430
            blocks=[Image(str(p),width=figure_width,height=figure_width*h/w),Spacer(1,7)]
            j=i+1
            while j<len(lines) and not lines[j].strip():j+=1
            if j<len(lines) and lines[j].strip().startswith('Figure '):
                blocks.append(Paragraph(inline(lines[j].strip()),styles['MainText']))
                i=j+1
            else:i+=1
            story.append(KeepTogether(blocks));continue
        if line.startswith('#'):
            depth=len(line)-len(line.lstrip('#'));label=line[depth:].strip()
            if label.startswith('4.4 Visualization'):
                story.append(PageBreak())
            style={1:'MainTitle',2:'SubTitle',3:'SectionTitle',4:'SubSectionTitle'}.get(depth,'SubSectionTitle')
            story.append(Paragraph(inline(label),styles[style]));i+=1;continue
        paragraph=[line];i+=1
        while i<len(lines) and lines[i].strip() and not lines[i].strip().startswith(('#','|','$$','![')) and not re.match(r'^\d+\.\s',lines[i].strip()):
            paragraph.append(lines[i].strip());i+=1
        story.append(Paragraph(inline(' '.join(paragraph)),styles['MainText']))
    doc=SimpleDocTemplate(str(OUT),pagesize=(612,792),rightMargin=72,leftMargin=72,topMargin=48,bottomMargin=54,
        title='SEFI:ELEMENTS - A Falsifiable Geometric Reconstruction of the Isotope Chart',author='Jason D. Dutton')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    shutil.copy2(OUT,ROOT/OUT.name)
    print(str(OUT))

if __name__=='__main__':main()
