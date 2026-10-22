from pathlib import Path
import csv,math
from reportlab.graphics.shapes import Drawing,Line,PolyLine,String,Rect
from reportlab.graphics import renderSVG
from reportlab.lib import colors
from phase_model import critical_detuning
P=Path(__file__).resolve().parent;rows=list(csv.DictReader((P/'results/sampled_trajectories.csv').open()))
d=Drawing(1000,430);navy=colors.HexColor('#17314d');blue=colors.HexColor('#267ab8');red=colors.HexColor('#b74b37');gray=colors.HexColor('#dce3e9')
d.add(Rect(0,0,1000,430,fillColor=colors.white,strokeColor=None))
d.add(String(40,399,'Classical phase stability with recoverable actuator degradation',fontSize=18,fillColor=navy))
d.add(String(40,375,'Specified engineering model; no hardware data or SEFI-specific advantage established',fontSize=11,fillColor=navy))
def panel(x,y,w,h,xvals,yvals,xmax,ymax,title,xlabel,ylabel):
    for t in xvals:
        xx=x+t/xmax*w;d.add(Line(xx,y,xx,y+h,strokeColor=gray));d.add(String(xx,y-18,str(t),textAnchor='middle',fontSize=10))
    for v in yvals:
        yy=y+v/ymax*h;d.add(Line(x,yy,x+w,yy,strokeColor=gray));d.add(String(x-9,yy-3,str(v),textAnchor='end',fontSize=10))
    d.add(Line(x,y,x+w,y,strokeColor=navy));d.add(Line(x,y,x,y+h,strokeColor=navy))
    d.add(String(x,y+h+23,title,fontSize=12,fillColor=navy));d.add(String(x,y+h+6,ylabel,fontSize=10));d.add(String(x+w/2,y-38,xlabel,fontSize=11,textAnchor='middle'))
def plot(points,x,y,w,h,xmax,ymax,color):d.add(PolyLine([(x+a/xmax*w,y+b/ymax*h) for a,b in points],strokeColor=color,strokeWidth=2,fillColor=None))
subset=[r for r in rows if float(r['detuning'])==1.3];maxphi=max(float(r['phase']) for r in subset);top=max(10,math.ceil(maxphi/10)*10)
panel(65,105,375,210,[0,30,60,90,120],list(range(0,top+1,10)),120,top,'A. Detuning 1.3 rad/s','Time (s)','Unwrapped relative phase (rad)')
for model,col in [('constant_gain',blue),('gain_degradation',red)]:plot([(float(r['t']),float(r['phase'])) for r in subset if r['model']==model],65,105,375,210,120,top,col)
panel(575,105,375,210,[0,1,2,3,4],[0,.5,1,1.5,2],4,2.1,'B. Derived stationary-lock boundary','Activity coefficient alpha (s)','Critical detuning (rad/s)')
plot([(i/50,critical_detuning(2,i/50,1)) for i in range(201)],575,105,375,210,4,2.1,red)
plot([(0,2),(4,2)],575,105,375,210,4,2.1,blue)
for x,col,label in [(70,blue,'Ideal constant gain'),(280,red,'Gain-degradation model')]:
    d.add(Line(x,44,x+26,44,strokeColor=col,strokeWidth=2));d.add(String(x+33,40,label,fontSize=11))
d.add(String(575,42,'K0 = 2 rad/s; beta = 1/s',fontSize=11,fillColor=navy))
(P/'figures').mkdir(exist_ok=True);renderSVG.drawToFile(d,str(P/'figures/phase_stability.svg'));print('Scientific SVG figure generated')
