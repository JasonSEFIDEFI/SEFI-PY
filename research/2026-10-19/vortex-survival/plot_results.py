"""Create a static research summary from recorded runs; no simulation rerun."""
from pathlib import Path
import json
from PIL import Image,ImageDraw,ImageFont
ROOT=Path(__file__).resolve().parent
im=Image.new('RGB',(1500,1040),'white');d=ImageDraw.Draw(im)
font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',20)
small=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',17)
bold=ImageFont.truetype('C:/Windows/Fonts/arialbd.ttf',23)
palette=['#175a7a','#b85b29','#7b4790','#427c45','#b8354f']
def read(name):return json.loads((ROOT/name).read_text())
def plot(box,title,curves,ymin,ymax):
    x0,y0,w,h=box
    values=[v for _,points in curves for _,v in points]
    ymin=min(ymin,min(values));ymax=max(ymax,max(values))
    pad=(ymax-ymin)*.06;ymin-=pad;ymax+=pad
    d.text((x0,y0-40),title,font=bold,fill='#24384b')
    X=lambda t:x0+w*t/24
    Y=lambda y:y0+h-h*(y-ymin)/(ymax-ymin)
    for j in range(5):
        v=ymin+(ymax-ymin)*j/4;yy=Y(v)
        label=f'{v:.4f}' if ymax-ymin<.1 else f'{v:.3g}'
        d.line((x0,yy,x0+w,yy),fill='#dce2e7');d.text((x0-72,yy-10),label,font=small,fill='#34495b')
    for t in [0,8,16,24]:d.text((X(t)-8,y0+h+10),str(t),font=small,fill='#34495b')
    for i,(label,points) in enumerate(curves):
        color=palette[i];d.line([(X(t),Y(y)) for t,y in points],fill=color,width=3)
        d.text((x0+(i%2)*290,y0+h+43+(i//2)*25),label,font=small,fill=color)

curves=[]
for name in ['control','small','medium','strong','phase']:
    s=read(name+'.json')['samples'];base=s[0]['carrier_rms_radius']
    curves.append((name,[(x['time'],x['carrier_rms_radius']/base) for x in s]))
plot((90,100,590,260),'Carrier RMS radius / initial value',curves,.995,1.005)
curves=[]
for m in [1,2,3,4]:
    s=read(f'angular_m{m}.json')['samples'];curves.append((f'angular mode {m}',[(x['time'],x['norm']) for x in s]))
plot((850,100,570,260),'Selected linear angular disturbance gain',curves,.8,1.8)
curves=[]
for name in ['refined','wide']:
    s=read(f'balance_{name}.json')['samples'];q0=s[0]['phi_charge']
    curves.append((name+' charge change',[(x['time'],x['phi_charge']-q0) for x in s]))
plot((90,650,590,240),'Original-field excess charge change',curves,-3,.3)
curves=[]
for name in ['refined','wide']:
    s=read(f'balance_{name}.json')['samples'];E=s[0]['energy']
    curves.append((name+' energy error',[(x['time'],1e10*(x['energy']-E)/abs(E)) for x in s]))
plot((850,650,570,240),'Relative rotating-frame energy error x 1e10',curves,-3,1)
d.text((90,20),'Finite-time moving-ring study | No independent Authorship dynamics',font=bold,fill='#24384b')
d.text((90,1005),'Time is dimensionless. Fixed boundaries; linear angular tests are not nonlinear 3D stability or matter identification.',font=small,fill='#34495b')
im.save(ROOT/'survival-summary.png')
