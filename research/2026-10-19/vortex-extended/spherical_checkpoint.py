"""Summarize only completed spherical searches, spectra and evolutions."""
from vacuum_spherical import *
from run_integrity import eligible
from PIL import Image,ImageDraw,ImageFont

tags=['vacuum_spherical_h0p4','vacuum_spherical_h0p2','vacuum_spherical_h0p1_L60','vacuum_spherical_h0p2_L80','vacuum_spherical_h0p1_L80']
profiles=[json.loads((ROOT/(x+'.json')).read_text()) for x in tags]
assert all(p['converged'] and p['reason']=='residual_tolerance' for p in profiles)
spectra=[json.loads(p.read_text()) for p in sorted(ROOT.glob('*_curvatures.json'))]
assert all(p['reason']=='completed' for p in spectra)
evolutions=[json.loads(p.read_text()) for p in sorted(ROOT.glob('radial_*T100.json'))]
assert len(evolutions)==4 and all(eligible(p) for p in evolutions)
fine=profiles[2];large=profiles[4]
energy=np.array([p['energy'] for p in profiles[:3]])
coarse=np.load(ROOT/'radial_amplitude_T100.npz')['state']
half=np.load(ROOT/'radial_amplitude_halfdt_T100.npz')['state']
summary=dict(scope='Separate lambda=12 spherical vacuum branch; no independent identity/Authorship dynamics or observed-particle identification.',
    profiles=profiles,grid_difference_ratio=float((energy[0]-energy[1])/(energy[1]-energy[2])),
    finest_domain_relative_energy_change=abs(fine['energy']-large['energy'])/large['energy'],
    binding_fraction=1-fine['energy_per_charge']/fine['free_phi_mass'],
    max_temporal_final_state_difference=float(abs(coarse-half).max()),
    minimum_computed_curvature=min(z[k]['eigenvalues'][0] for p in spectra for z in p['sectors'] for k in ['real_amplitudes','phi_phase','sigma_phase']),
    max_eigenpair_residual=max(z[k]['max_residual'] for p in spectra for z in p['sectors'] for k in ['real_amplitudes','phi_phase','sigma_phase']),
    evolutions=[{k:v for k,v in p.items() if k not in ['samples','scope']}|dict(final=p['samples'][-1],maximum_profile_change=max(s['amplitude_profile_relative_change'] for s in p['samples'])) for p in evolutions])
(ROOT/'SPHERICAL_CHECKPOINT.json').write_text(json.dumps(summary,indent=2)+'\n')

image=Image.new('RGB',(1400,1050),'white');draw=ImageDraw.Draw(image)
font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',20);small=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',16);title=ImageFont.truetype('C:/Windows/Fonts/arialbd.ttf',28)
draw.text((40,20),'Spherical vacuum candidate: convergence and finite-time checks',fill='#142337',font=title)
draw.text((40,60),'Separate lambda = 12 branch; Q = 100000. Dimensionless quantities. No particle identification.',fill='#333333',font=font)
colors=['#136b9a','#c96520','#36804c','#874ba2']
def panel(box,heading,xlabel,ylabel,series,xrange,yrange):
    x0,y0,x1,y1=box;draw.text((x0,y0),heading,fill='black',font=font)
    l=x0+85;r=x1-20;t=y0+50;b=y1-55
    def point(x,y):return (l+(x-xrange[0])/(xrange[1]-xrange[0])*(r-l),b-(y-yrange[0])/(yrange[1]-yrange[0])*(b-t))
    for frac in np.linspace(0,1,5):
        x=xrange[0]+frac*(xrange[1]-xrange[0]);y=yrange[0]+frac*(yrange[1]-yrange[0])
        xx,_=point(x,yrange[0]);_,yy=point(xrange[0],y)
        draw.line((l,yy,r,yy),fill='#e4e4e4');draw.text((x0+2,yy-8),f'{y:.4g}',font=small,fill='black');draw.text((xx-12,b+8),f'{x:.3g}',font=small,fill='black')
    draw.line((l,t,l,b,r,b),fill='black',width=2)
    draw.text((l+100,b+30),xlabel,font=small,fill='black');draw.text((l,t-23),ylabel,font=small,fill='black')
    for j,(label,xs,ys) in enumerate(series):
        pts=[point(float(x),float(y)) for x,y in zip(xs,ys)]
        draw.line(pts,fill=colors[j],width=3)
        for xx,yy in pts[::max(1,len(pts)//8)]:draw.ellipse((xx-3,yy-3,xx+3,yy+3),fill=colors[j])
        draw.text((r-170,t+8+j*22),label,font=small,fill=colors[j])
data=np.load(ROOT/'vacuum_spherical_h0p1_L80.npz')
panel((35,115,685,535),'A. Localized field amplitudes','radius','field amplitude',[(n,data['r'],data[k]) for n,k in [('Phi amplitude','f'),('Sigma amplitude','s')]],(0,80),(0,1))
panel((720,115,1370,535),'B. Three-grid energy convergence','grid spacing squared','1e6 x (energy / charge - 1.37736)',[('L = 60',[p['h']**2 for p in profiles[:3]],[1e6*(p['energy_per_charge']-1.37736) for p in profiles[:3]])],(0,.17),(0,9))
sp=next(p for p in spectra if p['source']=='vacuum_spherical_h0p1_L80')
panel((35,575,685,995),'C. Lowest energy curvatures','angular order l','curvature eigenvalue',[(label,range(5),[z[key]['eigenvalues'][0] for z in sp['sectors']]) for label,key in [('real amplitudes','real_amplitudes'),('Phi phase','phi_phase'),('Sigma phase','sigma_phase')]],(0,4),(0,.125))
series=[]
for p in evolutions:
    if p['epsilon'] and p['dt']==.02:series.append((p['kind'],[s['time'] for s in p['samples']],[s['rms_phi_radius'] for s in p['samples']]))
allr=[y for _,_,ys in series for y in ys];pad=(max(allr)-min(allr))*.15
panel((720,575,1370,995),'D. Perturbed physical radial evolution','physical time','RMS Phi radius',series,(0,100),(min(allr)-pad,max(allr)+pad))
draw.text((40,1015),'Curvatures are finite-grid energy tests. Evolution retains radial symmetry and fixed outer boundaries.',font=small,fill='#333333')
image.save(ROOT/'spherical-checkpoint.png')
print(json.dumps({k:v for k,v in summary.items() if k not in ['profiles','evolutions']},indent=2))
