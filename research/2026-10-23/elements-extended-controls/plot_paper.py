from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parent

def main():
    data=list(csv.DictReader((ROOT/'results/evaluated_ground_states.csv').open()))
    data=[r for r in data if int(r['Z'])>0 and int(r['N'])>0 and int(r['A'])>=16 and r['mass_estimated']=='False']
    metrics=list(csv.DictReader((ROOT/'results/mass_metrics.csv').open()))
    plt.rcParams.update({'font.size':10})
    fig,ax=plt.subplots(figsize=(6.5,4.3))
    dots=ax.scatter([int(r['N']) for r in data],[int(r['Z']) for r in data],
                    c=[float(r['binding_per_A_keV'])/1000 for r in data],s=10,cmap='viridis')
    ax.set(xlabel='Neutron number N',ylabel='Proton number Z',title='Measured AME2020 mass subset')
    fig.colorbar(dots,ax=ax,label='Binding energy per nucleon (MeV)')
    fig.tight_layout();fig.savefig(ROOT/'results/isotope_chart.png',dpi=200);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.5,4.3));x=np.arange(3)
    pairs=[('SEMF','Liquid-drop baseline'),('SEMF_known_shell_features','Baseline + known shells'),
           ('SEFI_geometric_prototype','Geometric prototype'),('conventional_matched_kernel','Matched conventional kernel')]
    for j,(method,label) in enumerate(pairs):
        values=[next(float(r['rmse']) for r in metrics if r['split']==s and r['method']==method)
                for s in ['random_4200','neutron_frontier','proton_extrapolation']]
        ax.bar(x+(j-1.5)*.2,values,.2,label=label)
    ax.set(xticks=x,xticklabels=['Random','Neutron frontier','Proton extrapolation'],
           ylabel='Binding-energy RMSE (MeV)',title='Held-out errors: equivalent methods tie')
    ax.legend(fontsize=8,loc='upper left');fig.tight_layout()
    fig.savefig(ROOT/'results/prediction_errors.png',dpi=200);plt.close(fig)

if __name__=='__main__':main()
