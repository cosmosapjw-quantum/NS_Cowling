#!/usr/bin/env python3
"""Plot actual common-density SPG candidate rows, without interpolation/resampling."""
from __future__ import annotations
import argparse,gzip,hashlib,json,platform
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

LABELS={
    'AME2408_candidate':'AME2408 candidate (baseline)',
    'geometric_area__nine':'Geometric area, 9 nuclei',
    'geometric_area__eight_measured':'Geometric area, 8 measured nuclei',
    'fixed_symmetric_area__nine':'Fixed symmetric area, 9 nuclei',
    'fixed_symmetric_area__eight_measured':'Fixed symmetric area, 8 measured nuclei',
}
COLORS=['#202020','#0072B2','#56B4E9','#D55E00','#CC79A7']
STYLES=['-','-','--','-','--']
CAPTION=('Actual common-density reoptimized p = 3 candidate prescriptions.\n'
         'Lines connect computed rows; no interpolation or resampling.\n'
         'Not a confidence band or recovery of the original SPG surface fit.')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path,help='common_grid129 directory containing model/label data')
    ap.add_argument('output',type=Path,help='output prefix, without extension')
    ap.add_argument('--model',choices=[f'SPG_M{i}' for i in range(1,6)],default='SPG_M4')
    args=ap.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    paths={label:args.input/args.model/label/'common_grid_attempts.json.gz' for label in LABELS}
    series={}
    for label,path in paths.items():
        with gzip.open(path,'rt') as stream:rows=json.load(stream)
        if not all(row['best']['status']=='STATIONARY_CANDIDATE' for row in rows):
            raise RuntimeError('plot refuses to bridge nonstationary rows: '+label)
        series[label]={
            'n':np.array([row['n_B'] for row in rows]),
            'P':np.array([row['best']['pressure_MeV_fm3'] for row in rows]),
            'Ye':np.array([row['best']['composition']['Ye'] for row in rows]),
            'Z':np.array([row['best']['composition']['Zav'] for row in rows]),
        }
    baseline=series['AME2408_candidate']
    for label,data in series.items():
        if not np.array_equal(data['n'],baseline['n']):
            raise RuntimeError('common density grid is not byte-value identical: '+label)
        if not all(np.isfinite(x).all() for x in data.values()):
            raise RuntimeError('nonfinite row: '+label)
        if not (np.all(data['n']>0) and np.all(np.diff(data['n'])>0) and np.all(data['P']>0)):
            raise RuntimeError('nonpositive pressure or unordered density: '+label)
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':10,
                         'axes.labelsize':11,'axes.titlesize':12,
                         'pdf.fonttype':42,'ps.fonttype':42}):
        fig,axes=plt.subplots(2,2,figsize=(10.8,8.2))
        for (label,data),color,style in zip(series.items(),COLORS,STYLES):
            width=2.1 if label=='AME2408_candidate' else 1.65
            kwargs=dict(label=LABELS[label],color=color,linestyle=style,linewidth=width)
            axes[0,0].plot(data['n'],data['P'],**kwargs)
            axes[0,1].plot(data['n'],data['Ye'],**kwargs)
            axes[1,0].plot(data['n'],data['Z'],**kwargs)
            axes[1,1].plot(data['n'],100*(data['P']/baseline['P']-1),**kwargs)
        axes[0,0].set_yscale('log')
        axes[0,0].set_ylabel(r'$P\;[\mathrm{MeV\,fm}^{-3}]$')
        axes[0,1].set_ylabel(r'$Y_e$')
        axes[1,0].set_ylabel(r'$Z$')
        axes[1,1].set_ylabel(r'$100(P/P_{\rm baseline}-1)\;[\%]$')
        titles=['(a) Pressure','(b) Electron fraction','(c) Nuclear charge','(d) Pressure change']
        for ax,title in zip(axes.flat,titles):
            ax.set_xscale('log');ax.set_xlabel(r'$n_B\;[\mathrm{fm}^{-3}]$')
            ax.set_title(title,loc='left');ax.grid(True,which='major',alpha=.25,linewidth=.6)
            ax.set_xlim(baseline['n'][0],baseline['n'][-1])
        handles,labels=axes[0,0].get_legend_handles_labels()
        fig.suptitle(args.model.replace('_',' ')+'  |  Reoptimized surface-prescription response',
                     fontsize=15,y=.985)
        fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.948),
                   ncol=2,frameon=False,fontsize=10,columnspacing=2.5)
        fig.subplots_adjust(left=.09,right=.98,bottom=.16,top=.82,wspace=.25,hspace=.35)
        fig.text(.5,.025,CAPTION,ha='center',va='bottom',fontsize=9,color='#444444')
        for suffix in ('.png','.pdf'):
            path=args.output.with_suffix(suffix)
            fig.savefig(path,dpi=180,facecolor='white',
                        metadata={'Creator':'NS_Cowling reproducible scientific plotter'} if suffix=='.pdf' else None)
        plt.close(fig)
    result={'schema':'spg-actual-common-density-plot-v1','model':args.model,
            'caption':CAPTION,'row_count_per_prescription':len(baseline['n']),
            'common_support_fm3':[float(baseline['n'][0]),float(baseline['n'][-1])],
            'source_calibrated':False,'confidence_band':False,'resampled_or_fitted':False,
            'source_sha256':{str(path):sha(path) for path in paths.values()},
            'plotter_sha256':sha(__file__),
            'software':{'python':platform.python_version(),'numpy':np.__version__,
                        'matplotlib':matplotlib.__version__},
            'output_sha256':{str(args.output.with_suffix(suffix)):sha(args.output.with_suffix(suffix))
                             for suffix in ('.png','.pdf')}}
    args.output.with_suffix('.receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('model','row_count_per_prescription','common_support_fm3','output_sha256')},indent=2))
if __name__=='__main__':main()
