"""Plot only frozen numeric ANM evidence; no structure coordinates needed."""
import os
for n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[n]='2'
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def rows(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def main():
    p=argparse.ArgumentParser();p.add_argument('--release',type=Path,required=True);a=p.parse_args();r=a.release
    out=r/'figures';out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':140,'savefig.dpi':180})
    registry=[]
    def save(fig,name,sources,note):
        fig.savefig(out/(name+'.png'),bbox_inches='tight');fig.savefig(out/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
        registry.append({'figure':name,'source_csv_or_npz':';'.join(sources),'caption':note})
    q=rows(r/'network_numerical_QC.csv');labels=[x['structure_id'].replace('M_PUB_ALL','Mouse').replace('H_AF','Human').replace('_full_chainA',' full chain').replace('_domain',' domain')+' '+str(int(float(x['cutoff_A'])))+' A' for x in q]
    fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained');x=np.arange(len(q))
    ax[0].barh(x,[float(z['lowest_positive_eigenvalue']) for z in q],color=['#247C95','#72B7C3','#9AABB1']*2)
    ax[0].set_xscale('log');ax[0].set_yticks(x,labels);ax[0].invert_yaxis();ax[0].set_xlabel('Lowest positive eigenvalue (relative units)');ax[0].axvline(1e-6,color='#A33F3F',ls='--',lw=1);ax[0].set_title('A. Lowest positive eigenvalue')
    ax[1].barh(x,[float(z['eigenpair_max_relative_residual'])*1e17 for z in q],color='#247C95');ax[1].set_xlim(0,7);ax[1].set_xticks([0,2,4,6]);ax[1].set_yticks(x,labels);ax[1].invert_yaxis();ax[1].set_xlabel(r'Relative residual ($\times 10^{-17}$)');ax[1].set_title('B. Eigenpair residuals')
    fig.suptitle('ANM numerical checks',fontsize=13)
    save(fig,'F1_network_and_solver_QC',['network_numerical_QC.csv'],'One connected component and six rigid-body zero modes per network. Relative eigenvalues are not physical frequencies.')
    m=rows(r/'low20_relative_MSF.csv');boundary=rows(r/'boundary_sensitivity.csv');cut=rows(r/'cutoff_sensitivity.csv')
    fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for j,(sid,name) in enumerate([('H_AF','Human'),('M_PUB_ALL','Mouse')]):
        d={c:[z for z in m if z['structure_id']==sid+'_domain' and float(z['cutoff_A'])==c] for c in [15.,13.]}
        for c,color in [(15.,'#175D79'),(13.,'#A6632C')]:axes[j,0].plot([int(z['residue_number']) for z in d[c]],[float(z['MSF_low20_normalized_mean1']) for z in d[c]],color=color,lw=1.4,label=f'{int(c)} A cutoff')
        axes[j,0].set_title(name+' domain: cutoff sensitivity');axes[j,0].set_ylabel('Low20 MSF / node-set mean')
        full={int(z['residue_number']):float(z['MSF_low20_normalized_mean1']) for z in m if z['structure_id']==sid+'_full_chainA'}
        xx=[int(z['residue_number']) for z in d[15.]];yy=np.array([full[i] for i in xx]);yy/=yy.mean()
        axes[j,1].plot(xx,[float(z['MSF_low20_normalized_mean1']) for z in d[15.]],label='Domain network',color='#175D79',lw=1.4)
        axes[j,1].plot(xx,yy,label='Full chain, restricted to domain / mean',color='#AA574E',lw=1.2)
        b=next(z for z in boundary if z['domain_scope']==sid+'_domain');axes[j,1].set_title(name+f": boundary (rho={float(b['domain_vs_full_restricted_MSF_Spearman']):.3f})");axes[j,1].set_ylabel('Restricted-profile MSF / mean')
        for a2 in axes[j]:a2.set_xlabel('Source chain-A residue number');a2.grid(alpha=.12)
    handles=[axes[0,0].lines[0],axes[0,0].lines[1],axes[0,1].lines[1]]
    fig.legend(handles,['Domain, 15 A','Domain, 13 A','Full chain projected onto domain, 15 A'],loc='outside lower center',ncol=3,fontsize=9)
    fig.suptitle('ANM relative mobility and boundary sensitivity',fontsize=13)
    save(fig,'F2_relative_MSF_and_boundary',['low20_relative_MSF.csv','boundary_sensitivity.csv','cutoff_sensitivity.csv'],'Weighted low20 ANM covariance diagonal, normalized to mean one. Full-chain curves are renormalized on the common domain. No temperature or MD RMSF interpretation.')
    f=rows(r/'conformation_diagnostics.csv');fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for sid,col in [('H_AF_domain','#175D79'),('M_PUB_ALL_domain','#A6632C')]:
        for mode in [1,2]:
            for sign in [-1,1]:
                z=[v for v in f if v['structure_id']==sid and int(v['mode_1based'])==mode and int(v['direction'])==sign];z.sort(key=lambda v:float(v['declared_CA_RMS_displacement_A']))
                xx=[float(v['declared_CA_RMS_displacement_A']) for v in z];ls='-' if sign>0 else '--'
                lab=('Human' if sid.startswith('H') else 'Mouse')+f' mode {mode} '+('+' if sign>0 else '-')
                for a2,key in zip(axes.ravel()[:3],['max_CA_displacement_A','max_consecutive_CA_distance_delta_abs_A','max_middle_CA_pair_delta_abs_A']):a2.plot(xx,[float(v[key]) for v in z],ls=ls,color=col,alpha=.8,lw=1.2,marker='o' if mode==1 else 's',ms=3,label=lab)
    axes[0,0].set_ylabel('Maximum CA displacement (A)');axes[0,0].set_title('A. Local movement can exceed the imposed global RMS')
    axes[0,1].set_ylabel('Max adjacent-CA distance change (A)');axes[0,1].axhline(.2,color='#B84C4C',lw=1,ls=':');axes[0,1].set_title('B. Frozen diagnostic flag at 0.20 A')
    axes[1,0].set_ylabel('Max middle15 CA-pair distance change (A)');axes[1,0].set_title('C. Middle-site pair geometry changes remain small')
    am=[.25,.5,1.];xx=np.arange(3);w=.35
    for dx,sid,c,label in [(-w/2,'H_AF_domain','#175D79','Human'),(w/2,'M_PUB_ALL_domain','#A6632C','Mouse')]:
        yy=[sum(z['coarse_reference_diagnostic_pass']=='True' for z in f if z['structure_id']==sid and float(z['declared_CA_RMS_displacement_A'])==amp) for amp in am]
        axes[1,1].bar(xx+dx,yy,width=w,color=c,label=label)
    axes[1,1].set_xticks(xx,[str(v) for v in am]);axes[1,1].set_ylim(0,4.5);axes[1,1].set_yticks(range(5));axes[1,1].set_ylabel('Coarse diagnostic pass / 4 frames');axes[1,1].set_title('D. All frames retained, including failed flags');axes[1,1].legend()
    for a2 in axes.ravel():a2.set_xlabel('Imposed global CA RMS displacement (A)');a2.grid(axis='y',alpha=.12)
    for a2 in axes.ravel()[:3]:a2.set_xticks(am)
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=4,fontsize=8)
    fig.suptitle('ANM traversal geometry diagnostics',fontsize=13)
    save(fig,'F3_traversal_geometry_diagnostics',['conformation_diagnostics.csv'],'First two modes, both signs, all predeclared amplitudes retained. The 0.5 A display setting has 2/4 human and 4/4 mouse coarse diagnostic passes. None is full-atom docking-ready.')
    nrows=rows(r/'CA_node_index.csv');fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained');local=[]
    for j,(sid,name) in enumerate([('H_AF_domain','Human'),('M_PUB_ALL_domain','Mouse')]):
        z=np.load(r/f'{sid}_cutoff15_modes.npz');vv=z['eigenvectors'];nn=z['residue_numbers'];meta={int(x['residue_number']):x for x in nrows if x['structure_id']==sid}
        for k in [0,1]:
            a2=axes[j,k];amp=.5*np.sqrt(len(nn))*np.linalg.norm(vv[:,k].reshape(len(nn),3),axis=1);sites=[i for i,x in enumerate(nn) if meta[int(x)]['middle_site']=='True'];imax=int(np.argmax(amp))
            a2.plot(nn,amp,color='#175D79' if j==0 else '#A6632C',lw=1.3);a2.scatter(nn[sites],amp[sites],s=18,color='#D47922',label='Middle15 CA nodes',zorder=5)
            a2.set_title(name+f' mode {k+1}; max at residue {nn[imax]}');a2.set_xlabel('Source chain-A residue number');a2.set_ylabel('CA displacement magnitude (A)');a2.grid(alpha=.12)
            local.append({'structure_id':sid,'mode_1based':k+1,'imposed_global_CA_RMS_displacement_A':.5,'maximum_CA_displacement_A':float(amp[imax]),'maximum_displacement_residue':int(nn[imax]),'source_B_column_at_max':meta[int(nn[imax])]['source_B_column_not_mobility'],'middle15_fraction_of_squared_mode_norm':float(np.sum(vv[:,k].reshape(len(nn),3)[sites]**2))})
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=1,fontsize=9)
    fig.suptitle('Low-mode displacement localization',fontsize=13)
    save(fig,'F4_low_mode_displacement_localization',['CA_node_index.csv','H_AF_domain_cutoff15_modes.npz','M_PUB_ALL_domain_cutoff15_modes.npz'],'Magnitude of fixed 0.5 A CA RMS traversal, identical for +/- direction. Residue-wise localization does not identify channel gating or drug action.')
    with (r/'mode_localization.csv').open('w',encoding='utf-8',newline='') as h:w=csv.DictWriter(h,local[0].keys());w.writeheader();w.writerows(local)
    with (r/'FIGURE_SOURCE_INDEX.csv').open('w',encoding='utf-8',newline='') as h:w=csv.DictWriter(h,registry[0].keys());w.writeheader();w.writerows(registry)
    print('Rendered four actual scientific PNG/PDF figures with CSV/NPZ provenance')
if __name__=='__main__':main()
