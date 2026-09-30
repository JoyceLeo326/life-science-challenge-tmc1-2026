"""Auditable C-alpha ANM; external source coordinates; no docking or MD."""
import os
for _n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[_n]='2'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import argparse, ast, csv, hashlib, json, platform, time
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import scipy
from scipy.linalg import eigh, subspace_angles
from scipy.spatial import cKDTree, distance
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits, threadpool_info

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,j): p.parent.mkdir(parents=True,exist_ok=True); p.write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8')
def table(p,rows):
    if not rows:return
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,rows[0].keys());w.writeheader();w.writerows(rows)

def read_ca(p,lo,hi):
    rows=[]; keys=set()
    for l in p.read_text(encoding='utf-8').splitlines():
        if l.startswith('ENDMDL'):break
        if not l.startswith('ATOM  ') or l[12:16].strip()!='CA' or l[21]!='A' or l[16] not in [' ','A']:continue
        res=int(l[22:26]);icode=l[26]
        if not lo<=res<=hi:continue
        key=(res,icode)
        if key in keys:raise ValueError('Duplicate CA/alternate locations require explicit resolution: '+str(key))
        keys.add(key)
        rows.append({'residue_number':res,'insertion_code':icode.strip(),'residue_name':l[17:20].strip(),
            'source_B_column':float(l[60:66]),'line':l,'xyz':np.array([float(l[k:k+8]) for k in (30,38,46)])})
    rows.sort(key=lambda r:(r['residue_number'],r['insertion_code']))
    if len(rows)<100:raise ValueError('Insufficient domain CA nodes')
    xyz=np.array([r['xyz'] for r in rows]);assert np.isfinite(xyz).all()
    if cKDTree(xyz).query_pairs(0.1):raise ValueError('Duplicate/near-coincident CA nodes')
    return rows,xyz

def hessian(xyz,cutoff):
    n=len(xyz); h=np.zeros((3*n,3*n));pairs=np.array(sorted(cKDTree(xyz).query_pairs(cutoff)),dtype=int)
    for i,j in pairs:
        d=xyz[j]-xyz[i];b=np.outer(d,d)/np.dot(d,d)
        a=slice(3*i,3*i+3);c=slice(3*j,3*j+3)
        h[a,a]+=b;h[c,c]+=b;h[a,c]-=b;h[c,a]-=b
    adjacency=coo_matrix((np.ones(2*len(pairs)),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(n,n)).tocsr()
    nc,labels=connected_components(adjacency,directed=False)
    return h,pairs,nc,np.asarray(adjacency.sum(axis=1)).ravel()

def prody_reference_hessian(xyz,cutoff,path,expected):
    """Execute only the reviewed, hash-checked upstream buildHessian function.

    This does not claim a full ProDy installation. No upstream imports/parser,
    eigen-solver or extension module is executed. The wrappers implement the
    documented ndarray, constant-gamma, dense/non-KDTree branch only.
    """
    if sha(path)!=expected:raise ValueError('Upstream reference source SHA mismatch')
    tree=ast.parse(path.read_text(encoding='utf-8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='ANMBase')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='buildHessian')
    module=ast.Module(body=[fn],type_ignores=[])
    class Log:
        def __getattr__(self,name):return lambda *a,**k:None
    class Receiver:
        def _reset(self):pass
    def check(x):
        if not isinstance(x,np.ndarray) or x.ndim!=2 or x.shape[1]!=3:raise TypeError('Invalid CA ndarray')
    env={'np':np,'LOGGER':Log(),'checkCoords':check,'checkENMParameters':lambda c,g:(float(c),float(g),lambda d,i,j:float(g))}
    exec(compile(ast.fix_missing_locations(module),str(path),'exec'),env)
    obj=Receiver();env['buildHessian'](obj,xyz,cutoff,1.0,kdtree=False,sparse=False)
    return obj._hessian

def kabsch_rmsd(x,y):
    x=x-x.mean(0);y=y-y.mean(0);u,s,vt=np.linalg.svd(x.T@y);d=np.linalg.det(u@vt)
    rot=u@np.diag([1.,1.,d])@vt
    return float(np.sqrt(np.mean(np.sum((x@rot-y)**2,axis=1))))

def ca_pdb(p,nodes,xyz,remarks):
    p.parent.mkdir(parents=True,exist_ok=True)
    lines=['REMARK '+s for s in remarks]
    for row,x in zip(nodes,xyz):
        l=row['line'];lines.append(l[:30]+''.join(f'{v:8.3f}' for v in x)+l[54:])
    p.write_text('\n'.join(lines)+'\nEND\n',encoding='utf-8')

def analyze(cfg,paths,out):
    start=time.perf_counter();started_utc=datetime.now(timezone.utc).isoformat(); internal=out/'internal'; public=out/'release'
    internal.mkdir(parents=True,exist_ok=True);public.mkdir(parents=True,exist_ok=True)
    summaries=[]; modes_rows=[]; nodes_rows=[]; msf_rows=[]; frame_rows=[]; vector_rows=[]; sensitivity=[]; all_eigen=[]; scope_cache={}
    box=np.asarray(cfg['box']['center']);half=np.asarray(cfg['box']['size'])/2
    for spec in cfg['structures']:
        sid=spec['id'];p=Path(paths[sid]);assert sha(p)==spec['source_sha256'],sid+' input hash mismatch'
        nodes,xyz=read_ca(p,*spec['domain_residue_envelope']);n=len(nodes)
        observed={r['residue_number']:i for i,r in enumerate(nodes)}
        sites=spec['middle_residues'];missing=[r for r in sites if r not in observed]
        if missing:raise ValueError(sid+' incomplete required middle sites: '+str(missing))
        si=np.array([observed[r] for r in sites]);site_ref=distance.pdist(xyz[si])
        serial=np.array([r['residue_number'] for r in nodes]);consecutive=np.where(np.diff(serial)==1)[0]
        adjacent_ref=np.linalg.norm(xyz[consecutive+1]-xyz[consecutive],axis=1)
        nonadjacent=np.abs(serial[:,None]-serial[None,:])>1
        upper=np.triu(nonadjacent,1);ref_dist=distance.squareform(distance.pdist(xyz))
        ref_clashes=int(np.sum((ref_dist<cfg['coarse_diagnostics']['nonadjacent_CA_distance_A'])&upper))
        degree_base=None; primary=None;secondary=None
        ca_pdb(internal/sid/'static_domain_CA.pdb',nodes,xyz,['REFERENCE STATIC SOURCE; chain A TM envelope only','Not a full atomistic docking receptor'])
        for cut in spec.get('cutoffs_A',cfg['cutoffs_A']):
            t=time.perf_counter();h,pairs,nc,degree=hessian(xyz,float(cut))
            hr=prody_reference_hessian(xyz,float(cut),Path(paths['prody_reference_source']),cfg['reference_source_sha256'])
            reference_error=float(np.max(np.abs(h-hr)));del hr
            if reference_error>1e-10:raise ValueError('Independent/upstream Hessians disagree')
            vals,vec=eigh(h,driver='evd',overwrite_a=False,check_finite=False)
            tol=cfg['eigenvalue_zero_tolerance'];zero=int(np.sum(np.abs(vals)<=tol));negative=int(np.sum(vals < -tol))
            keep=np.where(vals>tol)[0][:cfg['retained_nonzero_modes']];v=vec[:,keep];lam=vals[keep]
            # Fix only the mathematically arbitrary sign for reproducible files.
            for k in range(v.shape[1]):
                if v[np.argmax(np.abs(v[:,k])),k]<0:v[:,k]*=-1
            residual=np.linalg.norm(h@v-v*lam[None,:],axis=0)/np.linalg.norm(h,ord='fro')
            orth=float(np.max(np.abs(v.T@v-np.eye(len(keep)))))
            rigid=[];centered=xyz-xyz.mean(0)
            for axis in np.eye(3):rigid.append(np.tile(axis,(n,1)).ravel())
            for axis in np.eye(3):rigid.append(np.cross(np.tile(axis,(n,1)),centered).ravel())
            q,_=np.linalg.qr(np.array(rigid).T)
            rigid_res=float(np.max(np.linalg.norm(h@q,axis=0)/np.linalg.norm(h,ord='fro')))
            rigid_overlap=float(np.max(np.abs(q.T@v)))
            valid=nc==1 and zero==6 and negative==0 and float(residual.max())<1e-10 and orth<1e-10 and rigid_res<1e-10
            summary={'structure_id':sid,'cutoff_A':cut,'nodes':n,'Hessian_dimension':3*n,'Hessian_memory_MiB':round(h.nbytes/2**20,3),
                'graph_edges':len(pairs),'connected_components':nc,'min_degree':int(degree.min()),'max_degree':int(degree.max()),
                'zero_mode_count':zero,'negative_eigenvalues':negative,'eigenvalue_zero_tolerance':tol,
                'lowest_positive_eigenvalue':float(lam[0]),'largest_eigenvalue':float(vals[-1]),'retained_nonzero_modes':len(keep),
                'symmetry_max_abs':float(np.max(np.abs(h-h.T))),'upstream_Hessian_max_abs_error':reference_error,
                'eigenpair_max_relative_residual':float(residual.max()),'orthogonality_max_abs_error':orth,
                'rigid_body_max_relative_residual':rigid_res,'nonzero_modes_rigid_basis_max_abs_overlap':rigid_overlap,
                'ANM_numerical_QC_pass':valid,'elapsed_seconds':round(time.perf_counter()-t,3),'source_sha256':spec['source_sha256']}
            summaries.append(summary);print(sid,cut,'nodes',n,'components',nc,'zero',zero,'QC',valid,flush=True)
            for idx,value in enumerate(vals):all_eigen.append({'structure_id':sid,'cutoff_A':cut,'eigenvalue_index_0based':idx,'eigenvalue_relative_units':float(value),'zero_by_frozen_tolerance':abs(float(value))<=tol})
            for k,e in enumerate(lam):modes_rows.append({'structure_id':sid,'cutoff_A':cut,'nonzero_mode_1based':k+1,'eigenvalue_relative_units':float(e),'inverse_eigenvalue_relative_units':float(1/e),'relative_residual':float(residual[k])})
            msf=((v.reshape(n,3,-1)**2)/lam[None,None,:]).sum(axis=(1,2));normalized=msf/msf.mean()
            for i,row in enumerate(nodes):
                msf_rows.append({'structure_id':sid,'cutoff_A':cut,'chain':'A','residue_number':row['residue_number'],'MSF_low20_relative_units':float(msf[i]),'MSF_low20_normalized_mean1':float(normalized[i]),'middle_site':row['residue_number'] in sites,'source_B_column_not_mobility':row['source_B_column']})
            np.savez_compressed(public/f'{sid}_cutoff{cut:g}_modes.npz',eigenvalues=lam,eigenvectors=v,residue_numbers=serial)
            if cut==cfg['primary_cutoff_A']:
                primary=(v,lam,normalized,valid);degree_base=degree
            else:secondary=(v,lam,normalized,valid)
            del h,vec
        for i,r in enumerate(nodes):nodes_rows.append({'structure_id':sid,'node_index_0based':i,'chain':'A','residue_number':r['residue_number'],'insertion_code':r['insertion_code'],'residue_name':r['residue_name'],'source_B_column_not_mobility':r['source_B_column'],'degree_primary15A':int(degree_base[i]),'middle_site':r['residue_number'] in sites})
        v,lam,norm,valid=primary;scope_cache[sid]={'serial':serial,'vectors':v,'MSF':norm,'valid':valid}
        if secondary is not None:
            s,sl,snorm,svalid=secondary
            rmsip=float(np.sqrt(np.sum((v.T@s)**2)/v.shape[1]));angles=np.degrees(subspace_angles(v,s))
            sensitivity.append({'structure_id':sid,'cutoffs_compared':'15A vs 13A','low20_RMSIP':rmsip,'largest_principal_angle_deg':float(angles.max()),'MSF_low20_normalized_Spearman':float(spearmanr(norm,snorm).statistic),'both_networks_numerical_QC_pass':valid and svalid})
        if not valid or not spec.get('generate_traversals',True):continue
        # No stochastic sample, candidate scores or docking results are used.
        for mode in cfg['conformation_modes_1based']:
            unit=v[:,mode-1].reshape(n,3)
            for amp in cfg['CA_RMS_displacement_levels_A']:
                for sign in [-1,1]:
                    delta=sign*amp*np.sqrt(n)*unit;new=xyz+delta
                    fid=f'{sid}_mode{mode}_{"minus" if sign<0 else "plus"}_rms{amp:g}A'
                    actual=float(np.sqrt(np.mean(np.sum(delta**2,axis=1))))
                    maximum=float(np.linalg.norm(delta,axis=1).max())
                    adjacency=np.linalg.norm(new[consecutive+1]-new[consecutive],axis=1)
                    adj_change=float(np.max(np.abs(adjacency-adjacent_ref)))
                    site_change=distance.pdist(new[si])-site_ref
                    max_site=float(np.max(np.abs(site_change)))
                    all_dist=distance.squareform(distance.pdist(new));clashes=int(np.sum((all_dist<cfg['coarse_diagnostics']['nonadjacent_CA_distance_A'])&upper))
                    margin=float((half-np.abs(new[si]-box)).min())
                    diagnostic=adj_change<=cfg['coarse_diagnostics']['max_adjacent_CA_distance_change_A'] and max_site<=cfg['coarse_diagnostics']['max_middle_CA_pair_distance_change_A'] and clashes<=ref_clashes and margin>=0
                    frame_rows.append({'structure_id':sid,'frame_id':fid,'mode_1based':mode,'direction':sign,'declared_CA_RMS_displacement_A':amp,
                        'actual_CA_RMS_displacement_A':actual,'CA_best_fit_RMSD_A':kabsch_rmsd(new,xyz),'max_CA_displacement_A':maximum,
                        'middle_CA_RMS_displacement_A':float(np.sqrt(np.mean(np.sum(delta[si]**2,axis=1)))),
                        'max_middle_CA_displacement_A':float(np.linalg.norm(delta[si],axis=1).max()),'middle_CA_pair_delta_RMS_A':float(np.sqrt(np.mean(site_change**2))),
                        'max_middle_CA_pair_delta_abs_A':max_site,'max_consecutive_CA_distance_delta_abs_A':adj_change,
                        'min_consecutive_CA_distance_A':float(adjacency.min()),'max_consecutive_CA_distance_A':float(adjacency.max()),
                        'nonadjacent_CA_pairs_under3A':clashes,'reference_nonadjacent_CA_pairs_under3A':ref_clashes,'middle15_CA_min_fixed_box_margin_A':margin,
                        'coarse_reference_diagnostic_pass':diagnostic,'illustrative_primary_amplitude':amp==cfg['primary_visualization_amplitude_A'],
                        'eligible_for_full_atom_docking':False,'reason':'C-alpha traversal only; full atom reconstruction/relaxation and independent geometry gate not performed'})
                    for i,row in enumerate(nodes):vector_rows.append({'frame_id':fid,'structure_id':sid,'residue_number':row['residue_number'],'dx_A':float(delta[i,0]),'dy_A':float(delta[i,1]),'dz_A':float(delta[i,2])})
                    ca_pdb(internal/sid/(fid+'.pdb'),nodes,new,['DETERMINISTIC LINEAR ANM TRAVERSAL; not MD/time/thermal sampling','CA domain only; no full-atom docking eligibility','Input reference preserved separately; imposed CA RMS '+str(amp)+' A'])
        np.savez_compressed(internal/sid/'reference_and_low20.npz',reference_CA=xyz,site_indices=si,eigenvectors=v,eigenvalues=lam,residue_numbers=serial)
    boundary=[]
    for pair in cfg.get('boundary_comparisons',[]):
        a=scope_cache[pair['domain']];b=scope_cache[pair['full_chain']]
        ix=np.array([int(np.where(b['serial']==r)[0][0]) for r in a['serial']])
        restricted=b['vectors'].reshape(len(b['serial']),3,-1)[ix].reshape(3*len(ix),-1)
        rank=int(np.linalg.matrix_rank(restricted));basis,_=np.linalg.qr(restricted)
        boundary.append({'domain_scope':pair['domain'],'full_chain_scope':pair['full_chain'],'common_CA_nodes':len(ix),
            'restricted_full_low20_rank':rank,'domain_vs_restricted_full_low20_RMSIP':float(np.sqrt(np.sum((a['vectors'].T@basis)**2)/a['vectors'].shape[1])),
            'domain_vs_full_restricted_MSF_Spearman':float(spearmanr(a['MSF'],b['MSF'][ix]).statistic),
            'both_networks_numerical_QC_pass':a['valid'] and b['valid'],'interpretation':'Boundary-condition sensitivity; not evidence of physiological fluctuation amplitude'})
    table(public/'network_numerical_QC.csv',summaries);table(public/'all_eigenvalues.csv',all_eigen);table(public/'low20_modes.csv',modes_rows)
    table(public/'CA_node_index.csv',nodes_rows);table(public/'low20_relative_MSF.csv',msf_rows)
    table(public/'cutoff_sensitivity.csv',sensitivity);table(public/'conformation_diagnostics.csv',frame_rows);table(public/'conformation_displacements.csv',vector_rows)
    table(public/'boundary_sensitivity.csv',boundary)
    receipt={'started_utc':started_utc,'completed_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':round(time.perf_counter()-start,3),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        'implementation':'transparent SciPy ANM, independently matched to reviewed hash-checked ProDy v2.6.1 buildHessian source; full ProDy package not installed',
        'network_count':len(summaries),'numerically_valid_networks':sum(r['ANM_numerical_QC_pass'] for r in summaries),'CA_traversal_count':len(frame_rows),
        'primary_illustrative_frames':sum(r['illustrative_primary_amplitude'] for r in frame_rows),'full_atom_docking_eligible_frames':0,
        'maximum_threads_per_BLAS_library':2,'threadpool_actual':threadpool_info(),'random_sampling_used':False,'MD_performed':False,'docking_performed':False}
    dump(public/'ACTUAL_RUN_RECEIPT.json',receipt);print(json.dumps(receipt,ensure_ascii=False),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--protocol',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    cfg=json.loads(a.protocol.read_text(encoding='utf-8'));paths=json.loads(a.inputs.read_text(encoding='utf-8'))
    with threadpool_limits(limits=2):analyze(cfg,paths,a.out.resolve())
if __name__=='__main__':main()
