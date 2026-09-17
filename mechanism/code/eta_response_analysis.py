"""Conditional eta response interpretation. Cached inputs read-only; no plotting."""
from pathlib import Path
from datetime import datetime, timezone
from itertools import combinations
import hashlib, json, re, os, argparse
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, eye, load_npz
from scipy.sparse.linalg import eigs, expm_multiply
from scipy.sparse.csgraph import connected_components
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits
threadpool_limits(1)

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('DMWA_OUTPUT_ROOT', str(PACKAGE)))
OUT = ROOT/'data/eta'
REPORT = ROOT/'reports'
INPUT = Path(os.environ.get('DMWA_INPUT_DIR', str(PACKAGE/'data/eta')))
for p in [OUT, REPORT]: p.mkdir(parents=True, exist_ok=True)
GAINS = (.5,.8,.9)
TIMES = np.arange(1601)*.05
PHASES = ['Pre','Early','Mid','Late','Post']
SLOPE_TOL = .0005
SOLVER_TOL = 1e-7
sources = []

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def track(path,role):
    sources.append(dict(role=role,path=Path(path).name,sha256=sha(path),bytes=Path(path).stat().st_size))

def save(data,name):
    if not isinstance(data,pd.DataFrame): data=pd.DataFrame(data)
    data.to_csv(OUT/f'{name}.csv',index=False)

def leading(W):
    S=csr_matrix(W)
    if S.shape[0]<=8:
        val,vec=np.linalg.eig(S.toarray()); j=np.argmax(val.real); value=val[j]; v=vec[:,j]
    else:
        val,vec=eigs(S,k=1,which='LR',v0=np.ones(S.shape[0]),tol=1e-13,maxiter=10000)
        value=val[0];v=vec[:,0]
    residual=np.linalg.norm(S@v-value*v)/(np.linalg.norm(v)*max(abs(value),1e-15))
    assert abs(value.imag)<1e-9
    return float(value.real),float(residual)

def parts(W):
    W=np.asarray(W,float)
    if not np.isfinite(W).all() or (W<0).any(): raise ValueError('Nonfinite or negative W')
    R=float(W.sum()/len(W));rho,res=leading(W)
    if R<=0 or rho<=0: raise ValueError('R or rho not positive')
    return dict(m=len(W),R=R,rho=rho,eta=R/rho,eigen_residual=res,
                components=connected_components(csr_matrix(W),directed=True,connection='strong')[0])

def make_scaffold(edges,n=5):
    adj=np.zeros((n,n),int)
    for i,j in edges:adj[i-1,j-1]=adj[j-1,i-1]=1
    groups=[tuple((i-1,j-1)) for i,j in edges]
    groups += [e for e in combinations(range(n),3) if all(adj[i,j] for i,j in combinations(e,2))]
    H=np.zeros((len(groups),n),int)
    for i,e in enumerate(groups):H[i,list(e)]=1
    overlap=H@H.T; np.fill_diagonal(overlap,0)
    sizes=H.sum(1); degree=(overlap>0).sum(1)
    den=degree[:,None]+degree[None,:]
    bias=np.divide(degree[None,:],den,out=np.full(den.shape,.5),where=den>0)
    W=overlap/sizes[:,None]*np.where(sizes[:,None]==sizes[None,:],bias,1)
    return H,W

def slopes(logmass):
    out={}
    for name,lo,hi in [('20_40',20,40),('40_80',40,80)]:
        keep=(TIMES>=lo)&(TIMES<=hi)
        out[name]=float(np.polyfit(TIMES[keep],logmass[keep],1)[0])
    return out

def responses(W,p,metadata,example=False):
    A=csr_matrix(W/p['R']);m=len(W);rhoA=p['rho']/p['R'];x0=np.ones(m)/m
    rate_rows=[]; curves=[]
    for g in GAINS:
        s=g*rhoA-1
        # Spectral shift stabilizes relative numerical accuracy after long decay.
        K=g*(A.T.tocsr()-rhoA*eye(m,format='csr'))
        Y=expm_multiply(K,x0,start=0,stop=80,num=len(TIMES),endpoint=True,traceA=-g*rhoA*m)
        ode=solve_ivp(lambda t,y:K@y,(0,80),x0,method='DOP853',t_eval=TIMES,rtol=1e-11,atol=1e-13)
        if not ode.success: raise RuntimeError(ode.message)
        Z=ode.y.T
        err=float(np.max(np.linalg.norm(Y-Z,axis=1)/np.maximum(np.linalg.norm(Y,axis=1),1e-30)))
        logmass=np.log(Y.sum(1))+s*TIMES
        logmass_ode=np.log(Z.sum(1))+s*TIMES
        sl=slopes(logmass);sl_ode=slopes(logmass_ode)
        slope_err=max(abs(sl[k]-s) for k in sl)
        solver_slope_err=max(abs(sl[k]-sl_ode[k]) for k in sl)
        direct_err=np.nan
        if example:
            J=-eye(m,format='csr')+g*A.T
            direct=solve_ivp(lambda t,x:J@x,(0,80),x0,method='DOP853',t_eval=TIMES,rtol=1e-11,atol=1e-13)
            X=Y*np.exp(s*TIMES[:,None]);mask=np.linalg.norm(X,axis=1)>1e-8
            direct_err=float(np.max(np.linalg.norm(X[mask]-direct.y.T[mask],axis=1)/np.linalg.norm(X[mask],axis=1)))
        initial_slope=float(np.sum((-eye(m,format='csr')+g*A.T)@x0))
        row=dict(**metadata,gain=g,eta=p['eta'],R=p['R'],rho=p['rho'],m=m,
                 spectral_rate=s,critical_gain=p['eta'],recovery_time=-1/s if s<0 else np.nan,
                 slope_20_40=sl['20_40'],slope_40_80=sl['40_80'],
                 slope_interval_difference=sl['40_80']-sl['20_40'],slope_abs_error=slope_err,
                 resolved=slope_err<=SLOPE_TOL,stable=s<0,
                 solver_relative_error=err,solver_slope_error=solver_slope_err,
                 solver_pass=err<=SOLVER_TOL and solver_slope_err<=SOLVER_TOL,
                 direct_unshifted_relative_error=direct_err,initial_mass_slope=initial_slope,
                 initial_slope_error=abs(initial_slope-(g-1)),
                 peak_total_activity=float(np.exp(logmass.max())),
                 final_total_activity=float(np.exp(logmass[-1])))
        rate_rows.append(row)
        curves.append(pd.DataFrame(dict(**metadata,gain=g,time=TIMES,
                                       total_activity=np.exp(logmass),log_total_activity=logmass,
                                       ode_log_total_activity=logmass_ode)))
    return rate_rows,pd.concat(curves,ignore_index=True)

def invariances(W,p,metadata):
    rows=[]
    for label,M in [('transpose',W.T),('scale_0.5',.5*W),('scale_2',2*W)]:
        q=parts(M)
        rows.append(dict(**metadata,transformation=label,eta=q['eta'],eta_error=abs(q['eta']-p['eta']),
                         raw_critical_gain=1/q['rho'],normalized_critical_gain=q['eta'],
                         eigen_residual=q['eigen_residual']))
    rows.append(dict(**metadata,transformation='critical_gain',eta=p['eta'],eta_error=abs(p['eta']*p['rho']/p['R']-1),
                     raw_critical_gain=1/p['rho'],normalized_critical_gain=p['eta'],eigen_residual=p['eigen_residual']))
    return rows

def examples():
    specifications={'Chain-attached':[(1,2),(1,3),(2,3),(3,4),(4,5)],
                    'Hub-attached':[(1,2),(1,3),(2,3),(1,4),(1,5)]}
    pars=[];curves=[];rates=[];bounds=[];inv=[];hetero=[];edgelist=[];group_rows=[]
    for model,edges in specifications.items():
        H,W=make_scaffold(edges);p=parts(W);pars.append(dict(model=model,**p,normalized_R=(W/p['R']).sum()/len(W)))
        slug=model.lower().replace('-','_')
        pd.DataFrame(H).to_csv(OUT/f'{slug}_incidence.csv',index=False)
        pd.DataFrame(W).to_csv(OUT/f'{slug}_W.csv',index=False)
        pd.DataFrame(W/p['R']).to_csv(OUT/f'{slug}_W_normalized.csv',index=False)
        for i,j in edges:edgelist.append(dict(model=model,source=i,target=j))
        for k,h in enumerate(H,1):group_rows.append(dict(model=model,group=k,regions='-'.join(str(i+1) for i in np.flatnonzero(h))))
        r,c=responses(W,p,dict(model=model),example=True);rates.extend(r);curves.append(c)
        for g in np.arange(241)*.005:
            direct=float(np.linalg.eigvals(-np.eye(len(W))+g*W.T/p['R']).real.max())
            bounds.append(dict(model=model,gain=g,spectral_rate=g/p['eta']-1,direct_eigen_rate=direct,eta=p['eta']))
        inv.extend(invariances(W,p,dict(model=model)))
        for name,d in [('Ascending',np.linspace(.5,1.5,len(W))),('Descending',np.linspace(1.5,.5,len(W)))]:
            next_generation=(W.T/p['R'])/d[:,None]
            rr,_=leading(next_generation)
            critical=1/rr
            hetero.append(dict(model=model,leakage_assignment=name,mean_leakage=d.mean(),critical_gain=critical,
                               eta=p['eta'],threshold_minus_eta=critical-p['eta'],
                               direct_rate_at_eta=np.linalg.eigvals(-np.diag(d)+p['eta']*W.T/p['R']).real.max()))
    save(pars,'example_parameters');save(pd.concat(curves),'example_trajectories');save(rates,'example_rates')
    save(bounds,'example_gain_curve');save(edgelist,'example_edges');save(group_rows,'example_groups')
    save(inv,'example_invariances');save(hetero,'heterogeneous_leakage_audit')
    print('Examples',pd.DataFrame(pars).to_string(index=False),flush=True)

def patients():
    """Re-run the 120 published, de-identified derived matrices (no raw SEEG)."""
    selection=INPUT/'patient_selected_windows.csv'
    track(selection,'frozen seizure-stage midpoint selection')
    selected=pd.read_csv(selection)
    assert len(selected)==120 and not selected.duplicated(['seizure_id','phase']).any()
    assert selected.seizure_id.nunique()==24 and selected.patient_id.nunique()==14
    save(selected,'patient_selected_windows')
    save(pd.DataFrame(columns=['seizure_id','phase','reason']),'patient_missing_stages')
    all_metrics=[];all_rates=[];all_inv=[];errors=[]
    trajectory_file=OUT/'patient_trajectories.csv';first=True
    for sid in sorted(selected.seizure_id.unique()):
        for _,r in selected[selected.seizure_id==sid].iterrows():
            meta=dict(patient_id=r.patient_id,seizure_id=sid,phase=r.phase,window_idx=int(r.window_idx))
            path=INPUT/'matrices'/f'{sid}_{r.phase}_W.npz'
            track(path,'derived DMWA matrix')
            W=load_npz(path).toarray()
            try:
                p=parts(W)
                all_metrics.append(dict(**meta,**p,cached_R=r.R,cached_rho=r.rho,cached_eta=r.eta,
                    R_error=abs(p['R']-r.R),rho_error=abs(p['rho']-r.rho),eta_error=abs(p['eta']-r.eta)))
                rr,cc=responses(W,p,meta);all_rates.extend(rr)
                cc.to_csv(trajectory_file,index=False,mode='w' if first else 'a',header=first);first=False
                all_inv.extend(invariances(W,p,meta))
            except Exception as e:
                errors.append(dict(**meta,error=repr(e)))
                print('FAILURE',meta,repr(e),flush=True)
        save(all_metrics,'patient_matrix_metrics');save(all_rates,'patient_rates');save(all_inv,'patient_invariances')
        save(pd.DataFrame(errors,columns=['patient_id','seizure_id','phase','window_idx','error']),'patient_failures')
        print('Completed',sid,'matrices',len(all_metrics),'responses',len(all_rates),flush=True)
    if errors: raise RuntimeError(f'{len(errors)} matrix analyses failed; see patient_failures.csv')


def report():
    rates=pd.read_csv(OUT/'patient_rates.csv');metrics=pd.read_csv(OUT/'patient_matrix_metrics.csv')
    ex=pd.read_csv(OUT/'example_rates.csv');inv=pd.read_csv(OUT/'patient_invariances.csv')
    epar=pd.read_csv(OUT/'example_parameters.csv');hetero=pd.read_csv(OUT/'heterogeneous_leakage_audit.csv')
    errors=pd.read_csv(OUT/'patient_failures.csv')
    summary=rates.groupby('gain').agg(n=('eta','size'),n_resolved=('resolved','sum'),n_solver_pass=('solver_pass','sum'),
        n_stable=('stable','sum'),max_slope_error=('slope_abs_error','max'),max_solver_relative_error=('solver_relative_error','max'),
        min_eta=('eta','min'),max_eta=('eta','max')).reset_index()
    save(summary,'patient_response_summary')
    qa=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),patients=int(metrics.patient_id.nunique()),
        seizures=int(metrics.seizure_id.nunique()),selected_matrices=int(len(pd.read_csv(OUT/'patient_selected_windows.csv'))),
        successful_matrices=len(metrics),failed_matrices=len(errors),responses=len(rates),
        reducible_matrices=int((metrics.components>1).sum()),unresolved_responses=int((~rates.resolved).sum()),
        numerical_solver_failures=int((~rates.solver_pass).sum()),max_solver_relative_error=float(rates.solver_relative_error.max()),
        max_solver_slope_error=float(rates.solver_slope_error.max()),max_asymptotic_slope_error=float(rates.slope_abs_error.max()),
        max_cached_eta_error=float(metrics.eta_error.max()),max_eigen_residual=float(metrics.eigen_residual.max()),
        max_invariance_error=float(inv.eta_error.max()),max_initial_slope_error=float(rates.initial_slope_error.max()),
        example_max_direct_ODE_error=float(ex.direct_unshifted_relative_error.max()),
        example_eta=dict(zip(epar.model,epar.eta)),heterogeneous_leakage_thresholds=hetero.to_dict('records'),
        raw_signals_reprocessed=False,physiological_validation=False,identity_requires_inferential_p_value=False,
        frozen_plan_sha256=sha(PACKAGE/'plans/eta_analysis_plan.md'),analysis_code_sha256=sha(Path(__file__)))
    (REPORT/'eta_numerical_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
    save(sources,'source_manifest')
    print(json.dumps(qa,indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--examples-only',action='store_true',help='Run the two fully synthetic scaffolds only.')
    args=parser.parse_args()
    track(PACKAGE/'plans/eta_analysis_plan.md','original frozen analysis plan')
    examples()
    if not args.examples_only:
        patients();report()
