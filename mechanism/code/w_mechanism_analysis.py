"""Conditional W-model closure and response analysis; no plotting in Python."""
from pathlib import Path
import json,itertools,time,hashlib,os
import numpy as np
import pandas as pd
from scipy import linalg
from scipy.sparse.linalg import expm_multiply
from threadpoolctl import threadpool_limits
threadpool_limits(1)
PACKAGE=Path(__file__).resolve().parents[1]
ROOT=Path(os.environ.get('DMWA_OUTPUT_ROOT', str(PACKAGE)))
OUT=ROOT/'data'/'w';OUT.mkdir(parents=True,exist_ok=True)
REPORT=ROOT/'reports';REPORT.mkdir(exist_ok=True)
TIMES=np.linspace(0,10,101)
G=.5;DECAY=1.;TOL=1e-10

def scaffold(seed,family,n=12):
    rng=np.random.default_rng(seed);A=np.zeros((n,n),bool)
    for u in range(n):A[u,(u+1)%n]=A[(u+1)%n,u]=True
    for u in range(n):
        for v in range(u+1,n):
            if family=='Random':p=.20
            elif family=='Core-periphery':p=.70 if u<4 and v<4 else (.25 if (u<4)!=(v<4) else .08)
            else:p=.65 if u//4==v//4 else .05
            if rng.random()<p:A[u,v]=A[v,u]=True
    return A

def toy_graph(name):
    if name=='Triangle':n=3;edges=[(0,1),(0,2),(1,2)]
    elif name=='Triangle + branches':n=5;edges=[(0,1),(0,2),(1,2),(0,3),(0,4)]
    else:n=4;edges=[(0,1),(0,2),(1,2),(1,3),(2,3)]
    A=np.zeros((n,n),bool)
    for u,v in edges:A[u,v]=A[v,u]=True
    return A

def build(A):
    n=len(A);edges=list(zip(*np.where(np.triu(A,1))))
    triangles=[x for x in itertools.combinations(range(n),3) if A[x[0],x[1]] and A[x[0],x[2]] and A[x[1],x[2]]]
    groups=[tuple(x) for x in edges]+triangles;m=len(groups)
    H=np.zeros((m,n))
    for i,e in enumerate(groups):H[i,list(e)]=1
    size=H.sum(1);O=H@H.T;np.fill_diagonal(O,0)
    degree=(O>0).sum(1);den=degree[:,None]+degree[None,:]
    bias=np.divide(degree[None,:],den,out=np.full_like(O,.5),where=den>0)
    cov=O/size[:,None];same=size[:,None]==size[None,:]
    Ws={'Coverage':cov,'Constant bias':cov*np.where(same,.5,1.),'Full DMWA':cov*np.where(same,bias,1.)}
    return H,size,degree,groups,Ws

def projection(H,size,degree,kind):
    n=H.shape[1];m=len(H)
    if kind=='Node':keys=[0]*m
    elif kind=='Node + order':keys=[int(s) for s in size]
    else:keys=[(int(s),int(d)) for s,d in zip(size,degree)]
    unique=sorted(set(keys));P=np.zeros((n*len(unique),m));Ms=[]
    for i,key in enumerate(unique):
        ix=[j for j,k in enumerate(keys) if k==key]
        P[i*n:(i+1)*n,ix]=H[ix].T/size[ix]
        Ms.append((H[ix].T/size[ix])@H[ix])
    merge=np.tile(np.eye(n),(1,len(unique)))
    return P,merge,unique,Ms

def class_generator(keys,Ms,kappa,variant):
    n=len(Ms[0]);q=len(keys);B=np.zeros((n*q,n*q))
    for b,target in enumerate(keys):
        for a,source in enumerate(keys):
            if variant=='Coverage':coef=1.
            elif variant=='Constant bias':coef=.5 if source==target else 1.
            else:
                sa,da=source;sb,db=target
                coef=1. if sa!=sb else (db/(da+db) if da+db else .5)
            B[b*n:(b+1)*n,a*n:(a+1)*n]=kappa*coef*Ms[b]
        diagonal_penalty=1. if variant=='Coverage' else .5
        B[b*n:(b+1)*n,b*n:(b+1)*n]-=(DECAY+kappa*diagonal_penalty)*np.eye(n)
    return B

def integrate(J,Z):
    return expm_multiply(J,Z,start=0,stop=10,num=101,endpoint=True,traceA=np.trace(J))

def intg(y):return np.trapezoid(y,TIMES,axis=0)

def same_node_pair(H,groups,C,seed,is_toy):
    m=len(H)
    if is_toy:
        z1=np.zeros(m);z2=np.zeros(m);z1[groups.index((0,1,2))]=1
        for edge in [(0,1),(0,2),(1,2)]:z2[groups.index(edge)]=1/3
        return np.c_[z1,z2],'Fixed triangle versus its three edges'
    rng=np.random.default_rng(seed+1000000)
    v=rng.normal(size=m);d=v-linalg.pinv(C,rtol=1e-12)@(C@v)
    if np.max(abs(d))<1e-10:return None,'Projection injective or null direction numerically absent'
    d*=.8/m/np.max(abs(d));z0=np.full(m,1/m)
    return np.c_[z0+d,z0-d],'Geometry-only fixed random null direction'

def summarize(values):
    v=np.asarray(values,float);v=v[np.isfinite(v)]
    return dict(n=len(v),median=float(np.median(v)),q25=float(np.quantile(v,.25)),q75=float(np.quantile(v,.75)),
                mean=float(v.mean()),minimum=float(v.min()),maximum=float(v.max()))

closure=[];response=[];separation=[];witness=[];meta=[];initial=[];operators=[];toytime=[];constraints=[]
topology=[];groupdata=[];failures=[]
start=time.time()
sets=[('toy_'+str(i),name,None,toy_graph(name),True) for i,name in enumerate(['Triangle','Triangle + branches','Double triangle'])]
families=['Random','Core-periphery','Modular']
sets += [(f'network_{i+1:03d}',families[i%3],930000+i,scaffold(930000+i,families[i%3]),False) for i in range(100)]
for ii,(network,family,seed,A,is_toy) in enumerate(sets):
    H,size,degree,groups,Ws=build(A);n=len(A);m=len(H);C=H.T/size
    Pmap={k:projection(H,size,degree,k) for k in ['Node','Node + order','Node + order + degree']}
    Z,pair_type=same_node_pair(H,groups,C,seed,is_toy)
    rC=np.linalg.matrix_rank(C,tol=1e-10)
    md=dict(network=network,family=family,is_toy=is_toy,seed=seed,n=n,m=m,n_pairs=int((size==2).sum()),n_triangles=int((size==3).sum()),
            node_rank=rC,node_nullity=m-rC,has_pair=Z is not None,pair_type=pair_type)
    for kind,(P,merge,keys,Ms) in Pmap.items():
        md[kind+'_dimension']=len(P);md[kind+'_active_dimension']=int((np.linalg.norm(P,axis=1)>0).sum());md[kind+'_rank']=np.linalg.matrix_rank(P,tol=1e-10)
    meta.append(md)
    for u,v in zip(*np.where(np.triu(A,1))):topology.append(dict(network=network,u=int(u),v=int(v)))
    for j,grp in enumerate(groups):groupdata.append(dict(network=network,group=j,members=';'.join(map(str,grp)),size=int(size[j]),degree=int(degree[j])))
    if Z is not None:
        e=float(np.max(abs(C@Z[:,0]-C@Z[:,1])));minimum=float(Z.min());masserr=float(abs(Z[:,0].sum()-Z[:,1].sum()))
        assert e<1e-10 and minimum>=-1e-12 and masserr<1e-10
        witness.append(dict(network=network,same_node_error=e,min_group_state=minimum,mass_difference=masserr,group_state_L1=float(np.abs(Z[:,0]-Z[:,1]).sum())))
        if is_toy:
            for state in range(2):
                for j in range(m):initial.append(dict(network=network,state=state+1,group=j,members=''.join(chr(65+v) for v in groups[j]),mass=Z[j,state]))
    # Accessible member-mean states: C^T x lies in range(C^T), orthogonal to ker(C).
    if Z is not None:
        d=Z[:,0]-Z[:,1];N=np.eye(m)-linalg.pinv(C,rtol=1e-12)@C
        constraints.append(dict(network=network,nullpair_distance_from_member_mean_subspace=float(np.linalg.norm(N@d)/np.linalg.norm(d)),
                                initial_node_error=float(np.linalg.norm(C@d))))
    for variant,W in Ws.items():
        R=float(W.sum()/m);rho=float(max(abs(linalg.eigvals(W))));eta=R/rho if rho>0 else np.nan;kappa=G/R
        J=-DECAY*np.eye(m)+kappa*W.T
        alpha=float(np.max(linalg.eigvals(J).real))
        operators.append(dict(network=network,variant=variant,R=R,rho=rho,eta=eta,g=G,decay=DECAY,spectral_abscissa=alpha,stable=alpha<0))
        truth=integrate(J,Z) if Z is not None else None
        ytrue=np.einsum('nm,tms->tns',C,truth) if Z is not None else None
        if ytrue is not None:
            assert np.isfinite(ytrue).all()
            totals=ytrue.sum(1);profiles=ytrue/np.maximum(totals[:,None,:],1e-300)
            tv=.5*np.abs(profiles[:,:,0]-profiles[:,:,1]).sum(1)
            total_frac=float(intg(abs(totals[:,0]-totals[:,1]))/intg(totals.mean(1)))
            mean_tv=float(intg(tv)/10)
            sep=dict(network=network,family=family,is_toy=is_toy,variant=variant,total_fraction=total_frac,spatial_TV_mean=mean_tv,
                     spatial_TV_t1=float(tv[10]),spatial_TV_t10=float(tv[-1]),spatial_TV_max=float(tv.max()),
                     total_ratio_t1=float(totals[10,0]/totals[10,1]),total_ratio_t10=float(totals[-1,0]/totals[-1,1]))
            separation.append(sep)
            if is_toy:
                for ti,t in enumerate(TIMES):
                    for state in range(2):
                        for node in range(n):toytime.append(dict(network=network,family=family,variant=variant,time=t,state=state+1,node=chr(65+node),value=ytrue[ti,node,state],total=totals[ti,state],fraction=profiles[ti,node,state],TV=tv[ti]))
        for kind,(P,merge,keys,Ms) in Pmap.items():
            pinv=linalg.pinv(P,rtol=1e-12);null=np.eye(m)-pinv@P
            singular=linalg.svdvals(P);nonzero=singular[singular>singular[0]*1e-12]
            condition=float(nonzero[0]/nonzero[-1]) if len(nonzero) else np.inf
            denom=linalg.norm(P@W.T,'fro');res=linalg.norm(P@W.T@null,'fro')/denom if denom>1e-15 else 0.
            B=P@J@pinv
            exact_expected=(variant=='Coverage' and kind=='Node') or (variant=='Constant bias' and kind=='Node + order') or (variant=='Full DMWA' and kind=='Node + order + degree')
            analytic_error=np.nan
            if exact_expected:
                if variant=='Coverage':Ba=-DECAY*np.eye(n)+kappa*(C@H-np.eye(n))
                else:Ba=class_generator(keys,Ms,kappa,variant)
                analytic_error=float(linalg.norm(P@J-Ba@P,'fro')/max(linalg.norm(P@J,'fro'),1e-15))
                assert res<TOL and analytic_error<TOL,(network,variant,kind,res,analytic_error)
                B=Ba
            closure.append(dict(network=network,family=family,is_toy=is_toy,variant=variant,projection=kind,residual=float(res),
                                exact_expected=exact_expected,analytic_error=analytic_error,dimension=len(P),rank=np.linalg.matrix_rank(P,tol=1e-10),
                                observed_closed=bool(res<TOL),nonzero_spectrum_condition=condition))
            # Fixed finite-time readout tests: all variants at node level and all contexts for full W.
            if Z is not None and (kind=='Node' or variant=='Full DMWA' or exact_expected):
                yt=integrate(B,P@Z);yn=np.einsum('nk,tks->tns',merge,yt)
                finite=np.isfinite(yn).all()
                if finite:
                    sqerr=np.sum((yn-ytrue)**2,axis=(1,2));sqnorm=np.sum(ytrue**2,axis=(1,2))
                    rel=float(np.sqrt(intg(sqerr)/intg(sqnorm)))
                    maxerr=float(np.max(abs(yn-ytrue)))
                    neg=float((yn<-1e-10).mean())
                    contextneg=float((yt<-1e-10).mean())
                    contextmin=float(yt.min())
                else:rel=np.nan;maxerr=np.nan;neg=np.nan;contextneg=np.nan;contextmin=np.nan;failures.append(dict(network=network,variant=variant,projection=kind,failure='nonfinite projected response'))
                response.append(dict(network=network,family=family,is_toy=is_toy,variant=variant,projection=kind,relative_L2_error=rel,
                                     max_absolute_error=maxerr,negative_fraction=neg,finite=bool(finite),exact_expected=exact_expected,
                                     context_negative_fraction=contextneg,minimum_context_state=contextmin,
                                     reduced_spectral_abscissa=float(np.max(linalg.eigvals(B).real))))
                if exact_expected:assert finite and rel<1e-8,(network,variant,kind,rel)
    if (ii+1)%10==0:print(f'{ii+1}/{len(sets)} scaffolds, {time.time()-start:.1f}s',flush=True)

frames={'closure':closure,'response_errors':response,'same_node_separation':separation,'witness_checks':witness,'scaffold_manifest':meta,
        'operators':operators,'toy_initial_states':initial,'toy_trajectories':toytime,'admissibility_limits':constraints,'scaffold_edges':topology,'hyperedges':groupdata}
for name,rows in frames.items():pd.DataFrame(rows).to_csv(OUT/(name+'.csv'),index=False)
pd.DataFrame(failures,columns=['network','variant','projection','failure']).to_csv(OUT/'failures.csv',index=False)
cd=pd.DataFrame(closure);rd=pd.DataFrame(response);sd=pd.DataFrame(separation)
cs=[]
for (variant,proj),q in cd[~cd.is_toy].groupby(['variant','projection']):
    cs.append(dict(variant=variant,projection=proj,**summarize(q.residual),n_exact=int((q.residual<TOL).sum())))
pd.DataFrame(cs).to_csv(OUT/'closure_summary.csv',index=False)
rs=[]
for (variant,proj),q in rd[~rd.is_toy].groupby(['variant','projection']):rs.append(dict(variant=variant,projection=proj,**summarize(q.relative_L2_error)))
pd.DataFrame(rs).to_csv(OUT/'response_summary.csv',index=False)
ss=[]
for variant,q in sd[~sd.is_toy].groupby('variant'):
    for metric in ['total_fraction','spatial_TV_mean']:ss.append(dict(variant=variant,metric=metric,**summarize(q[metric])))
pd.DataFrame(ss).to_csv(OUT/'separation_summary.csv',index=False)
od=pd.DataFrame(operators);md=pd.DataFrame(meta)
acceptance=dict(n_formal=100,n_toys=3,n_no_hidden_pair=int((~md.has_pair).sum()),n_failures=len(failures),
    exact_identity_max_residual=float(cd[cd.exact_expected].residual.max()),
    exact_formula_max_error=float(cd[cd.exact_expected].analytic_error.max()),
    exact_response_max_relative_error=float(rd[rd.exact_expected].relative_L2_error.max()),
    n_unstable_operators=int((~od.stable).sum()),minimum_eta=float(od.eta.min()),
    max_same_node_initial_error=float(pd.DataFrame(witness).same_node_error.max()),
    maximum_nonzero_projection_condition=float(cd.nonzero_spectrum_condition.max()),
    maximum_node_negative_fraction=float(rd.negative_fraction.max()),maximum_context_negative_fraction=float(rd.context_negative_fraction.max()),
    seconds=time.time()-start,
    plan_sha256=hashlib.sha256((PACKAGE/'plans/W_frozen_plan.md').read_bytes()).hexdigest(),
    interpretation='Conditional W-defined model; independent group states; no independent neural-dynamics validation')
(OUT/'acceptance.json').write_text(json.dumps(acceptance,indent=2),encoding='utf8')
print(json.dumps(acceptance,indent=2),flush=True)
