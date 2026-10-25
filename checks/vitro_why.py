"""VITRO_WHY  Why the re-extracted in vitro data differ from the stored traces (RESULTS.md
section 42). Needs ../zigzag_scratch/orig_vitro (checks/vitro_orig_verify.m output: new
tracking read with the original code). Part 1: per movie, which new traces appear in
the stored data, and how stored traces relate to new ones. Part 2: brightness and net
change for stored / new tracking + original code / new tracking + fixed reader, and
without the movies absent from the stored data. Usage: python checks/vitro_why.py"""
import numpy as np, scipy.io as sio, glob, os
Z='../zigzag_scratch/orig_vitro/'
def L(p):
    m=sio.loadmat(p,squeeze_me=True); k=[k for k in m if not k.startswith('__')][0]; return [np.ravel(x).astype(float) for x in m[k][1:]]
pairs={'protein_RNA_1_1':'vitro_1000ngul_0p1mgml','protein_RNA_10_1':'vitro_1000ng_10ug','protein_RNA_100_1':'vitro_1000ngul_1mgml'}
for s,f in pairs.items():
    S=L(f'in_vitro_experiments/{s}_data.mat'); N=L(Z+f+'/data.mat')
    # movie of each new trace
    mv=[]
    for i,x in enumerate(sorted(glob.glob(Z+f+'/*_xy.mat'))):
        xy=sio.loadmat(x,squeeze_me=True); mv+= [i+1]*sum(np.size(t)-1>=24 for t in xy['x'])
    mv=np.array(mv)
    # for each new trace: share of its values found in the stored data
    allS=np.concatenate(S)
    cov=np.array([np.isin(n,allS).mean() for n in N])
    print(f'== {s} <- {f}: stored {len(S)} traces ({sum(len(t) for t in S)} points), new {len(N)} traces ({sum(len(t) for t in N)} points)')
    print('  new traces whose values are (mostly) in the stored data:', (cov>0.5).sum(), '| not at all:', (cov==0).sum())
    for m in np.unique(mv):
        sel=mv==m
        print(f'   movie {m:2d}: new traces {sel.sum():3d}, of which in stored {int((cov[sel]>0.5).sum()):3d}; new median level {np.median([np.median(N[i]) for i in np.where(sel)[0]]):6.0f}')
    # stored traces: relation to their best new trace
    rel={'identical':0,'stored is a piece of a longer new trace':0,'stored is longer (new trace is a piece)':0,'other':0}
    for t in S:
        best=max(range(len(N)),key=lambda j: np.isin(t,N[j]).mean()); n=N[best]; sh=np.isin(t,n).mean()
        if len(n)==len(t) and np.array_equal(n,t): rel['identical']+=1
        elif sh>0.95 and len(n)>len(t): rel['stored is a piece of a longer new trace']+=1
        elif np.isin(n,t).mean()>0.95 and len(t)>len(n): rel['stored is longer (new trace is a piece)']+=1
        else: rel['other']+=1
    print('  ',rel)

from scipy.stats import theilslopes
Z='../zigzag_scratch/orig_vitro/'
def L(p):
    m=sio.loadmat(p,squeeze_me=True); k=[k for k in m if not k.startswith('__')][0]; return [np.ravel(x).astype(float) for x in m[k][1:]]
def slope(y):
    ok=np.isfinite(y)&(y>0); t=np.arange(len(y)); return theilslopes(np.log(y[ok]),t[ok])[0]*360*100
pairs={'protein_RNA_1_1':'vitro_1000ngul_0p1mgml','protein_RNA_10_1':'vitro_1000ng_10ug','protein_RNA_100_1':'vitro_1000ngul_1mgml'}
print('version | traces >=180 | median length (all) | brightness raw | net change %/h (granule mean)')
for s,f in pairs.items():
    V={'stored':(L(f'in_vitro_experiments/{s}_data.mat'),L(f'in_vitro_experiments/{s}_background.mat'),200),
       'new tracking, original code':(L(Z+f+'/data.mat'),L(Z+f+'/background.mat'),200),
       'new tracking, fixed reader (R1)':(L(f'in_vitro_reextracted/after_R1_redchannel/{f}_data.mat'),L(f'in_vitro_reextracted/after_R1_redchannel/{f}_background.mat'),200)}
    print('==',s)
    for k,(S,B,off) in V.items():
        keep=[i for i,x in enumerate(S) if len(x)>=180]
        br=np.median([np.median(S[i]) for i in keep]); net=np.mean([slope((S[i]-B[i])/(B[i]-off)) for i in keep])
        print(f'  {k:32s} | {len(keep):3d} | {np.median([len(x) for x in S]):5.0f} | {br:6.0f} | {net:6.1f}')
print('--- new tracking, original code, without the movies absent from the stored data')
import glob
for s,f,drop in [('protein_RNA_10_1','vitro_1000ng_10ug',[5]),('protein_RNA_100_1','vitro_1000ngul_1mgml',[5,7])]:
    S=L(Z+f+'/data.mat'); B=L(Z+f+'/background.mat'); mv=[]
    for i,x in enumerate(sorted(glob.glob(Z+f+'/*_xy.mat'))):
        xy=sio.loadmat(x,squeeze_me=True); mv+=[i+1]*sum(np.size(t)-1>=24 for t in xy['x'])
    keep=[i for i,x in enumerate(S) if len(x)>=180 and mv[i] not in drop]
    allk=[i for i,x in enumerate(S) if len(x)>=180]
    print(f'  {s}: dropped movies {drop} had {len(allk)-len(keep)} of {len(allk)} long traces; left {len(keep)} | brightness {np.median([np.median(S[i]) for i in keep]):.0f} | net {np.mean([slope((S[i]-B[i])/(B[i]-200)) for i in keep]):.1f}')
