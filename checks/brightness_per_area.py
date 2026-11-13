"""BRIGHTNESS_PER_AREA  Per granule (>= 180 frames, as burstMeasures.m), medians over frames of: spot mean
(brightness of one granule pixel, = Table I brightness), background (one cell / dilute-phase pixel),
above = spot mean - background (granule light per pixel above the background), area = spot sum / above
(pixels), total = spot sum. Table: median over granules (and over movie medians for "above");
p: shuffle whole movies (10000x), difference of movie medians. RESULTS.md section 83.
Usage: py -3 checks/brightness_per_area.py -> analysis_outputs/granule_points/brightness_per_area.csv"""
import scipy.io as sio, numpy as np, pandas as pd, sys
sys.stdout.reconfigure(encoding='utf-8')
def L(p):
    m=sio.loadmat(p,squeeze_me=True); k=[k for k in m if not k.startswith('__')][0]; return m[k]
rows=[]; per=[]
for root,conds in [('in_vivo_experiments',['4PP7','4PP7-FUS','30PP7','30PP7-FUS','5Qb','10Qb']),('in_vitro_experiments',['protein_RNA_1_1','protein_RNA_10_1','protein_RNA_100_1'])]:
  for c in conds:
    d=f'{root}/{c}/'; S=L(d+'data.mat')[1:]; B=L(d+'background.mat')[1:]; G=L(d+'spotsum.mat')[1:]; mov=np.ravel(L(d+'movies.mat'))
    for s,b,g,m in zip(S,B,G,mov):
      if np.ravel(s).size<180: continue
      s,b,g=[np.ravel(v).astype(float) for v in (s,b,g)]
      a=s-b; ok=a>0
      per.append(dict(cond=c,movie=m,spot=np.median(s),bg=np.median(b),above=np.median(a),area=np.median(g[ok]/a[ok]) if ok.any() else np.nan,total=np.median(g)))
P=pd.DataFrame(per); P.to_csv('analysis_outputs/granule_points/brightness_per_area.csv',index=False)
cols=['spot','bg','above','area','total']
gm=P.groupby('cond',sort=False)[cols].median()
mm=P.groupby(['cond','movie'],sort=False)[cols].median().groupby('cond',sort=False).median()
n=P.groupby('cond',sort=False).agg(granules=('spot','size'),movies=('movie','nunique'))
out=n.join(gm.round(0).astype(int)).join(mm[['above']].round(0).astype(int).rename(columns={'above':'above_movies'}))
print(out.to_string())
# shuffle movies p for above, PP7 vs FUS
rng=np.random.default_rng(1)
def p(c1,c2,col='above'):
  A=P[P.cond==c1].groupby('movie')[col].median().values; Bv=P[P.cond==c2].groupby('movie')[col].median().values
  obs=np.median(A)-np.median(Bv); z=np.r_[A,Bv]; k=len(A); c=0
  for _ in range(10000):
    rng.shuffle(z); c+=abs(np.median(z[:k])-np.median(z[k:]))>=abs(obs)-1e-9
  return round(obs), (c+1)/10001
for c1,c2 in [('4PP7','4PP7-FUS'),('30PP7','30PP7-FUS'),('5Qb','10Qb'),('protein_RNA_1_1','protein_RNA_100_1')]:
  for col in ['spot','above','area','total']: print(c1,c2,col,p(c1,c2,col))
