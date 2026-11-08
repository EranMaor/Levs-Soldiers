import numpy as np, itertools
from py_app_calc import app_calc
import channel_sim as cs

def alphabet(w):
    L=len(w); vecs={}
    for e in itertools.product(range(4),repeat=L):
        v=np.zeros(4)
        for l in range(L): v[e[l]]+=w[l]
        vecs.setdefault(tuple(np.round(v,9)),v)
    return np.array(list(vecs.values()))

def table2(ps):
    T=np.array([[1-ps,ps/2,0,ps/2],[ps/5,1-ps,ps/5,3*ps/5],[0,ps/2,1-ps,ps/2],[ps/5,3*ps/5,ps/5,1-ps]])
    return T.T   # [y,z]

rng=np.random.default_rng(0)
configs=[
 dict(w=[1/3,2/3], n=10, r=6, pi=0.03, pd=0.05, ps=0.05, Dmin=-4, Dmax=3, asym=True),
 dict(w=[1/3,2/3], n=12, r=8, pi=0.02, pd=0.04, ps=0.08, Dmin=-3, Dmax=3, asym=False),
 dict(w=[1/8,1/8,3/8,3/8], n=8, r=5, pi=0.05, pd=0.05, ps=0.04, Dmin=-3, Dmax=2, asym=True),
 dict(w=[1/3,2/3], n=9, r=4, pi=0.0, pd=0.0, ps=0.06, Dmin=0, Dmax=0, asym=False),
]
out=open('cases.txt','w')
for ci,cf in enumerate(configs):
    TH=alphabet(cf['w']); m=len(TH)
    for attempt in range(200):
        c=rng.integers(0,m,cf['n'])
        X=[[int(rng.choice(4,p=TH[c[i]])) for i in range(cf['n'])] for _ in range(cf['r'])]
        ps=cf['ps']
        p_sub = table2(ps) if cf['asym'] else (np.eye(4)*(1-ps)+(np.ones((4,4))-np.eye(4))*ps/3)
        np.random.seed(ci*1000+attempt)
        ch=cs.Channel(cf['pi'],cf['pd'],cf['Dmin'],cf['Dmax'],p_sub)
        Y=ch.simulate(X)
        if all(cf['Dmin']<=len(y)-cf['n']<=cf['Dmax'] for y in Y): break
    else: raise SystemExit("no valid draw")
    P=app_calc(Y,TH,np.full((cf['n'],m),1/m),cf['pi'],cf['pd'],cf['Dmin'],cf['Dmax'],p_sub)
    out.write(f"CASE {ci} {cf['n']} {cf['r']} {m} {cf['pi']} {cf['pd']} {cf['Dmin']} {cf['Dmax']}\n")
    out.write("THETA\n")
    for row in TH: out.write(' '.join(repr(float(v)) for v in row)+'\n')
    out.write("PSUB\n")
    for row in p_sub: out.write(' '.join(repr(float(v)) for v in row)+'\n')
    out.write("Y\n")
    for y in Y: out.write(' '.join(str(int(v)) for v in y)+'\n')
    out.write("POST\n")
    for row in P: out.write(' '.join(repr(float(v)) for v in row)+'\n')
out.close()
print("ok", [ (cf['n'],cf['r']) for cf in configs])
