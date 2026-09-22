import math
import numpy as np

#the channel sim
class Channel:
    def probsim(self, p):
        if np.random.random()<p:
            return True
        else:
            return False
    def __init__(self, p_ins, p_del, Dmin, Dmax, p_sub=np.eye(4)):
        self.p_ins = p_ins
        self.p_del = p_del
        self.Dmin = Dmin
        self.Dmax = Dmax
        self.p_sub= p_sub
    def drift(self, x):
        d=0
        drifted_x=[]
        for ind in x:
            if d > self.Dmin:
                if  self.probsim(self.p_del):
                    d+=-1
                    continue
                elif d < self.Dmax and self.probsim((self.p_ins)/(1-self.p_del)):
                    d+=1
                    drifted_x.append(ind)
            elif self.probsim(self.p_ins):
                d+=1
                drifted_x.append(ind)
            drifted_x.append(ind)
        return(drifted_x)
    def simulate(self, X):
        DX = [self.drift(x) for x in X]
        Y=[]
        for strand in DX:
            current=[]
            for ind in strand:
                nucleo = np.random.choice([0,1,2,3],p=self.p_sub[:, ind])
                current.append(nucleo)
            Y.append(current)
        return Y


#probabilites of insersions, deletions, substitutions given in the paper.
CH1_P_INS = 0.017         
CH1_P_DEL = 0.040
CH1_P_S = 0.043            



def build_SUB(p_s=CH1_P_S): #this is the substitution probablity matrix
    
    return [
        [1.0 - p_s, p_s / 2.0, 0.0, p_s / 2.0],            
        [p_s / 5.0, 1.0 - p_s, p_s / 5.0, 3.0 * p_s / 5.0],  
        [0.0, p_s / 2.0, 1.0 - p_s, p_s / 2.0],            
        [p_s / 5.0, 3.0 * p_s / 5.0, p_s / 5.0, 1.0 - p_s],  
    ]


def build_TRANS(D, p_ins=CH1_P_INS, p_del=CH1_P_DEL): #transition matrix from drift
    
    ND = len(D)
    D_MIN, D_MAX = min(D), max(D)
    T = [[0.0] * ND for _ in range(ND)]
    for ai, di in enumerate(D):
        for bi, dj in enumerate(D):
            if dj == di + 1 and di < D_MAX:
                T[ai][bi] = p_ins
            elif dj == di - 1 and di > D_MIN:
                T[ai][bi] = p_del
            elif dj == di:
                if di == D_MIN:
                    T[ai][bi] = 1.0 - p_ins
                elif di == D_MAX:
                    T[ai][bi] = 1.0 - p_del
                else:
                    T[ai][bi] = 1.0 - p_ins - p_del
            else:
                T[ai][bi] = 0.0
    return T


def D_from_bounds(Dmin, Dmax):
    #list of drift states we're allowing
    return list(range(int(Dmin), int(Dmax) + 1))


def default_D(Y, n, p_ins=CH1_P_INS, p_del=CH1_P_DEL, slack=1):

    nets = [len(y) - n for y in Y]
    lo = min([0, min(nets)] + [-(2.0 * n * p_del + slack) - 1])
    hi = max([0, max(nets)] + [2.0 * n * p_ins + slack + 1])
    return D_from_bounds(math.floor(lo), math.ceil(hi))


#basic functions to be used later
def norm(v):
    #normalizes list of probablities so it sums to 1
    s = 0.0
    for x in v:
        s += x
    if s > 0.0:
        return [x / s for x in v]
    return [1.0 / len(v)] * len(v)


def one_hot_in_D(D, value):
    #creates prob distribution where prob of one thing is 1 and rest is 0
    if value in D:
        out = [0.0] * len(D)
        out[D.index(value)] = 1.0
        return out
    return None


def emission(y, j0, j1, P_sub): #prob of observed DNA bases between two drift positions

    if j0 > j1:
        return [1.0, 1.0, 1.0, 1.0]          
    if j0 < 0 or j1 >= len(y):
        return [0.0, 0.0, 0.0, 0.0]          
    out = [1.0, 1.0, 1.0, 1.0]
    for j in range(j0, j1 + 1):
        row = P_sub[y[j]]                    
        for x in range(len(out)):
            out[x] *= row[x]
    return out

def _sub_cols(P_sub):

    m = len(P_sub)
    return np.array([[P_sub[z][y] for z in range(m)] for y in range(m)])

def _sub_for_emission(P_sub_row_oriented):

    return _sub_cols(P_sub_row_oriented)


def shared_recipes(THETA, n):

    return [THETA for _ in range(n)]


def _recipes_for(recipes, n, m):

    if recipes is None:
        raise ValueError("recipes is required")
    if isinstance(recipes, list) and len(recipes) == n and isinstance(recipes[0], list) \
            and isinstance(recipes[0][0], list):
        return recipes
    return shared_recipes(recipes, n)

def chi_x_to_theta(MxChi, recipes_i, m):
    out = [0.0] * m
    for r in range(m):
        row = recipes_i[r]
        s = 0.0
        for x in range(m):
            s += row[x] * MxChi[x]
        out[r] = s
    return norm(out)


def chi_theta_to_x(MthetaChi, recipes_i, m):
    nb = len(recipes_i[0])                     
    out = [0.0] * nb
    for x in range(nb):
        s = 0.0
        for r in range(m):
            s += recipes_i[r][x] * MthetaChi[r]
        out[x] = s
    return norm(out)


def theta_excluding(MxTheta, t, m):
    r = len(MxTheta)
    out = [1.0] * m
    for t0 in range(r):
        if t0 == t:
            continue
        for q in range(m):
            out[q] *= MxTheta[t0][q]
    return norm(out)


def theta_all(MxTheta, m):
    out = [1.0] * m
    for t in range(len(MxTheta)):
        for q in range(m):
            out[q] *= MxTheta[t][q]
    return norm(out)


def phi_to_dnext(y, i, F_i, MxPhi, D, P_sub, TRANS_local):
    ND = len(D)
    out = [0.0] * ND
    for a in range(ND):
        if F_i[a] == 0.0:
            continue
        for b in range(ND):
            if TRANS_local[a][b] == 0.0:
                continue
            em = emission(y, i + D[a], i + D[b], P_sub)
            wsum = 0.0
            for x in range(4):
                wsum += em[x] * MxPhi[x]
            out[b] += F_i[a] * TRANS_local[a][b] * wsum
    return norm(out)


def phi_to_dprev(y, i, B_next, MxPhi, D, P_sub, TRANS_local):
    ND = len(D)
    out = [0.0] * ND
    for a in range(ND):
        for b in range(ND):
            if TRANS_local[a][b] == 0.0 or B_next[b] == 0.0:
                continue
            em = emission(y, i + D[a], i + D[b], P_sub)
            wsum = 0.0
            for x in range(4):
                wsum += em[x] * MxPhi[x]
            out[a] += B_next[b] * TRANS_local[a][b] * wsum
    return norm(out)



#inbound
def inbound(recipes, Y, D=None, *, P_sub=None, P_trans=None, m=4, #uses the onehot prob function since we know drift at start is 0 so forward message starts with that
            codebook=None, paper_line15=True):
    n = len(recipes)
    R = len(Y)
    if D is None:
        D = default_D(Y, n)
    D = list(D)
    ND = len(D)
    if m is None:
        m = len(recipes[0])
    if P_sub is None:
        P_sub = _sub_for_emission(build_SUB(CH1_P_S))
    if P_trans is None:
        P_trans = build_TRANS(D)
    if codebook is None:
        codebook = list(range(m))
    F = [[None] * (n + 1) for _ in range(R)]
    B = [[None] * (n + 1) for _ in range(R)]

    prior = [1.0 / m] * m

    out_of_D = []
    for t in range(R):

        F[t][0] = one_hot_in_D(D, 0)

        n_prime = len(Y[t])
        declared_printed = n - n_prime
        declared_eq9 = n_prime - n
        declared = declared_printed if paper_line15 else declared_eq9

        oh = one_hot_in_D(D, declared)
        if oh is None:
            out_of_D.append((t, declared, len(Y[t])))
            oh = [0.0] * ND
        B[t][n] = oh

    ctx = {
        "n": n, "R": R, "D": D, "ND": ND, "m": m,
        "P_sub": P_sub, "P_trans": P_trans, "codebook": codebook,
        "prior": prior, "out_of_D": out_of_D,
        "paper_line15": paper_line15,
    }
    return F, B, ctx

def forward(recipes, Y, F, B, D=None, *, ctx=None, P_sub=None, P_trans=None, #at each position uses current drift info to estimate the DNA base then convert this info to info about the symbol. and also combines info from other strands
            m=4, codebook=None):

    ctx = _resolve_ctx(ctx, D, P_sub, P_trans, m, codebook, recipes, Y)
    D, m, P_sub, P_trans = ctx["D"], ctx["m"], ctx["P_sub"], ctx["P_trans"]
    n, R = ctx["n"], ctx["R"]
    U = [1.0 / len(D)] * len(D)
    recipes = _recipes_for(recipes, n, m)

    for i in range(n):
        MxTheta = []
        for t in range(R):
            MxChi = phi_to_x_p(Y[t], i, F[t][i], U, D, P_sub, P_trans)
            MxTheta.append(chi_x_to_theta(MxChi, recipes[i], m))

        for t in range(R):
            MthetaChi = theta_excluding(MxTheta, t, m)
            mphi = chi_theta_to_x(MthetaChi, recipes[i], m)
            F[t][i + 1] = phi_to_dnext(Y[t], i, F[t][i], mphi, D, P_sub, P_trans)
    return F

def backward(recipes, Y, F, B, D=None, *, ctx=None, P_sub=None, P_trans=None, #same as forward except starts from end and also uses info from later observations.
             m=4, codebook=None):
    ctx = _resolve_ctx(ctx, D, P_sub, P_trans, m, codebook, recipes, Y)
    D, m, P_sub, P_trans = ctx["D"], ctx["m"], ctx["P_sub"], ctx["P_trans"]
    n, R = ctx["n"], ctx["R"]
    U = [1.0 / len(D)] * len(D)
    recipes = _recipes_for(recipes, n, m)

    for i in range(n - 1, -1, -1):
        MxTheta = []
        for t in range(R):
            MxChi = phi_to_x_p(Y[t], i, U, B[t][i + 1], D, P_sub, P_trans)
            MxTheta.append(chi_x_to_theta(MxChi, recipes[i], m))

        for t in range(R):
            MthetaChi = theta_excluding(MxTheta, t, m)
            mphi = chi_theta_to_x(MthetaChi, recipes[i], m)
            B[t][i] = phi_to_dprev(Y[t], i, B[t][i + 1], mphi, D, P_sub, P_trans)

    return B



def outbound(recipes, Y, F, B, D=None, *, ctx=None, P_sub=None, P_trans=None, #at each c_i calcultes prob distribution for each symbol
             m=4, codebook=None):

    ctx = _resolve_ctx(ctx, D, P_sub, P_trans, m, codebook, recipes, Y)
    D, m, P_sub, P_trans = ctx["D"], ctx["m"], ctx["P_sub"], ctx["P_trans"]
    n, R, cb = ctx["n"], ctx["R"], ctx["codebook"]
    recipes = _recipes_for(recipes, n, m)

    posteriors = []
    for i in range(n):
        MxTheta = []
        for t in range(R):
            MxChi = phi_to_x_p(Y[t], i, F[t][i], B[t][i + 1], D, P_sub, P_trans)
            MxTheta.append(chi_x_to_theta(MxChi, recipes[i], m))

        full = theta_all(MxTheta, m)
        ev = [full[cb[v]] for v in range(m)]
        posteriors.append(norm(ev))

    return posteriors


def phi_to_x_p(y, i, F_i, B_next, D, P_sub, P_trans):
    """phi_to_x with the transition matrix threaded in (no module state)."""
    ND = len(D)
    acc = [0.0] * 4
    for a in range(ND):
        fa = F_i[a]
        if fa == 0.0:
            continue
        for b in range(ND):
            w = fa * P_trans[a][b] * B_next[b]
            if w == 0.0:
                continue
            em = emission(y, i + D[a], i + D[b], P_sub)
            for x in range(4):
                acc[x] += w * em[x]
    return norm(acc)

def _resolve_ctx(ctx, D, P_sub, P_trans, m, codebook, recipes, Y):
    """Build or repair the context dictionary the three later stages share."""
    n = len(recipes)
    if ctx is None:
        if D is None:
            D = default_D(Y, n)
        D = list(D)
        if P_sub is None:
            P_sub = _sub_for_emission(build_SUB(CH1_P_S))
        if P_trans is None:
            P_trans = build_TRANS(D)
        if codebook is None:
            codebook = list(range(m))
        ctx = {"n": n, "R": len(Y), "D": D, "ND": len(D), "m": m,
               "P_sub": P_sub, "P_trans": P_trans, "codebook": codebook,
               "prior": [1.0 / m] * m, "out_of_D": [],
               "paper_line15": None}
    if m is None:
        m = ctx.get("m", 4)
        ctx["m"] = m
    return ctx


def run_pipeline(recipes, Y, D=None, *, P_sub=None, P_trans=None, m=4,
                 codebook=None, paper_line15=True):
    """The four stages end to end.  Returns (posteriors, ctx)."""
    F, B, ctx = inbound(recipes, Y, D, P_sub=P_sub, P_trans=P_trans, m=m,
                        codebook=codebook, paper_line15=paper_line15)
    forward(recipes, Y, F, B, ctx=ctx)
    backward(recipes, Y, F, B, ctx=ctx)
    post = outbound(recipes, Y, F, B, ctx=ctx)
    return post, ctx


def hard_decisions(posteriors):
    """argmax over symbols at every position."""
    out = []
    for p in posteriors:
        best, bi = -1.0, 0
        for v in range(len(p)):
            if p[v] > best:
                best, bi = p[v], v
        out.append(bi)
    return out


#a demo to try

def _demo():
    n, r, m = 12, 4, 4
    P_I, P_D, P_S = CH1_P_INS, CH1_P_DEL, CH1_P_S
    Dmin, Dmax = -1, +1
    rng = np.random.RandomState(20260921)

    print("DNA decoding demo")

    A = [0.85, 0.10, 0.03, 0.02]
    T = [0.05, 0.80, 0.10, 0.05]
    G = [0.02, 0.03, 0.90, 0.05]
    C = [0.10, 0.05, 0.05, 0.80]
    THETA = [A, T, G, C]
    codebook = list(range(m))
    recipes = shared_recipes(THETA, n)

    #creating the true codeword written strands
    c_true = [int(v) for v in rng.randint(0, m, size=n)]
    X = []
    for t in range(r):
        strand = []
        for i in range(n):
            theta = THETA[codebook[c_true[i]]]
            strand.append(int(rng.choice([0, 1, 2, 3], p=theta)))
        X.append(strand)

    sub_cols = _sub_cols(build_SUB(P_S))
    ch = Channel(P_I, P_D, Dmin, Dmax, sub_cols)

    np.random.seed(20260921)
    Y = ch.simulate(X)

    print("True codeword:", " ".join(str(v) for v in c_true))

    print("Reads:")
    for t in range(r):
        print(" ", "".join("ATGC"[b] for b in Y[t]))

    D = D_from_bounds(Dmin, Dmax)

    print("\nResult:")
    post, ctx = run_pipeline(
        recipes, Y, D, m=m, codebook=codebook,
        paper_line15=True
    )

    dec = hard_decisions(post)
    wrong = sum(dec[i] != c_true[i] for i in range(n))

    print("Decoded:", " ".join(str(v) for v in dec))
    print("Wrong:", wrong, "/", n)

    return 0


if __name__ == "__main__":
    raise SystemExit(_demo())
