import numpy as np
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


def channel_params(channel, alpha, n):
    # Channels CH-1 to CH-5 of Table 1, with alpha = error scale.
    # p_sub[y, x] = p(read y | written x), so each column sums to 1.
    # Drift bounds follow Sec. IV-A: Dmin = -ceil(2 n p_del), Dmax = ceil(2 n p_ins).
    if channel == 1:     # DNA channel, Nanopore asymmetric substitutions (Table 2)
        p_ins, p_del, p_s = 0.017 * alpha, 0.040 * alpha, 0.043 * alpha
        table2 = np.array([            # rows = written, columns = read (as printed in Table 2)
            [1 - p_s,   p_s / 2,     0,         p_s / 2],
            [p_s / 5,   1 - p_s,     p_s / 5,   3 * p_s / 5],
            [0,         p_s / 2,     1 - p_s,   p_s / 2],
            [p_s / 5,   3 * p_s / 5, p_s / 5,   1 - p_s],
        ])
        p_sub = table2.T
    elif channel == 2:   # deletion only
        p_ins, p_del, p_s = 0.0, 0.1 * alpha, 0.0
        p_sub = None
    elif channel == 3:   # insertion only
        p_ins, p_del, p_s = 0.1 * alpha, 0.0, 0.0
        p_sub = None
    elif channel == 4:   # substitution only
        p_ins, p_del, p_s = 0.0, 0.0, 0.1 * alpha
        p_sub = None
    elif channel == 5:   # insertion, deletion and substitution, equal shares
        p_ins = p_del = p_s = 0.1 * alpha / 3
        p_sub = None
    else:
        raise ValueError("channel must be 1, 2, 3, 4 or 5")

    if p_sub is None:    # 4-ary symmetric substitution (Eq. 44)
        p_sub = (1 - p_s) * np.eye(4) + (p_s / 3) * (np.ones((4, 4)) - np.eye(4))

    return {
        "p_ins": p_ins,
        "p_del": p_del,
        "p_s": p_s,
        "p_sub": p_sub,
        "Dmin": -int(np.ceil(2 * n * p_del)),
        "Dmax": int(np.ceil(2 * n * p_ins)),
    }

        
