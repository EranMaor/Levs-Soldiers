def app_calc(Y, THETA, a_priori_probs, p_ins, p_del, Dmin, Dmax, p_sub=np.eye(4)):
    # Y = list of all received dna strands
    # Theta = list of all possible composite symbol vectors
    r = len(Y)
    n = len(a_priori_probs)
    m = len(THETA)

    def p_d_given_next(d, d_next):
        if d_next == d + 1:
            return p_ins if d < Dmax else 0.0
        elif d_next == d - 1:
            return p_del if d > Dmin else 0.0
        elif d_next == d:
            if d == Dmin:
                return 1 - p_ins
            elif d == Dmax:
                return 1 - p_del
            else:
                return 1 - p_ins - p_del
        else:
            return 0.0
   
    def phi_prob(x, d_i, d_i1, y_t, i):
        condprob = p_d_given_next(d_i, d_i1)

        start = i + d_i
        end = i + d_i1

        for j in range(start, end + 1):
            if j < 0 or j >= len(y_t):   # out of scope of sequence, so wont work
                return 0.0
            
            condprob *= p_sub[y_t[j], x]
        return condprob

    # inbound calc a_priori_probs is M_psi_theta
    init_dprob = np.zeros((n+1,Dmax+1-Dmin))
    init_dprob[0,-Dmin] = 1.0
    dprob = [init_dprob for x in range(r)]. #init Mdphi
    # forward calc
    for i in range(n):
        for t in range(r):
            #for each x, sum over all di, and di1
            for x in range(4):
                for d_i in range(Dmin, Dmax+1):
                    for d_i1 in [d_i-1,d_i,d_i+1]:
                        pass
            





