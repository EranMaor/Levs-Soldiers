import numpy as np
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

    def message_v_to_f(v_node, f_node, F):
        return np.prod(np.array[[message_f_to_v(v_node, f) for f in F[~f_node]]])

    def message_f_to_v(Sa, phi_prob, v_node, f_node, V):
        return np.sum(np.array(phi_prob*np.prod(np.array(message_v_to_f))))


    # INBOUND CALC 
    # a_priori_probs is M_psi_theta


    # FORWARD CALC
    init_dprob = np.zeros((n+1,Dmax+1-Dmin))
    init_dprob[0,-Dmin] = 1.0
    M_d_phi = [init_dprob.copy() for _ in range(r)]

    for i in range(n):
        for t in range(r):
            M_x_chi = np.zeros(4)
            #for each x, sum over all di, and di1
            for x in range(4):
                for d_i in range(Dmin, Dmax+1):
                    for d_i1 in [d_i-1,d_i,d_i+1]:
                        phi_prob(x, d_i, d_i1, Y[t], i)


    #BACKWARD CALC
    M_d_phi = []
    for t in range(r):
        d_nprob = np.zeros((n+1,Dmax+1-Dmin))
        n_prime_t = len(Y[t])
        d_nprob[n, (n - n_prime_t) - Dmin] = 1.0 if Dmin <= (n - n_prime_t) <= Dmax else 0.0
        M_d_phi.append(d_nprob)

    for i in range(n-1,0,-1):
        for t in range(r):
            pass


    #OUTBOUND CALC
    for i in range(n):
        for t in range(r):

            





