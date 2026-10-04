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


    # INBOUND CALC 
    # a_priori_probs is M_psi_theta


    # FORWARD CALC
    init_dprob = np.zeros((n+1,Dmax+1-Dmin))
    init_dprob[0,-Dmin] = 1.0
    M_d_phi = [init_dprob.copy() for _ in range(r)]

    for i in range(n):
        M_x_chi = np.zeros((4,r))
        for t in range(r):

            #for each x, sum over all di, and di1
            for x in range(4):
                s=0
                for d_i in range(Dmin, Dmax+1):
                    prob_di = M_d_phi[t][i,d_i-Dmin]
                    if prob_di == 0:
                        continue
                    for d_i1 in [d_i-1,d_i,d_i+1]:
                        s+=prob_di * phi_prob(x, d_i, d_i1, Y[t], i)
                M_x_chi[x,t]=s
        M_theta_chi =  THETA @ M_x_chi
        cum_belief = np.prod(M_theta_chi,axis=1)
        cum_belief_mat = np.tile(cum_belief[:, None], (1, r))
        M_chi_theta = cum_belief_mat / M_theta_chi
        M_x_phi = THETA.T @ M_chi_theta

        for t in range(r):

            for d_i1 in range(Dmin, Dmax+1):
                s=0
                for d_i in [d_i1-1,d_i1,d_i1+1]:
                    prob_di = M_d_phi[t][i,d_i-Dmin]
                    if prob_di == 0:
                        continue
                    for x in range(4):
                        prob_x = M_x_phi[x,t]
                        s+=prob_di * phi_prob(x, d_i, d_i1, Y[t], i) * prob_x
                M_d_phi[t][i+1,d_i1-Dmin] = s
            M_d_phi[t][i+1,:] = M_d_phi[t][i+1,:] * (1/((M_d_phi[t][i+1,:]).sum))


        
        



    M_d_phi_back=[]
    init_dprob = np.zeros((n+1,Dmax+1-Dmin))
    for t in range(r):
        M_d_phi_back.append(init_dprob.copy())
        M_d_phi_back[t][n, len(Y[t])-n-Dmin] = 1


    for i in reversed(range(n)):
        M_x_chi = np.zeros((4,r))
        for t in range(r):

            #for each x, sum over all di, and di1
            for x in range(4):
                s=0
                for d_i1 in range(Dmin, Dmax+1):
                    prob_di = M_d_phi_back[t][i,d_i1-Dmin]
                    if prob_di == 0:
                        continue
                    for d_i in [d_i1-1,d_i1,d_i1+1]:
                        s+=prob_di * phi_prob(x, d_i, d_i1, Y[t], i)
                M_x_chi[x,t]=s
        M_theta_chi =  THETA @ M_x_chi
        cum_belief = np.prod(M_theta_chi,axis=1)
        cum_belief_mat = np.tile(cum_belief[:, None], (1, r))
        M_chi_theta = cum_belief_mat / M_theta_chi
        M_x_phi = THETA.T @ M_chi_theta

        for t in range(r):

            for d_i in range(Dmin, Dmax+1):
                s=0
                for d_i1 in [d_i-1,d_i,d_i+1]:
                    prob_di = M_d_phi_back[t][i+1,d_i1-Dmin]
                    if prob_di == 0:
                        continue
                    for x in range(4):
                        prob_x = M_x_phi[x,t]
                        s+=prob_di * phi_prob(x, d_i, d_i1, Y[t], i) * prob_x
                M_d_phi_back[t][i+1,d_i1-Dmin] = s
            M_d_phi_back[t][i+1,:] = M_d_phi_back[t][i,:] * (1/((M_d_phi_back[t][i+1,:]).sum()))

    #OUTBOUND CALC


            





