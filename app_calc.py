import numpy as np
def app_calc(Y, THETA, a_priori_probs, p_ins, p_del, Dmin, Dmax, p_sub=np.eye(4)):
    # Y = list of all received dna strands
    # Theta = m x 4 matrix, one composite symbol vector per row, list of all possible composite symbols
    # a_priori_probs = n x m prior over the composite symbols at each position
    # p_sub[y, x] = p(read y | written x)
    # Arrays are indexed [read t, position i, drift state d - Dmin, base x], so every loop over
    # t, d and x in Algorithm 1 is a numpy operation and only the loop over i stays.
    # In the np.einsum index strings: t = read, a = drift state, x = base.
    THETA = np.asarray(THETA, dtype=float)
    a_priori_probs = np.asarray(a_priori_probs, dtype=float)
    p_sub = np.asarray(p_sub, dtype=float)
    r = len(Y)
    n = len(a_priori_probs)
    m = len(THETA)
    D = np.arange(Dmin, Dmax + 1)
    ND = len(D)

    # p(d_next | d) for the three possible moves, as vectors over d
    p_stay = 1 - p_ins * (D < Dmax) - p_del * (D > Dmin) # staying in the same state
    p_up = p_ins * (D < Dmax) # moving up to the next state (insertion)
    p_down = p_del * (D > Dmin) # moving down to the previous state (deletion)

    # phi(x, d_i, d_i1) for every read, position and drift state.
    # stay: one read symbol at i + d_i,  insertion: two symbols,  deletion: none (no x dependence)
    PAD = max(-Dmin, Dmax) + 2
    width = max(max(len(read) for read in Y), n) + 2 * PAD + 2
    read_probs = np.zeros((r, width, 4))   # read_probs[t, j + PAD, x] = p_sub[y_t[j], x], zero outside the strand
    for t in range(r):
        read = np.asarray(Y[t], dtype=int)
        read_probs[t, PAD:PAD + len(read)] = p_sub[read]
    read_index = np.arange(n)[:, None] + D[None, :] + PAD            # n x ND, index i + d (padded)
    phi_stay = (p_stay[None, :, None] * read_probs[:, read_index]).transpose(1, 0, 2, 3)   # n x r x ND x 4
    phi_ins = (p_up[None, :, None] * read_probs[:, read_index]
               * read_probs[:, read_index + 1]).transpose(1, 0, 2, 3)
    phi_del = p_down                                                                     # ND

    # a read whose length implies a tail drift outside [Dmin, Dmax] can't be explained, so it is ignored
    # *can't happen with our current implementation
    tail = np.array([len(read) for read in Y]) - n
    valid = (tail >= Dmin) & (tail <= Dmax)

    def shift_down(drift_probs):   # out[:, d] = drift_probs[:, d - 1]
        out = np.zeros_like(drift_probs)
        out[:, 1:] = drift_probs[:, :-1]
        return out

    def shift_up(drift_probs):     # out[:, d] = drift_probs[:, d + 1]
        out = np.zeros_like(drift_probs)
        out[:, :-1] = drift_probs[:, 1:]
        return out

    def chi_to_theta(M_x_chi):
        # columns are rescaled to max 1 to avoid underflow, which doesn't change any normalised result
        M_chi_theta = THETA @ M_x_chi
        peak = M_chi_theta.max(axis=0)
        usable = valid & (peak > 0)
        M_chi_theta[:, ~usable] = 1.0
        M_chi_theta[:, usable] /= peak[usable]
        return M_chi_theta

    # needed because we vectorised
    def theta_to_x(M_chi_theta, i):
        # leave-one-out product over reads (prefix * suffix products, no division)
        ones = np.ones((m, 1))
        prefix = np.cumprod(np.hstack([ones, M_chi_theta[:, :-1]]), axis=1)
        suffix = np.cumprod(np.hstack([M_chi_theta[:, 1:], ones])[:, ::-1], axis=1)[:, ::-1]
        M_theta_chi = prefix * suffix * a_priori_probs[i][:, None]
        peak = M_theta_chi.max(axis=0)
        peak[peak == 0] = 1.0
        M_theta_chi /= peak
        return THETA.T @ M_theta_chi        # M_x_phi, 4 x r

    def gains(i, M_x_phi):
        # gain = how well each move (stay / insertion / deletion) from drift d fits the read, given the letter guess M_x_phi
        # sum over x of phi * M_x_phi for each read and starting drift state, r x ND each
        letter_probs = M_x_phi.T[:, None, :]
        gain_stay = (phi_stay[i] * letter_probs).sum(axis=2)
        gain_ins = (phi_ins[i] * letter_probs).sum(axis=2)
        gain_del = M_x_phi.sum(axis=0)[:, None] * phi_del[None, :]
        return gain_stay, gain_ins, gain_del

    def normalise(drift_probs):
        total = drift_probs.sum(axis=1, keepdims=True)
        total[total == 0] = 1.0
        return drift_probs / total

    # INBOUND CALC
    # a_priori_probs is M_psi_theta

    # FORWARD CALC
    M_d_phi = np.zeros((r, n + 1, ND))
    M_d_phi[:, 0, -Dmin] = 1.0

    for i in range(n):
        forward_drift = M_d_phi[:, i]
        M_x_chi = np.einsum('ta,tax->xt', forward_drift,
                            phi_stay[i] + phi_ins[i] + phi_del[None, :, None])
        M_x_phi = theta_to_x(chi_to_theta(M_x_chi), i)

        gain_stay, gain_ins, gain_del = gains(i, M_x_phi)
        M_d_phi[:, i + 1] = normalise(forward_drift * gain_stay
                                      + shift_down(forward_drift * gain_ins)
                                      + shift_up(forward_drift * gain_del))

    # BACKWARD CALC
    M_d_phi_back = np.zeros((r, n + 1, ND))
    M_d_phi_back[valid, n, tail[valid] - Dmin] = 1.0
    M_d_phi_back[~valid, n] = 1.0 / ND

    for i in reversed(range(n)):
        backward_drift = M_d_phi_back[:, i + 1]
        backward_drift_up, backward_drift_down = shift_up(backward_drift), shift_down(backward_drift)
        M_x_chi = (np.einsum('ta,tax->xt', backward_drift, phi_stay[i])
                   + np.einsum('ta,tax->xt', backward_drift_up, phi_ins[i])
                   + (backward_drift_down * phi_del[None, :]).sum(axis=1)[None, :])
        M_x_phi = theta_to_x(chi_to_theta(M_x_chi), i)

        gain_stay, gain_ins, gain_del = gains(i, M_x_phi)
        M_d_phi_back[:, i] = normalise(backward_drift * gain_stay
                                       + backward_drift_up * gain_ins
                                       + backward_drift_down * gain_del)

    # OUTBOUND CALC
    APP = np.zeros((n, m))
    for i in range(n):
        forward_drift = M_d_phi[:, i]
        backward_drift = M_d_phi_back[:, i + 1]
        backward_drift_up, backward_drift_down = shift_up(backward_drift), shift_down(backward_drift)
        M_x_chi = np.einsum('ta,tax->xt', forward_drift,
                            phi_stay[i] * backward_drift[:, :, None]
                            + phi_ins[i] * backward_drift_up[:, :, None]
                            + (phi_del[None, :] * backward_drift_down)[:, :, None])
        M_chi_theta = chi_to_theta(M_x_chi)

        # product over all reads times the prior, in the log domain
        with np.errstate(divide='ignore'):
            log_psi = np.log(M_chi_theta).sum(axis=1) + np.log(a_priori_probs[i])
        peak = log_psi.max()
        if np.isfinite(peak):
            M_theta_psi = np.exp(log_psi - peak)
            APP[i] = M_theta_psi / M_theta_psi.sum()

    return APP
