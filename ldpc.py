"""Soft LDPC over GF(256): FFT-based sum-product decoder (Kaneko & Shibata).

Ports gf256.m, gf_helpers.m, make_ldpc.m, ldpc_encode.m, ldpc_decode.m, ldpc_check.m.
Symbols are integers 0..255. Addition is XOR. The primitive polynomial is
x^8 + x^4 + x^3 + x^2 + 1 (0x11D).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# GF(256) tables
# ---------------------------------------------------------------------------

def _build_gf256():
    prim = 0x11D
    exp = np.zeros(512, dtype=np.int32)
    log = np.zeros(256, dtype=np.int32)
    x = 1
    for k in range(255):
        exp[k] = x
        log[x] = k
        x <<= 1
        if x >= 256:
            x ^= prim
    for k in range(255, 512):
        exp[k] = exp[k - 255]
    mul = np.zeros((256, 256), dtype=np.int32)
    for aa in range(1, 256):
        mul[aa, 1:] = exp[log[aa] + log[1:]]
    inv = np.zeros(256, dtype=np.int32)
    inv[1:] = exp[255 - log[1:]]
    return {"EXP": exp, "LOG": log, "MUL": mul, "INV": inv}


GF = _build_gf256()


def gf_mul(A, B):
    """Elementwise GF(256) product with numpy broadcasting."""
    A = np.asarray(A, dtype=np.int32)
    B = np.asarray(B, dtype=np.int32)
    A2, B2 = np.broadcast_arrays(A, B)
    out = np.zeros(A2.shape, dtype=np.int32)
    nz = (A2 != 0) & (B2 != 0)
    out[nz] = GF["EXP"][GF["LOG"][A2[nz]] + GF["LOG"][B2[nz]]]
    return out


def gf_xor_cols(M):
    """XOR-reduce columns of M (R x C) -> R x 1."""
    M = np.asarray(M, dtype=np.int32).copy()
    while M.shape[1] > 1:
        c = M.shape[1]
        h = c // 2
        M = np.hstack([M[:, :h] ^ M[:, h:2 * h], M[:, 2 * h:]])
    return M


def _hadamard256():
    H = np.array([[1]], dtype=np.float64)
    for _ in range(8):
        H = np.block([[H, H], [H, -H]])
    return H


HADAMARD = _hadamard256()


# ---------------------------------------------------------------------------
# Code construction
# ---------------------------------------------------------------------------

def make_ldpc(N, M, dv=2, seed=1, cache_dir="results"):
    """Random LDPC code over GF(256): length N, M checks, column weight dv.

    Check-node degrees are as equal as possible. Non-zero edge coefficients are
    uniform in 1..255. A systematic encoder is derived by Gauss-Jordan over GF(256).
    The code (without decoder index tables) is cached under cache_dir.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"ldpc_N{N}_M{M}_dv{dv}_s{seed}.npz"
    if cache.is_file():
        data = np.load(cache, allow_pickle=False)
        code = {k: data[k] for k in data.files}
        return _add_indices(code)

    rng = np.random.default_rng(seed)
    E = N * dv
    base, extra = divmod(E, M)
    deg = np.full(M, base, dtype=np.int32)
    if extra:
        deg[rng.permutation(M)[:extra]] = base + 1

    sockets = np.repeat(np.arange(M, dtype=np.int32), deg)
    edge_chk = rng.permutation(sockets)
    edge_var = np.arange(E, dtype=np.int32) // dv

    # remove duplicate (variable, check) pairs by random swaps
    dup = np.array([0])  # non-empty so the loop runs at least once
    for _ in range(5000):
        key = edge_var * M + edge_chk
        _, first = np.unique(key, return_index=True)
        keep = np.zeros(E, dtype=bool)
        keep[first] = True
        dup = np.flatnonzero(~keep)
        if dup.size == 0:
            break
        swap_with = rng.integers(0, E, size=dup.size)
        edge_chk[dup], edge_chk[swap_with] = edge_chk[swap_with], edge_chk[dup].copy()
    if dup.size:
        raise RuntimeError("make_ldpc: could not remove duplicate edges")

    coef = rng.integers(1, 256, size=E, dtype=np.int32)

    # Gauss-Jordan over GF(256) -> systematic encoder
    A = np.zeros((M, N), dtype=np.int32)
    A[edge_chk, edge_var] = coef
    piv = []
    p = 0
    for c in range(N):
        if p == M:
            break
        rws = np.flatnonzero(A[p:M, c] != 0)
        if rws.size == 0:
            continue
        r0 = p + int(rws[0])
        if r0 != p:
            A[[p, r0]] = A[[r0, p]]
        A[p] = gf_mul(GF["INV"][A[p, c]], A[p])
        rows = np.flatnonzero(A[:, c] != 0)
        rows = rows[rows != p]
        if rows.size:
            f = A[rows, c]
            A[rows] = A[rows] ^ gf_mul(f[:, None], A[p][None, :])
        piv.append(c)
        p += 1

    piv = np.asarray(piv, dtype=np.int32)
    rk = len(piv)
    nonpiv = np.setdiff1d(np.arange(N, dtype=np.int32), piv, assume_unique=False)
    code = {
        "N": np.int32(N),
        "M": np.int32(M),
        "dv": np.int32(dv),
        "E": np.int32(E),
        "seed": np.int32(seed),
        "edge_var": edge_var,
        "edge_chk": edge_chk,
        "coef": coef,
        "rank": np.int32(rk),
        "K": np.int32(N - rk),
        "rate": np.float64((N - rk) / N),
        "piv": piv,
        "nonpiv": nonpiv,
        "P": A[:rk, nonpiv].astype(np.int32),
        "deg": deg,
    }
    np.savez(cache, **code)
    return _add_indices(code)


def _add_indices(code):
    """Index tables used by the FFT-SPA decoder (cheap to rebuild)."""
    E = int(code["E"])
    M = int(code["M"])
    deg = np.asarray(code["deg"], dtype=np.int32)
    dcmax = int(deg.max())
    edge_chk = np.asarray(code["edge_chk"], dtype=np.int32)
    coef = np.asarray(code["coef"], dtype=np.int32)

    chk_edges = np.full((M, dcmax), E, dtype=np.int32)  # E = pad / dummy edge
    edge_slot = np.zeros(E, dtype=np.int32)             # C-order index into M x dcmax
    for j in range(M):
        eds = np.flatnonzero(edge_chk == j)
        chk_edges[j, :eds.size] = eds
        edge_slot[eds] = j * dcmax + np.arange(eds.size, dtype=np.int32)

    code = dict(code)
    code["chk_edges"] = chk_edges
    code["edge_slot"] = edge_slot
    code["dcmax"] = np.int32(dcmax)
    code["hinv"] = GF["INV"][coef]
    return code


# ---------------------------------------------------------------------------
# Encode / syndrome check
# ---------------------------------------------------------------------------

def ldpc_encode(code, info):
    """Systematic encoding. info: length-K symbols (0..255). Returns length-N codeword."""
    info = np.asarray(info, dtype=np.int32).ravel()
    c = np.zeros(int(code["N"]), dtype=np.int32)
    c[code["nonpiv"]] = info
    prodm = gf_mul(code["P"], info[None, :])          # rank x K
    par = gf_xor_cols(prodm).ravel()
    c[code["piv"]] = par
    return c


def ldpc_check(code, c):
    """True if all parity checks are satisfied by symbol vector c (0..255)."""
    c = np.asarray(c, dtype=np.int32).ravel()
    pr = gf_mul(code["coef"], c[code["edge_var"]])    # E
    pv = np.concatenate([pr, [0]])                     # E+1, pad = 0
    s = gf_xor_cols(pv[code["chk_edges"]])             # M x 1
    return bool(np.all(s == 0))


# ---------------------------------------------------------------------------
# Soft decoder (FFT-SPA)
# ---------------------------------------------------------------------------

def ldpc_decode(code, gamma, max_iter=20):
    """FFT-based sum-product decoder over GF(256).

    gamma : N x 256 channel likelihoods (rows need not be normalised)
    returns (chat, ok, iters, ext) where
      chat : length-N hard decisions (0..255)
      ok   : True if syndrome is zero
      iters: iterations used
      ext  : N x 256 extrinsic message (product of check->variable messages),
             i.e. what the outer code tells the inner decoder (used as its prior)
    """
    q = 256
    N = int(code["N"])
    E = int(code["E"])
    dv = int(code["dv"])
    M = int(code["M"])
    dcmax = int(code["dcmax"])
    Hd = HADAMARD
    hinv = code["hinv"]
    coef = np.asarray(code["coef"], dtype=np.int32)
    chk_edges = code["chk_edges"]
    edge_slot = code["edge_slot"]
    u_vals = np.arange(q, dtype=np.int32)

    gamma = np.asarray(gamma, dtype=np.float64)
    gamma = np.maximum(gamma, 1e-12)
    gamma = gamma / gamma.sum(axis=1, keepdims=True)
    Gam3 = gamma[None, :, :]                          # 1 x N x q

    cv = np.full((E, q), 1.0 / q)                      # check -> variable
    chat = np.zeros(N, dtype=np.int32)
    ok = False
    iters = 0

    # permutation tables: c = hinv * u, and u = h * c
    c_from_u = GF["MUL"][hinv[:, None], u_vals[None, :]]       # E x q
    u_from_c = GF["MUL"][coef[:, None], u_vals[None, :]]       # E x q
    edge_rows = np.arange(E)[:, None]

    for it in range(1, max_iter + 1):
        iters = it
        # ---- variable-node update ----------------------------------------
        CV = cv.reshape(N, dv, q).transpose(1, 0, 2)            # dv x N x q
        VC = np.empty_like(CV)
        for k in range(dv):
            Pk = Gam3.copy()
            for k2 in range(dv):
                if k2 != k:
                    Pk = Pk * CV[k2]
            VC[k] = Pk[0]
        vc = VC.transpose(1, 0, 2).reshape(E, q)
        vc = vc / vc.sum(axis=1, keepdims=True)

        # ---- check-node update (FFT-SPA) ---------------------------------
        vcx = np.vstack([vc, np.concatenate([[1.0], np.zeros(q - 1)])[None, :]])
        PU = vcx[edge_rows, c_from_u]                          # E x q
        T = PU @ Hd
        Tx = np.vstack([T, np.ones((1, q))])                   # E+1 x q
        TT = Tx[chk_edges]                                     # M x dcmax x q

        pre = np.ones_like(TT)
        suf = np.ones_like(TT)
        for k in range(1, dcmax):
            pre[:, k, :] = pre[:, k - 1, :] * TT[:, k - 1, :]
        for k in range(dcmax - 2, -1, -1):
            suf[:, k, :] = suf[:, k + 1, :] * TT[:, k + 1, :]
        LOO = pre * suf
        Rm = (LOO.reshape(M * dcmax, q) @ Hd) / q
        PUo = Rm[edge_slot]
        cvn = PUo[edge_rows, u_from_c]
        cvn = np.maximum(cvn, 1e-13)
        cv = cvn / cvn.sum(axis=1, keepdims=True)

        # ---- hard decision + syndrome ------------------------------------
        CV = cv.reshape(N, dv, q).transpose(1, 0, 2)
        B = Gam3.copy()
        for k in range(dv):
            B = B * CV[k]
        B = B.reshape(N, q)
        chat = np.argmax(B, axis=1).astype(np.int32)
        if ldpc_check(code, chat):
            ok = True
            break

    CV = cv.reshape(N, dv, q).transpose(1, 0, 2)
    Pe = np.ones((1, N, q))
    for k in range(dv):
        Pe = Pe * CV[k]
    ext = Pe.reshape(N, q)
    ext = ext / ext.sum(axis=1, keepdims=True)
    return chat, ok, iters, ext


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # GF tables: inverses and distributivity
    for a in range(1, 256):
        assert GF["MUL"][a, GF["INV"][a]] == 1
    a, b, c = 37, 201, 94
    assert GF["MUL"][a, b ^ c] == (GF["MUL"][a, b] ^ GF["MUL"][a, c])
    print("[ok] GF(256) tables: inverses, distributivity")

    code = make_ldpc(256, 128, 2, seed=1, cache_dir="results_selftest")
    rng = np.random.default_rng(0)
    for _ in range(3):
        cw = ldpc_encode(code, rng.integers(0, 256, size=int(code["K"])))
        assert ldpc_check(code, cw)
    print("[ok] systematic LDPC encoder produces valid codewords")
    print(f"     N={code['N']} M={code['M']} K={code['K']} rate={float(code['rate']):.3f} rank={code['rank']}")

    # soft decode of q-ary symmetric errors
    N = 256
    pe = 0.12
    nerr = 0
    for k in range(4):
        cw = ldpc_encode(code, rng.integers(0, 256, size=int(code["K"])))
        rx = cw.copy()
        flip = rng.random(N) < pe
        nflip = int(flip.sum())
        if nflip:
            rx[flip] = (rx[flip] + rng.integers(1, 256, size=nflip)) % 256
        gam = np.full((N, 256), pe / 255.0)
        gam[np.arange(N), rx] = 1.0 - pe
        chat, ok, iters, _ = ldpc_decode(code, gam, max_iter=30)
        nerr += int(not np.array_equal(chat, cw))
        print(f"     trial {k}: ok={ok} iters={iters} match={np.array_equal(chat, cw)}")
    assert nerr == 0
    print("[ok] FFT-SPA corrects 12% q-ary symmetric symbol errors")
    print("All LDPC self-tests passed.")
