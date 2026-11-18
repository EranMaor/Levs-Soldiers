"""Mutual-information curves and concatenated LDPC + APP decoding.

Ports estimate_mi.m, sim_mi_point.m, mi_curve.m, sim_concat_block.m, get_cfg.m.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import time

import numpy as np

from app_calc import app_calc
from channel_sim import Channel, channel_params
from composite_alphabet import create_alphabet, get_weights, composite_encode
from ldpc import make_ldpc, ldpc_encode, ldpc_decode


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

@dataclass
class MIConfig:
    n_words: int = 60
    alphas: np.ndarray = field(default_factory=lambda: np.arange(0.0, 1.41, 0.2))
    r_list: tuple = (30, 60, 100)
    r_sel: tuple = (50, 100)


@dataclass
class BLERConfig:
    Ms: tuple = (640, 768, 896)
    r_list: tuple = (50, 100)
    r14_list: tuple = (50, 75, 100)
    alphas: np.ndarray = field(default_factory=lambda: np.arange(0.2, 1.21, 0.2))
    target_errors: int = 5
    max_blocks: int = 20
    stop_bler: float = 0.5


@dataclass
class Config:
    mode: str = "quick"
    cache_dir: str = "results"
    n: int = 128
    N: int = 2048
    dv: int = 2
    ldpc_seed: int = 1
    max_rounds: int = 20
    max_iter_outer: int = 20
    mi: MIConfig = field(default_factory=MIConfig)
    bler: BLERConfig = field(default_factory=BLERConfig)


def get_cfg(mode: str = "quick") -> Config:
    """Simulation settings. mode = 'smoke' | 'quick' | 'full'."""
    cfg = Config(mode=mode)
    if mode == "quick":
        pass  # defaults above match MATLAB quick
    elif mode == "smoke":
        cfg.mi = MIConfig(n_words=4, alphas=np.array([0.0, 1.0]), r_list=(30,), r_sel=(50,))
        cfg.bler = BLERConfig(Ms=(640,), r_list=(50,), r14_list=(50,),
                              alphas=np.array([0.3]), target_errors=1, max_blocks=1)
    elif mode == "full":
        cfg.mi = MIConfig(
            n_words=400,
            alphas=np.arange(0.0, 1.41, 0.1),
            r_list=tuple(range(30, 101, 10)),
            r_sel=(50, 100),
        )
        cfg.bler = BLERConfig(
            Ms=(512, 640, 768, 896),
            r_list=(50, 100),
            r14_list=tuple(range(50, 101, 10)),
            alphas=np.arange(0.2, 1.41, 0.05),
            target_errors=100,
            max_blocks=2000,
        )
    else:
        raise ValueError("mode must be 'smoke', 'quick' or 'full'")
    Path(cfg.cache_dir).mkdir(parents=True, exist_ok=True)
    return cfg


# ---------------------------------------------------------------------------
# Mutual information
# ---------------------------------------------------------------------------

def estimate_mi(ctrue, cdec, m, correct=True):
    """I(C;U) in bits between uniform symbols C and MAP decisions U (Eq. 47)."""
    ctrue = np.asarray(ctrue, dtype=int).ravel()
    cdec = np.asarray(cdec, dtype=int).ravel()
    J = np.zeros((m, m), dtype=np.float64)
    np.add.at(J, (ctrue, cdec), 1)
    N = J.sum()
    P = J / N
    pu = P.sum(axis=0)
    mask = P > 0
    Hcu = -np.sum(P[mask] * np.log2(P[mask]))
    mask_u = pu > 0
    Hu = -np.sum(pu[mask_u] * np.log2(pu[mask_u]))
    if correct:
        Hcu += (np.count_nonzero(P) - 1) / (2 * N * np.log(2))
        Hu += (np.count_nonzero(pu) - 1) / (2 * N * np.log(2))
    return float(np.log2(m) + Hu - Hcu)


def sim_mi_point(THETA, chan, alpha, r, n_words, n, rng=None):
    """Estimate I(C;U) (bits/symbol) for one operating point. No outer code."""
    if rng is None:
        rng = np.random.default_rng()
    cp = channel_params(chan, alpha, n)
    m = len(THETA)
    prior = np.full((n, m), 1.0 / m)
    c_true = np.empty((n_words, n), dtype=np.int32)
    c_dec = np.empty((n_words, n), dtype=np.int32)
    times = np.empty(n_words)

    for w in range(n_words):
        c = rng.integers(0, m, size=n)
        X = composite_encode(THETA, c, r, rng=rng)
        # Channel still uses the legacy global RNG for the IDS draw
        Y = Channel(cp["p_ins"], cp["p_del"], cp["Dmin"], cp["Dmax"], cp["p_sub"]).simulate(X)
        t0 = time.perf_counter()
        post = app_calc(Y, THETA, prior, cp["p_ins"], cp["p_del"],
                        cp["Dmin"], cp["Dmax"], cp["p_sub"])
        times[w] = time.perf_counter() - t0
        c_true[w] = c
        c_dec[w] = post.argmax(axis=1)

    I = estimate_mi(c_true, c_dec, m)
    return I, float(times.mean())


def mi_curve(cfg, chan, wname, r, alphas=None):
    """I(C;U) versus alpha for one (channel, weight vector, r). Caches each point."""
    if alphas is None:
        alphas = cfg.mi.alphas
    alphas = np.asarray(alphas, dtype=float)
    THETA = create_alphabet(get_weights(wname))
    cache = Path(cfg.cache_dir) / f"mi_ch{chan}_{wname}_r{r}_w{cfg.mi.n_words}.npz"

    saved_alphas = np.array([], dtype=float)
    saved_I = np.array([], dtype=float)
    if cache.is_file():
        data = np.load(cache)
        saved_alphas = data["alphas"]
        saved_I = data["I"]

    I = np.full(alphas.shape, np.nan)
    for k, alpha in enumerate(alphas):
        hit = np.where(np.abs(saved_alphas - alpha) < 1e-9)[0]
        if hit.size:
            I[k] = saved_I[hit[0]]
            continue
        t0 = time.perf_counter()
        I[k], _ = sim_mi_point(THETA, chan, float(alpha), r, cfg.mi.n_words, cfg.n)
        saved_alphas = np.append(saved_alphas, alpha)
        saved_I = np.append(saved_I, I[k])
        np.savez(cache, alphas=saved_alphas, I=saved_I)
        print(f"  CH-{chan} {wname} r={r:3d} alpha={alpha:.2f}  "
              f"I={I[k]:.3f} bit  ({time.perf_counter() - t0:.0f} s)")
    return I


# ---------------------------------------------------------------------------
# Concatenated block (LDPC encode + inner/outer iterations)
# ---------------------------------------------------------------------------

def sim_concat_block(code, THETA, cp, r, n, max_rounds=20, max_iter_outer=20, rng=None):
    """One block of the concatenated scheme (Fig. 1 / Table 4).

    random info -> LDPC encode over GF(256) -> N/n inner words -> composite
    synthesis (r strands each) -> IDS channel -> iterate:
      inner APP (app_calc) -> outer FFT-SPA -> extrinsic feedback as inner prior
    Stops when the outer syndrome is zero or after max_rounds rounds.
    """
    if rng is None:
        rng = np.random.default_rng()
    N = int(code["N"])
    m = len(THETA)
    W = N // n
    if N % n != 0:
        raise ValueError(f"N={N} must be divisible by n={n}")
    if m != 256:
        raise ValueError("concatenated scheme expects |Theta|=256 (w41) to match GF(256)")

    info = rng.integers(0, 256, size=int(code["K"]))
    c = ldpc_encode(code, info)

    Ys = []
    for w in range(W):
        symbols = c[w * n:(w + 1) * n]          # GF symbol s -> THETA row s
        X = composite_encode(THETA, symbols, r, rng=rng)
        Y = Channel(cp["p_ins"], cp["p_del"], cp["Dmin"], cp["Dmax"], cp["p_sub"]).simulate(X)
        Ys.append(Y)

    prior = np.full((N, m), 1.0 / m)
    tin = tout = 0.0
    rounds = 0
    chat = np.zeros(N, dtype=np.int32)

    for rd in range(1, max_rounds + 1):
        rounds = rd
        t0 = time.perf_counter()
        G = []
        for w in range(W):
            G.append(app_calc(
                Ys[w], THETA, prior[w * n:(w + 1) * n],
                cp["p_ins"], cp["p_del"], cp["Dmin"], cp["Dmax"], cp["p_sub"],
            ))
        gamma = np.vstack(G)                    # N x 256
        tin += time.perf_counter() - t0

        t0 = time.perf_counter()
        chat, ok, _, ext = ldpc_decode(code, gamma, max_iter=max_iter_outer)
        tout += time.perf_counter() - t0
        if ok:
            break
        prior = ext

    return {
        "err": not np.array_equal(chat, c),
        "rounds": rounds,
        "tin": tin,
        "tout": tout,
        "symerr": int(np.count_nonzero(chat != c)),
        "chat": chat,
        "c": c,
    }


def bler_curve(cfg, M, r, alphas=None):
    """Block error rate versus alpha (CH-1, w41, LDPC over GF(256)).

    Runs blocks until cfg.bler.target_errors errors or cfg.bler.max_blocks blocks
    per point. Cached per (M, r). Stops sweeping alpha once BLER exceeds stop_bler.
    Returns dict with alphas, nerr, nblk.
    """
    if alphas is None:
        alphas = np.asarray(cfg.bler.alphas, dtype=float)
    else:
        alphas = np.asarray(alphas, dtype=float)

    cache = Path(cfg.cache_dir) / (
        f"bler_M{M}_r{r}_e{cfg.bler.target_errors}_b{cfg.bler.max_blocks}.npz"
    )
    nerr = np.full(alphas.shape, np.nan)
    nblk = np.full(alphas.shape, np.nan)
    if cache.is_file():
        data = np.load(cache)
        if data["alphas"].shape == alphas.shape and np.allclose(data["alphas"], alphas):
            nerr = data["nerr"]
            nblk = data["nblk"]

    code = make_ldpc(cfg.N, M, cfg.dv, cfg.ldpc_seed, cfg.cache_dir)
    THETA = create_alphabet(get_weights("w41"))
    n = cfg.n
    rng = np.random.default_rng(cfg.ldpc_seed)

    for k, alpha in enumerate(alphas):
        if np.isfinite(nerr[k]):
            if nerr[k] / nblk[k] > cfg.bler.stop_bler:
                break
            continue
        cp = channel_params(1, float(alpha), n)
        err_count = 0
        blk_count = 0
        t0 = time.perf_counter()
        while err_count < cfg.bler.target_errors and blk_count < cfg.bler.max_blocks:
            res = sim_concat_block(
                code, THETA, cp, r, n,
                max_rounds=cfg.max_rounds,
                max_iter_outer=cfg.max_iter_outer,
                rng=rng,
            )
            err_count += int(res["err"])
            blk_count += 1
        nerr[k] = err_count
        nblk[k] = blk_count
        np.savez(cache, alphas=alphas, nerr=nerr, nblk=nblk)
        bler = err_count / blk_count
        print(f"  M={M} r={r:3d} alpha={alpha:.2f} : "
              f"BLER = {err_count}/{blk_count} = {bler:.3g}  "
              f"({time.perf_counter() - t0:.0f} s)")
        if bler > cfg.bler.stop_bler:
            break

    return {"alphas": alphas, "nerr": nerr, "nblk": nblk}


# ---------------------------------------------------------------------------
# Smoke / wiring test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    cfg = get_cfg("smoke")
    print("== MI smoke (CH-1, w41, r=30, alpha=0, small n) ==")
    THETA = create_alphabet(get_weights("w41"))
    # shorter word for a fast wiring check; paper uses n=128
    I, tword = sim_mi_point(THETA, chan=1, alpha=0.0, r=30, n_words=4, n=32,
                            rng=np.random.default_rng(0))
    print(f"  I={I:.3f} bit/symbol  ({tword:.2f} s/word)")

    print("== concat smoke (small LDPC, one block) ==")
    # N=256, n=64 -> 4 inner words; much faster than paper N=2048
    code = make_ldpc(256, 128, dv=2, seed=1, cache_dir=cfg.cache_dir)
    cp = channel_params(1, 0.3, n=64)
    res = sim_concat_block(code, THETA, cp, r=30, n=64,
                           max_rounds=10, max_iter_outer=20,
                           rng=np.random.default_rng(1))
    print(f"  err={res['err']} rounds={res['rounds']} "
          f"symerr={res['symerr']} tin={res['tin']:.1f}s tout={res['tout']:.1f}s")
    print("Done.")
