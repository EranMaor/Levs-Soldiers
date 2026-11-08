# MATLAB replication of Kaneko & Shibata, *A Posteriori Probability Calculation for Composite-Symbol DNA Channels With Synchronization Errors* (IEEE Access 2026)

Reproduces the paper's Figs. 8–15 with our own code: inner APP decoder (Algorithm 1),
IDS channel, superposition composite alphabets, and an outer **regular LDPC code over GF(2^8)
with an FFT-based sum-product decoder** (written from scratch, no toolboxes needed).
Runs in MATLAB (R2016b+) and was tested in Octave 8.

## Run
```matlab
selftest            % 1-2 min, checks every building block
mode = 'smoke'; run_all   % seconds, wiring test only
mode = 'quick'; run_all   % tens of minutes, rough curves (default)
mode = 'full';  run_all   % paper-scale sample sizes (hours/days; start a parpool first)
```
Individual figures: `fig8_mutual_information(get_cfg('quick'))`, `fig9_error_types`, `fig10_weight_vectors`,
`fig11_12_complexity`, `fig13_14_bler`, `fig15_overall_complexity`. Everything is cached in `./results`
(per operating point), so interrupted runs resume and Fig. 9/10 reuse Fig. 8 data. `parfor` is used for
words/blocks. Edit `get_cfg.m` to change sample sizes, alpha grids, r lists.

## Files
| file | role |
|---|---|
| `app_calc.m` | **Algorithm 1**: forward/backward/outbound on the segmented factor graph, vectorised over reads and drift states |
| `channel_params.m`, `channel_sim.m` | Table 1/2 channels CH-1..5, drift (Eqs. 7–9) + substitution |
| `create_alphabet.m`, `get_weights.m`, `composite_encode.m` | superposition alphabet (Eq. 45, Table 3), strand synthesis |
| `gf256.m`, `gf_helpers.m`, `make_ldpc.m`, `ldpc_encode.m`, `ldpc_decode.m`, `ldpc_check.m` | GF(256) arithmetic, random LDPC (N=2048, column weight 2), systematic encoder, FFT-SPA |
| `estimate_mi.m`, `sim_mi_point.m`, `mi_curve.m` | I(C;U) (Eq. 47) |
| `sim_concat_block.m`, `bler_curve.m` | concatenated scheme with inner/outer iterations (Table 4) |
| `fig*.m`, `run_all.m`, `get_cfg.m` | figures and settings |

## What was verified
* `app_calc.m` matches the reference Python implementation (`app_calc.py`) to **5e-16** on random cases,
  including the asymmetric Table-2 substitution matrix, m = 100, and the substitution-only (|D| = 1) case.
* Alphabets reproduce Table 3 exactly: |Theta| = 256 / 100 / 16, rates 8.00 / 6.64 / 4.00.
* Noiseless (alpha = 0) mutual information: **6.12 bit (r=50) and 6.74 bit (r=100)** vs ~6.1 / ~6.7 in Figs. 8/10.
* CH-4 (substitution only), alpha = 1.4, r = 100: converges to ~5.5 bit with enough samples (paper ~5.4–5.5).
* LDPC: encoder output satisfies all checks; FFT-SPA corrects q-ary symmetric errors; full-size codes
  (M = 512/640/768/896) have full rank and rates 12/16, 11/16, 10/16, 9/16 as in Table 4.
* Whole chain at full size decodes blocks correctly (inner APP ~80% of the run time, as in the paper).

## Modelling choices and known differences (please read)
1. **Synthesis model.** Each x_i^t is drawn i.i.d. from theta_i, independently across strands (the paper's
   "probabilistic synthesis", Eq. 12 "≃"). This reproduces the paper's alpha = 0 rate loss.
2. **MI estimator bias.** I(C;U) over a 256x256 confusion matrix is biased *upward* for small samples
   (about +0.3 bit at 4k symbols, still +0.05 at 40k). We apply a Miller–Madow correction, but 'quick' mode
   (60 words) is only a rough preview; use >= 400 words per point ('full').
3. **Prior handling.** The a-priori message from the outer decoder is used *inside* the inner decoder
   (theta -> chi messages) and excluded from its output, i.e. the output is the extrinsic message handed to
   the outer code. (The printed line 9 of Algorithm 1 omits the psi_i term; without it feedback would do nothing.)
4. **Tail drift sign.** Line 15 prints `1[d = n - n'_t]`; Eq. (9) implies d_n = n'_t - n. We use Eq. (9)
   (the printed version roughly doubled the symbol error rate in our tests of the Python code).
5. **Outer code.** Random column-weight-2 parity-check matrix with check degrees as equal as possible (exactly
   regular only when 2N/M is an integer, e.g. M = 512), random non-zero coefficients, fixed seed
   (`cfg.ldpcSeed`). Real systematic encoding via GF(256) Gauss–Jordan, not an all-zero-codeword shortcut.
   The outer decoder restarts from scratch in each inner–outer round (the paper does not specify warm starts).
   f (outer symbol -> composite vector) is the identity on THETA's row order.
6. **Reads whose length implies a drift outside [Dmin, Dmax]** are ignored by the decoder (never happens in
   practice with the paper's D = [-ceil(2 n pd), ceil(2 n pi)]).
7. **Complexity figures.** Absolute times/memory are not comparable with the paper's C/OpenMP numbers.
   Times are measured wall-clock; "memory" is the size of app_calc's main arrays (+ an analytic estimate for
   the FFT-SPA working set in Fig. 15). Scaling trends (linear in r, inner decoder dominates) are what to compare.
   Our memory does not depend on m (the paper's does, because it stores per-edge distributions over theta).
8. **Numerics.** Messages are max-normalised, leave-one-out products use prefix/suffix products (no division),
   outbound products are taken in the log domain; FFT-SPA probabilities are floored at 1e-12/1e-13.

## Rough run times (Octave, one core; MATLAB is typically faster)
~0.2–0.4 s per inner word (r = 50–100), ~6–8 s per concatenated block at the paper's operating points.
Full Fig. 8 (600 points x 400 words) is ~17 core-hours; paper-scale BLER (100 errors per point near 1e-4)
needs ~10^6 blocks and is only practical with a cluster — use `cfg.bler.maxBlocks` to cap it (zero-error
points are drawn as 'v' markers at 1/blocks).
