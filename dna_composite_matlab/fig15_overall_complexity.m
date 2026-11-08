function fig15_overall_complexity(cfg)
% FIG15_OVERALL_COMPLEXITY  Paper Fig. 15 / Table 5: time and memory of the whole concatenated
%   decoder (inner APP + outer FFT-SPA) for seven operating points (a)-(g).
%   Memory = array bytes of one inner word's messages (x #inner words decoded at once = 1 here)
%   plus an estimate of the FFT-SPA working set. Times in seconds per block.
if nargin < 1, cfg = get_cfg('quick'); end
rr = [50 50 50 100 100 100 100];   MM = [640 768 896 512 640 768 896];
aa = [0.25 0.60 0.70 0.65 0.75 0.90 1.05];
THETA = create_alphabet(get_weights('w41'));  n = cfg.n;
tin = zeros(1, 7);  tout = tin;  rounds = tin;  mem_in = tin;  mem_out = tin;
for k = 1:7
    code = make_ldpc(cfg.N, MM(k), cfg.dv, cfg.ldpcSeed, cfg.cacheDir);
    cp = channel_params(1, aa(k), n);
    for b = 1:cfg.cx.nBlocks
        res = sim_concat_block(code, THETA, cp, rr(k), n, cfg.maxRounds, cfg.maxIterOuter);
        tin(k) = tin(k) + res.tin / cfg.cx.nBlocks;  tout(k) = tout(k) + res.tout / cfg.cx.nBlocks;
        rounds(k) = rounds(k) + res.rounds / cfg.cx.nBlocks;
    end
    X = composite_encode(THETA, randi(256, 1, n), rr(k));
    Y = channel_sim(X, cp.pins, cp.pdel, cp.Dmin, cp.Dmax, cp.Psub);
    [~, st] = app_calc(Y, THETA, [], cp.pins, cp.pdel, cp.Dmin, cp.Dmax, cp.Psub, n);
    mem_in(k)  = st.bytes / 2^20;
    mem_out(k) = 8 * 256 * (8*code.E + 5*code.M*code.dcmax + 2*code.E) / 2^20;   % estimate
    fprintf('  (%c) r=%3d M=%d alpha=%.2f : inner %.1fs outer %.1fs  avg rounds %.2f  mem %.1f/%.1f MiB\n', ...
        'a'+k-1, rr(k), MM(k), aa(k), tin(k), tout(k), rounds(k), mem_in(k), mem_out(k));
end
save(fullfile(cfg.cacheDir, 'complexity_fig15.mat'), 'tin', 'tout', 'rounds', 'mem_in', 'mem_out', 'rr', 'MM', 'aa');
f = figure('Position', [100 100 900 420]);
subplot(1, 2, 1);  bar([tin; tout].', 'stacked');  title('computation time per block');  ylabel('time (sec)');
set(gca, 'XTickLabel', num2cell('abcdefg'));  legend({'inner', 'outer'}, 'Location', 'northwest');
subplot(1, 2, 2);  bar([mem_in; mem_out].', 'stacked');  title('memory requirement');  ylabel('size (MiB)');
set(gca, 'XTickLabel', num2cell('abcdefg'));  legend({'inner', 'outer'}, 'Location', 'northwest');
save_fig(f, cfg, 'fig15_overall_complexity');
end
