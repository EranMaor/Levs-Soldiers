function S = bler_curve(cfg, M, r, alphas)
% BLER_CURVE  Block error rate versus alpha for outer check length M and depth r (CH-1,
%   w41, random LDPC over GF(2^8)). Runs blocks until cfg.bler.targetErrors errors or
%   cfg.bler.maxBlocks blocks per point; cached per point; stops once BLER > stopBLER.
f = fullfile(cfg.cacheDir, sprintf('bler_M%d_r%d_e%d_b%d.mat', M, r, cfg.bler.targetErrors, cfg.bler.maxBlocks));
S.alphas = alphas(:).';  S.nerr = nan(size(S.alphas));  S.nblk = nan(size(S.alphas));
if exist(f, 'file'), S0 = load(f); if isequal(S0.alphas, S.alphas), S = S0; end, end
code  = make_ldpc(cfg.N, M, cfg.dv, cfg.ldpcSeed, cfg.cacheDir);
THETA = create_alphabet(get_weights('w41'));
n = cfg.n;  B = cfg.bler.batch;
for k = 1:numel(S.alphas)
    if ~isnan(S.nerr(k)), if S.nerr(k)/S.nblk(k) > cfg.bler.stopBLER, break; else, continue; end, end
    cp = channel_params(1, S.alphas(k), n);
    nerr = 0;  nblk = 0;  t0 = tic;
    while nerr < cfg.bler.targetErrors && nblk < cfg.bler.maxBlocks
        errs = zeros(B, 1);
        parfor b = 1:B
            res = sim_concat_block(code, THETA, cp, r, n, cfg.maxRounds, cfg.maxIterOuter);
            errs(b) = res.err;
        end
        nerr = nerr + sum(errs);  nblk = nblk + B;
    end
    S.nerr(k) = nerr;  S.nblk(k) = nblk;
    save(f, '-struct', 'S', '-v7');
    fprintf('  M=%d r=%3d alpha=%.2f : BLER = %d/%d = %.3g  (%.0f s)\n', M, r, S.alphas(k), nerr, nblk, nerr/nblk, toc(t0));
    if nerr/nblk > cfg.bler.stopBLER, break; end
end
end
