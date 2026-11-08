function I = mi_curve(cfg, chan, wname, r, alphas)
% MI_CURVE  I(C;U) versus alpha for one (channel, weight vector, r). Results are cached
%   per point in cfg.cacheDir so interrupted runs resume and figures can share data.
THETA = create_alphabet(get_weights(wname));
f = fullfile(cfg.cacheDir, sprintf('mi_ch%d_%s_r%d_w%d.mat', chan, wname, r, cfg.mi.nWords));
S.alphas = []; S.I = [];
if exist(f, 'file'), S = load(f); end
I = nan(size(alphas));
for k = 1:numel(alphas)
    j = find(abs(S.alphas - alphas(k)) < 1e-9, 1);
    if ~isempty(j), I(k) = S.I(j); continue; end
    t0 = tic;
    I(k) = sim_mi_point(THETA, chan, alphas(k), r, cfg.mi.nWords, cfg.n);
    S.alphas(end+1) = alphas(k);  S.I(end+1) = I(k);
    save(f, '-struct', 'S', '-v7');
    fprintf('  CH-%d %s r=%3d alpha=%.2f  I=%.3f bit  (%.0f s)\n', chan, wname, r, alphas(k), I(k), toc(t0));
end
end
