function fig10_weight_vectors(cfg)
% FIG10_WEIGHT_VECTORS  Paper Fig. 10: I(C;U) for w41, w42, w21 on CH-1, r = 50 and 100.
if nargin < 1, cfg = get_cfg('quick'); end
wn = {'w41', 'w42', 'w21'};
f = figure('Position', [100 100 900 420]);
for q = 1:numel(cfg.mi.rsel)
    r = cfg.mi.rsel(q);
    subplot(1, numel(cfg.mi.rsel), q);  hold on;
    for k = 1:3
        [c, mk] = plot_style(k);
        I = mi_curve(cfg, 1, wn{k}, r, cfg.mi.alphas);
        plot(cfg.mi.alphas, I, ['-' mk], 'Color', c, 'MarkerSize', 5);
    end
    axis([0 1.4 0 7]);  grid on;  title(sprintf('r=%d', r));
    xlabel('Error scale \alpha');  ylabel('Mutual information (bit)');
    legend({'w_{4,1}', 'w_{4,2}', 'w_{2,1}'}, 'Location', 'southwest');
end
save_fig(f, cfg, 'fig10_weight_vectors');
end
