function fig9_error_types(cfg)
% FIG9_ERROR_TYPES  Paper Fig. 9: effect of deletion (CH-2), insertion (CH-3) and
%   substitution (CH-4) errors on I(C;U), w41, r = 50 and 100.
if nargin < 1, cfg = get_cfg('quick'); end
names = {'Del', 'Ins', 'Sub'};  chans = [2 3 4];
f = figure('Position', [100 100 900 420]);
for q = 1:numel(cfg.mi.rsel)
    r = cfg.mi.rsel(q);
    subplot(1, numel(cfg.mi.rsel), q);  hold on;
    for k = 1:3
        [c, mk] = plot_style(k);
        I = mi_curve(cfg, chans(k), 'w41', r, cfg.mi.alphas);
        plot(cfg.mi.alphas, I, ['-' mk], 'Color', c, 'MarkerSize', 5);
    end
    axis([0 1.4 0 7]);  grid on;  title(sprintf('r=%d', r));
    xlabel('Error scale \alpha');  ylabel('Mutual information (bit)');
    legend(names, 'Location', 'southwest');
end
save_fig(f, cfg, 'fig9_error_types');
end
