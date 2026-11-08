function fig8_mutual_information(cfg)
% FIG8_MUTUAL_INFORMATION  Paper Fig. 8: I(C;U) versus error scale alpha for CH-1..CH-5,
%   weight vector w41, several sequencing depths r.
if nargin < 1, cfg = get_cfg('quick'); end
ttl = {'CH-1: DNA channel model', 'CH-2: Deletion channel', 'CH-3: Insertion channel', ...
       'CH-4: Substitution channel', 'CH-5: IDS channel'};
f = figure('Position', [100 100 850 1000]);
for ch = 1:5
    subplot(3, 2, ch);  hold on;
    for k = 1:numel(cfg.mi.rlist)
        [c, mk] = plot_style(k);
        I = mi_curve(cfg, ch, 'w41', cfg.mi.rlist(k), cfg.mi.alphas);
        plot(cfg.mi.alphas, I, ['-' mk], 'Color', c, 'MarkerSize', 5);
    end
    axis([0 1.4 0 7]);  grid on;  title(ttl{ch});
    xlabel('Error scale \alpha');  ylabel('Mutual information (bit)');
end
subplot(3, 2, 6);  hold on;
for k = 1:numel(cfg.mi.rlist)
    [c, mk] = plot_style(k);  plot(nan, nan, ['-' mk], 'Color', c);
end
legend(arrayfun(@(r) sprintf('r=%d', r), cfg.mi.rlist, 'UniformOutput', false), 'Location', 'west');
axis off;
save_fig(f, cfg, 'fig8_mutual_information');
end
