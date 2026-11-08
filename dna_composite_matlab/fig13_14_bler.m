function fig13_14_bler(cfg)
% FIG13_14_BLER  Paper Figs. 13-14: block error rate of the concatenated scheme
%   (outer regular LDPC over GF(2^8), N = 2048; inner APP, w41, CH-1) versus alpha.
%   Fig. 13: M in cfg.bler.Ms, r in cfg.bler.rlist.  Fig. 14: M = 768, r in cfg.bler.r14list.
if nargin < 1, cfg = get_cfg('quick'); end
f = figure('Position', [100 100 700 520]);  hold on;
for i = 1:numel(cfg.bler.Ms)
    [c, ~] = plot_style(i);
    for j = 1:numel(cfg.bler.rlist)
        S = bler_curve(cfg, cfg.bler.Ms(i), cfg.bler.rlist(j), cfg.bler.alphas);
        plot_bler(S, c, ternary(cfg.bler.rlist(j) == max(cfg.bler.rlist), '-', '--'));
    end
end
grid on;  set(gca, 'YScale', 'log');  axis([0.2 1.4 1e-5 1]);
xlabel('Error scale \alpha');  ylabel('Block error rate (BLER)');
title('BLER vs error scale (solid: largest r, dashed: smaller r; colours: M)');
legend(arrayfun(@(M) sprintf('M=%d', M), cfg.bler.Ms, 'UniformOutput', false), 'Location', 'southeast');
save_fig(f, cfg, 'fig13_bler');

f = figure('Position', [100 100 700 520]);  hold on;
for j = 1:numel(cfg.bler.r14list)
    [c, ~] = plot_style(j);
    S = bler_curve(cfg, 768, cfg.bler.r14list(j), cfg.bler.alphas);
    plot_bler(S, c, '-');
end
grid on;  set(gca, 'YScale', 'log');  axis([0.2 1.4 1e-5 1]);
xlabel('Error scale \alpha');  ylabel('Block error rate (BLER)');  title('BLER vs error scale, M = 768');
legend(arrayfun(@(r) sprintf('r=%d', r), cfg.bler.r14list, 'UniformOutput', false), 'Location', 'southeast');
save_fig(f, cfg, 'fig14_bler_vs_r');
end

function plot_bler(S, c, ls)
ok = ~isnan(S.nerr);
b  = S.nerr(ok) ./ S.nblk(ok);  a = S.alphas(ok);
pos = b > 0;
plot(a(pos), b(pos), ls, 'Color', c, 'Marker', 'o', 'MarkerSize', 4);          % measured BLER
if any(~pos)                                                                     % no error seen: upper bound 1/blocks
    nb = S.nblk(ok);
    plot(a(~pos), 1 ./ nb(~pos), 'v', 'Color', c, 'MarkerSize', 5);
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
