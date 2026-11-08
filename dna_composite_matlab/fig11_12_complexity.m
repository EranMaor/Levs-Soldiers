function fig11_12_complexity(cfg)
% FIG11_12_COMPLEXITY  Paper Figs. 11-12: computation time (s/word) and memory (MiB) of the
%   inner APP calculation versus sequencing depth r, CH-1, alpha = 1.0.
%   Absolute numbers differ from the paper (C/OpenMP, 12 threads); the trends are what
%   matter. Memory = size of the main message/factor arrays of app_calc (single word).
if nargin < 1, cfg = get_cfg('quick'); end
combos = {'w41', 128; 'w41', 96; 'w41', 64; 'w42', 128; 'w21', 128};
lab = {'w_{4,1}, n=128', 'w_{4,1}, n=96', 'w_{4,1}, n=64', 'w_{4,2}, n=128', 'w_{2,1}, n=128'};
rl = cfg.cx.rlist;
T = nan(size(combos, 1), numel(rl));  B = T;
for c = 1:size(combos, 1)
    THETA = create_alphabet(get_weights(combos{c, 1}));
    for k = 1:numel(rl)
        [~, T(c, k), by] = sim_mi_point(THETA, 1, 1.0, rl(k), cfg.cx.nWords, combos{c, 2});
        B(c, k) = by / 2^20;
        fprintf('  %s n=%d r=%3d : %.3f s/word, %.1f MiB\n', combos{c, 1}, combos{c, 2}, rl(k), T(c, k), B(c, k));
    end
end
save(fullfile(cfg.cacheDir, 'complexity_fig11_12.mat'), 'T', 'B', 'rl', 'combos');
f1 = figure('Position', [100 100 560 420]);  hold on;
for c = 1:size(T, 1), [col, mk] = plot_style(c);  plot(rl, T(c, :), ['-' mk], 'Color', col); end
grid on;  xlabel('r');  ylabel('time (sec)');  title('Computation time (CH-1, \alpha = 1.0)');  legend(lab, 'Location', 'northwest');
save_fig(f1, cfg, 'fig11_computation_time');
f2 = figure('Position', [100 100 560 420]);  hold on;
for c = 1:size(B, 1), [col, mk] = plot_style(c);  plot(rl, B(c, :), ['-' mk], 'Color', col); end
grid on;  xlabel('r');  ylabel('memory size (MiB)');  title('Memory requirement (CH-1, \alpha = 1.0)');  legend(lab, 'Location', 'northwest');
save_fig(f2, cfg, 'fig12_memory');
end
