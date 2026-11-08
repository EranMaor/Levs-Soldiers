function [c, mk] = plot_style(k)
% PLOT_STYLE  k-th colour / marker (cycles).
C = [0.58 0.00 0.83; 0.00 0.62 0.45; 0.34 0.71 0.91; 0.90 0.62 0.00; 0.94 0.89 0.26; 0.00 0.45 0.70; 0.84 0.37 0.00; 0.00 0.00 0.00];
M = {'+', 'x', '*', 's', 'o', 'd', '^', 'v'};
k = mod(k-1, size(C,1)) + 1;
c = C(k, :);  mk = M{k};
end
