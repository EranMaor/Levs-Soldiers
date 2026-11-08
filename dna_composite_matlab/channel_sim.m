function Y = channel_sim(X, pins, pdel, Dmin, Dmax, Psub)
% CHANNEL_SIM  Insertion/deletion/substitution channel of Sec. II-A for all strands at once.
%   X    : r x n matrix (0..3), one strand per row
%   Y    : 1 x r cell array of received strands (row vectors, variable length)
%   Per symbol: deletion w.p. pdel (if d > Dmin), insertion (duplicate) w.p. pins
%   (if d < Dmax), as in Eqs. (7)-(9); then 4-ary memoryless substitution Psub(y+1,z+1).
[r, n] = size(X);
d   = zeros(r, 1);
out = zeros(r, 2*n);
len = zeros(r, 1);
rows = (1:r).';
for i = 1:n
    u       = rand(r, 1);
    candel  = d > Dmin;
    candins = d < Dmax;
    isdel   = candel & (u < pdel);
    isins   = candins & ~isdel & (u < candel*pdel + pins);
    cnt     = 1 - double(isdel) + double(isins);          % symbols emitted: 0, 1 or 2
    d       = d - double(isdel) + double(isins);
    x       = X(:, i);
    w1 = cnt >= 1;  w2 = cnt == 2;
    out(sub2ind([r, 2*n], rows(w1), len(w1) + 1)) = x(w1);
    out(sub2ind([r, 2*n], rows(w2), len(w2) + 2)) = x(w2);
    len = len + cnt;
end
mask = bsxfun(@le, 1:2*n, len);                           % r x 2n, true where a symbol was written
z    = out(mask);
cumP = cumsum(Psub, 1);
U    = rand(1, numel(z));
c    = cumP(:, z.' + 1);                                  % 4 x K
ys   = sum(bsxfun(@gt, U, c(1:3, :)), 1).';
yall = zeros(r, 2*n);  yall(mask) = ys;
Y = cell(1, r);
for t = 1:r, Y{t} = yall(t, 1:len(t)); end
end
