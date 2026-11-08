function THETA = create_alphabet(w)
% CREATE_ALPHABET  Superposition composite alphabet (paper Eq. 45).
%   w     : layer weights, sum(w) = 1
%   THETA : m x 4 matrix, one row per DISTINCT composite vector (duplicates merge,
%           so e.g. w = [1/8 1/8 3/8 3/8] gives m = 100, as in Table 3).
w = w(:).';
if abs(sum(w) - 1) > 1e-9, error('weights must sum to 1 (got %.12f)', sum(w)); end
L = numel(w);
combos = zeros(4^L, L);
for l = 1:L
    combos(:, l) = repmat(kron((0:3).', ones(4^(L-l), 1)), 4^(l-1), 1);
end
V = zeros(4^L, 4);
for l = 1:L
    V = V + w(l) * (bsxfun(@eq, combos(:, l), 0:3));
end
[~, ia] = unique(round(V * 1e9), 'rows', 'first');
THETA = V(sort(ia), :);
end
