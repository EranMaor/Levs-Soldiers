function X = composite_encode(THETA, cidx, r)
% COMPOSITE_ENCODE  Synthesis of r strands from composite symbols.
%   cidx : 1 x n row indices into THETA (outer symbol -> composite vector, f = identity)
%   X    : r x n matrix with entries in 0..3; x_i^t ~ Categorical(theta_i), independent
%          across strands (the "probabilistic synthesis" of the paper, Eq. 12).
n   = numel(cidx);
cum = cumsum(THETA(cidx, :), 2);            % n x 4
U   = rand(r, n);
X   = double(bsxfun(@gt, U, cum(:, 1).')) + double(bsxfun(@gt, U, cum(:, 2).')) ...
    + double(bsxfun(@gt, U, cum(:, 3).'));
end
