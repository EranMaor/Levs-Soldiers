function [post, stats] = app_calc(Y, THETA, prior, pins, pdel, Dmin, Dmax, Psub, n)
% APP_CALC  Inner decoder: Algorithm 1 of Kaneko & Shibata (IEEE Access 2026).
%   Forward-backward (BCJR-type) message passing on the segmented factor graph,
%   vectorised over the r reads and the drift states.
%
%   Y      : 1 x r cell array of received strands (values 0..3, variable length)
%   THETA  : m x 4 composite alphabet (row k = composite vector of outer symbol k)
%   prior  : n x m a-priori probabilities of the outer symbols ([] = uniform)
%   pins, pdel : insertion / deletion probabilities;  Dmin <= 0 <= Dmax drift bounds
%   Psub   : 4 x 4,  Psub(y+1, z+1) = p(y|z)
%   n      : strand length at the channel input
%   post   : n x m, post(i,k) ~ Pr(c_i = k | Y)  (rows sum to 1).  The prior is used
%            INSIDE the message passing (theta -> chi messages) but is excluded from
%            the output, i.e. post is the extrinsic message handed to the outer code.
r  = numel(Y);   m = size(THETA, 1);
D  = Dmin:Dmax;  ND = numel(D);
if isempty(prior), prior = ones(n, m) / m; end

% ---- emission tables ------------------------------------------------------
PAD  = max(abs(Dmin), Dmax) + 3;
L    = cellfun(@numel, Y);
Wd   = max(max(L), n) + 2*PAD + 2;
Q    = zeros(4, Wd, r);                       % Q(x, j+PAD, t) = p(y_t(j) | z = x)
for t = 1:r
    y = Y{t}(:).';
    Q(:, PAD + (1:numel(y)), t) = Psub(y + 1, :).';
end
tail  = L - n;                                % tail drift d_n^t = n'_t - n  (Eq. 9)
valid = (tail >= Dmin) & (tail <= Dmax);      % reads outside the drift range are ignored
idx0  = bsxfun(@plus, 1:n, D(:)) + PAD;       % ND x n : y-index i + d (1-based)
Q0 = reshape(Q(:, idx0(:),     :), [4, ND, n, r]);
Q1 = reshape(Q(:, idx0(:) + 1, :), [4, ND, n, r]);
clear Q
stayv = 1 - pins*(D < Dmax) - pdel*(D > Dmin);   % Eq. (8)
insv  = pins*(D < Dmax);
delv  = pdel*(D > Dmin);
% W0: no drift change (1 output symbol), W1: insertion (2 output symbols)
W0 = permute(bsxfun(@times, reshape(stayv, [1, ND]), Q0),      [1 2 4 3]);   % 4 x ND x r x n
W1 = permute(bsxfun(@times, reshape(insv,  [1, ND]), Q0.*Q1),  [1 2 4 3]);
clear Q0 Q1
del3 = reshape(delv, [1, ND, 1]);             % deletion: no output symbol, no x-dependence

% ---- forward calculation ------------------------------------------------------
alpha = zeros(ND, r, n+1);
alpha(1 - Dmin, :, 1) = 1;                    % initial drift d_0 = 0
for i = 1:n
    a   = alpha(:, :, i);
    W0i = W0(:, :, :, i);   W1i = W1(:, :, :, i);
    Mx  = reshape(sum((W0i + W1i + del3) .* reshape(a, [1, ND, r]), 2), [4, r]);
    mxphi = mix(Mx, THETA, valid, prior(i, :));
    [G0, G1, Gm] = gains(W0i, W1i, delv, mxphi, ND, r);
    na = a .* G0;
    na(2:end, :)   = na(2:end, :)   + a(1:end-1, :) .* G1(1:end-1, :);
    na(1:end-1, :) = na(1:end-1, :) + a(2:end, :)   .* Gm(2:end, :);
    s = sum(na, 1);  s(s == 0) = 1;
    alpha(:, :, i+1) = bsxfun(@rdivide, na, s);
end

% ---- backward calculation ------------------------------------------------------
beta = zeros(ND, r, n+1);
for t = 1:r
    if valid(t), beta(tail(t) - Dmin + 1, t, n+1) = 1;   % tail drift pinned by |y_t|
    else,        beta(:, t, n+1) = 1/ND; end
end
for i = n:-1:1
    b  = beta(:, :, i+1);
    bp = [b(2:end, :); zeros(1, r)];          % bp(k) = b(k+1)
    bm = [zeros(1, r); b(1:end-1, :)];        % bm(k) = b(k-1)
    W0i = W0(:, :, :, i);   W1i = W1(:, :, :, i);
    Mx  = reshape(sum(W0i .* reshape(b,  [1, ND, r]) + W1i .* reshape(bp, [1, ND, r]) ...
                      + del3 .* reshape(bm, [1, ND, r]), 2), [4, r]);
    mxphi = mix(Mx, THETA, valid, prior(i, :));
    [G0, G1, Gm] = gains(W0i, W1i, delv, mxphi, ND, r);
    nb = b .* G0 + bp .* G1 + bm .* Gm;
    s = sum(nb, 1);  s(s == 0) = 1;
    beta(:, :, i) = bsxfun(@rdivide, nb, s);
end

% ---- outbound calculation ------------------------------------------------------
post = zeros(n, m);
for i = 1:n
    a  = alpha(:, :, i);
    b  = beta(:, :, i+1);
    bp = [b(2:end, :); zeros(1, r)];
    bm = [zeros(1, r); b(1:end-1, :)];
    W0i = W0(:, :, :, i);   W1i = W1(:, :, :, i);
    Mx  = reshape(sum(reshape(a, [1, ND, r]) .* (W0i .* reshape(b,  [1, ND, r]) ...
                      + W1i .* reshape(bp, [1, ND, r]) + del3 .* reshape(bm, [1, ND, r])), 2), [4, r]);
    Mchi = THETA * Mx;                                    % m x r
    bad  = ~valid | ~(max(Mchi, [], 1) > 0);
    Mchi(:, bad) = 1;
    Mchi = bsxfun(@rdivide, Mchi, max(Mchi, [], 1));
    lg = sum(log(Mchi), 2);                               % product over ALL reads (log domain)
    mx = max(lg);
    if ~isfinite(mx), p = ones(m, 1); else, p = exp(lg - mx); end
    post(i, :) = (p / sum(p)).';
end
if nargout > 1
    stats.bytes = 8 * (numel(alpha) + numel(beta) + numel(W0) + numel(W1));
    stats.nvalid = nnz(valid);
end
end

% =========================================================================
function mxphi = mix(Mx, THETA, valid, priorRow)
% x -> chi -> theta -> chi -> x : combine the evidence of all OTHER reads (Alg. 1, line 9).
m = size(THETA, 1);
Mchi = THETA * Mx;                                        % m x r  (chi -> theta messages)
bad  = ~valid | ~(max(Mchi, [], 1) > 0);                  % ignore unusable reads
Mchi(:, bad) = 1;
Mchi = bsxfun(@rdivide, Mchi, max(Mchi, [], 1));
loo  = loo_prod(Mchi);                                    % theta -> chi (leave-one-out product)
loo  = bsxfun(@times, loo, priorRow(:));
z = max(loo, [], 1);  z(z == 0) = 1;
loo = bsxfun(@rdivide, loo, z);
mxphi = THETA.' * loo;                                    % chi -> x messages, 4 x r
s = sum(mxphi, 1);  s(s == 0) = 1;
mxphi = bsxfun(@rdivide, mxphi, s);
end

function loo = loo_prod(M)
% product over all columns except the own one (prefix/suffix products, no division)
[m, r] = size(M);
if r == 1, loo = ones(m, 1); return; end
pre = cumprod([ones(m, 1), M(:, 1:r-1)], 2);
cp  = cumprod(fliplr(M(:, 2:r)), 2);
suf = [fliplr(cp), ones(m, 1)];
loo = pre .* suf;
end

function [G0, G1, Gm] = gains(W0i, W1i, delv, mxphi, ND, r)
mx3 = reshape(mxphi, [4, 1, r]);
G0 = reshape(sum(W0i .* mx3, 1), [ND, r]);
G1 = reshape(sum(W1i .* mx3, 1), [ND, r]);
Gm = delv(:) * sum(mxphi, 1);
end
