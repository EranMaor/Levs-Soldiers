function [chat, ok, iters, ext] = ldpc_decode(code, gamma, maxIter)
% LDPC_DECODE  FFT-based sum-product decoder over GF(256) (Walsh-Hadamard transform).
%   gamma : N x 256 channel likelihoods (rows need not be normalised)
%   chat  : N x 1 hard decisions (0..255);  ok: syndrome satisfied
%   ext   : N x 256 extrinsic message (product of check->variable messages), i.e.
%           what the outer code "tells" the inner decoder (used as its prior).
q = 256;  N = code.N;  E = code.E;  dv = code.dv;  M = code.M;  dcmax = code.dcmax;
persistent Hd
if isempty(Hd), Hd = 1; for b = 1:8, Hd = [Hd Hd; Hd -Hd]; end, end
gamma = max(gamma, 1e-12);
gamma = bsxfun(@rdivide, gamma, sum(gamma, 2));
Gam3  = reshape(gamma, [1, N, q]);
cv    = ones(E, q) / q;                 % check -> variable messages
chat  = zeros(N, 1);  ok = false;  iters = 0;
for it = 1:maxIter
    iters = it;
    % ---- variable-node update ---------------------------------------------
    CV = reshape(cv, [dv, N, q]);
    VC = zeros(dv, N, q);
    for k = 1:dv
        Pk = Gam3;
        for k2 = [1:k-1, k+1:dv], Pk = Pk .* CV(k2, :, :); end
        VC(k, :, :) = Pk;
    end
    vc = reshape(VC, [E, q]);
    vc = bsxfun(@rdivide, vc, sum(vc, 2));
    % ---- check-node update (FFT-SPA) ---------------------------------------
    vcx = [vc; [1, zeros(1, q-1)]];
    PU  = vcx(code.Lidx1);                          % messages about u = h*c
    T   = PU * Hd;                                  % Walsh-Hadamard transform
    Tx  = [T; ones(1, q)];
    TT  = reshape(Tx(code.chkEdges(:), :), [M, dcmax, q]);
    pre = ones(size(TT));  suf = ones(size(TT));
    for k = 2:dcmax,          pre(:, k, :) = pre(:, k-1, :) .* TT(:, k-1, :); end
    for k = dcmax-1:-1:1,     suf(:, k, :) = suf(:, k+1, :) .* TT(:, k+1, :); end
    LOO = pre .* suf;
    Rm  = (reshape(LOO, [M*dcmax, q]) * Hd) / q;    % inverse transform
    PUo = Rm(code.posOfEdge, :);
    cvn = PUo(code.Lidx2);                          % back to c-domain
    cvn = max(cvn, 1e-13);
    cv  = bsxfun(@rdivide, cvn, sum(cvn, 2));
    % ---- hard decision + syndrome ------------------------------------------
    CV = reshape(cv, [dv, N, q]);
    B = Gam3;  for k = 1:dv, B = B .* CV(k, :, :); end
    B = reshape(B, [N, q]);
    [~, ix] = max(B, [], 2);
    chat = ix - 1;
    if ldpc_check(code, chat), ok = true; break; end
end
CV = reshape(cv, [dv, N, q]);
Pe = ones(1, N, q);  for k = 1:dv, Pe = Pe .* CV(k, :, :); end
ext = reshape(Pe, [N, q]);
ext = bsxfun(@rdivide, ext, sum(ext, 2));
end
