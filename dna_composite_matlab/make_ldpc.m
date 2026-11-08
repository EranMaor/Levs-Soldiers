function code = make_ldpc(N, M, dv, seed, cacheDir)
% MAKE_LDPC  Random LDPC code over GF(256): length N, M checks, column weight dv.
%   Check-node degrees are as equal as possible ("regular" when dv*N/M is an
%   integer). Non-zero coefficients are uniform in 1..255. A systematic
%   encoder is derived by Gauss-Jordan elimination over GF(256) (cached on disk).
if nargin < 5, cacheDir = 'results'; end
if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
cache = fullfile(cacheDir, sprintf('ldpc_N%d_M%d_dv%d_s%d.mat', N, M, dv, seed));
if exist(cache, 'file')
    S = load(cache);  code = add_indices(S.code);  return;
end
G = gf256();
set_seed(seed);
E = N*dv;
base = floor(E/M);  extra = E - base*M;
deg = base*ones(1, M);  deg(randperm(M, extra)) = base + 1;
sockets = zeros(1, E);  p0 = 0;
for j = 1:M, sockets(p0 + (1:deg(j))) = j;  p0 = p0 + deg(j); end
edgeChk = sockets(randperm(E)).';
edgeVar = ceil((1:E).'/dv);               % edge e = (v-1)*dv + k
% remove duplicate (variable, check) pairs by random swaps
for iter = 1:5000
    key = (edgeVar-1)*M + edgeChk;
    [~, ia] = unique(key, 'first');
    dup = setdiff((1:E).', ia);
    if isempty(dup), break; end
    for e = dup.'
        f = randi(E);
        tmp = edgeChk(e); edgeChk(e) = edgeChk(f); edgeChk(f) = tmp;
    end
end
if ~isempty(dup), error('make_ldpc: could not remove duplicate edges'); end
coef = randi(255, E, 1);

% ---- Gauss-Jordan elimination over GF(256) -> systematic encoder ----------
A = zeros(M, N);
A(sub2ind([M, N], edgeChk, edgeVar)) = coef;
piv = zeros(1, M);  p = 0;
for c = 1:N
    if p == M, break; end
    rws = find(A(p+1:M, c) ~= 0);
    if isempty(rws), continue; end
    r0 = p + rws(1);  p = p + 1;
    if r0 ~= p, tmp = A(p, :); A(p, :) = A(r0, :); A(r0, :) = tmp; end
    A(p, :) = gf_helpers('mul', G.INV(A(p, c)+1), A(p, :));    % pivot -> 1
    rows = find(A(:, c) ~= 0);  rows(rows == p) = [];
    if ~isempty(rows)
        f = A(rows, c);
        A(rows, :) = bitxor(A(rows, :), gf_helpers('mul', f, A(p, :)));
    end
    piv(p) = c;
end
rk = p;
piv = piv(1:rk);
nonpiv = setdiff(1:N, piv);
code.N = N; code.M = M; code.dv = dv; code.E = E; code.seed = seed;
code.edgeVar = edgeVar; code.edgeChk = edgeChk; code.coef = coef;
code.rank = rk; code.K = N - rk; code.rate = (N - rk)/N;
code.piv = piv; code.nonpiv = nonpiv; code.P = A(1:rk, nonpiv);
code.deg = deg;
save(cache, 'code', '-v7');
code = add_indices(code);
end

function code = add_indices(code)
% Index tables used by the FFT-SPA decoder (cheap to rebuild, not cached).
G = gf256();  q = 256;
E = code.E; M = code.M;
dcmax = max(code.deg);
chkEdges = (E+1) * ones(M, dcmax);
posOfEdge = zeros(E, 1);
for j = 1:M
    eds = find(code.edgeChk == j);
    chkEdges(j, 1:numel(eds)) = eds.';
    posOfEdge(eds) = j + M*((1:numel(eds)).' - 1);
end
hinv = G.INV(code.coef + 1).';  hinv = hinv(:);
e = (1:E).';
cidx = G.MUL(hinv + 1, :);                      % E x q : c = hinv * u
code.Lidx1 = sub2ind([E+1, q], repmat(e, 1, q), cidx + 1);
uidx = G.MUL(code.coef + 1, :);                 % E x q : u = h * c
code.Lidx2 = sub2ind([E, q], repmat(e, 1, q), uidx + 1);
code.chkEdges = chkEdges;  code.posOfEdge = posOfEdge;  code.dcmax = dcmax;
end
