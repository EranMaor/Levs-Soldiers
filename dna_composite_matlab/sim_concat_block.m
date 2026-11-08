function res = sim_concat_block(code, THETA, cp, r, n, maxRounds, maxIterOuter)
% SIM_CONCAT_BLOCK  One block of the concatenated scheme of Fig. 1 / Table 4.
%   random information -> systematic LDPC encoding over GF(256) -> N/n inner words ->
%   composite-symbol synthesis (r strands each) -> IDS channel -> iterative decoding:
%   inner APP (app_calc) -> outer FFT-SPA -> extrinsic feedback as inner prior -> ...
%   Stops when the outer syndrome is zero or after maxRounds inner-outer iterations.
N = code.N;  m = size(THETA, 1);  W = N / n;
info = randi([0 255], 1, code.K);
c    = ldpc_encode(code, info);
Ys   = cell(W, 1);
for w = 1:W
    X = composite_encode(THETA, c((w-1)*n + (1:n)) + 1, r);   % f: GF symbol s -> THETA row s+1
    Ys{w} = channel_sim(X, cp.pins, cp.pdel, cp.Dmin, cp.Dmax, cp.Psub);
end
prior = ones(N, m) / m;
tin = 0;  tout = 0;  rounds = 0;  chat = zeros(N, 1);
for rd = 1:maxRounds
    rounds = rd;
    t0 = tic;
    G = cell(W, 1);
    parfor w = 1:W
        G{w} = app_calc(Ys{w}, THETA, prior((w-1)*n + (1:n), :), cp.pins, cp.pdel, cp.Dmin, cp.Dmax, cp.Psub, n);
    end
    gamma = vertcat(G{:});
    tin = tin + toc(t0);
    t0 = tic;
    [chat, ok, ~, ext] = ldpc_decode(code, gamma, maxIterOuter);
    tout = tout + toc(t0);
    if ok, break; end
    prior = ext;
end
res.err = ~isequal(chat(:).', c);
res.rounds = rounds;  res.tin = tin;  res.tout = tout;
res.symerr = nnz(chat(:).' ~= c);
end
