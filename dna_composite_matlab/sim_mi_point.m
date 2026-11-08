function [I, tword, bytes] = sim_mi_point(THETA, chan, alpha, r, nWords, n)
% SIM_MI_POINT  Estimate I(C;U) (bits/symbol) for one operating point.
%   Symbols c_i are uniform over the m composite symbols, no outer code; U is the MAP
%   decision from the inner APP (Eq. 46). Also returns seconds per word and the size
%   (bytes) of the main message arrays of app_calc.
cp = channel_params(chan, alpha, n);
m  = size(THETA, 1);
ct = zeros(nWords, n);  cd = zeros(nWords, n);  tw = zeros(nWords, 1);  by = zeros(nWords, 1);
parfor w = 1:nWords
    c  = randi(m, 1, n);
    X  = composite_encode(THETA, c, r);
    Y  = channel_sim(X, cp.pins, cp.pdel, cp.Dmin, cp.Dmax, cp.Psub);
    t0 = tic;
    [post, st] = app_calc(Y, THETA, [], cp.pins, cp.pdel, cp.Dmin, cp.Dmax, cp.Psub, n);
    tw(w) = toc(t0);  by(w) = st.bytes;
    [~, u] = max(post, [], 2);
    ct(w, :) = c;  cd(w, :) = u.';
end
I = estimate_mi(ct(:), cd(:), m);
tword = mean(tw);  bytes = mean(by);
end
