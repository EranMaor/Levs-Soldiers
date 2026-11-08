function I = estimate_mi(ctrue, cdec, m, correct)
% ESTIMATE_MI  I(C;U) in bits between uniform symbols C and MAP decisions U (Eq. 47),
%   from the pooled m x m confusion matrix of (c_i, u_i). H(C) = log2(m) exactly.
%   correct = true (default) applies a Miller-Madow bias correction to the entropy
%   estimates (the plug-in estimate is biased upward for sparse 256x256 histograms).
if nargin < 4, correct = true; end
J  = accumarray([ctrue(:), cdec(:)], 1, [m m]);
N  = sum(J(:));  P = J / N;
pu = sum(P, 1);
Hcu = -sum(P(P > 0) .* log2(P(P > 0)));
Hu  = -sum(pu(pu > 0) .* log2(pu(pu > 0)));
if correct
    Hcu = Hcu + (nnz(P)  - 1) / (2*N*log(2));
    Hu  = Hu  + (nnz(pu) - 1) / (2*N*log(2));
end
I = log2(m) + Hu - Hcu;
end
