function cp = channel_params(ch, alpha, n)
% CHANNEL_PARAMS  Channels CH-1..CH-5 of Table 1 (alpha = error scaling parameter).
%   cp.Psub(y+1, z+1) = p(y|z)  (each column sums to 1).
%   Drift bounds follow Sec. IV-A: Dmin = -ceil(2 n pd), Dmax = ceil(2 n pi).
a = alpha;
switch ch
    case 1   % DNA channel, Nanopore asymmetric substitutions (Table 2)
        pins = 0.017*a;  pdel = 0.040*a;  ps = 0.043*a;
        T = [1-ps,  ps/2,  0,     ps/2;        % rows z, columns y  (Table 2)
             ps/5,  1-ps,  ps/5,  3*ps/5;
             0,     ps/2,  1-ps,  ps/2;
             ps/5,  3*ps/5, ps/5, 1-ps];
        Psub = T.';
    case 2, pins = 0;       pdel = 0.1*a;   ps = 0;       Psub = eye(4);   % deletion
    case 3, pins = 0.1*a;   pdel = 0;       ps = 0;       Psub = eye(4);   % insertion
    case 4, pins = 0;       pdel = 0;       ps = 0.1*a;   Psub = [];       % substitution
    case 5, pins = 0.1*a/3; pdel = 0.1*a/3; ps = 0.1*a/3; Psub = [];       % IDS
    otherwise, error('channel must be 1..5');
end
if isempty(Psub)      % 4-ary symmetric (Eq. 44)
    Psub = (1 - ps) * eye(4) + (ps/3) * (ones(4) - eye(4));
end
cp.pins = pins;  cp.pdel = pdel;  cp.ps = ps;  cp.Psub = Psub;
cp.Dmin = -ceil(2*n*pdel);  cp.Dmax = ceil(2*n*pins);
cp.Dmin = cp.Dmin + 0;      % avoid -0
end
