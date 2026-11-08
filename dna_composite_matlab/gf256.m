function G = gf256()
% GF256  Log/antilog and multiplication tables for GF(2^8), primitive poly x^8+x^4+x^3+x^2+1.
%   Symbols are integers 0..255 (bit-vector representation, addition = XOR).
%   G.EXP(k+1) = alpha^k (k = 0..511),  G.LOG(a+1) = log_alpha(a) (a >= 1)
%   G.MUL(a+1,b+1) = a*b,               G.INV(a+1)  = 1/a   (a >= 1)
persistent GG
if ~isempty(GG), G = GG; return; end
prim = 285;
EXP = zeros(1, 512);  LOG = zeros(1, 256);
x = 1;
for k = 0:254
    EXP(k+1) = x;  LOG(x+1) = k;
    x = bitshift(x, 1);
    if x >= 256, x = bitxor(x, prim); end
end
for k = 255:511, EXP(k+1) = EXP(k-255+1); end
[A, B] = ndgrid(0:255, 0:255);
MUL = zeros(256);
nz  = (A > 0) & (B > 0);
MUL(nz) = EXP(LOG(A(nz)+1) + LOG(B(nz)+1) + 1);
INV = zeros(1, 256);
INV(2:256) = EXP(255 - LOG(2:256) + 1);
GG = struct('EXP', EXP, 'LOG', LOG, 'MUL', MUL, 'INV', INV);
G = GG;
end
