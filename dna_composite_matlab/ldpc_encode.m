function c = ldpc_encode(code, info)
% LDPC_ENCODE  Systematic encoding. info: 1 x K symbols (0..255). Returns 1 x N codeword.
G = gf256();  %#ok<NASGU>
c = zeros(1, code.N);
c(code.nonpiv) = info;
prodm = gf_helpers('mul', code.P, info(:).');            % rank x K
par = gf_helpers('xorcols', prodm);                      % rank x 1
c(code.piv) = par.';
end
