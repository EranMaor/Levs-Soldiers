function ok = ldpc_check(code, c)
% LDPC_CHECK  True if all parity checks are satisfied by symbol vector c (0..255).
c = c(:);
pr = gf_helpers('mul', code.coef, c(code.edgeVar));       % E x 1
pv = [pr; 0];
s  = gf_helpers('xorcols', pv(code.chkEdges));            % M x 1
ok = all(s == 0);
end
