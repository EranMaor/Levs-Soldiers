% SELFTEST  Fast regression checks of every building block (about 1-2 minutes).
set_seed(1);  pass = @(s) fprintf('  [ok] %s\n', s);

% --- GF(256) / LDPC ----------------------------------------------------------
G = gf256();
assert(all(arrayfun(@(a) G.MUL(a+1, G.INV(a+1)+1), 1:255) == 1));
a = 37; b = 201; c = 94;
assert(G.MUL(a+1, bitxor(b, c)+1) == bitxor(G.MUL(a+1, b+1), G.MUL(a+1, c+1)));
pass('GF(256) tables: inverses, distributivity');
code = make_ldpc(256, 128, 2, 1, 'results_selftest');
for k = 1:3
    cw = ldpc_encode(code, randi([0 255], 1, code.K));
    assert(ldpc_check(code, cw));
end
pass('systematic LDPC encoder produces valid codewords');
N = 256;  pe = 0.12;  nerr = 0;
for k = 1:4
    cw = ldpc_encode(code, randi([0 255], 1, code.K));  rx = cw;
    fl = rand(1, N) < pe;  rx(fl) = mod(rx(fl) + randi(255, 1, nnz(fl)), 256);
    gam = (pe/255) * ones(N, 256);  gam(sub2ind([N 256], (1:N)', rx(:)+1)) = 1 - pe;
    chat = ldpc_decode(code, gam, 30);  nerr = nerr + ~isequal(chat(:).', cw);
end
assert(nerr == 0);  pass('FFT-SPA corrects 12% q-ary symmetric symbol errors');

% --- alphabet and channel ---------------------------------------------------------
assert(size(create_alphabet(get_weights('w41')), 1) == 256);
assert(size(create_alphabet(get_weights('w42')), 1) == 100);
assert(size(create_alphabet(get_weights('w21')), 1) == 16);
pass('composite alphabets reproduce Table 3 sizes (256 / 100 / 16)');
cp = channel_params(1, 1.0, 128);
Y  = channel_sim(randi([0 3], 3000, 128), cp.pins, cp.pdel, cp.Dmin, cp.Dmax, cp.Psub);
L  = cellfun(@numel, Y);
assert(abs(mean(L) - (128 + 128*(cp.pins - cp.pdel))) < 0.5);
assert(all(L - 128 >= cp.Dmin) && all(L - 128 <= cp.Dmax));
pass('channel: mean output length and drift bounds');

% --- inner decoder: noiseless limit vs paper (Fig. 8/10: ~6.1 bit at r = 50) ----------
T = create_alphabet(get_weights('w41'));
I = sim_mi_point(T, 1, 0.0, 50, 40, 128);
fprintf('       I(C;U) at alpha=0, r=50: %.2f bit (paper ~6.1)\n', I);
assert(abs(I - 6.1) < 0.3);  pass('inner APP decoder reproduces the paper''s rate loss at alpha = 0');

% --- one concatenated block ---------------------------------------------------------
codeL = make_ldpc(2048, 768, 2, 1, 'results_selftest');
res = sim_concat_block(codeL, T, channel_params(1, 0.3, 128), 50, 128, 20, 20);
assert(~res.err);  pass('concatenated block (inner APP + outer FFT-SPA) decodes at alpha = 0.3');
disp('All self-tests passed.');
