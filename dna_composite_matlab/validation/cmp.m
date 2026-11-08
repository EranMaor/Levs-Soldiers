addpath('..');
fid = fopen('cases.txt'); worst = 0;
while true
    hdr = fgetl(fid); if ~ischar(hdr), break; end
    v = sscanf(hdr(6:end), '%f').';  ci=v(1); n=v(2); r=v(3); m=v(4); pi_=v(5); pd_=v(6); Dmin=v(7); Dmax=v(8);
    fgetl(fid); TH = zeros(m,4);  for k=1:m, TH(k,:) = sscanf(fgetl(fid), '%f').'; end
    fgetl(fid); Ps = zeros(4,4);  for k=1:4, Ps(k,:) = sscanf(fgetl(fid), '%f').'; end
    fgetl(fid); Y = cell(1,r);    for t=1:r, Y{t} = sscanf(fgetl(fid), '%f').'; end
    fgetl(fid); P = zeros(n,m);   for i=1:n, P(i,:) = sscanf(fgetl(fid), '%f').'; end
    post = app_calc(Y, TH, [], pi_, pd_, Dmin, Dmax, Ps, n);
    d = max(abs(post(:) - P(:)));  worst = max(worst, d);
    [~,a1] = max(post,[],2); [~,a2] = max(P,[],2);
    fprintf('case %d (n=%d r=%d m=%d ND=%d): max|MATLAB - python| = %.3g   argmax agree: %d\n', ci, n, r, m, Dmax-Dmin+1, d, isequal(a1,a2));
end
fclose(fid);
fprintf('worst difference over all cases: %.3g\n', worst);
