% Test byindex.proxyfun.

run ../../startup.m

rng(42);

% Shared proxy surface used by all scenarios
ntheta = 10;  nphi = 10;  r_sep = 2.0;
[ct, wt] = polytens.lege.pts(ntheta);
theta = acos(ct(:));  phi = (0:nphi-1).' * (2*pi/nphi);
[TH, PH] = meshgrid(theta, phi);  TH = TH(:);  PH = PH(:);
[~, WT]  = meshgrid(wt(:), ones(nphi,1));  WT = WT(:);
pw  = WT * (2*pi/nphi);
pr  = r_sep * [sin(TH).*cos(PH), sin(TH).*sin(PH), cos(TH)].';
pn  = pr ./ vecnorm(pr);
pin = @(dr) vecnorm(dr) < r_sep;

%% Now run the tests

test_single_surfer(pr, pn, pw, pin);
test_with_transpose(pr, pn, pw, pin);
test_transmission_kernel(pr, pn, pw, pin);


function test_single_surfer(pr, pn, pw, pin)
% Single surfer SLP: proxy-based ID gives low-rank approximation of src->targ block.

eps = 1e-10;  tol = 1e-9;

S    = geometries.sphere(2, 10, [0;0;0], 6, 1);
kern = kernel3d('l', 's');
ctr  = [2; 0; 0];  r_src = 0.5;  r_targ = 2.0;  l = 2*r_src;

srcids  = find(vecnorm(S.r - ctr) < r_src).';
targids = find(vecnorm(S.r - ctr) > r_targ).';
novers  = {NaN*zeros(S.npatches,1)};
src2targ = byindex.kernbyindex(targids, srcids, S, kern, eps, novers);

[Kpxy, ~] = byindex.proxyfun(srcids, targids, l, ctr, S, kern, [], pr, pn, pw, pin, false);
check_lr('single surfer', src2targ, Kpxy, tol);

end


function test_with_transpose(pr, pn, pw, pin)
% SLP forward + DLP transpose, ifaddtrans=true: skeleton compresses both blocks.
% Also checks row counts: ifaddtrans=false -> npxy rows, true -> 2*npxy.

eps = 1e-10;  tol = 1e-9;

S     = geometries.sphere(2, 10, [0;0;0], 6, 1);
kern  = kernel3d('l', 's');
kern2 = kernel3d('l', 'd');
ctr   = [2; 0; 0];  r_src = 0.5;  r_targ = 2.0;  l = 2*r_src;

srcids  = find(vecnorm(S.r - ctr) < r_src).';
targids = find(vecnorm(S.r - ctr) > r_targ).';
novers  = {NaN*zeros(S.npatches,1)};

src2targ    = byindex.kernbyindex(targids, srcids, S, kern,  eps, novers);
targ2src_dlp = byindex.kernbyindex(srcids, targids, S, kern2, eps, novers);

[Kpxy2, ~] = byindex.proxyfun(srcids, targids, l, ctr, S, kern, kern2, pr, pn, pw, pin, true);

[~, R, jpvt] = qr(Kpxy2, 'vector');
d = abs(diag(R));  k = max(sum(d > tol * d(1)), 1);
skel = jpvt(1:k);

A       = src2targ(:, skel);
err_fwd = norm(A*(A\src2targ) - src2targ, 'fro') / norm(src2targ, 'fro');
assert(err_fwd < tol, 'with-transpose forward: err %.2e > tol %.2e', err_fwd, tol);

B           = targ2src_dlp(skel, :);
targ2src_lr = B.' * (B.' \ targ2src_dlp.');
err_trans   = norm(targ2src_lr.' - targ2src_dlp, 'fro') / norm(targ2src_dlp, 'fro');
assert(err_trans < tol, 'with-transpose trans: err %.2e > tol %.2e', err_trans, tol);

[Kpxy2_notrans, ~] = byindex.proxyfun(srcids, targids, l, ctr, S, kern, kern2, pr, pn, pw, pin, false);
npxy = size(pr, 2);
assert(size(Kpxy2_notrans, 1) == npxy,   'ifaddtrans=false should have %d rows, got %d', npxy,   size(Kpxy2_notrans,1));
assert(size(Kpxy2,         1) == 2*npxy, 'ifaddtrans=true should have %d rows, got %d',  2*npxy, size(Kpxy2,        1));

end


function test_transmission_kernel(pr, pn, pw, pin)
% Two surfers, transmission kernel (mixed opdims): proxy ID gives low-rank block.

eps = 1e-10;  tol = 1e-9;

zk3 = 1.3 + 0.1i;
S1 = geometries.sphere(2, 8, [0;0;0], 4, 1);
S2 = geometries.sphere(2, 8, [8;0;0], 4, 1);
kern = [kernel3d('h', 'trans_rep', zk3), kernel3d('h', 's',      zk3); ...
        kernel3d('h', 'trans',     zk3), kernel3d('h', 's2trans', zk3)];

ctr  = [2; 0; 0];  r_src = 0.5;  r_targ = 2.0;  l = 2*r_src;

s1_pts = find(vecnorm(S1.r - ctr) < r_src).';
s2_pts = find(vecnorm(S2.r - ctr) > r_targ).';
srcids  = reshape([2*s1_pts-1; 2*s1_pts], 1, []);
targids = S1.npts + reshape([2*s2_pts-1; 2*s2_pts], 1, []);

surfers = [S1, S2];
novers  = {NaN*zeros(S1.npatches,1), NaN*zeros(S1.npatches,1); ...
           NaN*zeros(S2.npatches,1), NaN*zeros(S2.npatches,1)};
src2targ = byindex.kernbyindex(targids, srcids, surfers, kern, eps, novers);

[Kpxy, ~] = byindex.proxyfun(srcids, targids, l, ctr, surfers, kern, [], pr, pn, pw, pin, false);
check_lr('transmission kernel', src2targ, Kpxy, tol);

end


function check_lr(label, src2targ, Kpxy, tol)
% Helper: verify Kpxy gives a low-rank factorization of src2targ.
[~, R, jpvt] = qr(Kpxy, 'vector');
d    = abs(diag(R));
k    = max(sum(d > tol * d(1)), 1);
skel = jpvt(1:k);
A    = src2targ(:, skel);
err  = norm(A*(A\src2targ) - src2targ, 'fro') / norm(src2targ, 'fro');
assert(err < tol, '%s: err %.2e > tol %.2e', label, err, tol);
end
