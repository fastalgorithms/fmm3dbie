% Test byindex.kernbyindexeval and byindex.proxyfuneval.

run ../../startup.m

rng(42);

%% Now run the tests

test_single_surfer_laplace();
test_two_surfer_helmholtz();
test_proxyfuneval_lowrank();


function test_single_surfer_laplace()
% Single surfer, Laplace SLP: kernbyindexeval matches surferkernevalmat.

eps = 1e-10;  tol = 1e-9;

S    = geometries.sphere(2, 1, [0;0;0], 4, 1);
kern = kernel3d('l', 's');

rng(1);
ntarg = 50;
targs.r = 3 * (2*rand(3,ntarg) - 1);
i = sort(randperm(ntarg,  20)).';
j = sort(randperm(S.npts, 30)).';

check_eval('single-surfer Laplace', i, j, S, kern, targs, eps, tol, 0);

end


function test_two_surfer_helmholtz()
% Two surfers, Helmholtz: kernbyindexeval matches surferkernevalmat.

eps = 1e-10;  tol = 1e-9;

zk = 1.3 + 0.1i;
S1 = geometries.sphere(2, 1, [0;0;0], 3, 1);
S2 = geometries.sphere(2, 1, [6;0;0], 3, 1);
kern = [kernel3d('h', 's', zk), kernel3d('h', 'd', zk)];

rng(2);
ntarg = 40;
targs.r = [3*(2*rand(3,ntarg/2)-1), [6;0;0] + 3*(2*rand(3,ntarg/2)-1)];
i = sort(randperm(ntarg,          15)).';
j = sort(randperm(S1.npts+S2.npts, 25)).';

check_eval('two-surfer Helmholtz', i, j, [S1,S2], kern, targs, eps, tol, 1);

end


function test_proxyfuneval_lowrank()
% proxyfuneval gives a low-rank factorization of an off-surface src->targ block.

eps = 1e-10;  tol = 1e-9;

S    = geometries.sphere(2, 1, [0;0;0], 6, 1);
kern = kernel3d('l', 's');
ctr  = [2; 0; 0];  r_src = 0.5;  r_targ = 2.0;  l = 2*r_src;

srcids = find(vecnorm(S.r - ctr) < r_src).';

rng(3);
ntarg = 200;
th = acos(2*rand(1,ntarg)-1);  ph = 2*pi*rand(1,ntarg);
targs.r = ctr + (r_targ + rand(1,ntarg)) .* [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];

novers    = {NaN*zeros(S.npatches,1)};
src2targ  = byindex.kernbyindexeval((1:ntarg).', srcids, S, kern, targs, eps, novers);

ntheta = 20;  nphi = 40;  r_sep = 2.0;
[ct, wt] = polytens.lege.pts(ntheta);
theta = acos(ct(:));  phi = (0:nphi-1).' * (2*pi/nphi);
[TH, PH] = meshgrid(theta, phi);  TH = TH(:);  PH = PH(:);
[~, WT]  = meshgrid(wt(:), ones(nphi,1));  WT = WT(:);
pw = WT * (2*pi/nphi);
pr = r_sep * [sin(TH).*cos(PH), sin(TH).*sin(PH), cos(TH)].';
pn = pr ./ vecnorm(pr);
pin = @(dr) vecnorm(dr) < r_sep;

[Kpxy, ~] = byindex.proxyfuneval(srcids, (1:ntarg).', l, ctr, S, kern, ...
    targs.r, pr, pn, pw, pin);

[~, R, jpvt] = qr(Kpxy, 'vector');
d = abs(diag(R));  k = max(sum(d > tol * d(1)), 1);
skel = jpvt(1:k);
A    = src2targ(:, skel);
err  = norm(A*(A\src2targ) - src2targ, 'fro') / norm(src2targ, 'fro');
assert(err < tol, 'proxyfuneval low-rank: err %.2e > tol %.2e', err, tol);

end


function check_eval(label, i, j, surfers, kern, targs, eps, tol, ifoversamp)
% Helper: compare kernbyindexeval against full surferkernevalmat.
opts = []; opts.ifoversamp = ifoversamp;
opts_ns = opts; opts_ns.nonsmoothonly = true;
[Qns, novers] = surferkernevalmat(surfers, kern, targs, eps, opts_ns);
Afull = surferkernevalmat(surfers, kern, targs, eps, opts);
Atest = byindex.kernbyindexeval(i, j, surfers, kern, targs, eps, novers, Qns);
err   = norm(Atest - Afull(i,j), 'fro') / norm(Afull(i,j), 'fro');
assert(err < tol, '%s (ifoversamp=%d): err %.2e > tol %.2e', label, ifoversamp, err, tol);
end
