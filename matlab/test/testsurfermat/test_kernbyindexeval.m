%TEST_KERNBYINDEXEVAL  Verify byindex.kernbyindexeval and byindex.proxyfuneval.
%
% Scenario 1: single surfer, Laplace SLP. kernbyindexeval matches
%   surferkernevalmat, tested with ifoversamp=1 and ifoversamp=0.
% Scenario 2: two surfers, Helmholtz. Same, with ifoversamp=1 and 0.
% Scenario 3: proxyfuneval gives a low-rank factorization of an off-surface
%   src->targ block.

run ../../startup.m

rng(42);
eps = 1e-10;
tol = 1e-9;

fprintf('=== test_kernbyindexeval ===\n\n');

function check_eval(label, i, j, surfers, kern, targs, eps, tol, ifoversamp)
    opts = []; opts.ifoversamp = ifoversamp;
    opts_ns = opts; opts_ns.nonsmoothonly = true;
    [Qns, novers] = surferkernevalmat(surfers, kern, targs, eps, opts_ns);
    Afull = surferkernevalmat(surfers, kern, targs, eps, opts);
    Aref  = Afull(i, j);
    Atest = byindex.kernbyindexeval(i, j, surfers, kern, targs, eps, novers, Qns);
    err = norm(Atest - Aref, 'fro') / norm(Aref, 'fro');
    fprintf('%s (ifoversamp=%d): rel err = %.2e\n', label, ifoversamp, err);
    assert(err < tol, 'FAIL: %s err %.2e > tol %.2e', label, err, tol);
end

%% Scenario 1: single surfer, Laplace SLP
fprintf('Scenario 1: single surfer, Laplace SLP ...\n');

S    = geometries.sphere(2, 2, [0;0;0], 4, 1);
kern = kernel3d('l', 's');

rng(1);
ntarg = 50;
targs.r = 3 * (2*rand(3,ntarg) - 1);

i1 = sort(randperm(ntarg,   20)).';
j1 = sort(randperm(S.npts,  30)).';

check_eval('Scenario 1', i1, j1, S, kern, targs, eps, tol, 0);

%% Scenario 2: two surfers, Helmholtz
fprintf('Scenario 2: two surfers, Helmholtz ...\n');

zk = 1.3 + 0.1i;
S1 = geometries.sphere(2, 2, [0;0;0], 3, 1);
S2 = geometries.sphere(2, 2, [6;0;0], 3, 1);
kern2 = [kernel3d('h', 's', zk), kernel3d('h', 'd', zk)];

rng(2);
ntarg2 = 40;
targs2.r = [3*(2*rand(3,ntarg2/2)-1), [6;0;0] + 3*(2*rand(3,ntarg2/2)-1)];

i2 = sort(randperm(ntarg2,         15)).';
j2 = sort(randperm(S1.npts+S2.npts, 25)).';

check_eval('Scenario 2', i2, j2, [S1,S2], kern2, targs2, eps, tol, 1);

%% Scenario 3: proxyfuneval low-rank test
fprintf('Scenario 3: proxyfuneval low-rank factorization ...\n');

S3    = geometries.sphere(2, 2, [0;0;0], 6, 1);
kern3 = kernel3d('l', 's');
ctr3  = [2; 0; 0]; r_src = 0.5; r_targ = 2.0; l3 = 2*r_src;

srcids = find(vecnorm(S3.r - ctr3) < r_src).';

rng(3);
ntarg3 = 200;
th = acos(2*rand(1,ntarg3)-1); ph = 2*pi*rand(1,ntarg3);
targs3.r = ctr3 + (r_targ + rand(1,ntarg3)) .* [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];

novers3 = {NaN*zeros(S3.npatches,1)};
src2targ3 = byindex.kernbyindexeval((1:ntarg3).', srcids, S3, kern3, targs3, eps, novers3);

ntheta = 20; nphi = 40; r_sep = 2.0;
[ct, wt] = polytens.lege.pts(ntheta);
theta = acos(ct(:)); phi = (0:nphi-1).' * (2*pi/nphi);
[TH, PH] = meshgrid(theta, phi); TH = TH(:); PH = PH(:);
[~, WT] = meshgrid(wt(:), ones(nphi,1)); WT = WT(:);
pw3 = WT * (2*pi/nphi);
pr3 = r_sep * [sin(TH).*cos(PH), sin(TH).*sin(PH), cos(TH)].';
pn3 = pr3 ./ vecnorm(pr3);
pin3 = @(dr) vecnorm(dr) < r_sep;

[Kpxy3, ~] = byindex.proxyfuneval(srcids, (1:ntarg3).', l3, ctr3, S3, kern3, ...
    targs3.r, pr3, pn3, pw3, pin3);

[~, R, jpvt] = qr(Kpxy3, 'vector');
d = abs(diag(R)); k = max(sum(d > tol * d(1)), 1);
skel = jpvt(1:k);
A = src2targ3(:, skel);
err3 = norm(A*(A\src2targ3) - src2targ3, 'fro') / norm(src2targ3, 'fro');
fprintf('Scenario 3: rank = %d,  err = %.2e\n', k, err3);
assert(err3 < tol, 'FAIL: Scenario 3 err %.2e > tol %.2e', err3, tol);

fprintf('\nAll tests passed.\n');
