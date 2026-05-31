%TEST_PROXYFUN  Test byindex.proxyfun.
% Scenario 1: single surfer, ifaddtrans=false. Proxy-based ID should give
%   low-rank approximation of the true src->targ block.
% Scenario 2: single surfer, kern=SLP forward, kern2=DLP transpose, ifaddtrans=true.
%   The stacked [forward; trans] block finds a skeleton that compresses both.
%   Negative test: ifaddtrans=false fails to compress the DLP targ->src block.
% Scenario 3: two surfers, transmission kernel (mixed opdims), ifaddtrans=false.

run ../../startup.m

rng(42);
eps = 1e-10;
tol = 1e-9;

fprintf('=== test_proxyfun ===\n\n');

novers1 = {@(S) NaN*zeros(S.npatches,1)};

% proxy surface (unit-scale sphere, GL x equispaced)
ntheta = 10; nphi = 10; r_sep = 2.0;
[ct, wt] = polytens.lege.pts(ntheta);
theta = acos(ct(:)); phi = (0:nphi-1).' * (2*pi/nphi);
[TH, PH] = meshgrid(theta, phi); TH = TH(:); PH = PH(:);
[~, WT] = meshgrid(wt(:), ones(nphi,1)); WT = WT(:);
pw  = WT * (2*pi/nphi);
pr  = r_sep * [sin(TH).*cos(PH), sin(TH).*sin(PH), cos(TH)].';
pn  = pr ./ vecnorm(pr);
pin = @(dr) vecnorm(dr) < r_sep;

function check_lr(label, src2targ, Kpxy, tol)
    [~, R, jpvt] = qr(Kpxy, 'vector');
    d = abs(diag(R));
    k = max(sum(d > tol * d(1)), 1);
    skel = jpvt(1:k);
    A = src2targ(:, skel);
    src2targ_lr = A * (A \ src2targ);
    err = norm(src2targ_lr - src2targ, 'fro') / norm(src2targ, 'fro');
    fprintf('%s: rank = %d,  err = %.2e\n', label, k, err);
    assert(err < tol, 'FAIL: %s err %.2e > tol %.2e', label, err, tol);
end

%% Scenario 1: single surfer, ifaddtrans=false
fprintf('Scenario 1: single surfer, ifaddtrans=false ...\n');

S    = geometries.sphere(2, 10, [0;0;0], 6, 1);
kern = kernel3d('l', 's');
ctr  = [2; 0; 0]; r_src = 0.5; r_targ = 2.0; l = 2*r_src;

srcids  = find(vecnorm(S.r - ctr) < r_src).';
targids = find(vecnorm(S.r - ctr) > r_targ).';
novers  = {NaN*zeros(S.npatches,1)};
src2targ = byindex.kernbyindex(targids, srcids, S, kern, eps, novers);

[Kpxy, ~] = byindex.proxyfun(srcids, targids, l, ctr, S, kern, [], ...
    pr, pn, pw, pin, false);
check_lr('Scenario 1', src2targ, Kpxy, tol);

%% Scenario 2: single surfer, kern=SLP forward, kern2=DLP transpose, ifaddtrans=true
fprintf('Scenario 2: single surfer, kern=SLP, kern2=DLP, ifaddtrans=true ...\n');

kern2 = kernel3d('l', 'd');
targ2src_dlp = byindex.kernbyindex(srcids, targids, S, kern2, eps, novers);

[Kpxy2, ~] = byindex.proxyfun(srcids, targids, l, ctr, S, kern, kern2, ...
    pr, pn, pw, pin, true);

[~, R, jpvt] = qr(Kpxy2, 'vector');
d = abs(diag(R)); k = max(sum(d > tol * d(1)), 1);
skel2 = jpvt(1:k);

A = src2targ(:, skel2);
err_fwd = norm(A*(A\src2targ) - src2targ, 'fro') / norm(src2targ, 'fro');
fprintf('Scenario 2 forward: rank = %d,  err = %.2e\n', k, err_fwd);
assert(err_fwd < tol, 'FAIL: Scenario 2 forward err %.2e > tol %.2e', err_fwd, tol);

B = targ2src_dlp(skel2, :);
targ2src_lr = B.' * (B.' \ targ2src_dlp.');
err_trans = norm(targ2src_lr.' - targ2src_dlp, 'fro') / norm(targ2src_dlp, 'fro');
fprintf('Scenario 2 trans:   rank = %d,  err = %.2e\n', k, err_trans);
assert(err_trans < tol, 'FAIL: Scenario 2 trans err %.2e > tol %.2e', err_trans, tol);

% Structural test: ifaddtrans=false gives npxy rows, ifaddtrans=true gives 2*npxy.
[Kpxy2_notrans, ~] = byindex.proxyfun(srcids, targids, l, ctr, S, kern, kern2, ...
    pr, pn, pw, pin, false);
npxy = size(pr, 2);
assert(size(Kpxy2_notrans, 1) == npxy,   'FAIL: ifaddtrans=false should have %d rows, got %d', npxy,   size(Kpxy2_notrans,1));
assert(size(Kpxy2,         1) == 2*npxy, 'FAIL: ifaddtrans=true  should have %d rows, got %d', 2*npxy, size(Kpxy2,        1));
fprintf('Scenario 2 structural: ifaddtrans=false rows=%d, ifaddtrans=true rows=%d  [correct]\n', size(Kpxy2_notrans,1), size(Kpxy2,1));

%% Scenario 3: two surfers, transmission kernel (mixed opdims), ifaddtrans=false
% S1 carries the [2,2]-opdim transmission representation
% S2 carries a Helmholtz single layer
fprintf('Scenario 3: two surfers, transmission kernel, ifaddtrans=false ...\n');

zk3 = 1.3 + 0.1i;
S1 = geometries.sphere(2, 8, [0;0;0], 4, 1);
S2 = geometries.sphere(2, 8, [8;0;0], 4, 1);
kern3 = [kernel3d('h', 'trans_rep', zk3), kernel3d('h', 's',      zk3); ...
         kernel3d('h', 'trans',     zk3), kernel3d('h', 's2trans', zk3)];

ctr3 = [2; 0; 0]; r_src3 = 0.5; r_targ3 = 2.0; l3 = 2*r_src3;

% Column layout: [1..2*S1.npts | 2*S1.npts+1..2*S1.npts+S2.npts]  (opdim 2,1)
% Row layout:    [1..S1.npts   | S1.npts+1..S1.npts+2*S2.npts]    (opdim 1,2)
s1_pts = find(vecnorm(S1.r - ctr3) < r_src3).';
s2_pts = find(vecnorm(S2.r - ctr3) > r_targ3).';
srcids3  = reshape([2*s1_pts-1; 2*s1_pts], 1, []);           % opdim-2 cols on S1
targids3 = S1.npts + reshape([2*s2_pts-1; 2*s2_pts], 1, []); % opdim-2 rows on S2

surfers3 = [S1, S2];
novers3  = {NaN*zeros(S1.npatches,1), NaN*zeros(S1.npatches,1); ...
            NaN*zeros(S2.npatches,1), NaN*zeros(S2.npatches,1)};
src2targ3 = byindex.kernbyindex(targids3, srcids3, surfers3, kern3, eps, novers3);

[Kpxy3, ~] = byindex.proxyfun(srcids3, targids3, l3, ctr3, surfers3, kern3, [], ...
    pr, pn, pw, pin, false);
check_lr('Scenario 3', src2targ3, Kpxy3, tol);

fprintf('\nAll tests passed.\n');
