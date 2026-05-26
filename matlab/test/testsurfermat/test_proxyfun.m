%TEST_PROXYFUN  Test byindex.proxyfun by checking it gives a low-rank
% factorization of a well-separated source-to-target block via an ID.

run ../../startup.m

rng(42);
zk = 1.3 + 0.1i;
eps = 1e-10;
tol = 1e-9;

fprintf('=== test_proxyfun ===\n\n');

% Big sphere surfer, radius 1, centered at origin
S = geometries.sphere(2, 10, [0;0;0], 6, 1);
kern = kernel3d('h', 's', zk);
kern = kernel3d('l', 's');

% Source ids: points inside ball of radius r_src
r_src  = 0.3;
r_sep  = 1;   % proxy sphere radius (unit-scale, will be rescaled by l)
r_targ = 2;   % targets outside this ball

ctr = [2; 0; 0];
l   = 2*r_src;    % box half-length — proxy points live at r_sep/l in unit scale

srcids  = find(vecnorm(S.r-ctr) < r_src).';
targids = find(vecnorm(S.r-ctr) > r_targ).';

assert(~isempty(srcids),  'No source points found in inner ball');
assert(~isempty(targids), 'No target points found outside outer ball');

% True src->targ block
src2targ = byindex.kernbyindex(targids, srcids, S, kern, eps);

% Proxy surface: random points on a sphere of radius r_sep (unit-scale: r_sep/l)
nproxy = 400;
pr = randn(3, nproxy);
pr = (r_sep / l) .* pr ./ vecnorm(pr);   % unit-scale proxy points
pn = pr ./ vecnorm(pr);                   % outward normals
pw = ones(nproxy, 1)/nproxy*4*pi;                     % equal weights

% pin: logical mask — is a (unit-scale recentered) point inside the proxy sphere?
pin = @(dr) vecnorm(dr) < (r_sep / l);

% slf = column indices for srcids (scalar kernel, opdims=[1,1])
slf = srcids;
% nbr = row indices for targids — proxyfun will filter to those outside proxy
nbr = targids;

ifaddtrans = false;

[Kpxy, nbr_filt] = byindex.proxyfun(slf, nbr, l, ctr, S, kern, [], ...
    pr, pn, pw, pin, ifaddtrans);

% % Interpolative decomposition of src2proxy via pivoted QR
% [~, R, jpvt] = qr(Kpxy.', 'vector');
% d = abs(diag(R));
% k = max(sum(d > tol * d(1)), 1);
% skel = jpvt(1:k);
% 
% 
% % Check that src2targ columns live in span of skeleton columns
% src2tarag_skel = src2targ(:, skel);
% coeffs = src2targ_skel \ src2targ;
% src2targ_lr = src2targ_skel * coeffs;

[skel,rem,T] = id(Kpxy,tol);
errpxy = norm(Kpxy(:,rem) - Kpxy(:,skel)*T,'fro')/norm(Kpxy,'fro')
err = norm(src2targ(:,rem) - src2targ(:,skel)*T,'fro')/norm(src2targ,'fro')

assert(err < tol, 'FAIL: low-rank error %.2e exceeds tol %.2e', err, tol);
% err = norm(src2targ - src2targ_lr, 'fro') / norm(src2targ, 'fro');
% fprintf('Rank of proxy approximation: %d  (out of %d sources, %d targets)\n', ...
%     k, length(srcids), length(targids));
% fprintf('Low-rank factorization error: %.2e\n', err);
% assert(err < tol, 'FAIL: low-rank error %.2e exceeds tol %.2e', err, tol);
% fprintf('PASS\n');

