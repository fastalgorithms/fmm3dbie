%TEST_PROXYFUN  Test byindex.proxyfun by checking it gives a low-rank
% factorization of a well-separated source-to-target block via an ID.

run ../../startup.m

rng(42);
zk = 1.3 + 0.1i;
eps = 1e-10;
tol = 1e-9;

fprintf('=== test_proxyfun ===\n\n');

% Single sphere surfer. Sources are points near ctr on one side of the sphere,
% targets are points far from ctr on the other side.
S    = geometries.sphere(2, 10, [0;0;0], 6, 1);
kern = kernel3d('l', 's');

ctr    = [2; 0; 0];   % near one pole of the sphere
r_src  = 0.5;         % sources within this radius of ctr
r_targ = 2.0;         % targets beyond this radius from ctr
l      = 2*r_src;     % box half-length = 0.6
r_sep  = 2.0;         % proxy sphere unit-scale radius
                       % physical proxy radius = r_sep*l = 1.2
                       % source-to-proxy gap = 1.2 - 0.3 = 0.9  (ratio 4:1)
                       % target-to-proxy gap = 2.0 - 1.2 = 0.8

srcids  = find(vecnorm(S.r - ctr) < r_src).';
targids = find(vecnorm(S.r - ctr) > r_targ).';

assert(~isempty(srcids),  'No source points found in inner ball');
assert(~isempty(targids), 'No target points found outside outer ball');

% True src->targ block
src2targ = byindex.kernbyindex(targids, srcids, S, kern, eps);

% Proxy surface: GL in cos(theta) x equispaced in phi on sphere of radius r_sep
ntheta = 20;
nphi   = 20;
[ct, wt] = polytens.lege.pts(ntheta);   % ntheta GL nodes on [-1,1] for cos(theta)
theta = acos(ct(:));                     % (ntheta,1)
phi   = (0:nphi-1).' * (2*pi/nphi);     % (nphi,1), equispaced
wphi  = 2*pi/nphi;                       % scalar weight per phi point

% Tensor product: (ntheta*nphi) proxy points
[TH, PH] = meshgrid(theta, phi);         % (nphi,ntheta) each
TH = TH(:); PH = PH(:);

% Weights: w_theta * w_phi, with sin(theta) already in GL weight (integrating
% over cos(theta) absorbs the Jacobian: int f dOmega = int_{-1}^{1} int_0^{2pi} f dphi d(cos theta))
[~, WT] = meshgrid(wt(:), ones(nphi,1));
WT = WT(:);                              % (ntheta*nphi,1), GL weights
pw = WT * wphi;                          % (ntheta*nphi,1), solid-angle weights

% pr in unit scale: proxyfun does pxy = pr*l + ctr internally,
% so pass points at radius r_sep (no /l factor here)
pr = r_sep * [sin(TH).*cos(PH), sin(TH).*sin(PH), cos(TH)].';  % (3,nproxy)
pn = pr ./ vecnorm(pr);                   % outward normals

% pin: unit-scale test — a point dr (already divided by l) is inside
% the proxy sphere if its unit-scale radius < r_sep
pin = @(dr) vecnorm(dr) < r_sep;

slf = srcids;
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