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

% --- scatter plot: sources, targets, proxy sphere ---
figure(1); clf;
scatter3(S.r(1,srcids),  S.r(2,srcids),  S.r(3,srcids),  20, 'r', 'filled'); hold on;
scatter3(S.r(1,targids), S.r(2,targids), S.r(3,targids), 10, 'b', 'filled');
scatter3(ctr(1), ctr(2), ctr(3), 100, 'k', 'p', 'filled');  % box center
legend('sources','targets','ctr');
axis equal; grid on;
title(sprintf('r\\_src=%.2f  r\\_sep*l=%.2f  r\\_targ=%.2f', r_src, r_sep*l, r_targ));
drawnow;

% Proxy surface: GL in cos(theta) x equispaced in phi on sphere of radius r_sep
ntheta = 20;
nphi   = 40;
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

% Add physical proxy points to scatter plot
pxy_phys = pr * l + ctr;
figure(1);
scatter3(pxy_phys(1,:), pxy_phys(2,:), pxy_phys(3,:), 5, 'g', 'filled');
legend('sources','targets','ctr','proxy');
drawnow;

% pin: unit-scale test — a point dr (already divided by l) is inside
% the proxy sphere if its unit-scale radius < r_sep
pin = @(dr) vecnorm(dr) < r_sep;

slf = srcids;
nbr = targids;

ifaddtrans = false;

% [Kpxy, nbr_filt] = byindex.proxyfun(slf, nbr, l, ctr, S, kern, [], ...
%     pr, pn, pw, pin, ifaddtrans);

srcs = []; srcs.r = S.r(:,slf);
Kpxy = kern.eval(srcs,struct('r',pr+ctr)).*S.wts(slf).';


% Kpxy = sqrt(pw).*Kpxy./sqrt(S.wts(slf)).';
% src2targ = sqrt(S.wts(targids)).*src2targ./sqrt(S.wts(slf)).';
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

% Color targets by per-row residual norm
resid   = src2targ(:,rem) - src2targ(:,skel)*T;  % ntarg x nrem
row_err = vecnorm(resid, 2, 2);                   % ntarg x 1

figure(2); clf;
scatter3(S.r(1,targids), S.r(2,targids), S.r(3,targids), 20, row_err, 'filled'); hold on;
scatter3(S.r(1,srcids),  S.r(2,srcids),  S.r(3,srcids),  30, 'r', 'filled');
scatter3(pxy_phys(1,:),  pxy_phys(2,:),  pxy_phys(3,:),   5, 'g', 'filled');
scatter3(ctr(1), ctr(2), ctr(3), 100, 'k', 'p', 'filled');
colorbar; colormap(hot);
legend('targets (colored by row residual)','sources','proxy','ctr');
axis equal; grid on;
title('per-target residual norm  ||resid(i,:)||');
drawnow;

assert(err < tol, 'FAIL: low-rank error %.2e exceeds tol %.2e', err, tol);