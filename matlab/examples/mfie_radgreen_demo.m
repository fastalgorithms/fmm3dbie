% MFIE for a perfect conductor with a radgreen kernel:
%
%    J/2 - n x curl S_k[J] = n x H_inc
%
% with n x curl S_k built from helm3d.green by kernel3d(g, 'mfie') and its
% fmm from hfmm3d by the 'fmm' option. The scattered field curl S_k[J] is
% compared with the NRCCIE solution from em3d.pec.

zk  = 1.5;
eps = 1e-7;
S = geometries.ellipsoid([1, 0.8, 1.2], [2, 2, 2], [0;0;0], 8);
N = S.npts;

% the fmm of the green, from which the fmms of the kernels are built (one
% hfmm3d call with a density for each component of J)
g = @(s,t) helm3d.green(zk, s, t);
fmm = @(eps, s, pg, t, pgt) hfmm3d(eps, zk, s, pg, t, pgt);
Kmfie = kernel3d(g, 'mfie', 'zk', zk, 'fmm', fmm);   % n x curl S_k, opdims [3 3]
Kcurl = kernel3d(g, 'curls', 'zk', zk, 'fmm', fmm);  % curl S_k,     opdims [3 3]

% incident plane wave
dir = [pi/3, 0];  pol = [1; 0.2];
[Einc, Hinc] = em3d.planewave(zk, dir, pol, S);
rhs = cross(S.n, Hinc);

% solve; densities are Cartesian, interleaved [Jx; Jy; Jz] per point
tic;
J = surfermatsolve(S, kernel3d.eye(0.5*eye(3)) - Kmfie, ...
                                      rhs(:), eps, 1e-10);
toc;
J = reshape(J, 3, N);
fprintf('max |n.J| / max |J|: %.1e\n', max(abs(sum(S.n.*J))) / max(abs(J(:))));

% scattered field at exterior targets, smooth quadrature
r = randn(3, 50);  targ = struct('r', 3*r./vecnorm(r));
H = reshape(surferkerneval(S, Kcurl, J(:), targ, eps), 3, []);

% compare with NRCCIE
tic;
dens = em3d.pec.solver(S, Einc, Hinc, eps, zk, 1.0, struct('eps_gmres', 1e-10));
[~, He] = em3d.pec.eval(S, dens, targ, eps, zk, 1.0);
fprintf('NRCCIE solve: %.1fs\n', toc);
fprintf('MFIE vs NRCCIE, scattered H at r = 3: %.2e\n', ...
        norm(H - He, 'fro')/norm(He, 'fro'));

figure(1); clf
plot(S, real(J(1,:))); colorbar; title('Re J_x')

% total H on the xy plane, at exterior points
npl = 200;
xx = linspace(-2.5, 2.5, npl);
[X, Y, Z] = meshgrid(xx, xx, 0);
pl = struct('r', [X(:).'; Y(:).'; Z(:).']);
ext = ~surferinterior(S, pl);
tic;
Ht = NaN(3, npl^2);
Ht(:,ext) = reshape(surferkerneval(S, Kcurl, J(:), struct('r', pl.r(:,ext)), eps), 3, []);
[~, Hi] = em3d.planewave(zk, dir, pol, struct('r', pl.r(:,ext)));
Ht(:,ext) = Ht(:,ext) + Hi;
fprintf('H on %d exterior plane points: %.1fs\n', nnz(ext), toc);

figure(2); clf
pcolor(X, Y, reshape(real(Ht(3,:)), npl, npl)); shading interp; hold on
plot(S, real(J(1,:)))
colorbar; axis equal tight
title('$\Re H_z$', 'Interpreter', 'latex')
set(gca, 'fontsize', 14)

figure(3); clf
pcolor(X, Y, reshape(vecnorm(Ht, 2, 1), npl, npl)); shading interp; hold on
plot(S, abs(J(1,:)))
colorbar; axis equal tight
title('$|H|$', 'Interpreter', 'latex')
set(gca, 'fontsize', 14)
