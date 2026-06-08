% stok3d_mixed_bc.m  --  Stokes mixed boundary condition demo
%
% S1: traction-free (u = S[sigma1],            BIE: 1/2 sigma1 + S'[sigma1] + L = rhs1)
% S2: no-slip       (u = aS[sigma2] + D[sigma2], BIE: 1/2 sigma2 + (aS+D)    = rhs2)
%
% Incident flow: uniform u_inf = [1;0;0].

alpha = 1;
eps   = 1e-8;

S1 = geometries.ellipsoid([1,1,1.5], 3*[1,1,1], [0;0;0], 6);
S2 = geometries.ellipsoid([1,1.5,1], 3*[1,1,1], [3;0;0], 6);
srfrs = [S1, S2];

kerns(2,2) = kernel3d();
kerns(1,1) = kernel3d('stokes', 'strac');
kerns(1,2) = kernel3d('stokes', 'ctrac', [alpha;1]);
kerns(2,1) = kernel3d('stokes', 's');
kerns(2,2) = kernel3d('stokes', 'c',     [alpha;1]);

%% S1 surface moments (for L stabilizer on traction-free body)

wi1       = S1.wts;
area1     = sum(wi1);
centroid1 = (S1.r * wi1) / area1;
dx1       = S1.r - centroid1;
rmoi_inv1 = inv(trace(dx1*(dx1.*wi1')')*eye(3) - dx1*(dx1.*wi1')');
icomps1   = [1; S1.npatches+1];

%% Build system

fprintf('Building quadrature corrections...\n'); tic;
opts_build = struct; opts_build.corrections = 1; opts_build.unif_nover = 1;
n1 = 3*S1.npts;  n2 = 3*S2.npts;  npts_tot = n1 + n2;
[Qcor, Sovers] = surfermat(srfrs, kerns, eps, opts_build);
Qcor = Qcor + 0.5*speye(npts_tot);
fprintf('  t = %.2f s\n', toc);

matvec = @(x) stok3d.mixed_matvec(x, srfrs, alpha, eps, Qcor, Sovers, ...
    area1, centroid1, rmoi_inv1);

%% Solve

u_inf = [1;0;0];
rhs = [zeros(n1,1); repmat(-u_inf, S2.npts, 1)];

fprintf('GMRES...\n'); tic;
[sol, flag, relres, iter] = gmres(matvec, rhs, [], eps, 200);
fprintf('  t = %.2f s, flag=%d, relres=%.2e, iter=%d\n', toc, flag, relres, iter(end));

%% Plot

nplot = 40;
xx = linspace(-2, 6, nplot) + 1e-2;
yy = linspace(-3, 3, nplot) + 3e-2;
[XX, YY] = meshgrid(xx, yy);
targs = []; targs.r = [XX(:).'; YY(:).'; zeros(1, numel(XX))];

uscat = surferkerneval(srfrs, kerns(2,:), sol, targs, eps);
utot  = reshape(uscat, 3, []) + u_inf;

u_on_S1 = reshape(surferkerneval(srfrs, kerns(2,:), sol, S1, eps), 3, S1.npts) + u_inf;
u_surf_norm = [vecnorm(u_on_S1), zeros(1, S2.npts)];

figure(1); clf
UX = reshape(utot(1,:), size(XX));
UY = reshape(utot(2,:), size(XX));
h = pcolor(XX, YY, reshape(vecnorm(utot), size(XX))); h.EdgeColor = 'none';
hold on
hs = streamslice(XX, YY, UX, UY);
set(hs, 'Color', 'k', 'LineWidth', 1.5);
plot(srfrs, u_surf_norm); colorbar
title('Stokes mixed BC: |u|, z=0'); xlabel('x'); ylabel('y')
