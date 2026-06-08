% stok3d_mixed_bc.m  --  Stokes mixed boundary condition demo
%
% S1: traction-free (u = S[sigma1],           BIE: 1/2 sigma1 + S'[sigma1] + L = rhs1)
% S2: no-slip       (u = aS[sigma2]+D[sigma2], BIE: 1/2 sigma2 + (aS+D)    = rhs2)
%
% Analytic tests on single bodies before the coupled solve:
%   Test 1: S1 alone, traction BC. Source inside S2 region -> extinction inside S1.
%   Test 2: S2 alone, no-slip BC. Source inside S2 -> extinction outside S2.

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

kern_s  = kernel3d('stokes', 's');
kern_c  = kernel3d('stokes', 'c',     [alpha;1]);
kern_sp = kernel3d('stokes', 'strac');
kern_cp = kernel3d('stokes', 'ctrac', [alpha;1]);

ntest = 8;
th = linspace(0,2*pi,ntest+1); th = th(1:end-1);

%% S1 surface moments (needed for L stabilizer)

wi1 = S1.wts;
area1     = sum(wi1);
centroid1 = (S1.r * wi1) / area1;
dx1 = S1.r - centroid1;
rmoi_inv1 = inv(trace(dx1*(dx1.*wi1')')*eye(3) - dx1*(dx1.*wi1')');
icomps1 = [1; S1.npatches+1];

%% Test 1: S1 alone, traction BC
% Source outside S1, solve zero exterior traction BC.
% Verify via force balance: int f_tot dS = 0 (since f_tot = 0 pointwise).

opts1 = struct; opts1.corrections = 1; opts1.unif_nover = 1;
[Qsp1, Sovers1] = surfermat(S1, kern_sp, eps, opts1);
Qsp1 = Qsp1 + 0.5*speye(3*S1.npts);
mv1 = @(x) surfermatapply(S1, kern_sp, x, eps, Sovers1, Qsp1) + ...
    reshape(stok3d.traction.apply_mob_L(reshape(x,3,[]), S1, icomps1, area1, centroid1, rmoi_inv1), [], 1);

src1 = struct('r', [3;0;0], 'd', [1;0;0]);
t_inc1 = reshape(kern_sp.eval(src1, S1) * src1.d, 3, S1.npts);
[sol1, ~, rr1] = gmres(mv1, -t_inc1(:), [], eps, 200);
fprintf('Test 1 (traction BC), gmres relres = %.2e\n', rr1);

%% Test 2: S2 alone, no-slip BC
% Source inside S2 -> solve zero velocity -> u_tot = 0 outside S2 (exterior extinction).

opts2 = struct; opts2.corrections = 1;
[Qc2, Sovers2] = surfermat(S2, kern_c, eps, opts2);
Qc2 = Qc2 + 0.5*speye(3*S2.npts);
mv2 = @(x) surfermatapply(S2, kern_c, x, eps, Sovers2, Qc2);

src2 = struct('r', [3;0;0.3], 'd', [0;1;0]);
u_inc2 = reshape(kern_s.eval(src2, S2) * src2.d, 3, S2.npts);
[sol2, ~, rr2] = gmres(mv2, -u_inc2(:), [], eps, 200);
fprintf('Test 2 (no-slip BC),  gmres relres = %.2e\n', rr2);

% Exterior extinction: source inside -> u_tot = 0 outside
ttest2_out = struct('r', [3+2*cos(th); 2*sin(th); zeros(1,ntest)]);
u_inc2_out  = reshape(kern_s.eval(src2, ttest2_out) * src2.d, 3, ntest);
u_scat2_out = reshape(surferkerneval(S2, kern_c, sol2, ttest2_out, eps), 3, ntest);
fprintf('Test 2: max |u_tot| outside S2 = %.2e  (should be ~0)\n', ...
    max(vecnorm(u_inc2_out + u_scat2_out)));

%% Block-level verification: embed single-body solutions into full system

fprintf('\nBuilding quadrature corrections...\n'); tic;
opts_build = struct; opts_build.corrections = 1; opts_build.unif_nover = 1;
n1 = 3*S1.npts;  n2 = 3*S2.npts;  npts_tot = n1 + n2;
[Qcor, Sovers] = surfermat(srfrs, kerns, eps, opts_build);
Qcor = Qcor + 0.5*speye(npts_tot);
fprintf('  t = %.2f s\n', toc);

matvec = @(x) stok3d.mixed_matvec(x, srfrs, alpha, eps, Qcor, Sovers, ...
    area1, centroid1, rmoi_inv1);

% For each density alone, compare matvec blocks against surferkerneval directly.
% (1,1) block: (1/2 I + S' + L)[sigma1] at S1
% (2,1) block: S[sigma1] at S2
% (1,2) block: (aS+D)'[sigma2] at S1
% (2,2) block: (1/2 I + aS+D)[sigma2] at S2

x1 = [sol1; zeros(n2,1)];
mv_x1 = matvec(x1);
A11_mv = mv_x1(1:n1);
A21_mv = mv_x1(n1+1:end);
A11_ref = mv1(sol1);   % single-body matvec
A21_ref = reshape(surferkerneval(S1, kern_s, sol1, S2, eps), [], 1);
fprintf('(1,1) block: matvec vs single-body mv  = %.2e\n', norm(A11_mv - A11_ref)/norm(A11_ref));
fprintf('(2,1) block: matvec vs surferkerneval  = %.2e\n', norm(A21_mv - A21_ref)/norm(A21_ref));

x2 = [zeros(n1,1); sol2];
mv_x2 = matvec(x2);
A12_mv = mv_x2(1:n1);
A22_mv = mv_x2(n1+1:end);
A12_ref = reshape(surferkerneval(S2, kern_cp, sol2, S1, eps), [], 1);
A22_ref = mv2(sol2);   % single-body matvec
fprintf('(1,2) block: matvec vs surferkerneval  = %.2e\n', norm(A12_mv - A12_ref)/norm(A12_ref));
fprintf('(2,2) block: matvec vs single-body mv  = %.2e\n', norm(A22_mv - A22_ref)/norm(A22_ref));

%% Main solve: uniform incident flow

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
