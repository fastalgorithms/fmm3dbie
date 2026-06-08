% stok3d_mobility.m  --  Stokes mobility problem demo
%
% Representation: u = S[sigma] + S[sigma_0]
%   sigma_0 = F_i/|Gamma_i| + tau_i^{-1} T_i x (x - x_{c,i})
%
% BIE: 1/2 sigma + S'[sigma] + L[sigma] = -(1/2 sigma_0 + S'[sigma_0])
%
% Analytic check (ncomp=1): isolated unit sphere, U=F/(6pi), Omega=T/(8pi).
% For ncomp>1: verify no-slip residual on the surface.

ncomp  = 2;
shifts = [0, 4; 0, 0; 0, 0];

surfers = cell(1, ncomp);
for ic = 1:ncomp
    surfers{ic} = geometries.sphere(1, 3, shifts(:,ic).', 6, 11);
end
S = merge([surfers{1}, surfers{2}]);
fprintf('npts = %d\n', S.npts);

icomps = zeros(ncomp+1, 1);
icomps(1) = 1;
for ic = 1:ncomp
    icomps(ic+1) = icomps(ic) + surfers{ic}.npatches;
end

rng(42);
forces  = randn(3, ncomp);
torques = randn(3, ncomp);
trans_vels_ex = forces  / (6*pi);
rot_vels_ex   = torques / (8*pi);

eps = 1e-8;

%% Surface moments per component

npts = S.npts;
pts  = @(ic) S.ixyzs(icomps(ic)) : S.ixyzs(icomps(ic+1))-1;

area     = zeros(1, ncomp);
centroid = zeros(3, ncomp);
rmoi_inv = zeros(3, 3, ncomp);

for ic = 1:ncomp
    xpts = S.r(:, pts(ic));
    wi   = S.wts(pts(ic));

    area(ic)        = sum(wi);
    centroid(:, ic) = (xpts * wi) / area(ic);

    dx   = xpts - centroid(:, ic);
    dxdx = dx * (dx .* wi')';
    rmoi_inv(:,:,ic) = inv(trace(dxdx)*eye(3) - dxdx);
end

%% sigma_0

sigma0 = zeros(3, npts);
for ic = 1:ncomp
    dx = S.r(:, pts(ic)) - centroid(:,ic);
    omega0 = rmoi_inv(:,:,ic) * torques(:,ic);
    sigma0(:, pts(ic)) = forces(:,ic)/area(ic) + cross(omega0*ones(1,numel(pts(ic))), dx);
end

%% Build S' quadrature correction (reused for RHS and matvec)

kerns_strac = kernel3d('stokes', 'strac');
fprintf('Building strac quadrature correction...\n'); tic;
opts_build = struct; opts_build.corrections = 1;
[Qstrac, Sovers_strac] = surfermat(S, kerns_strac, eps, opts_build);
Qstrac = Qstrac - 0.5*speye(3*npts);
fprintf('  t = %.2f s\n', toc);

%% RHS: -(1/2 sigma_0 + S'[sigma_0])

fprintf('RHS eval...\n'); tic;
Sprime_sigma0 = reshape(surfermatapply(S, kerns_strac, sigma0(:), eps, Sovers_strac, Qstrac), 3, npts);
rhs_vec = -Sprime_sigma0(:);
fprintf('  t = %.2f s\n', toc);

%% GMRES: (1/2 I + S'[.] + L[.]) sigma = rhs

matvec = @(x) matvec_mob(x, S, kerns_strac, eps, Sovers_strac, Qstrac, ...
    icomps, area, centroid, rmoi_inv);
fprintf('GMRES...\n'); tic;
[sol_vec, flag, relres, iter] = gmres(matvec, rhs_vec, [], eps, 200);
fprintf('  t = %.2f s, flag=%d, relres=%.2e, iter=%d\n', toc, flag, relres, iter(end));

sigma     = reshape(sol_vec, 3, npts);
sigma_tot = sigma + sigma0;

%% Velocity on surface and rigid-body extraction

kerns_s = kernel3d('stokes', 's');
fprintf('Evaluating u on surface...\n'); tic;
uvel = reshape(surferkerneval(S, kerns_s, sigma_tot(:), S, eps), 3, npts);
fprintf('  t = %.2f s\n', toc);

trans_vels = zeros(3, ncomp);
rot_vels   = zeros(3, ncomp);

for ic = 1:ncomp
    wi = S.wts(pts(ic));
    ui = uvel(:, pts(ic));
    dx = S.r(:, pts(ic)) - centroid(:,ic);

    trans_vels(:,ic) = (ui * wi) / area(ic);
    udiff = ui - trans_vels(:,ic);
    rot_vels(:,ic) = rmoi_inv(:,:,ic) * (cross(dx, udiff, 1) * wi);
end

%% No-slip residual (always) + analytic check (ncomp=1)

fprintf('\n component   no-slip rms\n');
for ic = 1:ncomp
    dx  = S.r(:, pts(ic)) - centroid(:,ic);
    res = uvel(:, pts(ic)) - trans_vels(:,ic) - cross(rot_vels(:,ic)*ones(1,numel(pts(ic))), dx);
    fprintf('    %d            %.2e\n', ic, rms(res(:)));
end

if ncomp == 1
    fprintf('\n component  |U err|/|U|  |Om err|/|Om|\n');
    for ic = 1:ncomp
        et = norm(trans_vels(:,ic) - trans_vels_ex(:,ic)) / norm(trans_vels_ex(:,ic));
        er = norm(rot_vels(:,ic)   - rot_vels_ex(:,ic))   / norm(rot_vels_ex(:,ic));
        fprintf('    %d         %.2e       %.2e\n', ic, et, er);
    end
end

%% Visualization

nplot = 80;
[XX, YY] = meshgrid(linspace(-3,7,nplot)+1e-2, linspace(-3,3,nplot)+3e-2);
targs_plot = []; targs_plot.r = [XX(:).'; YY(:).'; zeros(1,numel(XX))];

Dlap    = lap3d.eval(S, 'double', ones(npts,1), targs_plot, 1e-3);
outside = abs(Dlap) < abs(Dlap+1);
targs_out = []; targs_out.r = targs_plot.r(:, outside);

fprintf('\nPlot eval...\n'); tic;
uplot = reshape(surferkerneval(S, kerns_s, sigma_tot(:), targs_out, eps), 3, []);
fprintf('  t = %.2f s\n', toc);

ufull = nan(3, numel(XX));
ufull(:, outside) = uplot;
UX = reshape(ufull(1,:), size(XX));
UY = reshape(ufull(2,:), size(XX));

figure(2); clf
h = pcolor(XX, YY, reshape(vecnorm(ufull), size(XX))); h.EdgeColor = 'none';
hold on
hs = streamslice(XX, YY, UX, UY);
set(hs, 'Color', 'k', 'LineWidth', 1.5);
plot(S, zeros(npts,1)); colorbar
title('Stokes mobility: streamlines, z=0'); xlabel('x'); ylabel('y')

%% -----------------------------------------------------------------------

function v = matvec_mob(x, S, kern, eps, Sovers, Qstrac, ...
    icomps, area, centroid, rmoi_inv)
% Apply (1/2 I + S'[.] + L[.]) to x.  1/2 I is baked into Qstrac.
sig = reshape(x, 3, S.npts);
Sp  = surfermatapply(S, kern, x, eps, Sovers, Qstrac);
L   = stok3d.traction.apply_mob_L(sig, S, icomps, area, centroid, rmoi_inv);
v   = Sp(:) + L(:);
end
