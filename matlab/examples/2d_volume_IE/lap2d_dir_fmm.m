%
%  lap2d_dir_fmm
%
%  Solves the variable-coefficient Laplace equation
%
%    -\Delta u + V(r) u = f    in \Omega
%
%  with Dirichlet boundary conditions u = g on \partial\Omega,
%  using a coupled volume and boundary integral equation formulation.
%
%  This script uses the FMM to apply the v2v, v2b, and b2v operators.
%

%% Geometry and problem definition

S = geometries.disk([],[],[4 4 4],6);

nch = 4*4;
cparams = []; cparams.ta = pi/nch; cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t), nch, cparams);
chnkr = sort(chnkr);

V = eval_gauss(S.r);

eps = 1e-8;

%% Quadrature corrections

% Volume to volume
start = tic;
[Av2v_cor, nover] = lap2d.get_quad_cor_sub(S, eps);
fprintf('%5.2e s : time to compute v2v quadrature correction\n', toc(start))

% Boundary to volume (DLP: double-layer from boundary to volume)
l2d_d = kernel('l', 'd');
start = tic;
opts = []; opts.corrections = true;
Ab2v_cor = chunkerkernevalmat(chnkr, l2d_d, S.r(1:2,:), opts);
fprintf('%5.2e s : time to compute b2v quadrature correction\n', toc(start))

% Volume to boundary (SLP: single-layer from volume to boundary)
targinfo = [];
targinfo.r = [chnkr.r(:,:); 0*chnkr.r(1,:)];
start = tic;
[Av2b_cor, nover_v2b] = lap2d.get_quad_cor_v2b_dir(S, targinfo, eps);
fprintf('%5.2e s : time to compute v2b quadrature correction\n', toc(start))

% Boundary to boundary (dense, small): -1/2 I + D
l2d_dp = kernel('l', 'd');
start = tic;
lhs_22 = -0.5*eye(chnkr.npt) + chunkermat(chnkr, l2d_dp);
fprintf('%5.2e s : time to assemble b2b matrix\n', toc(start))

%% FMM-accelerated apply functions

v2v_apply = @(mu)  lap2d.apply_v2v(S, mu, Av2v_cor, nover, eps);
b2v_apply = @(rho) lap2d.apply_b2v_dir(S, chnkr, rho, Ab2v_cor, eps);
v2b_apply = @(mu)  lap2d.apply_v2b_dir(S, targinfo, mu, Av2b_cor, nover_v2b, eps);

%% Block operator and right-hand side

nv = S.npts;
nb = chnkr.npt;

lhs_11  = @(mu)  -mu(:) + V(:) .* v2v_apply(mu);
lhs_12  = @(rho)  V(:) .* b2v_apply(rho);
lhs_21  = @(mu)   v2b_apply(mu);
lhs_22f = @(rho)  lhs_22 * rho(:);

lhs = @(x) [lhs_11(x(1:nv))  + lhs_12(x(nv+1:end)); ...
            lhs_21(x(1:nv))  + lhs_22f(x(nv+1:end))];

% Manufactured solution: u = sin(x) sin(y)
rhs_vol = (-2 + V(:).') .* sin(S.r(1,:)) .* sin(S.r(2,:));
rhs_bc  = sin(chnkr.r(1,:)) .* sin(chnkr.r(2,:));

rhs = [rhs_vol rhs_bc].';

%% Solve

start = tic;
sol = gmres(lhs, rhs, [], 1e-10, 200);
fprintf('%5.2e s : time for FMM-accelerated GMRES\n', toc(start))

%% Evaluate solution and compute error

mu  = sol(1:nv);
rho = sol(nv+1:end);
u = v2v_apply(mu) + chunkerkerneval(chnkr, l2d_d, rho, S.r(1:2,:));
u = real(u);

ref_u = (sin(S.r(1,:)) .* sin(S.r(2,:))).';
err = abs(u - ref_u(:)) / max(abs(u));
fprintf('max relative error: %5.2e\n', max(err))

figure(1); clf
scatter(S.r(1,:), S.r(2,:), 8, log10(err));
title('log_{10} relative error'); colorbar

%%
function val = eval_gauss(targ)
    val = exp(-10*targ(1,:).^2 - 10*targ(2,:).^2);
    val = val(:);
end
