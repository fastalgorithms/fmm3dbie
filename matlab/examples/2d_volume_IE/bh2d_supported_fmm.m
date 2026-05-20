%
%  bh2d_supported_fmm
%
%  Solves the variable-coefficient biharmonic equation
%
%    \Delta^2 u + V(r) u = f    in \Omega
%
%  with simply-supported plate boundary conditions
%    u = g,  M_n[u] = h    on \partial\Omega,
%
%  where M_n is the bending moment operator depending on Poisson ratio nu.
%
%  The v2v block uses the dense matrix bh2d.v2v_matgen. The b2v, v2b, and
%  b2b blocks are also assembled as dense matrices and the full system is
%  solved with GMRES.
%

%% Geometry and problem definition

S = geometries.disk([],[],[4 4 4],8);

cparams = []; cparams.maxchunklen = 0.5;
chnkr = chunkerfunc(@(t) ellipse(t), cparams);
chnkr = sort(chnkr);
chnkr = makedatarows(chnkr, 2);

V = eval_gauss(S.r);

zk = 0;
nu = 0.1;
eps = 1e-8;

%% Curvature info (including derivatives for supported plate kernel)

kappa = signed_curvature(chnkr);
kp    = arclengthder(chnkr, kappa);
kpp   = arclengthder(chnkr, kp);

chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;

%% Quadrature and matrix assembly

% Volume to volume
start = tic;
A = bh2d.v2v_matgen(S, zk, eps);
v2v = eye(S.npts) + V .* A;
fprintf('%5.2e s : time to assemble v2v matrix\n', toc(start))

% Boundary to volume
fkern_b2v = @(s,t) chnk.flex2d.kern(0, s, t, 'supported_plate_eval', nu);
targetinfo_b2v = [];
targetinfo_b2v.r = S.r(1:2,:);
targetinfo_b2v.n = S.n(1:2,:);

start = tic;
b2v = chunkerkernevalmat(chnkr, fkern_b2v, targetinfo_b2v);
b2v = V .* b2v;
fprintf('%5.2e s : time to assemble b2v matrix\n', toc(start))

% Volume to boundary
targinfo = [];
targinfo.r = [chnkr.r(:,:); 0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:); 0*chnkr.n(1,:)];

start = tic;
v2b_dir   = bh2d.v2b_matgen_dir(S, zk, targinfo, eps);
v2b_supp2 = bh2d.v2b_matgen_supp2(S, zk, nu, targinfo, eps);
v2b = zeros(2*chnkr.npt, S.npts);
v2b(1:2:end,:) = v2b_dir;
v2b(2:2:end,:) = v2b_supp2;
fprintf('%5.2e s : time to assemble v2b matrix\n', toc(start))

% Boundary to boundary
fkern_log    = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_log', nu);
fkern_smooth = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_smooth', nu);

kappa_vec = signed_curvature(chnkr);
kappa_vec = kappa_vec(:);

opts_log    = []; opts_log.sing = 'log';
opts_smooth = []; opts_smooth.quad = 'native'; opts_smooth.sing = 'smooth';

start = tic;
M  = chunkermat(chnkr, fkern_log, opts_log);
M2 = chunkermat(chnkr, fkern_smooth, opts_smooth);

c0 = (nu - 1)*(nu + 3)*(2*nu - 1) / (2*(3 - nu));
M(2:2:end,1:2:end) = M(2:2:end,1:2:end) + M2 ...
    - c0 .* kappa_vec.^2 .* eye(chnkr.npt) ...
    - max(zk^2)/2 * eye(chnkr.npt);
b2b = M + 0.5*eye(2*chnkr.npt);
fprintf('%5.2e s : time to assemble b2b matrix\n', toc(start))

%% Block operator and right-hand side

nv = S.npts;
nb = 2*chnkr.npt;

lhs = [v2v, b2v;
       v2b, b2b];

% Manufactured solution: u = sin(x) sin(y)
rhs_vol = (4 + V(:).') .* sin(S.r(1,:)) .* sin(S.r(2,:));

rhs_bc = zeros(nb, 1);
% Dirichlet: u = sin(x)sin(y)
rhs_bc(1:2:end) = sin(chnkr.r(1,:)) .* sin(chnkr.r(2,:));
% Bending moment M_n[u] = -(1+nu)sin(x)sin(y) + 2(1-nu)cos(x)cos(y) n_x n_y
rhs_bc(2:2:end) = -(1+nu)*sin(chnkr.r(1,:)).*sin(chnkr.r(2,:)) + ...
    2*(1-nu)*cos(chnkr.r(1,:)).*cos(chnkr.r(2,:)) ...
    .*chnkr.n(1,:).*chnkr.n(2,:);

rhs = [rhs_vol, rhs_bc.'].';

%% Solve

start = tic;
sol = gmres(lhs, rhs, [], 1e-10, 100);
fprintf('%5.2e s : time for dense GMRES\n', toc(start))

%% Evaluate solution and compute error

mu  = sol(1:nv);
rho = sol(nv+1:end);

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_eval', nu);
u = A * mu + chunkerkerneval(chnkr, ikern, rho, S.r(1:2,:));

ref_u = (sin(S.r(1,:)) .* sin(S.r(2,:))).';
err = abs(u - ref_u(:)) / max(abs(u));
fprintf('max relative error: %5.2e\n', max(err))

figure; clf
scatter(S.r(1,:), S.r(2,:), 8, log10(err));
title('log_{10} relative error'); colorbar

%%
function val = eval_gauss(targ)
    val = exp(-10*targ(1,:).^2 - 10*targ(2,:).^2);
    val = val(:);
end
