%
%  bh2d_free_fmm
%
%  Solves the variable-coefficient biharmonic equation
%
%    \Delta^2 u + V(r) u = f    in \Omega
%
%  with free plate boundary conditions
%    M_n[u] = g,  V_n[u] = h    on \partial\Omega,
%
%  where M_n is the bending moment operator and V_n is the Kirchhoff
%  shear force operator, both depending on Poisson ratio nu.
%
%  The v2v block uses an FMM-accelerated operator (bh2d.apply_v2v via
%  lfmm2d). The b2v, v2b, and b2b blocks are assembled as dense matrices.
%

%% Geometry and problem definition

S = geometries.disk([],[],[4 4 4],6);

cparams = []; cparams.maxchunklen = 0.5;
chnkr = chunkerfunc(@(t) ellipse(t), cparams);
chnkr = sort(chnkr);
chnkr = makedatarows(chnkr, 2);

V = eval_gauss(S.r);

zk = 0;
nu = 0.3;
eps = 1e-8;

%% Curvature and tangent info

kappa = signed_curvature(chnkr);

dx = chnkr.d(1,:);
dy = chnkr.d(2,:);
ds = sqrt(dx.*dx + dy.*dy);
taux = dx ./ ds;
tauy = dy ./ ds;

%% Quadrature and matrix assembly

% Volume to volume: FMM-accelerated
start = tic;
[Av2v_cor, nover] = bh2d.get_quad_cor_sub(S, zk, eps);
fprintf('%5.2e s : time to compute v2v quadrature correction\n', toc(start))

v2v_apply = @(mu) bh2d.apply_v2v(S, zk, mu, Av2v_cor, nover, eps);
lhs_11 = @(mu) mu(:) + V(:) .* v2v_apply(mu);

% Boundary to volume
fkern_b2v = @(s,t) chnk.flex2d.kern(0, s, t, 'free_plate_eval', nu);
targetinfo_b2v = [];
targetinfo_b2v.r = S.r(1:2,:);
targetinfo_b2v.n = S.n(1:2,:);

% Hilbert transform needed to recover the correct representation
double  = @(s,t) chnk.lap2d.kern(s, t, 'd');
hilbert = @(s,t) chnk.lap2d.kern(s, t, 'hilb');
opts_pv = []; opts_pv.sing = 'pv';
H = chunkermat(chnkr, hilbert, opts_pv);

% densmap maps 2-component boundary density (eta, zeta) to the
% 3-component representation (eta, -H*eta, zeta) used by the kernel
npt = chnkr.npt;
densmap = zeros(3*npt, 2*npt);
densmap(1:3:end, 1:2:end) = eye(npt);
densmap(3:3:end, 2:2:end) = eye(npt);
densmap(2:3:end, 1:2:end) = -H;

start = tic;
b2v = chunkerkernevalmat(chnkr, fkern_b2v, targetinfo_b2v);
b2v = V .* (b2v * densmap);
fprintf('%5.2e s : time to assemble b2v matrix\n', toc(start))

% Volume to boundary
targinfo = [];
targinfo.r  = chnkr.r;
targinfo.n  = [chnkr.n(:,:); 0*chnkr.n(1,:)];
targinfo.kappa = kappa(:);
targinfo.du = [taux; tauy; zeros(size(taux))];

[sxyz, patch_inds, uvsloc, dists, flags] = get_closest_pts(S, targinfo);

targinfo.nu = nu;
start = tic;
v2b_supp = bh2d.matgen(S, zk, 'supp2', targinfo, eps, patch_inds, uvsloc);
v2b_free = bh2d.matgen(S, zk, 'free2', targinfo, eps, patch_inds, uvsloc);
v2b = zeros(2*chnkr.npt, S.npts);
v2b(1:2:end,:) = v2b_supp;
v2b(2:2:end,:) = v2b_free;
fprintf('%5.2e s : time to assemble v2b matrix\n', toc(start))

% Boundary to boundary
fkern_b2b = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate', nu);
opts_log = []; opts_log.sing = 'log';

start = tic;
sysmat1 = chunkermat(chnkr, fkern_b2b, opts_log);
opts_d  = []; opts_d.sing = 'log';
D = chunkermat(chnkr, double, opts_d);

b2b = zeros(2*chnkr.npt);
b2b(1:2:end,1:2:end) = sysmat1(1:4:end,1:2:end) - sysmat1(3:4:end,1:2:end)*H ...
                      + 2*((1+nu)/2)^2 * D * D;
b2b(2:2:end,1:2:end) = sysmat1(2:4:end,1:2:end) - sysmat1(4:4:end,1:2:end)*H;
b2b(1:2:end,2:2:end) = sysmat1(1:4:end,2:2:end) + sysmat1(3:4:end,2:2:end);
b2b(2:2:end,2:2:end) = sysmat1(2:4:end,2:2:end) + sysmat1(4:4:end,2:2:end);

Djump = -[-1/2 + (1/8)*(1+nu)^2, 0; 0, 1/2];
Djump = kron(eye(chnkr.npt), Djump);
b2b = Djump + b2b;
fprintf('%5.2e s : time to assemble b2b matrix\n', toc(start))

%% Block operator and right-hand side

nv = S.npts;
nb = 2*chnkr.npt;

lhs = @(x) [lhs_11(x(1:nv))   + b2v*x(nv+1:end); ...
            v2b*x(1:nv)        + b2b*x(nv+1:end)];

% Manufactured solution: u = sin(x) sin(y)
rhs_vol = (4 + V(:).') .* sin(S.r(1,:)) .* sin(S.r(2,:));

sinx = sin(chnkr.r(1,:));  cosx = cos(chnkr.r(1,:));
siny = sin(chnkr.r(2,:));  cosy = cos(chnkr.r(2,:));
nx = chnkr.n(1,:);  ny = chnkr.n(2,:);
taux = taux(:).';  tauy = tauy(:).';  kappa = kappa(:).';

% Bending moment M_n[u]
rhs_bc = zeros(nb, 1);
rhs_bc(1:2:end) = -(1+nu)*siny.*siny + ...
    2*(1-nu)*cosx.*cosy.*chnkr.n(1,:).*chnkr.n(2,:);

% Kirchhoff shear V_n[u]
gsxx  = (-sinx.*siny);
gsxy  = ( cosx.*cosy);
gsyy  = (-sinx.*siny);
gsxxx = (-cosx.*siny);
gsxxy = (-sinx.*cosy);
gsxyy = (-cosx.*siny);
gsyyy = (-sinx.*cosy);

rhs_bc(2:2:end) = ...
    gsxxx.*(nx.*nx.*nx) + gsxxy.*(3*nx.*nx.*ny) + ...
    gsxyy.*(3*nx.*ny.*ny) + gsyyy.*(ny.*ny.*ny) + ...
    (2-nu) .* ( gsxxx.*(taux.*taux.*nx) + ...
                gsxxy.*(taux.*taux.*ny + 2*taux.*tauy.*nx) + ...
                gsxyy.*(2*taux.*tauy.*ny + tauy.*tauy.*nx) + ...
                gsyyy.*(tauy.*tauy.*ny) ) + ...
    (1-nu) .* kappa .* ( ...
        gsxx.*(taux.*taux) + gsxy.*(2*taux.*tauy) + gsyy.*(tauy.*tauy) - ...
        gsxx.*(nx.*nx)     - gsxy.*(2*nx.*ny)     - gsyy.*(ny.*ny) );

rhs = [rhs_vol, rhs_bc.'].';

%% Solve

start = tic;
sol = gmres(lhs, rhs, [], 1e-10, 200);
fprintf('%5.2e s : time for FMM-accelerated GMRES\n', toc(start))

%% Evaluate solution and compute error

mu  = sol(1:nv);
rho = sol(nv+1:end);

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu);
u = v2v_apply(mu) + chunkerkerneval(chnkr, ikern, densmap*rho, S.r(1:2,:));

ref_u = (sin(S.r(1,:)) .* sin(S.r(2,:))).';
err = abs(u - ref_u(:)) / max(abs(u));
fprintf('max relative error: %5.2e\n', max(err))

figure(1); clf
plot(S, log10(S.patch_max(err)));
title('log_{10} relative error'); colorbar

%%
function val = eval_gauss(targ)
    val = exp(-10*targ(1,:).^2 - 10*targ(2,:).^2);
    val = val(:);
end
