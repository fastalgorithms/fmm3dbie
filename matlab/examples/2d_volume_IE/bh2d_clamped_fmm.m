%
%  bh2d_clamped_fmm
%
%  Solves the variable-coefficient biharmonic equation
%
%    \Delta^2 u + V(r) u = f    in \Omega
%
%  with clamped plate boundary conditions
%    u = g,  du/dn = h    on \partial\Omega,
%
%  using a coupled volume and boundary integral equation formulation.
%
%  The v2v block uses an FMM-accelerated operator (bh2d.apply_v2v via
%  lfmm2d). The b2v, v2b, and b2b blocks are assembled as dense matrices.
%

%% Geometry and problem definition

S = geometries.disk([],[],[4 4 4],8);

cparams = []; cparams.maxchunklen = 0.5;
chnkr = chunkerfunc(@(t) ellipse(t), cparams);
chnkr = sort(chnkr);

V = eval_gauss(S.r);

zk = 0;
eps = 1e-8;

%% Quadrature

% Volume to volume: FMM-accelerated
start = tic;
[Av2v_cor, nover] = bh2d.get_quad_cor_sub(S, zk, eps);
fprintf('%5.2e s : time to compute v2v quadrature correction\n', toc(start))

v2v_apply = @(mu) bh2d.apply_v2v(S, zk, mu, Av2v_cor, nover, eps);
lhs_11 = @(mu) mu(:) + V(:) .* v2v_apply(mu);

% Boundary to volume (dense)
fkern_b2v = @(s,t) chnk.flex2d.kern(0, s, t, 'clamped_plate_eval');
targetinfo_b2v = [];
targetinfo_b2v.r = S.r(1:2,:);
targetinfo_b2v.n = S.n(1:2,:);
start = tic;
b2v = V(:) .* chunkerkernevalmat(chnkr, fkern_b2v, targetinfo_b2v);
fprintf('%5.2e s : time to assemble b2v matrix\n', toc(start))

% Volume to boundary (dense)
targinfo = [];
targinfo.r = [chnkr.r(:,:); 0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:); 0*chnkr.n(1,:)];
start = tic;
v2b_dir = bh2d.v2b_matgen_dir(S, zk, targinfo, eps);
v2b_neu = bh2d.v2b_matgen_neu(S, zk, targinfo, eps);
v2b = zeros(2*chnkr.npt, S.npts);
v2b(1:2:end,:) = v2b_dir;
v2b(2:2:end,:) = v2b_neu;
fprintf('%5.2e s : time to assemble v2b matrix\n', toc(start))

% Boundary to boundary (dense)
fkern_b2b = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');
kappa = signed_curvature(chnkr);
kappa = kappa(:);
opts = []; opts.sing = 'log';
start = tic;
b2b = chunkermat(chnkr, fkern_b2b, opts);
b2b = b2b + 0.5*eye(2*chnkr.npt);
b2b(2:2:end,1:2:end) = b2b(2:2:end,1:2:end) - kappa.*eye(chnkr.npt);
fprintf('%5.2e s : time to assemble b2b matrix\n', toc(start))

%% Block operator and right-hand side

nv = S.npts;
nb = 2*chnkr.npt;

lhs = @(x) [lhs_11(x(1:nv))    + b2v*x(nv+1:end); ...
            v2b*x(1:nv)         + b2b*x(nv+1:end)];

% Manufactured solution: u = sin(x) sin(y)
rhs_vol = (4 + V(:).') .* sin(S.r(1,:)) .* sin(S.r(2,:));
rhs_bc  = zeros(nb, 1);
rhs_bc(1:2:end) = sin(chnkr.r(1,:))  .* sin(chnkr.r(2,:));
rhs_bc(2:2:end) = cos(chnkr.r(1,:))  .* sin(chnkr.r(2,:)) .* chnkr.n(1,:) ...
                + sin(chnkr.r(1,:))  .* cos(chnkr.r(2,:)) .* chnkr.n(2,:);

rhs = [rhs_vol rhs_bc'].';

%% Solve

start = tic;
sol = gmres(lhs, rhs, [], 1e-10, 200);
fprintf('%5.2e s : time for FMM-accelerated GMRES\n', toc(start))

%% Evaluate solution and compute error

mu  = sol(1:nv);
rho = sol(nv+1:end);
fkern_eval = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
u = v2v_apply(mu) + chunkerkerneval(chnkr, fkern_eval, rho, S.r(1:2,:));

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
