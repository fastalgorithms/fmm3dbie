%
%  helm2d_neu_fmm
%
%  Solves the variable-coefficient Helmholtz equation
%
%    -\Delta u - zk^2 u + V(r) u = f    in \Omega
%
%  with Neumann boundary conditions du/dn = g on \partial\Omega,
%  using a coupled volume and boundary integral equation formulation.
%
%  This script uses the FMM to apply the v2v, v2b, and b2v operators
%

%% Geometry and problem definition

zk = pi;

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
[v2v_cor, nover] = helm2d.get_quad_cor_sub(S, zk, eps);
fprintf('%5.2e s : time to compute v2v quadrature correction\n', toc(start))

% Boundary to volume
h2d_s = kernel('h', 's', zk);
start = tic;
opts = []; opts.corrections = true;
Ab2v_cor = chunkerkernevalmat(chnkr, h2d_s, S.r(1:2,:), opts);
fprintf('%5.2e s : time to compute b2v quadrature correction\n', toc(start))

% Volume to boundary
start = tic;
[Av2b_cor, nover_v2b] = helm2d.get_quad_cor_v2b_neu(S, zk, chnkr, eps);
fprintf('%5.2e s : time to compute v2b quadrature correction\n', toc(start))

% Boundary to boundary (dense, small)
h2d_sp = kernel('h', 'sp', zk);
start = tic;
lhs_22 = 0.5*eye(chnkr.npt) + chunkermat(chnkr, h2d_sp);
fprintf('%5.2e s : time to assemble b2b matrix\n', toc(start))

%% FMM-accelerated apply functions

v2v_apply  = @(mu)  helm2d.apply_v2v(S, zk, mu, v2v_cor, nover, eps);
b2v_apply  = @(rho) helm2d.apply_b2v_neu(S, zk, chnkr, rho, Ab2v_cor, eps);
v2b_apply  = @(mu)  helm2d.apply_v2b_neu(S, zk, chnkr, mu, Av2b_cor, nover_v2b, eps);

%% Block operator and right-hand side

nv = S.npts;
nb = chnkr.npt;

lhs_11 = @(mu)  -mu(:) + V(:) .* v2v_apply(mu);
lhs_12 = @(rho)  V(:) .* b2v_apply(rho);
lhs_21 = @(mu)   v2b_apply(mu);
lhs_22f = @(rho) lhs_22 * rho(:);

lhs = @(x) [lhs_11(x(1:nv))     + lhs_12(x(nv+1:end)); ...
            lhs_21(x(1:nv))     + lhs_22f(x(nv+1:end))];

% Manufactured solution: u = sin(x) sin(y)
rhs_vol = (-2 + zk^2 + V(:).')  .* sin(S.r(1,:))    .* sin(S.r(2,:));
rhs_bc  =  cos(chnkr.r(1,:))   .* sin(chnkr.r(2,:)) .* chnkr.n(1,:) ...
         + sin(chnkr.r(1,:))   .* cos(chnkr.r(2,:)) .* chnkr.n(2,:);

rhs = [rhs_vol rhs_bc].';

%% Solve

start = tic;
sol = gmres(lhs, rhs, [], 1e-10, 200);
fprintf('%5.2e s : time for FMM-accelerated GMRES\n', toc(start))

%% Evaluate solution and compute error

mu  = sol(1:nv);
rho = sol(nv+1:end);
u = v2v_apply(mu) + chunkerkerneval(chnkr, h2d_s, rho, S.r(1:2,:));
u = real(u);

ref_u = (sin(S.r(1,:)) .* sin(S.r(2,:))).';
err = abs(u - ref_u(:)) / max(abs(u));
fprintf('max relative error: %5.2e\n', max(err))

figure(1); clf
plot_nodes(S,log10(err));
title('log_{10} relative error'); colorbar

%%
function val = eval_gauss(targ)
    val = exp(-10*targ(1,:).^2 - 10*targ(2,:).^2);
    val = val(:);
end
