%
%  helm2d_neu_fmm_test
%
%  Test for the FMM-accelerated variable-coefficient Helmholtz Neumann
%  solver. Solves
%
%    -\Delta u - zk^2 u + V(r) u = f    in \Omega
%
%  with Neumann boundary conditions du/dn = g on \partial\Omega,
%  using the manufactured solution u = sin(x) sin(y), and asserts that the
%  max relative error is below 1e-7.
%

%% Geometry and problem definition

zk = pi;

S = geometries.disk([],[],[3 3 3],6);

nch = 4;
cparams = []; cparams.ta = pi/nch; cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t), nch, cparams);
chnkr = sort(chnkr);

V = eval_gauss(S.r);

eps = 1e-8;

%% Quadrature corrections

% Volume to volume
[v2v_cor, nover] = helm2d.get_quad_cor_sub(S, zk, eps);

% Boundary to volume
h2d_s = kernel('h', 's', zk);
opts = []; opts.corrections = true;
Ab2v_cor = chunkerkernevalmat(chnkr, h2d_s, S.r(1:2,:), opts);

% Volume to boundary
[Av2b_cor, nover_v2b] = helm2d.get_quad_cor_v2b_neu(S, zk, chnkr, eps);

% Boundary to boundary (dense, small)
h2d_sp = kernel('h', 'sp', zk);
lhs_22 = 0.5*eye(chnkr.npt) + chunkermat(chnkr, h2d_sp);

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

lhs = @(x) [lhs_11(x(1:nv))  + lhs_12(x(nv+1:end)); ...
            lhs_21(x(1:nv))  + lhs_22f(x(nv+1:end))];

% Manufactured solution: u = sin(x) sin(y)
rhs_vol = (-2 + zk^2 + V(:).')  .* sin(S.r(1,:))    .* sin(S.r(2,:));
rhs_bc  =  cos(chnkr.r(1,:))   .* sin(chnkr.r(2,:)) .* chnkr.n(1,:) ...
         + sin(chnkr.r(1,:))   .* cos(chnkr.r(2,:)) .* chnkr.n(2,:);

rhs = [rhs_vol rhs_bc].';

%% Solve

sol = gmres(lhs, rhs, [], 1e-10, 200);

%% Evaluate solution and check error

mu  = sol(1:nv);
rho = sol(nv+1:end);
u = v2v_apply(mu) + chunkerkerneval(chnkr, h2d_s, rho, S.r(1:2,:));
u = real(u);

ref_u = (sin(S.r(1,:)) .* sin(S.r(2,:))).';
err = max(abs(u - ref_u(:))) / max(abs(ref_u));
fprintf('helm2d_neu_fmm_test: max relative error = %5.2e\n', err)

assert(err < 1e-7, ...
    sprintf('helm2d_neu_fmm_test FAILED: max relative error %5.2e >= 1e-7', err));

disp('helm2d_neu_fmm_test PASSED')

%%
function val = eval_gauss(targ)
    val = exp(-10*targ(1,:).^2 - 10*targ(2,:).^2);
    val = val(:);
end
