%
%  helm2d_dir
%
%  Solves the variable-coefficient Helmholtz equation
%
%    -\Delta u - zk^2 u + V(r) u = f    in \Omega
%
%  with Dirichlet boundary conditions u = g on \partial\Omega,
%  using a coupled volume and boundary integral equation formulation.
%

%% Geometry and problem definition

zk = sqrt(2);

S = geometries.disk([],[],[4 4 4],8);

nch = 4*4;
cparams = []; cparams.ta = pi/nch; cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t), nch, cparams);
chnkr = sort(chnkr);

V = eval_gauss(S.r);

eps = 1e-9;

%% Quadrature

start = tic;
A = helm2d.slp_matgen(S, zk, eps);
l11 = -eye(S.npts) + V.*A;

h2d_d = kernel.helm2d('d', zk);
l12 = V.*chunkerkernevalmat(chnkr, h2d_d, S.r(1:2,:));

l21 = helm2d.v2b_dir(zk, S, chnkr, eps);

l22 = -0.5*eye(chnkr.npt) + chunkermat(chnkr, h2d_d);
fprintf('%5.2e s : time to assemble b2b matrix\n', toc(start))

%% Right hand side and solve

% Manufactured solution: u = sin(x) sin(y)
lhs = [l11, l12; l21, l22];

rhs_vol = (-2 + zk^2 + V(:).').*sin(S.r(1,:)).*sin(S.r(2,:));
rhs_bc  = sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));

rhs = [rhs_vol rhs_bc].';

start = tic;
sol = gmres(lhs, rhs, [], 1e-10, 200);
fprintf('%5.2e s : time for dense gmres\n', toc(start))

% Evaluate solution and compute error
mu  = sol(1:S.npts);
rho = sol(S.npts+1:end);
u = A*mu + chunkerkerneval(chnkr, h2d_d, rho, S.r(1:2,:));
u = real(u);

ref_u = (sin(S.r(1,:)).*sin(S.r(2,:))).';
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
