%
%  helm2d_neu
%
%  Solves the variable-coefficient Helmholtz equation
%
%    -\Delta u - zk^2 u + V(r) u = f    in \Omega
%
%  with Neumann boundary conditions du/dn = g on \partial\Omega,
%  using a coupled volume and boundary integral equation formulation.
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

%% Quadrature

% Volume to volume
start = tic;
A = helm2d.slp_matgen(S, zk, eps);
l11 = -eye(S.npts) + V.*A;
fprintf('%5.2e s : time to assemble v2v matrix\n', toc(start))

% Boundary to volume
h2d_s = kernel.helm2d('s', zk);
start = tic;
l12 = V.*chunkerkernevalmat(chnkr, h2d_s, S.r(1:2,:));
fprintf('%5.2e s : time to assemble b2v matrix\n', toc(start))

% Volume to boundary
start = tic;
l21 = helm2d.v2b_neu(zk, S, chnkr, eps);
fprintf('%5.2e s : time to assemble v2b matrix\n', toc(start))

% Boundary to boundary
h2d_sp = kernel.helm2d('sp', zk);
start = tic;
l22 = 0.5*eye(chnkr.npt) + chunkermat(chnkr, h2d_sp);
fprintf('%5.2e s : time to assemble b2b matrix\n', toc(start))

%% Right hand side and solve

% Manufactured solution: u = sin(x) sin(y)
lhs = [l11, l12; l21, l22];

rhs_vol = (-2 + zk^2 + V(:).').*sin(S.r(1,:)).*sin(S.r(2,:));
rhs_bc  = cos(chnkr.r(1,:)).*sin(chnkr.r(2,:)).*chnkr.n(1,:) ...
        + sin(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(2,:);

rhs = [rhs_vol rhs_bc].';

start = tic;
sol = gmres(lhs, rhs, [], 1e-10, 200);
fprintf('%5.2e s : time for dense gmres\n', toc(start))

% Evaluate solution and compute error
mu  = sol(1:S.npts);
rho = sol(S.npts+1:end);
u = A*mu + chunkerkerneval(chnkr, h2d_s, rho, S.r(1:2,:));
u = real(u);

ref_u = (sin(S.r(1,:)).*sin(S.r(2,:))).';
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
