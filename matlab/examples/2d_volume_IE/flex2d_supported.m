%
%  flex2d_supported
%
%  Solves the variable-coefficient flexural plate equation
%
%    a \Delta^2 u - b \Delta u - c u + V(r) u = f    in \Omega
%
%  with simply supported boundary conditions on \partial\Omega,
%  using a coupled volume and boundary integral equation formulation.
%  The volume domain \Omega is a disk discretized as a surfer object S,
%  and the boundary \partial\Omega is discretized as a chunkie object chnkr.
%  The manufactured solution is u = sin(x) sin(y).
%
%  Requires chunkie


%% Geometry and problem definition

a = -0.5;
b = 0.7;
c = 1/pi;
nu = 0.3;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));
zk = [zk1, zk2];

S = geometries.disk([],[],[4 4 4],8);

cparams = []; cparams.maxchunklen = 4/max(abs(zk));
chnkr = chunkerfunc(@(t) ellipse(t), cparams);
chnkr = sort(chnkr);
chnkr = makedatarows(chnkr, 2);

kappa = signed_curvature(chnkr);
kp  = arclengthder(chnkr, kappa);
kpp = arclengthder(chnkr, kp);
chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;

V = eval_gauss(S.r);

eps = 1e-8;

%% Quadrature

% Volume to volume
start = tic;
A = flex2d.v2v_matgen(S, zk, eps);
l11 = a*eye(S.npts) + V.*A;
fprintf('%5.2e s : time to assemble v2v matrix\n', toc(start))

% Boundary to volume
fkern = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_eval', nu);
targetinfo = [];
targetinfo.r = S.r(1:2,:);
targetinfo.n = S.n(1:2,:);

opts = []; opts.eps = eps;
start = tic;
b2v = chunkerkernevalmat(chnkr, fkern, targetinfo,opts);
l12 = V.*b2v;
fprintf('%5.2e s : time to assemble b2v matrix\n', toc(start))

% Volume to boundary
start = tic;
v2b_dir  = flex2d.v2b_matgen_dir(S, zk, chnkr, eps);
v2b_supp = flex2d.v2b_matgen_supp2(S, zk, nu, chnkr, eps);
l21 = zeros(2*chnkr.npt, S.npts);
l21(1:2:end,:) = v2b_dir;
l21(2:2:end,:) = v2b_supp;
fprintf('%5.2e s : time to assemble v2b matrix\n', toc(start))

% Boundary to boundary
fkern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_log', nu);
fkern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_smooth', nu);

opts  = []; opts.sing  = 'log';
opts2 = []; opts2.quad = 'native'; opts2.sing = 'smooth';

start = tic;
M  = chunkermat(chnkr, fkern1, opts);
M2 = chunkermat(chnkr, fkern2, opts2);

c0 = (nu - 1)*(nu + 3)*(2*nu - 1) / (2*(3 - nu));
M(2:2:end,1:2:end) = M(2:2:end,1:2:end) + M2 - c0.*kappa(:).^2.*eye(chnkr.npt) + b/(2*a)*eye(chnkr.npt);
l22 = M + 0.5*eye(2*chnkr.npt);
fprintf('%5.2e s : time to assemble b2b matrix\n', toc(start))

%% Right hand side and solve

% Manufactured solution: u = sin(x) sin(y)
lhs = [l11, l12;
       l21, l22];

rhs_vol = (4*a + 2*b - c + V(:).').*sin(S.r(1,:)).*sin(S.r(2,:));
rhs_bc = zeros(2*chnkr.npt, 1);
rhs_bc(1:2:end) = sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));
rhs_bc(2:2:end) = -(1+nu)*sin(chnkr.r(1,:)).*sin(chnkr.r(2,:)) + ...
                   2*(1-nu)*cos(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(1,:).*chnkr.n(2,:);

rhs = zeros(S.npts + 2*chnkr.npt, 1);
rhs(1:S.npts) = rhs_vol;
rhs(S.npts+1:end) = rhs_bc;

start = tic;
sol = gmres(lhs, rhs, [], 1e-10, 100);
fprintf('%5.2e s : time for dense gmres\n', toc(start))

% Evaluate solution and compute error
ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_eval', nu);
u = A*sol(1:S.npts) + b2v* sol(S.npts+1:end);
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
