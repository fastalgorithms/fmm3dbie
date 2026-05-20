%
%  flex2d_clamped
%
%  Solves the variable-coefficient flexural plate equation
%
%    a \Delta^2 u - b \Delta u - c u + V(r) u = f    in \Omega
%
%  with clamped boundary conditions u = g_0, du/dn = g_1 on \partial\Omega,
%  using a coupled volume and boundary integral equation formulation.
%  The volume domain \Omega is a disk discretized as a surfer object S,
%  and the boundary \partial\Omega is discretized as a chunkie object chnkr.
%  The manufactured solution is u = sin(x) sin(y).
%
%  Requires chunkie

%% Geometry and problem definition

a = 1;
b = 0.7;
c = 1/pi;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));
zk = [zk1 zk2];

S = geometries.disk([],[],[4 4 4],6);

cparams = []; cparams.maxchunklen = 4/max(abs(zk));
chnkr = chunkerfunc(@(t) ellipse(t), cparams);
chnkr = sort(chnkr);

V = eval_gauss(S.r);

eps = 1e-8;

%% Quadrature

% Volume to volume
start = tic;
A = flex2d.v2v_matgen(S, zk, eps);
l11 = a*eye(S.npts) + V.*A;
fprintf('%5.2e s : time to assemble v2v matrix\n', toc(start))

% Boundary to volume
fkern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
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
getnearquad = @(varargin) flex2d.getnearquad(varargin{1:17}, [], zk, varargin{18}, 'clamped');
rhskern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_bcs');
[flex_v2b_cor, norderup] = getnearquad_kern(S, rhskern, eps, getnearquad, chnkr);
[S_over, xinterp] = oversample(S, S.norders + norderup);
l21 = flex_v2b_cor + (rhskern(S_over, chnkr).*S_over.wts(:).')*xinterp;
fprintf('%5.2e s : time to assemble v2b matrix\n', toc(start))

% Boundary to boundary
fkern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');
kappa = signed_curvature(chnkr);
kappa = kappa(:);

opts = []; opts.sing = 'log';

start = tic;
b2b = chunkermat(chnkr, fkern, opts);
l22 = b2b + 0.5*eye(2*chnkr.npt);
l22(2:2:end,1:2:end) = l22(2:2:end,1:2:end) - kappa.*eye(chnkr.npt);
fprintf('%5.2e s : time to assemble b2b matrix\n', toc(start))

%% Right hand side and solve

% Manufactured solution: u = sin(x) sin(y)
% rhs = (a*(-4) - b*(-2) + c - V) * u = (-4a + 2b + c - V) * sin(x) sin(y)
lhs = [l11, l12;
       l21, l22];

rhs_vol = (-4*a - 2*b + c - V(:).').*sin(S.r(1,:)).*sin(S.r(2,:));
rhs_bc = zeros(2*chnkr.npt, 1);
rhs_bc(1:2:end) = -sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));
rhs_bc(2:2:end) = -cos(chnkr.r(1,:)).*sin(chnkr.r(2,:)).*chnkr.n(1,:) ...
                  - sin(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(2,:);

rhs = zeros(S.npts + 2*chnkr.npt, 1);
rhs(1:S.npts) = rhs_vol;
rhs(S.npts+1:end) = rhs_bc;

start = tic;
sol = gmres(lhs, rhs, [], 1e-10, 100);
fprintf('%5.2e s : time for dense gmres\n', toc(start))

% Evaluate solution and compute error
ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
u = A*sol(1:S.npts) + b2v*sol(S.npts+1:end);
u = real(u);

ref_u = (sin(S.r(1,:)).*sin(S.r(2,:))).';
err = abs(u + ref_u(:)) / max(abs(u));
fprintf('max relative error: %5.2e\n', max(err))

figure; clf
scatter(S.r(1,:), S.r(2,:), 8, log10(err));
title('log_{10} relative error'); colorbar

%%
function val = eval_gauss(targ)
    val = exp(-10*targ(1,:).^2 - 10*targ(2,:).^2);
    val = val(:);
end

