%
%  flex2d_free
%
%  Solves the variable-coefficient flexural plate equation
%
%    a \Delta^2 u - b \Delta u - c u + V(r) u = f    in \Omega
%
%  with free boundary conditions on \partial\Omega,
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
nu = 0.3;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));
zk = [zk1, zk2];
zk = 0;
S = geometries.disk([],[],[4 4 4],6);

cparams = []; cparams.maxchunklen = 4/max(abs(zk));
chnkr = chunkerfunc(@(t) ellipse(t), cparams);
chnkr = sort(chnkr);

V = eval_gauss(S.r);

eps = 1e-10;

%% Quadrature

% Volume to volume
start = tic;
A = flex2d.v2v_matgen(S, zk, eps);
l11 = a*eye(S.npts) + V.*A;
fprintf('%5.2e s : time to assemble v2v matrix\n', toc(start))

% Volume to boundary
targinfo = [];
targinfo.r = chnkr.r(:,:);
targinfo.n = chnkr.n(:,:);
targinfo.d = chnkr.d(:,:);
% targinfo.d = targinfo.d./vecnorm(targinfo.d);
targinfo.d2 = chnkr.d2(:,:);
kappa = signed_curvature(chnkr);
targinfo.kappa = kappa(:);
[sxyz, patch_inds, uvsloc, dists, flags] = get_closest_pts(S, targinfo);

start = tic;
v2b_supp = flex2d.v2b_matgen_supp2(S, zk, nu, targinfo, eps, patch_inds, uvsloc);
v2b_free = flex2d.v2b_matgen_free2(S, zk, nu, targinfo, eps, patch_inds, uvsloc);
l21 = zeros(2*chnkr.npt, S.npts);
l21(1:2:end,:) = v2b_supp;
l21(2:2:end,:) = v2b_free;
fprintf('%5.2e s : time to assemble v2b matrix\n', toc(start))

% Boundary to boundary
fkern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate', nu);
double  = @(s,t) chnk.lap2d.kern(s, t, 'd');
hilbert = @(s,t) chnk.lap2d.kern(s, t, 'hilb');

opts  = []; opts.sing  = 'log';
opts2 = []; opts2.sing = 'pv';

start = tic;
sysmat1 = chunkermat(chnkr, fkern1, opts);
D = chunkermat(chnkr, double,  opts);
H = chunkermat(chnkr, hilbert, opts2);

sysmat = zeros(2*chnkr.npt);
sysmat(1:2:end,1:2:end) = sysmat1(1:4:end,1:2:end) - sysmat1(3:4:end,1:2:end)*H + 2*((1+nu)/2)^2*D*D;
sysmat(2:2:end,1:2:end) = sysmat1(2:4:end,1:2:end) - sysmat1(4:4:end,1:2:end)*H;
sysmat(1:2:end,2:2:end) = sysmat1(1:4:end,2:2:end) + sysmat1(3:4:end,2:2:end);
sysmat(2:2:end,2:2:end) = sysmat1(2:4:end,2:2:end) + sysmat1(4:4:end,2:2:end);

Djump = -[-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];
Djump = kron(eye(chnkr.npt), Djump);
l22 = Djump + sysmat;
fprintf('%5.2e s : time to assemble b2b matrix\n', toc(start))

% Boundary to volume
fkern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu);
targetinfo = [];
targetinfo.r = S.r(1:2,:);
targetinfo.n = S.n(1:2,:);

% densmap converts 2-component boundary density (eta, zeta) to the
% 3-component representation (eta, -H*eta, zeta) used by the kernel
npt = chnkr.npt;
densmap = zeros(3*npt, 2*npt);
densmap(1:3:end, 1:2:end) = eye(npt);
densmap(3:3:end, 2:2:end) = eye(npt);
densmap(2:3:end, 1:2:end) = -H;

opts = []; opts.eps = eps;
start = tic;
b2v = chunkerkernevalmat(chnkr, fkern, targetinfo, opts);
l12 = V .* (b2v * densmap);
fprintf('%5.2e s : time to assemble b2v matrix\n', toc(start))

%% Right hand side and solve

% Manufactured solution: u = sin(x) sin(y)
lhs = [l11, l12;
       l21, l22];

rhs_vol = (4*a + 2*b - c + V(:).').*sin(S.r(1,:)).*sin(S.r(2,:));

nx = chnkr.n(1,:).';
ny = chnkr.n(2,:).';
dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';
ds = sqrt(dx.*dx + dy.*dy);
taux = dx./ds;
tauy = dy./ds;

hess = zeros(chnkr.npt, 3);
hess(:,1) = -sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));
hess(:,2) =  cos(chnkr.r(1,:)).*cos(chnkr.r(2,:));
hess(:,3) = -sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));

third = zeros(chnkr.npt, 4);
third(:,1) = -cos(chnkr.r(1,:)).*sin(chnkr.r(2,:));
third(:,2) = -sin(chnkr.r(1,:)).*cos(chnkr.r(2,:));
third(:,3) = -cos(chnkr.r(1,:)).*sin(chnkr.r(2,:));
third(:,4) = -sin(chnkr.r(1,:)).*cos(chnkr.r(2,:));

rhs_bc = zeros(2*chnkr.npt, 1);
rhs_bc(1:2:end) = (hess(:,1).*(nx.*nx) + hess(:,2).*(2*nx.*ny) + hess(:,3).*(ny.*ny)) + ...
                   nu.*(hess(:,1).*(taux.*taux) + hess(:,2).*(2*taux.*tauy) + hess(:,3).*(tauy.*tauy));
rhs_bc(2:2:end) = (third(:,1).*(nx.*nx.*nx) + third(:,2).*(3*nx.*nx.*ny) + ...
                   third(:,3).*(3*nx.*ny.*ny) + third(:,4).*(ny.*ny.*ny)) + ...
                  (2-nu).*(third(:,1).*(taux.*taux.*nx) + third(:,2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) + ...
                   third(:,3).*(2*taux.*tauy.*ny + tauy.*tauy.*nx) + third(:,4).*(tauy.*tauy.*ny)) + ...
                  (1-nu).*kappa(:).*(hess(:,1).*taux.*taux + hess(:,2).*(2*taux.*tauy) + hess(:,3).*tauy.*tauy ...
                   - hess(:,1).*nx.*nx - hess(:,2).*(2*nx.*ny) - hess(:,3).*ny.*ny);

rhs = zeros(S.npts + 2*chnkr.npt, 1);
rhs(1:S.npts) = rhs_vol;
rhs(S.npts+1:end) = rhs_bc;

start = tic;
sol = gmres(lhs, rhs, [], eps, 100);
fprintf('%5.2e s : time for dense gmres\n', toc(start))

% Evaluate solution and compute error
ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu);
u = A*sol(1:S.npts) + b2v*(densmap*sol(S.npts+1:end));
u = real(u);

ref_u = (sin(S.r(1,:)).*sin(S.r(2,:))).';
err = abs(u - ref_u(:)) / max(abs(u));
fprintf('max relative error: %5.2e\n', max(err))

figure; clf
plot(S, log10(S.patch_max(err)));
title('log_{10} relative error'); colorbar

%%
function val = eval_gauss(targ)
    val = exp(-10*targ(1,:).^2 - 10*targ(2,:).^2);
    val = val(:);
end
