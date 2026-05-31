% HELM3D_NEU  Sound-hard scattering from a ring of ellipsoids.
%

zk  = 5;
eps = 1e-6;

% reference obstacle
S0 = geometries.ellipsoid([1, 1.2, 0.8], [3,3,3], [0;0;0], 8);

% 6 obstacles in a ring with random rotations
nobj   = 6;
R_ring = 4;
rng(42);
theta   = (0:nobj-1)*(2*pi/nobj);
centres = R_ring * [cos(theta); sin(theta); zeros(1,nobj)];

Sarr = repmat(S0, 1, nobj);
for k = 1:nobj
    Sarr(k) = affine_transf(rotate(S0, 2*pi*rand(3,1)), eye(3), centres(:,k));
end

K  = kernel3d.helm3d('sp', zk);
n1   = S0.npts;
ntot = nobj * n1;

% self-correction once, reused for all diagonal blocks
tic;
opts_self = []; opts_self.corrections = 1; opts_self.selfquad = 1;
[cors_self, novers_self] = surfermat(S0, K, eps, opts_self);
fprintf('Self-correction: '); toc;

% off-diagonal corrections
% tic;
% opts_off = []; opts_off.corrections = 1; opts_off.selfquad = 0; opts_off.adaptive_correction = 1;
% [cors, novers] = surfermat(Sarr, K, eps, opts_off);
% fprintf('Off-diagonal corrections: '); toc;
cors = sparse(ntot,ntot);

for k = 1:nobj
    ii = (k-1)*n1 + (1:n1);
    cors(ii, ii) = cors_self;
end
cors = cors - 0.5*speye(ntot);

Smerge = merge(Sarr);
novers_merged = {repmat(novers_self{1,1}, nobj, 1)};

%% Get rhs and solve
[uinc, gradu_inc] = helm3d.planewave(zk, [1;0;0], Smerge);
dudn_inc = sum(Smerge.n(:,:) .* gradu_inc, 1).';
rhs = -dudn_inc;

sys_app = @(sig) surfermatapply(Smerge, K, sig, eps, novers_merged, cors);
tic; [sig, ~, relres] = gmres(sys_app, rhs, [], eps/10, 1000);
fprintf('GMRES relres=%.2e: ', relres); toc;

%%
% scattered field on z=0 slice
nplot  = 300;
extent = 7;
xx = linspace(-extent, extent, nplot);
[XX, YY] = meshgrid(xx, xx);
targs.r = [XX(:).'; YY(:).'; zeros(1, nplot^2)];

Ks = kernel3d.helm3d('s', zk);
tic; uscat = surferkerneval(Smerge, Ks, sig, targs, eps, []); fprintf('uscat: '); toc;
uscat = reshape(uscat, nplot, nplot);

uinc_slice = reshape(helm3d.planewave(zk, [1;0;0], targs), nplot, nplot);
utot = uinc_slice + uscat;

% Dirichlet trace on surface: u_tot = S_k[sigma] + u_inc
tic; usurf = surferkerneval(Smerge, Ks, sig, Smerge, eps) + uinc; fprintf('usurf: '); toc;

figure(1); clf
h = pcolor(XX, YY, real(utot));
h.EdgeColor = 'none';
colorbar;
title('$\Re(u^{\rm tot})$','Interpreter','latex');
set(gca,'fontsize',14); axis equal tight
hold on; plot(Smerge, real(usurf(:))); hold off
clim([-2,2])