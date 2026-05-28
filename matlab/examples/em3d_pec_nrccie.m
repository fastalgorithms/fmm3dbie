% EM3D_PEC_NRCCIE  PEC scattering from an ellipsoid via NRCCIE.

zk    = 7;
alpha = 0.5;
eps   = 1e-6;

S = geometries.ellipsoid([1,1,1.5],[2,2,2], [0;0;0], 8);
K = kernel3d.em3d('nrccie-bc', zk, alpha);

% corrections and solve
tic;
opts_self = []; opts_self.corrections = 1; opts_self.selfquad = 1;
[cors, novers] = surfermat(S, K, eps, opts_self);
cors = cors + 0.5*speye(size(cors));
fprintf('Corrections: '); toc;

[einc, hinc] = em3d.planewave(zk, [pi/3; pi/4], [1; 0], S);
rhs = nrccie_rhs(S, einc, hinc, alpha);

sys_app = @(dens) surfermatapply(S, K, dens, eps, novers, cors);
tic; sol = gmres(sys_app, rhs, [], eps/10, 1000); fprintf('GMRES: '); toc;

[Jvec, rho] = nrccie_to_cartesian(S, sol);
dens_eval = [Jvec; rho];

%%
% scattered field on z=0 slice
nplot = 200;
xx = linspace(-3, 3, nplot);
[XX, YY] = meshgrid(xx, xx);
targs.r = [XX(:).'; YY(:).'; zeros(1, nplot^2)];

Keval = kernel3d.em3d('nrccie-eval', zk);
tic;
Escat = surferkerneval(S, Keval, dens_eval(:), targs, eps, []);
toc;
Escat = reshape(Escat, 6, []);

[Einc, Hinc] = em3d.planewave(zk, [pi/3; pi/4], [1; 0], targs);
Etot = Einc + Escat(1:3, :);
Htot = Hinc + Escat(4:6, :);

% surface field from BCs: |E| = |rho|, |H| = |J|
Emag = abs(rho);

figure(1); clf
h = pcolor(XX, YY, reshape(vecnorm([Etot;Htot]), nplot, nplot));
h.EdgeColor = 'none';
colorbar; title('$|\mathbf{E}| + |\mathbf{H}|$ total','Interpreter','latex');
set(gca,'fontsize',14)
hold on; plot(S, Emag(:)); hold off

% field lines of Re(E) on the slice
Ex = reshape(real(Etot(1,:)), nplot, nplot);
Ey = reshape(real(Etot(2,:)), nplot, nplot);

figure(2); clf
h2 = pcolor(XX, YY, reshape(real(Etot(1,:)), nplot, nplot));
h2.EdgeColor = 'none';
colorbar;
title('$\Re E_x$ with field lines','Interpreter','latex');
set(gca,'fontsize',14)
hold on
sl = streamslice(XX, YY, Ex, Ey, 2);
set(sl, 'Color', 'k', 'linewidth', 1.5);
plot(S, Emag(:));
hold off
axis equal tight


% =========================================================================

function rhs = nrccie_rhs(S, einc, hinc, alpha)
% NRCCIE_RHS  NRCCIE right-hand side: (nxH - alpha*nxnxE, n.E) in orthonormal frame.
n = S.n(:,:);

nxH   = cross(n, hinc, 1);
ndotE = sum(n .* einc, 1);
nxnxE = ndotE .* n - einc;

zvec = nxH - alpha * nxnxE;

du = S.du(:,:);
ru = du ./ vecnorm(du);
rv = cross(n ./ vecnorm(n), ru, 1);

rhs1 = sum(ru .* zvec, 1);
rhs2 = sum(rv .* zvec, 1);
rhs3 = ndotE;

rhs = [rhs1; rhs2; rhs3];
rhs = rhs(:);
end


function [Jvec, rho] = nrccie_to_cartesian(S, sol)
% NRCCIE_TO_CARTESIAN  Convert orthonormal density to Cartesian J and rho.
npts = S.npts;
sol3 = reshape(sol, 3, npts);
j_ru = sol3(1,:);
j_rv = sol3(2,:);
rho  = sol3(3,:);

n  = S.n(:,:);
du = S.du(:,:);
ru = du ./ vecnorm(du);
rv = cross(n ./ vecnorm(n), ru, 1);

Jvec = bsxfun(@times, ru, j_ru) + bsxfun(@times, rv, j_rv);
end
