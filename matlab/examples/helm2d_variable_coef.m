% Adjoint Lippmann-Schwinger demo for variable-coefficient Helmholtz.
%
%   Delta u + zk^2 (1 + q) u = 0,   q supported in the unit disk D
%
% Using the representation: u = u_inc + V[sigma],  where
%
%   V[sigma](x) = int_D G(x,y) sigma(y) dA.

zk  = 4;
eps = 1e-8;

% flat disk in z = 0
S = geometries.disk([], [], [6 3 7], 8);

% contrast, vanishes to high order at r = 1
qfun = @(x,y) 2*max(1 - x.^2 - y.^2, 0).^4;

% incident plane wave
uinc = @(x,y) exp(1i*zk*x);

x = S.r(1,:).'; y = S.r(2,:).';
q = qfun(x,y);

% volume potential on the disk
K = kernel3d('helm2d', 's', zk);
V = surfermat(S, K, eps);

% solve
A     = eye(S.npts) - zk^2 * (q .* V);
rhs   = zk^2 * q .* uinc(x,y);
sigma = A \ rhs;

% evaluate field
usin = V*sigma;
utin = uinc(x,y) + usin;

nplot = 100;
[xx, yy] = meshgrid(linspace(-3, 3, nplot));
out  = xx.^2 + yy.^2 > 1;
targ = []; targ.r = [xx(out).'; yy(out).'; zeros(1, nnz(out))];

usout = nan(size(xx));
usout(out) = surferkerneval(S, K, sigma, targ, eps);
utout = uinc(xx, yy) + usout;

% plot
figure(1); clf
tiledlayout(1,2)
th = linspace(0, 2*pi);

nexttile
pcolor(xx, yy, real(usout)); hold on
plot(S, real(usin));
plot(cos(th), sin(th), 'k')
shading interp; view(2); axis equal tight; colorbar
title('Re u_{scat}')

nexttile
pcolor(xx, yy, real(utout)); hold on
plot(S, real(utin));
plot(cos(th), sin(th), 'k')
shading interp; view(2); axis equal tight; colorbar
title('Re u_{tot}')
