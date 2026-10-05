% Adjoint Lippmann-Schwinger demo for flexural wave scattering.
%
%   Delta^2 u - zk^4 (1 + q) u = 0,   q supported in the square D = [-1,1]^2
%
% Using the representation: u = u_inc + V[sigma], where
% 
%   V[sigma](x) = int_D G(x,y) sigma(y) dA.
%
% The system is solved with FMM-accelerated GMRES

zk  = 4;
eps = 1e-8;

% flat square [-1,1]^2 in z = 0
S = geometries.square(24, 8);

% contrast, vanishes to 8th order on the boundary of the square
qfun = @(x,y) 8*(max(1 - x.^2, 0).*max(1 - y.^2, 0)).^8;

% incident plane wave
uinc = @(x,y) exp(1i*zk*x);

x = S.r(1,:).'; y = S.r(2,:).';
q = qfun(x,y);

% volume potential on the square: near-field corrections and oversampling
K = kernel3d('flex2d', 's', zk);
opts = []; opts.corrections = true;
[cors, objover] = surfermat(S, K, eps, opts);

qt   = @(t) reshape(qfun(t.r(1,:), t.r(2,:)), 1, 1, []);
Kq   = -zk^4 * (qt * K);
corq = speye(S.npts) - zk^4 * q(:).'.* cors;
lhs  = @(dens) surfermatapply(S, Kq, dens, eps, objover, corq);

% solve
rhs   = zk^4 * q .* uinc(x,y);
sigma = gmres(lhs, rhs, [], eps, 200);

% evaluate field
usin = surfermatapply(S, K, sigma, eps, objover, cors);
utin = uinc(x,y) + usin;

nplot = 100;
[xx, yy] = meshgrid(linspace(-3, 3, nplot));
out  = max(abs(xx), abs(yy)) > 1;
targ = []; targ.r = [xx(out).'; yy(out).'; zeros(1, nnz(out))];

usout = nan(size(xx));
usout(out) = surferkerneval(S, K, sigma, targ, eps);
utout = uinc(xx, yy) + usout;

% plot
figure(1); clf
tiledlayout(1,2)
xb = [-1 1 1 -1 -1]; yb = [-1 -1 1 1 -1];

nexttile
pcolor(xx, yy, real(usout)); hold on
plot(S, real(usin));
plot(xb, yb, 'k')
shading interp; view(2); axis equal tight; colorbar
title('Re u_{scat}')

nexttile
pcolor(xx, yy, real(utout)); hold on
plot(S, real(utin));
plot(xb, yb, 'k')
shading interp; view(2); axis equal tight; colorbar
title('Re u_{tot}')
