% Adjoint Lippmann-Schwinger demo for flexural wave scattering from a
% plate with variable rigidity alpha.
%
%   Delta(alpha Delta u) - zk^4 u
%      - (1-nu)(alpha_xx u_yy - 2 alpha_xy u_xy + alpha_yy u_xx) = 0,
%
% with alpha - 1 supported in the square D = [-1,1]^2.
%
% Using the representation: u = u_inc + V[sigma], where
%
%   V[sigma](x) = int_D G(x,y) sigma(y) dA,   (Delta^2 - zk^4) G = delta,
%
% the density solves
%
%   alpha sigma + K[sigma] = -L u_inc
%
% with K the 'varcoef' kernel of kernel3d.flex2d and L the operator
% above. The system is solved with FMM-accelerated GMRES

zk  = 4;
nu  = 0.3;
eps = 1e-8;

% flat square [-1,1]^2 in z = 0
S = geometries.square(24, 8);

% incident plane wave
uinc = @(x,y) exp(1i*zk*x);

x = S.r(1,:).'; y = S.r(2,:).';

% plate coefficients at the nodes, and the target coefficients of K
p  = plate_pfun(S, zk);
cf = flex2d.plate_coefs(p, nu, zk);

% kernel and corrections
Kv = kernel3d('flex2d', 'varcoef', zk, nu, @(t) plate_pfun(t, zk));
opts = []; opts.corrections = true;
[corv, objover] = surfermat(S, Kv, eps, opts);
corv = spdiags(p.alpha(:), 0, S.npts, S.npts) + corv;
lhs  = @(dens) surfermatapply(S, Kv, dens, eps, objover, corv);

% right hand side -L u_inc. For the plane wave, (Delta^2 - zk^4) u_inc = 0,
% d_x Lap u_inc = -i zk^3 u_inc, Lap u_inc = u_xx = -zk^2 u_inc, and the
% other derivatives vanish
ui  = uinc(x,y);
rhs = -(-1i*zk^3*cf(1,:).' - zk^2*(cf(3,:).' + cf(5,:).') + cf(7,:).').*ui;

% solve
sigma = gmres(lhs, rhs, [], eps, 200);

% evaluate field
K = kernel3d('flex2d', 's', zk);
usin = surfermatapply(S, K, sigma, eps);
utin = ui + usin;

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


function p = plate_pfun(t, zk)
% plate coefficients alpha = 1 + amp*P(x)*P(y), P(s) = (1 - s^2)^8, which
% goes to 1 to 8th order on the boundary of the square, and beta = zk^4
amp = 8;
x = t.r(1,:); y = t.r(2,:);
wx = max(1 - x.^2, 0); wy = max(1 - y.^2, 0);
Px = wx.^8;  dPx = -16*x.*wx.^7;  d2Px = -16*wx.^7 + 224*x.^2.*wx.^6;
Py = wy.^8;  dPy = -16*y.*wy.^7;  d2Py = -16*wy.^7 + 224*y.^2.*wy.^6;

p = [];
p.alpha   = 1 + amp*Px.*Py;
p.dalpha  = amp*[dPx.*Py; Px.*dPy];
p.d2alpha = amp*[d2Px.*Py; dPx.*dPy; Px.*d2Py];
p.beta    = zk^4*ones(size(x));
end
