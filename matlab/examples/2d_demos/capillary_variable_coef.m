% Adjoint Lippmann-Schwinger demo for capillary-gravity wave scattering
% from a surface with variable tension b,
%
%   -div(b grad w) + gamma w - phi = 0,   w = d phi/dz on z = 0,
%
% with b - beta supported in the square D = [-1,1]^2, for two point
% sources.

beta  = 0.01;
gamma = 0.1;
eps   = 1e-8;

% flat square [-1,1]^2 in z = 0
S = geometries.square(5, 8);

% cap_varbeta kernel contains its identity term
pfun = @(t) tension_pfun(t, beta);
Kv = kernel3d.capillary('cap_varbeta', beta, gamma, 'pfun', pfun, ...
                        'rmax', 8);
K  = kernel3d.capillary('s', beta, gamma, 'gs', 'rmax', 8);
A  = surfermat(S, Kv, eps);
V  = surfermat(S, K, eps);

% two point sources outside D
src = []; src.r = [-2 1.5; 0.5 2];

% solve
sigma = A \ -Kv.eval(src, S);

% reciprocity, w_s at source 2 due to source 1 and vice versa
ws = surferkerneval(S, K, sigma, src, eps);
fprintf('reciprocity error %.2e\n', abs(ws(2,1) - ws(1,2))/abs(ws(2,1)));

% evaluate and plot the total fields
nplot = 100;
[xx, yy] = meshgrid(linspace(-3, 3, nplot));
out  = max(abs(xx), abs(yy)) > 1;
targ = []; targ.r = [xx(out).'; yy(out).'];
xb = [-1 1 1 -1 -1]; yb = [-1 -1 1 1 -1];

wtin  = K.eval(src, S) + V*sigma;
wtout = K.eval(src, targ) + surferkerneval(S, K, sigma, targ, eps);

figure(1); clf
tiledlayout(1,2)
for j = 1:2
    wplot = nan(size(xx));
    wplot(out) = wtout(:,j);

    nexttile
    pcolor(xx, yy, real(wplot)); hold on
    plot(S, real(wtin(:,j)));
    plot(xb, yb, 'k','LineWidth',2)
    plot(src.r(1,j), src.r(2,j), 'r*','LineWidth',2)
    shading interp; view(2); axis equal tight; colorbar
    clim(max(abs(real(wtin(:,j))))*[-1 1])
    title(sprintf('Re w_{tot}, source %d', j))
end


function p = tension_pfun(t, beta)
% surface tension b = beta (1 + amp*P(x)*P(y)), P(s) = (1 - s^2)^8, which
% goes to beta to 8th order on the boundary of the square
amp = 8;
x = t.r(1,:); y = t.r(2,:);
wx = max(1 - x.^2, 0); wy = max(1 - y.^2, 0);
Px = wx.^8;  dPx = -16*x.*wx.^7;
Py = wy.^8;  dPy = -16*y.*wy.^7;

p = [];
p.beta  = beta*(1 + amp*Px.*Py);
p.dbeta = beta*amp*[dPx.*Py; Px.*dPy];
end
