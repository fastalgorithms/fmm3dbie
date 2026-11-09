% Adjoint Lippmann-Schwinger demo for gravity wave scattering from a
% surface with variable gravity parameter h,
%
%   (2/h) w - phi = 0,   w = d phi/dz on z = 0,
%
% with h - g supported in the square D = [-1,1]^2, for two point
% sources.

g   = 8;
eps = 1e-8;

% flat square [-1,1]^2 in z = 0
S = geometries.square(5, 8);

% (2/h) delta + (2/h - 2/g) G_phi, the equation scaled by 1/g
c  = @(t) reshape(2./gravity_h(t, g), 1, 1, []);
Kp = kernel3d.gravity('s', g, 'gphi', 'rmax', 8);
Kv = (@(t) c(t) - 2/g)*Kp + kernel3d.eye(c);
K  = kernel3d.gravity('s', g, 'gs', 'rmax', 8);
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


function h = gravity_h(t, g)
% gravity parameter h = g/(1 + amp*P(x)*P(y)), P(s) = (1 - s^2)^8, which
% goes to g to 8th order on the boundary of the square
amp = 2;
x = t.r(1,:); y = t.r(2,:);
P = max(1 - x.^2, 0).^8.*max(1 - y.^2, 0).^8;

h = g./(1 + amp*P);
end
