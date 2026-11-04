% Test the surface wave kernels (capillary, iceflex, gravity): the fits
% against the root sums, finite depth corrections and the variable
% coefficient kernels.

% run ../startup.m

rng(42);

nu = 0.33;
rmax = 20;                  % for points()

cp = [1.3 -0.7];            % capillary beta, gamma
fp = [23 0.5];              % iceflex alpha, gamma
g = 1.7;

kerns = { ...
  'cap s gs',          kernel3d('capillary', 's',   cp(1), cp(2), 'gs', 'rmax', rmax); ...
  'cap s gphi',        kernel3d('capillary', 's',   cp(1), cp(2), 'gphi', 'rmax', rmax); ...
  'cap s s3d',         kernel3d('capillary', 's',   cp(1), cp(2), 's3d', 'rmax', rmax); ...
  'cap lap gphi',      kernel3d('capillary', 'lap', cp(1), cp(2), 'gphi', 'rmax', rmax); ...
  'cap lap gs+gphi',   kernel3d('capillary', 'lap', cp(1), cp(2), [0.3 0.5], 'rmax', rmax); ...
  'flex s gs',         kernel3d('iceflex', 's', fp(1), fp(2), 'gs', 'rmax', rmax); ...
  'flex s gphi',       kernel3d('iceflex', 's', fp(1), fp(2), 'gphi', 'rmax', rmax); ...
  'flex s s3d',        kernel3d('iceflex', 's', fp(1), fp(2), 's3d', 'rmax', rmax); ...
  'flex s gphi_bilap', kernel3d('iceflex', 's', fp(1), fp(2), 'gphi_bilap', 'rmax', rmax); ...
  'grav s gs',         kernel3d('gravity', 's',    g, 'gs', 'rmax', rmax); ...
  'grav s gphi',       kernel3d('gravity', 's',    g, 'gphi', 'rmax', rmax); ...
  'grav grad gs',      kernel3d('gravity', 'grad', g, 'gs', 'rmax', rmax); ...
};

% finite depth tables, one per family
d = 2;
fgws = {surfwave.fgw_setup(cp(1), cp(2), d, 3), ...
        surfwave.fgw_setup(fp(1), fp(2), d, 5), ...
        surfwave.fgw_setup(2/g, 0, d, 1)};

test_class(kerns, cp, fp, g);
test_eval(kerns, cp, fp, g);
test_floe(g, rmax, fgws{3});
test_finite_depth(cp, fp, g, rmax, d, fgws);
test_fgw_quadrature(cp, fp, g, d, fgws);
test_varcoef(cp, fp, g, nu, rmax);

fprintf('surfwave_kernelsTest: all checks passed.\n');


function check(name, err, tol)
fprintf('%-34s rel err %5.2e\n', name, err);
assert(err < tol, '%s: error %5.2e exceeds %5.2e', name, err, tol);
end


function [s, t] = points()
% points in the plane, targets with n, d, d2
s = []; s.r = [randn(2, 6); zeros(1, 6)];
th = 2*pi*rand(1, 5);
t = []; t.r = [3*cos(th); 3*sin(th); 0*th];
t.n = [cos(th); sin(th); 0*th];
t.d = [-3*sin(th); 3*cos(th); 0*th];
t.d2 = -t.r;
end


function test_class(kerns, cp, fp, g)
% zk the propagating root, 'rmax' required
zks = {'cap', @(z) cp(1)*z^3 + cp(2)*z - 1; 'flex', @(z) fp(1)*z^5 + fp(2)*z - 1; ...
       'grav', @(z) z - g/2};
for k = 1:size(kerns, 1)
    [name, K] = kerns{k,:};
    P = zks{strcmp(zks(:,1), strtok(name)), 2};
    assert(isreal(K.zk) && K.zk > 0 && abs(P(K.zk)) < 1e-12, ...
           '%s: zk = %g is not the propagating root', name, K.zk);
end

noargs = { ...
  'capillary', @() kernel3d.capillary('s', cp(1), cp(2)); ...
  'iceflex',   @() kernel3d.iceflex('s', fp(1), fp(2)); ...
  'gravity',   @() kernel3d.gravity('s', g); ...
};
for k = 1:size(noargs, 1)
    try
        noargs{k,2}();
        error('%s: no error without rmax', noargs{k,1});
    catch e
        assert(contains(e.message, 'rmax'), ...
               '%s without rmax: unexpected error ''%s''', noargs{k,1}, e.message);
    end
end
fprintf('%-34s ok\n', 'rmax required');

end


function test_eval(kerns, cp, fp, g)
% K.eval against the root sums
[rts_c, ejs_c] = surfwave.find_roots(cp(1), cp(2), 3);
[rts_f, ejs_f] = surfwave.find_roots(fp(1), fp(2), 5);
[s, t] = points();
ref = containers.Map();
xs = s.r(1:2,:);  xt = t.r(1:2,:);

lap3d = 1./(4*pi*sqrt((xt(1,:).' - xs(1,:)).^2 + (xt(2,:).' - xs(2,:)).^2));
rc = rts_c(:);  ec = ejs_c(:);
rf = rts_f(:);  ef = ejs_f(:);

[v, ~, hs] = surfwave.rootsum(rc, ec.*rc.^2/2, xs, xt);
ref('cap s gs') = v;
[v, ~, h] = surfwave.rootsum(rc, ec.*rc/4, xs, xt);
ref('cap s gphi') = v;
ref('cap lap gphi') = h(:,:,1) + h(:,:,3);
ref('cap lap gs+gphi') = 0.3*(hs(:,:,1) + hs(:,:,3)) + 0.5*(h(:,:,1) + h(:,:,3));
ref('cap s s3d') = lap3d + surfwave.rootsum(rc, ec/8, xs, xt);
ref('flex s gs') = surfwave.rootsum(rf, ef.*rf.^2/2, xs, xt, 1e-8);
[v, ~, ~, ~, f] = surfwave.rootsum(rf, ef.*rf/4, xs, xt, 1e-8);
ref('flex s gphi') = v;
ref('flex s gphi_bilap') = f(:,:,1) + 2*f(:,:,3) + f(:,:,5);
ref('flex s s3d') = lap3d + surfwave.rootsum(rf, ef/8, xs, xt, 1e-8);
% gravity: G_phi = rho/(2 pi r) + rho^2/4 K(rho r), rho = g/2, G_S = g G_phi
rho = g/2;
[v, gr] = surfwave.rootsum(rho, rho^2/4, xs, xt);
r = sqrt((xt(1,:).' - xs(1,:)).^2 + (xt(2,:).' - xs(2,:)).^2);
ref('grav s gphi') = v + rho./(2*pi*r);
ref('grav s gs') = g*ref('grav s gphi');
gr = gr - rho/(2*pi)*cat(3, xt(1,:).' - xs(1,:), xt(2,:).' - xs(2,:))./r.^3;
G = zeros(2*size(xt, 2), size(xs, 2));
G(1:2:end,:) = gr(:,:,1);  G(2:2:end,:) = gr(:,:,2);
ref('grav grad gs') = g*G;

for k = 1:size(kerns, 1)
    [name, K] = kerns{k,:};
    r = ref(name);
    A = K.eval(s, t);
    check(['eval ' name], norm(A - r, 'fro')/norm(r, 'fro'), 1e-10);
end
end


function test_floe(g, rmax, fgw)
% gravity 'floe', g (a Lap^2 + b - 1) G_phi, against the root sums
[s, t] = points();
xs = s.r(1:2,:);  xt = t.r(1:2,:);
rho = g/2;
[v, ~, ~, ~, f] = surfwave.rootsum(rho, rho^2/4, xs, xt, [], g);
bl = f(:,:,1) + 2*f(:,:,3) + f(:,:,5);
B = g*(0.7*bl + 0.1*v);
K = kernel3d.gravity('floe', g, 'coefs', [0.7 1.1], 'rmax', rmax);
check('grav floe eval', norm(K.eval(s, t) - B, 'fro')/norm(B, 'fro'), 1e-10);
check('grav floe diag', norm(K.diag(t) - 1.1*g)/(1.1*g), 1e-15);
% finite depth: plus the fgw corrections to G_phi and Lap^2 G_phi
B = B + g*(0.7*surfwave.fgw_green(fgw, 'bilap_gphi', xs, xt) ...
      + 0.1*surfwave.fgw_green(fgw, 'gphi', xs, xt));
K = kernel3d.gravity('floe', g, 'coefs', [0.7 1.1], 'fgw', fgw, 'rmax', rmax);
check('grav floe depth eval', norm(K.eval(s, t) - B, 'fro')/norm(B, 'fro'), 1e-10);
end


function test_finite_depth(cp, fp, g, rmax, d, fgws)
% iceflex, capillary and gravity 'depth' against the deep root sums plus SURFWAVE.FGW_GREEN
[rts, ejs] = surfwave.find_roots(fp(1), fp(2), 5);
rts = rts(:);  ejs = ejs(:);
[s, t] = points();
xs = s.r(1:2,:);  xt = t.r(1:2,:);
fgw = fgws{2};
gsdeep = surfwave.rootsum(rts, ejs.*rts.^2/2, xs, xt, 1e-8);
[gp, ~, ~, ~, f] = surfwave.rootsum(rts, ejs.*rts/4, xs, xt, 1e-8);
ref = {'gs',         gsdeep + surfwave.fgw_green(fgw, 'gs', xs, xt); ...
       'gphi',       gp + surfwave.fgw_green(fgw, 'gphi', xs, xt); ...
       'gphi_bilap', f(:,:,1) + 2*f(:,:,3) + f(:,:,5) ...
                     + surfwave.fgw_green(fgw, 'bilap_gphi', xs, xt)};
for k = 1:size(ref, 1)
    K = kernel3d.iceflex('s', fp(1), fp(2), ref{k,1}, 'fgw', fgw, 'rmax', rmax);
    check(['depth flex s ' ref{k,1}], ...
          norm(K.eval(s, t) - ref{k,2}, 'fro')/norm(ref{k,2}, 'fro'), 1e-10);
end
% zk is the finite depth propagating root
D = @(k) (fp(1)*k^4 + fp(2))*k*tanh(k*d) - 1;
check('depth zk root', abs(D(K.zk)), 1e-12);
% capillary
[rts, ejs] = surfwave.find_roots(cp(1), cp(2), 3);
rts = rts(:);  ejs = ejs(:);
fgw = fgws{1};
[v, gr] = surfwave.rootsum(rts, ejs.*rts/4, xs, xt);
[vd, grd] = surfwave.fgw_green(fgw, 'gphi', xs, xt);
gn = @(g) g(:,:,1).*t.n(1,:).' + g(:,:,2).*t.n(2,:).';
ref = {'s',  'gs',   surfwave.rootsum(rts, ejs.*rts.^2/2, xs, xt) ...
                     + surfwave.fgw_green(fgw, 'gs', xs, xt); ...
       'sp', 'gphi', gn(gr) + gn(grd)};
for k = 1:size(ref, 1)
    K = kernel3d.capillary(ref{k,1}, cp(1), cp(2), ref{k,2}, 'fgw', fgw, 'rmax', rmax);
    check(['depth cap ' ref{k,1} ' ' ref{k,2}], ...
          norm(K.eval(s, t) - ref{k,3}, 'fro')/norm(ref{k,3}, 'fro'), 1e-10);
end
D = @(k) (cp(1)*k^2 + cp(2))*k*tanh(k*d) - 1;
check('depth cap zk root', abs(D(K.zk)), 1e-12);
% gravity: G_phi is the fgw potential kernel for P(k) = 2k/g
fgw = fgws{3};
rho = g/2;
[v, gr] = surfwave.rootsum(rho, rho^2/4, xs, xt, [], g);
[vd, grd] = surfwave.fgw_green(fgw, 'gphi', xs, xt);
ref = {'s', 'gs', g*(v + vd);  'sp', 'gphi', gn(gr) + gn(grd)};
for k = 1:size(ref, 1)
    K = kernel3d.gravity(ref{k,1}, g, ref{k,2}, 'fgw', fgw, 'rmax', rmax);
    check(['depth grav ' ref{k,1} ' ' ref{k,2}], ...
          norm(K.eval(s, t) - ref{k,3}, 'fro')/norm(ref{k,3}, 'fro'), 1e-10);
end
D = @(k) 2*k/g*tanh(k*d) - 1;
check('depth grav zk root', abs(D(K.zk)), 1e-12);
% depth = inf is deep water
A = kernel3d.capillary('s', cp(1), cp(2), 'gs', 'depth', inf, 'rmax', rmax).eval(s, t);
B = kernel3d.capillary('s', cp(1), cp(2), 'gs', 'rmax', rmax).eval(s, t);
check('depth inf cap s gs', norm(A - B, 'fro')/norm(B, 'fro'), 1e-15);
end


function test_varcoef(cp, fp, g, nu, rmax)
% 'cap_varbeta' and 'flex_varalpha' eval and diag against the root sums
[s, t] = points();
xs = s.r(1:2,:);  xt = t.r(1:2,:);
x = xt(1,:).';  y = xt(2,:).';

[rts, ejs] = surfwave.find_roots(cp(1), cp(2), 3);
be = cp(1);  ga = cp(2);
pfun = @(t) struct('beta', be*(1 + 0.3*t.r(1,:).^2 + 0.2*t.r(2,:)), ...
                   'dbeta', be*[0.6*t.r(1,:); 0.2 + 0*t.r(1,:)]);
K = kernel3d.capillary('cap_varbeta', be, ga, 'pfun', pfun, 'rmax', rmax);
[gs, gr] = surfwave.rootsum(rts, ejs.*rts.^2/2, xs, xt);
gp = surfwave.rootsum(rts, ejs.*rts/4, xs, xt);
p = pfun(t);
B = (p.beta(:) - be).*(gp - ga/2*gs) ...
    - be/2*(p.dbeta(1,:).'.*gr(:,:,1) + p.dbeta(2,:).'.*gr(:,:,2));
check('cap_varbeta eval', norm(K.eval(s, t) - B, 'fro')/norm(B, 'fro'), 1e-10);
check('cap_varbeta diag', norm(K.diag(t) - p.beta(:))/norm(p.beta), 1e-15);

[rts, ejs] = surfwave.find_roots(fp(1), fp(2), 5);
al = fp(1);  ga = fp(2);
pfun = @(t) struct('alpha', al*(1 + 0.3*t.r(1,:).^2 + 0.2*t.r(1,:).*t.r(2,:)), ...
                   'dalpha', al*[0.6*t.r(1,:) + 0.2*t.r(2,:); 0.2*t.r(1,:)], ...
                   'd2alpha', al*[0.6 + 0*t.r(1,:); 0.2 + 0*t.r(1,:); 0*t.r(1,:)]);
K = kernel3d.iceflex('flex_varalpha', al, ga, 'nu', nu, 'pfun', pfun, ...
                     'rmax', rmax);
[gs, ~, h, d3] = surfwave.rootsum(rts, ejs.*rts.^2/2, xs, xt, 1e-8);
gp = surfwave.rootsum(rts, ejs.*rts/4, xs, xt, 1e-8);
p = pfun(t);
a = p.alpha(:);  da = p.dalpha.';  d2a = p.d2alpha.';
lap = h(:,:,1) + h(:,:,3);
glx = d3(:,:,1) + d3(:,:,3);  gly = d3(:,:,2) + d3(:,:,4);
B = (a - al).*(gp - ga/2*gs) + al/2*(2*(da(:,1).*glx + da(:,2).*gly) ...
    + (d2a(:,1) + d2a(:,3)).*lap - (1-nu)*(d2a(:,1).*h(:,:,3) ...
    - 2*d2a(:,2).*h(:,:,2) + d2a(:,3).*h(:,:,1)));
check('flex_varalpha eval', norm(K.eval(s, t) - B, 'fro')/norm(B, 'fro'), 1e-10);
check('flex_varalpha diag', norm(K.diag(t) - a)/norm(a), 1e-15);

% G_S = g (delta + G_phi): the delta in diag for 's' only
K = kernel3d.gravity('s', g, [0.3 0.5], 'rmax', rmax);
check('grav s G_S diag', norm(K.diag(t) - 0.3*g)/(0.3*g), 1e-15);
assert(isempty(kernel3d.gravity('s', g, 'gphi', 'rmax', rmax).diag) && ...
       isempty(kernel3d.gravity('sp', g, 'gs', 'rmax', rmax).diag), ...
       'grav: G_S delta only for type s');
end


function test_fgw_quadrature(cp, fp, g, d, fgws)
% SURFWAVE.FGW_GREEN against (1/2 pi) int (F - F_deep) J0(k r) k dk on a
% contour below the real poles, inside and beyond the table (R1)
[xg, wg] = lgwt(40);
cs = {[cp 3], [fp 5], [2/g 0 1]};
for q = 1:3
    a = cs{q}(1);  ga = cs{q}(2);  n = cs{q}(3);  S = fgws{q};
    r = [0.2 2 15 40 60 120]/S.kinf;
    P  = @(k) (a*k.^(n-1) + ga).*k;
    F1 = @(k) k.*tanh(k*d)./(P(k).*tanh(k*d) - 1) - k./(P(k) - 1);
    F2 = @(k) 1./(P(k).*tanh(k*d) - 1) - 1./(P(k) - 1);
    dl = min(0.02*S.kinf, 0.3/max(r));
    ka = 0.5*S.kinf;  kb = 1.5*S.k0*S.kinf + 1;
    bp = [linspace(0, ka, 200), linspace(ka, kb, 4000) - 1i*dl, ...
          linspace(kb, kb + 40/d, 8000)];
    c = (bp(1:end-1) + bp(2:end))/2;  h = (bp(2:end) - bp(1:end-1))/2;
    z = c + h.*xg;  w = h.*wg;  z = z(:);  w = w(:);
    J = besselj(0, z*r);
    G = [sum(w.*2.*F1(z).*z.*J, 1); sum(w.*F2(z).*z.*J, 1)]/(2*pi);
    t = [r; 0*r];
    B = [surfwave.fgw_green(S, 'gs', [0;0], t).'; ...
         surfwave.fgw_green(S, 'gphi', [0;0], t).'];
    check(sprintf('fgw n %d quadrature', n), ...
          max(abs(B(:) - G(:)))/max(abs(G(:))), 1e-10);
end
end

function [x, w] = lgwt(n)
% Gauss-Legendre nodes and weights on [-1, 1]
b = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));  w = 2*V(1, i).'.^2;
end
