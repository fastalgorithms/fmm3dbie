% Test the 2D kernels (lap2d, helm2d, flex2d) on a flat surfer.
%
% Every kernel is checked against the analytic volume potential on the
% unit disk D of the density
%
%   f(y) = J_n(kap r) e^{i n theta},    n = 2, kap = 1.3.
%
% Inside D the potential u = V[f] is a particular solution plus a regular
% homogeneous one, with the constants fixed by matching to the decaying
% exterior solution across r = 1:
%
%   Laplace     (G = -log r/(2 pi)):   u = f/kap^2 + C z^n
%   Helmholtz   (G = i/4 H_0(k r)):    u = -f/(k^2 - kap^2) + A J_n(k r) e^{i n theta}
%   flexural    (G = (G_k1 - G_k2)/(k1^2 - k2^2)):  the same difference of
%               the Helmholtz (or Laplace, if k2 = 0) potentials
%   biharmonic  (G = r^2 log r/(8 pi)): u = f/kap^4 + a z^n + b z^(n+1) conj(z)
%
% with z = x + i y. Derivatives of u follow from the ladder operators
% d/dx +- i d/dy, see DPM below. The Lippmann-Schwinger test uses a
% constant contrast on D, for which the solution is a single Bessel mode.
%
% Boundary targets are tagged with patch_id / uvs_targ.

test_lap2d_volume();
test_helm2d_volume();
test_flex2d_volume();
test_helm2d_lippmann_schwinger();
test_fmm();


function [S, tbd] = disk_and_boundary(norder)
S = geometries.disk([], [], [1,1,1], norder);
th = linspace(0, 2*pi, 41); th = th(1:end-1);
tbd = []; tbd.r = [cos(th); sin(th); zeros(size(th))];
tbd.n = [cos(th); sin(th); zeros(size(th))];
tbd = add_patch_ids(S, tbd);
end

function t = add_patch_ids(S, t)
% on-surface patch ids and uvs for boundary targets
[~, pid, uvs, ~, flags] = get_closest_pts(S, t);
assert(all(flags == 0), 'closest point search failed');
t.patch_id = pid(:);
t.uvs_targ = uvs;
end

function [n, kap] = density_mode()
n = 2; kap = 1.3;
end

function f = density(x, y)
% f = J_n(kap r) e^{i n theta} at the points (x,y), as a column
[n, kap] = density_mode();
z = x(:) + 1i*y(:);
f = besselj(n, kap*abs(z)).*exp(1i*n*angle(z));
end

function check(name, val, ex, tol)
err = norm(val(:) - ex(:), inf) / norm(ex(:), inf);
fprintf('%-40s rel err %5.2e\n', name, err);
assert(err < tol, '%s: error %5.2e exceeds %5.2e', name, err, tol);
end


function test_lap2d_volume()

tol = 1e-5;
eps = 1e-8;
[S, tbd] = disk_and_boundary(8);
[n, kap] = density_mode();

% real kernel, so use the real part of the density and of the potential
f = real(density(S.r(1,:), S.r(2,:)));
T = lap_terms(n, kap);
xb = tbd.r(1,:).'; yb = tbd.r(2,:).';
nb = tbd.n(1:2,:);

% on-surface (volume-to-volume)
A = surfermat(S, kernel3d('lap2d', 's'), eps);
uex = real(dxy(T, 0, 0, S.r(1,:).', S.r(2,:).'));
check('lap2d  s  (v2v)', A*f, uex, tol);

% normal derivative on the boundary curve
val = surferkerneval(S, kernel3d('lap2d', 'sp'), f, tbd, eps);
check('lap2d  sp (v2b)', val, real(dirder(T, xb, yb, {nb})), tol);

% gradient on the boundary curve, d/dx and d/dy interleaved
val = surferkerneval(S, kernel3d('lap2d', 'sg'), f, tbd, eps);
gex = [dxy(T, 1, 0, xb, yb), dxy(T, 0, 1, xb, yb)].';
check('lap2d  sg (v2b)', val, real(gex(:)), tol);

end


function test_helm2d_volume()

tol = 1e-5;
eps = 1e-8;
zk  = 2.3;
[S, tbd] = disk_and_boundary(8);
[n, kap] = density_mode();

f = density(S.r(1,:), S.r(2,:));
T = helm_terms(zk, n, kap);
xb = tbd.r(1,:).'; yb = tbd.r(2,:).';
nb = tbd.n(1:2,:);

A = surfermat(S, kernel3d('helm2d', 's', zk), eps);
uex = dxy(T, 0, 0, S.r(1,:).', S.r(2,:).');
check('helm2d s  (v2v)', A*f, uex, tol);

ub  = dxy(T, 0, 0, xb, yb);
dub = dirder(T, xb, yb, {nb});

val = surferkerneval(S, kernel3d('helm2d', 'sp', zk), f, tbd, eps);
check('helm2d sp (v2b)', val, dub, tol);

val = surferkerneval(S, kernel3d('helm2d', 'sg', zk), f, tbd, eps);
gex = [dxy(T, 1, 0, xb, yb), dxy(T, 0, 1, xb, yb)].';
check('helm2d sg (v2b)', val, gex(:), tol);

coefs = [1; 2];
val = surferkerneval(S, kernel3d('helm2d', 's2trans', zk, coefs), f, tbd, eps);
tex = [coefs(1)*ub, coefs(2)*dub].';
check('helm2d s2trans (v2b)', val, tex(:), tol);

end


function test_flex2d_volume()

tol = 1e-5;
eps = 1e-8;
nu  = 0.3;
[n, kap] = density_mode();

% flexural (scalar), general pair, pair with a zero, biharmonic
zks = {1.7, [1.7, 0.5i], [1.7, 0], 0};

[S, tbd] = disk_and_boundary(8);
f = density(S.r(1,:), S.r(2,:));

% the plate kernels also need the tangent d and its derivative d2
xb = tbd.r(1,:).'; yb = tbd.r(2,:).';
tbd.d  = [-yb.'; xb.'; zeros(size(xb.'))];
tbd.d2 = [-xb.'; -yb.'; zeros(size(xb.'))];
nb = tbd.n(1:2,:);
tb = tbd.d(1:2,:)./vecnorm(tbd.d(1:2,:));
curv = (tbd.d(1,:).*tbd.d2(2,:) - tbd.d2(1,:).*tbd.d(2,:)).' ...
    ./ vecnorm(tbd.d(1:2,:)).'.^3;

types = {'s', 'clamped_plate_bcs', 'supported_plate_bcs', 'free_plate_bcs'};

for iz = 1:numel(zks)
    zk = zks{iz};
    T  = flex_terms(zk, n, kap);
    zstr = mat2str(zk);

    % on-surface (volume-to-volume)
    A = surfermat(S, kernel3d('flex2d', 's', zk, nu), eps);
    uex = dxy(T, 0, 0, S.r(1,:).', S.r(2,:).');
    check(['flex2d s (v2v), zk = ' zstr], A*f, uex, tol);

    % boundary data from the derivatives of u on r = 1
    u    = dxy(T, 0, 0, xb, yb);
    un   = dirder(T, xb, yb, {nb});
    unn  = dirder(T, xb, yb, {nb, nb});
    utt  = dirder(T, xb, yb, {tb, tb});
    unnn = dirder(T, xb, yb, {nb, nb, nb});
    uttn = dirder(T, xb, yb, {tb, tb, nb});
    mom  = unn + nu*utt;
    shr  = unnn + (2-nu)*uttn + (1-nu)*curv.*(utt - unn);

    exs = {u, [u, un].', [u, mom].', [mom, shr].'};

    for it = 1:numel(types)
        K = kernel3d('flex2d', types{it}, zk, nu);
        val = surferkerneval(S, K, f, tbd, eps);
        if it == 1
            check(['flex2d s (v2b), zk = ' zstr], val, exs{it}, tol);
        else
            % check the two conditions separately, they scale differently
            ex = exs{it};
            check(['flex2d ' types{it} ' (1), zk = ' zstr], ...
                val(1:2:end), ex(1,:), tol);
            check(['flex2d ' types{it} ' (2), zk = ' zstr], ...
                val(2:2:end), ex(2,:), tol);
        end
    end
end

end


function test_helm2d_lippmann_schwinger()
% Adjoint Lippmann-Schwinger solve, as in helm2d_lippmann_schwinger_demo,
%
%   Delta u + zk^2 (1 + q) u = 0,   u = u_inc + V[sigma],
%   sigma - zk^2 q V[sigma] = zk^2 q u_inc,
%
% for a constant contrast q on the unit disk and the incident field
% u_inc = J_n(zk r) e^{i n theta}. Then u = a J_n(k1 r) e^{i n theta}
% inside, with k1 = zk sqrt(1+q), u_scat = b H_n(zk r) e^{i n theta}
% outside, and sigma = zk^2 q u.

tol = 1e-6;
eps = 1e-10;
zk  = 2.3;
q   = 0.5;
n   = 2;

[S, tbd] = disk_and_boundary(8);
x = S.r(1,:).'; y = S.r(2,:).';
r = sqrt(x.^2 + y.^2); ephi = exp(1i*n*atan2(y, x));
uinc = besselj(n, zk*r).*ephi;

K = kernel3d('helm2d', 's', zk);
V = surfermat(S, K, eps);
sigma = (eye(S.npts) - zk^2*q*V) \ (zk^2*q*uinc);

% transmission problem: u and du/dr continuous across r = 1
k1 = zk*sqrt(1 + q);
ab = [besselj(n, k1), -besselh(n, 1, zk); ...
      k1*dbesselj(n, k1), -zk*dbesselh(n, zk)] \ ...
     [besselj(n, zk); zk*dbesselj(n, zk)];

sigex = zk^2*q*ab(1)*besselj(n, k1*r).*ephi;
check('helm2d adjoint LS, sigma', sigma, sigex, tol);

% scattered field on the boundary curve
us = surferkerneval(S, K, sigma, tbd, eps);
usex = ab(2)*besselh(n, 1, zk)*exp(1i*n*atan2(tbd.r(2,:), tbd.r(1,:))).';
check('helm2d adjoint LS, u_scat on r = 1', us, usex, tol);

end


function test_fmm()
% fmm vs direct sum, skipped if fmm2d is not on the path

if ( exist(['fmm2d.' mexext], 'file') ~= 3 )
    fprintf('fmm2d not found, skipping fmm tests\n');
    return
end

tol = 1e-8;
S = geometries.disk([], [], [1,1,1], 4);
sigma = randn(S.npts, 1) + 1i*randn(S.npts, 1);
th = linspace(0, 2*pi, 31); th = th(1:end-1);
t = []; t.r = 1.5*[cos(th); sin(th); zeros(size(th))];
t.n = [cos(th); sin(th); zeros(size(th))];
t.d = [-sin(th); cos(th); zeros(size(th))];

kerns = {kernel3d('lap2d', 's'), kernel3d('lap2d', 'sp'), ...
         kernel3d('lap2d', 'sg'), ...
         kernel3d('helm2d', 's', 2.3), kernel3d('helm2d', 'sp', 2.3), ...
         kernel3d('helm2d', 'sg', 2.3), ...
         kernel3d('helm2d', 's2trans', 2.3, [1; 2]), ...
         kernel3d('flex2d', 's', 1.7), ...
         kernel3d('flex2d', 'clamped_plate_bcs', 1.7), ...
         kernel3d('flex2d', 'supported_plate_bcs', [1.7, 0.5i], 0.3), ...
         kernel3d('flex2d', 'clamped_plate_bcs', [1.7, 0]), ...
         kernel3d('flex2d', 's', 0), ...
         kernel3d('flex2d', 'clamped_plate_bcs', 0), ...
         kernel3d('flex2d', 'supported_plate_bcs', 0, 0.3)};

for i = 1:numel(kerns)
    K = kerns{i};
    assert(~isempty(K.fmm), 'missing fmm for %s %s', K.name, K.type);
    pf = K.fmm(1e-12, S, t, sigma);
    pd = K.eval(S, t)*sigma;
    err = norm(pf - pd) / norm(pd);
    fprintf('%-7s %-20s fmm err %5.2e\n', K.name, K.type, err);
    assert(err < tol);
end

assert(isempty(kernel3d('flex2d', 'free_plate_bcs', 1.7, 0.3).fmm));
assert(isempty(kernel3d('flex2d', 'free_plate_bcs', 0, 0.3).fmm));

end


%
%  Analytic volume potentials of f = J_n(kap r) e^{i n theta} on the unit
%  disk, valid for r <= 1.
%
%  A potential is stored as a list of terms T, one per row:
%     [1, coef, c, m]   coef * J_m(c r) e^{i m theta}
%     [0, coef, p, q]   coef * z^p conj(z)^q,   z = x + i y
%

function d = dbesselj(n, x)
d = (besselj(n-1, x) - besselj(n+1, x))/2;
end

function d = dbesselh(n, x)
d = (besselh(n-1, 1, x) - besselh(n+1, 1, x))/2;
end

function T = lap_terms(n, kap)
% -Delta u = f. Outside u = B r^(-n) e^{i n theta}, n >= 1.
P  = besselj(n, kap)/kap^2;
Pp = dbesselj(n, kap)/kap;
C  = -(Pp + n*P)/(2*n);
T  = [1, 1/kap^2, kap, n; ...
      0, C, n, 0];
end

function T = helm_terms(k, n, kap)
% (Delta + k^2) u = -f. Outside u = B H_n(k r) e^{i n theta}.
p  = -1/(k^2 - kap^2);
P  = p*besselj(n, kap);
Pp = p*kap*dbesselj(n, kap);
A  = (1i*pi/2)*(k*dbesselh(n, k)*P - Pp*besselh(n, 1, k));
T  = [1, p, kap, n; ...
      1, A, k, n];
end

function T = bh_terms(n, kap)
% Delta^2 u = f. Outside u = (B1 r^(-n) + B2 r^(2-n)) e^{i n theta},
% n >= 2; u, du/dr, Delta u and d/dr Delta u are continuous at r = 1.
Q  = -besselj(n, kap)/kap^2;
Qp = -dbesselj(n, kap)/kap;
b  = -(n*Q + Qp)/(8*n*(n+1));
B2 = (Q + 4*(n+1)*b)/(4*(1-n));
P  = besselj(n, kap)/kap^4;
Pp = dbesselj(n, kap)/kap^3;
a  = (2*B2 - n*P - Pp - (2*n+2)*b)/(2*n);
T  = [1, 1/kap^4, kap, n; ...
      0, a, n, 0; ...
      0, b, n+1, 1];
end

function T = flex_terms(zk, n, kap)
% wavenumbers chosen exactly as in FLEX2D.GREEN
zk = zk(:).';
if isscalar(zk)
    if abs(zk) < 1e-6
        T = bh_terms(n, kap);
        return
    end
    zk1 = zk; zk2 = 1i*zk;
elseif any(abs(zk) < 1e-6)
    zk1 = zk(abs(zk) >= 1e-6); zk2 = 0;
else
    zk1 = zk(1); zk2 = zk(2);
end
T1 = helm_terms(zk1, n, kap);
if zk2 == 0
    T2 = lap_terms(n, kap);
else
    T2 = helm_terms(zk2, n, kap);
end
T2(:,2) = -T2(:,2);
T = [T1; T2];
T(:,2) = T(:,2)/(zk1^2 - zk2^2);
end

function val = dpm(T, a, b, x, y)
%DPM  (d/dx + i d/dy)^a (d/dx - i d/dy)^b of the potential T, using
%   (d/dx +- i d/dy) J_m(c r) e^{i m theta} = -+ c J_{m+-1}(c r) e^{i (m+-1) theta}
%   (d/dx + i d/dy) = 2 d/dconj(z),  (d/dx - i d/dy) = 2 d/dz
z = x + 1i*y; r = abs(z); ephi = exp(1i*angle(z));
val = zeros(size(z));
for j = 1:size(T,1)
    coef = T(j,2);
    if T(j,1) == 1
        c = T(j,3); m = real(T(j,4)) + a - b;
        val = val + coef*(-c)^a*c^b*besselj(m, c*r).*ephi.^m;
    else
        p = real(T(j,3)); q = real(T(j,4));
        if ( b > p || a > q ), continue, end
        fac = 2^(a+b)*prod(p-b+1:p)*prod(q-a+1:q);
        val = val + coef*fac*z.^(p-b).*conj(z).^(q-a);
    end
end
end

function val = dxy(T, a, b, x, y)
%DXY  d^a/dx^a d^b/dy^b of the potential T at the points (x,y), using
%   d/dx = (D+ + D-)/2,  d/dy = (D+ - D-)/(2i),  D+- = d/dx +- i d/dy
P = 1;
for j = 1:a, P = conv(P, [1 1]);  end
for j = 1:b, P = conv(P, [1 -1]); end
m = a + b;
val = zeros(size(x));
for j = 0:m
    val = val + P(j+1)*dpm(T, m-j, j, x, y);
end
val = val/(2^a*(2i)^b);
end

function val = dirder(T, x, y, vecs)
%DIRDER  directional derivative of the potential T along the directions
% vecs{1}, ..., vecs{m}, each a (2,npts) array
m = numel(vecs);
val = zeros(size(x));
for ind = 0:2^m-1
    bits = bitget(ind, 1:m);       % 0: d/dx, 1: d/dy
    w = ones(size(x));
    for j = 1:m
        w = w.*vecs{j}(bits(j)+1,:).';
    end
    val = val + w.*dxy(T, m-sum(bits), sum(bits), x, y);
end
end
