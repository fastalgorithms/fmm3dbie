function S = fgw_setup(a, gam, d, n, opts)
%SURFWAVE.FGW_SETUP   Precompute the finite depth corrections for
%   SURFWAVE.FGW_GREEN.
%
%   S = SURFWAVE.FGW_SETUP(a,gam,d,n) for the free surface condition of
%   SURFWAVE.KERNEL and depth d. SURFWAVE.FGW_SETUP(a,gam,d,n,opts) sets
%   Ncheb (24), R1 (42) and verbose (false). S = [] for d = inf.
%
% See also SURFWAVE.FGW_GREEN, SURFWAVE.KERNEL

S = [];
if ( isinf(d) ), return; end
if nargin < 5, opts = struct(); end
Ncheb = getopt(opts, 'Ncheb', 24);
R1    = getopt(opts, 'R1', 42);
verbose = getopt(opts, 'verbose', false);
validateattributes_(a, gam, d);

% scaling
kinf = deep_root(a, gam, n);              % a k^n + gam k = 1
ap = a*kinf^n;  gp = gam*kinf;  dp = d*kinf;
ap = ap/(ap+gp); gp = 1 - ap;             % enforce ap + gp = 1 exactly
S = struct();
S.a = a; S.gam = gam; S.d = d; S.n = n; S.kinf = kinf;
S.ap = ap; S.gp = gp; S.dp = dp;
S.Ncheb = Ncheb; S.R1 = R1;

% roots (internal units)
Pf   = @(k) (ap*k.^(n-1) + gp).*k;
Pm1  = @(k) (k - 1).*(ap*polyval(ones(1, n), k) + gp);         % P(k) - 1 without cancellation near k = 1
dPf  = @(k) n*ap*k.^(n-1) + gp;
Ef   = @(k) exp(-2*k*dp);
Dt   = @(k) Pm1(k) - Ef(k).*(Pf(k) + 1);                       % = Delta*(1+E), entire
dDt  = @(k) dPf(k).*(1 - Ef(k)) + 2*dp*Ef(k).*(Pf(k) + 1);
S.k0 = real_root_k0(Dt, dDt);
qs = roots(dpoly(ap, gp, n));
for it = 1:5, qs = qs - Pm1(qs)./dPf(qs); end                   % polish
[~, i1] = min(abs(qs - 1)); qs(i1) = 1;                        % the deep-water real root is exactly 1
q1 = qs(real(qs) > 0 & imag(qs) > 0);                         % first quadrant (none for n = 1, 3)
S.qs = qs(:); S.q1 = q1(:);
S.k1 = complex_roots(Dt, dDt, q1, ap, gp, dp, n, S.k0);
S.kappa = imag_roots(ap, gp, dp, n, R1);

% residues, finite depth
resF = @(p, j) p .* (j==1) .* (p.*(1-Ef(p))./(1+Ef(p))) .* (1+Ef(p))./dDt(p) + ...
              p .* (j==2) .* (1+Ef(p))./dDt(p);
S.rho0 = [resF(S.k0, 1), resF(S.k0, 2)];
S.rho1 = [resF(S.k1, 1), resF(S.k1, 2)];
S.rhon = [real(resF(1i*S.kappa, 1)), real(resF(1i*S.kappa, 2))];   % evanescent residues (real)
% residues, deep water
resD = @(q, j) q .* ((j==1).*q + (j==2)) ./ dPf(q);
S.Ainf = [resD(1, 1), resD(1, 2)];                             % at kinf = 1
S.Aq1  = [resD(q1, 1), resD(q1, 2)];

% Bessel derivative matrices: d^n/dx^n C0(x) = sum_m DJ(n+1,m+1) C_m(x)
S.DJ = [1 0 0 0 0; 0 -1 0 0 0; -1/2 0 1/2 0 0; 0 3/4 0 -1/4 0; 3/8 0 -1/2 0 1/8];   % J, Y, H
S.DK = [1 0 0 0 0; 0 -1 0 0 0;  1/2 0 1/2 0 0; 0 -3/4 0 -1/4 0; 3/8 0  1/2 0 1/8];  % K

% deep water tail series, from the imaginary axis, Pt(y) = -i P(iy)
sg = (-1)^((n-1)/2);  h = (n+1)/2;
M = 80; c = zeros(1, M+1); c(1) = 1;
for m = 1:M
    v = -gp^2*c(m);
    if m >= h, v = v - 2*sg*ap*gp*c(m-h+1); end
    if m >= n, v = v - ap^2*c(m-n+1); end
    c(m+1) = v;
end
g1c = c;  g2c = gp*c + sg*ap*[zeros(1, h-1) c(1:end-h+1)];
mm = 0:M;
fac = 2.^(2*mm+1) .* gamma(mm + 1.5).^2 / pi^2;               % int_0^inf u^{2m+2} K0(u) du / pi^2
S.tailc = [g1c.*fac; g2c.*fac];                                 % tau_j ~ sum_m tailc(j,m+1) r^{-2m-3}

% tables
if verbose, fprintf('fgw_setup: kinf=%.6g  k0''=%.12g\n', kinf, S.k0); end
t0 = tic;
S.difftab = build_diff_table(S, Pf, Pm1, Ef, Dt, Ncheb, R1, verbose);
if verbose, fprintf('fgw_setup: difference-kernel table built in %.1f s (%d panels)\n', toc(t0), numel(S.difftab.edges)-1); end
end

function c = dpoly(a, b, m)
% coefficients of a x^m + b x - 1
c = zeros(1, m+1);
c(1) = c(1) + a;  c(m) = c(m) + b;  c(m+1) = -1;
end

function v = getopt(o, name, default)
if isfield(o, name), v = o.(name); else, v = default; end
end

function validateattributes_(a, gam, d)
if ~(isscalar(a) && isscalar(gam) && isscalar(d) && isreal(gam) && a > 0 && d > 0 && isfinite(gam))
    error('fgw_setup:args', 'a and d must be positive scalars and gam a real scalar');
end
end

function k = deep_root(a, gam, n)
% positive root of a k^n + gam k - 1: bisection + Newton
lo = 0; hi = 1;
while a*hi^n + gam*hi < 1, hi = 2*hi; end
for it = 1:200
    k = (lo+hi)/2;
    if a*k^n + gam*k < 1, lo = k; else, hi = k; end
    if hi - lo < 1e-15*hi, break; end
end
for it = 1:4, k = k - (a*k^n + gam*k - 1)/(n*a*k^(n-1) + gam); end
end

function k0 = real_root_k0(Dt, dDt)
lo = 1; hi = 1.5;
while Dt(hi) < 0, hi = 2*hi; end
for it = 1:200
    k0 = (lo+hi)/2;
    if Dt(k0) < 0, lo = k0; else, hi = k0; end
    if hi - lo < 4*eps*hi, break; end
end
for it = 1:3
    step = Dt(k0)/dDt(k0);
    if isfinite(step), k0 = k0 - step; end
end
k0 = max(k0, 1);
end

function k1 = complex_roots(Dt, dDt, q1, ap, gp, dp, n, k0)
% first-quadrant complex roots of Dt (none or one), counted by the argument
% principle, found by Newton from several seeds
u = roots(dpoly(ap*dp, gp*dp, (n+1)/2));                      % shallow water: k^2 = u
Kb = 3*max([abs(q1(:)); sqrt(abs(u)); k0]) + 2;
nw = round(winding_number(Dt, 0.05, Kb, 0.05, Kb));
k1 = zeros(0, 1);
if nw == 0, return; end
if nw > 1
    error('fgw_setup:roots', '%d first-quadrant roots (at most one supported)', nw);
end
seeds = [q1(:); sqrt(u(imag(u) > 0)); 1.2*exp(1i*pi/3); 2*exp(1i*pi/4); 0.7*exp(1i*1.2)];
for s = 1:numel(seeds)
    z = seeds(s);
    for it = 1:60
        dz = Dt(z)/dDt(z); z = z - dz;
        if abs(dz) < 1e-15*abs(z), break; end
    end
    if abs(dz) < 1e-15*abs(z) && real(z) > 1e-8 && imag(z) > 1e-8
        k1 = z;
        return
    end
end
error('fgw_setup:roots', 'no first-quadrant complex root found');
end

function nw = winding_number(f, x0, x1, y0, y1)
N = 6000;
t = (0:N-1)'/N;
z = [x0 + (x1-x0)*t + 1i*y0; x1 + 1i*(y0 + (y1-y0)*t); x1 - (x1-x0)*t + 1i*y1; x0 + 1i*(y1 - (y1-y0)*t)];
v = f(z);
ph = angle(v([2:end 1])./v);
if max(abs(ph)) > pi/2
    N2 = 60000; t = (0:N2-1)'/N2;
    z = [x0 + (x1-x0)*t + 1i*y0; x1 + 1i*(y0 + (y1-y0)*t); x1 - (x1-x0)*t + 1i*y1; x0 + 1i*(y1 - (y1-y0)*t)];
    v = f(z); ph = angle(v([2:end 1])./v);
end
nw = sum(ph)/(2*pi);
end

function kap = imag_roots(ap, gp, dp, n, R1)
% imaginary roots i*kappa of Delta:  Pt(kappa) sin(kappa d) + cos(kappa d) = 0,
% Pt = (sg ap kappa^(n-1) + gp) kappa, from sign changes on a fine grid
sg = (-1)^((n-1)/2);
h = @(x) (sg*ap*x.^(n-1) + gp).*x.*sin(x*dp) + cos(x*dp);
xmax = max(50/R1, 3.5*pi/dp);
x = linspace(0, xmax, ceil(64*xmax*dp/pi) + 2);
hx = h(x);
ic = find(sign(hx(1:end-1)) ~= sign(hx(2:end)));
kap = zeros(numel(ic), 1);
for q = 1:numel(ic)
    lo = x(ic(q)); hi = x(ic(q)+1); flo = hx(ic(q));
    for it = 1:200
        xm = (lo+hi)/2; fx = h(xm);
        if (fx > 0) == (flo > 0), lo = xm; flo = fx; else, hi = xm; end
        if hi - lo < 4*eps*hi, break; end
    end
    kap(q) = (lo+hi)/2;
end
end

function tab = build_diff_table(S, Pf, Pm1, Ef, Dt, N, R1, verbose)
dp = S.dp;
Dk = @(k) 2*Ef(k) ./ (Dt(k).*Pm1(k));                             % common factor: D1 = k*Dk, D2 = P*Dk
% panels in r, widths growing with r
thw = min([angle([S.k1; S.q1]) - 10*pi/180; 40*pi/180]);
w0 = 0.5*min(2*dp, 1);
edges = 0;
while edges(end) < R1
    w = min(2, max(w0, 0.5*(2*dp + edges(end)*tan(thw))));
    edges(end+1) = min(edges(end) + w, R1); %#ok<AGROW>
end
np = numel(edges) - 1;
xc = cos((2*(1:N)'-1)*pi/(2*N));                                 % Chebyshev nodes (1st kind), on [-1,1]
Tm = cos((0:N-1)'*acos(xc'));                                     % N x N: T_k(x_i)
C2 = (2/N)*Tm; C2(1,:) = C2(1,:)/2;                               % coefficient matrix
coef = zeros(N, np, 20);
for p = 1:np
    rr = (edges(p)+edges(p+1))/2 + (edges(p+1)-edges(p))/2*xc;   % N nodes in this panel
    [zD, wD] = dip_contour(S, edges(p+1));                        % contour adapted to the panel's largest r
    DD = Dk(zD); D1 = wD.*DD.*zD; D2 = wD.*DD.*Pf(zD);
    vals = zeros(N, 20);
    for i = 1:N
        vals(i, :) = diff_quad_dip(S, rr(i), zD, D1, D2);
    end
    for q = 1:20, coef(:, p, q) = C2*vals(:, q); end
    if verbose && mod(p, 10) == 0, fprintf('   panel %d/%d\n', p, np); end
end
tab = struct('edges', edges, 'coef', coef, 'N', N, 'coefT', permute(coef, [2 1 3]));
end

function [z, w] = dip_contour(S, rmax)
% contour along the real axis, dipping by delta below the real poles
dp = S.dp; k0 = S.k0;
delta = min([0.3; 0.4*imag([S.k1; S.q1]); 0.5*S.kappa(1); 2/max(rmax, 1e-3)]);
ka = 1 - 2*delta; kb = k0 + 2*delta;
Kmax = kb + 2 + 22/dp;
wosc = 0.8*2*pi/max(rmax, 1e-3);                                  % one oscillation period of J0(k r)
[xg, wg] = gauss_legendre(16);
z = []; w = [];
% real [0, ka]: graded near 0 (scale 1/(2d)), then limited by the oscillation period
bp = 0; t = min([0.25/dp, ka/8, wosc]);
while t < ka*(1-1e-12), bp(end+1) = t; t = min(2*t, t + wosc); end %#ok<AGROW>
bp(end+1) = ka;
[z, w] = add_segments(z, w, bp, xg, wg);
[z, w] = add_segments(z, w, [ka, ka - 1i*delta], xg, wg);         % down
nh = ceil((kb - ka)/min(delta, wosc));                            % across, panels <= delta (pole distance)
bp = ka + (0:nh)/nh*(kb - ka) - 1i*delta;
[z, w] = add_segments(z, w, bp, xg, wg);
[z, w] = add_segments(z, w, [kb - 1i*delta, kb], xg, wg);         % up
wid = min([1, 1/dp, wosc]);                                       % real [kb, Kmax]
nr = ceil((Kmax - kb)/wid);
bp = kb + (0:nr)/nr*(Kmax - kb);
[z, w] = add_segments(z, w, bp, xg, wg);
end

function [z, w] = add_segments(z, w, bp, xg, wg)
for s = 1:numel(bp)-1
    c = (bp(s)+bp(s+1))/2; h = (bp(s+1)-bp(s))/2;
    z = [z; c + h*xg]; w = [w; h*wg]; %#ok<AGROW>
end
end

function vals = diff_quad_dip(S, r, z, D1, D2)
% the Hankel transforms on the contour, radial (cols 1-10) and
% ((1/r) d/dr)^m (cols 11-20) derivatives
x = z*r;
J = bessel_orders(@besselj, x);                                   % M x 5 (complex): J_0..J_4
Jn = J*S.DJ.';                                                    % J0^{(n)}(kr), n = 0..4
if r == 0
    Jt = ones(numel(z), 1)*(1./(2.^(0:4).*factorial(0:4)));       % lim J_m(x)/x^m
else
    Jt = J./(x.^(0:4));                                           % J_m(x)/x^m
end
vals = zeros(1, 20);
for n = 0:4
    kp = z.^(n+1);
    vals(n+1)    = sum(D1.*kp.*Jn(:, n+1));
    vals(5+n+1)  = sum(D2.*kp.*Jn(:, n+1));
    kq = (-1)^n*z.^(2*n+1);
    vals(10+n+1) = sum(D1.*kq.*Jt(:, n+1));
    vals(15+n+1) = sum(D2.*kq.*Jt(:, n+1));
end
vals = vals/(2*pi);
end

function [x, w] = gauss_legendre(n)
b = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
J = diag(b, 1) + diag(b, -1);
[V, D] = eig(J);
[x, ix] = sort(diag(D));
w = 2*V(1, ix).'.^2;
x = x(:); w = w(:);
end

function B = bessel_orders(fun, x)
% B(:, m+1) = fun(m, x), m = 0..4
x = x(:);
B = zeros(numel(x), 5);
for m = 0:4
    B(:, m+1) = fun(m, x);
end
end
