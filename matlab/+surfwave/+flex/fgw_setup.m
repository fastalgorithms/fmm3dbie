function S = fgw_setup(a, gam, d, opts)
%FGW_SETUP  Precompute everything needed to evaluate the flexural-gravity (floating
%   ice sheet, finite depth) point-load kernels and their r-derivatives.
%
%   S = FGW_SETUP(a, gam, d)   with a, d > 0 and gam real (any consistent units: a ~ L^5, gam ~ L, d ~ L)
%   gam < 0 (plate inertia exceeds buoyancy, omega^2 > rho g/m) is supported as long as the root structure
%   is the same as for gam > 0 (checked below; an error is raised otherwise).
%
%   The kernels (in the Fourier variable k = |xi|, T = tanh(k d)):
%       F1(k) = k T / ((a k^4 + gam) k T - 1)      (plate deflection kernel)
%       F2(k) = 1   / ((a k^4 + gam) k T - 1)      (surface velocity-potential kernel)
%   deep-water counterparts F1inf = k/((a k^4+gam)k - 1), F2inf = 1/((a k^4+gam)k - 1),
%   and the transform convention
%       G(r) = (2 pi)^-2 \int_{R^2} F(|xi|) e^{i xi.x} dxi = (1/2pi) \int_0^inf F(k) J0(k r) k dk ,
%   with every real pole passed by the OUTGOING prescription (Im k -> 0+, i.e. contour below
%   the pole; time dependence e^{-i omega t}; far field ~ e^{i k0 r} / sqrt(r)).
%
%   FGW_EVAL(S, r, kind) then returns d^n/dr^n G_j(r), n = 0..4, j = 1,2 for
%       kind = 'diff' : transform of F_j - F_jinf   (the "difference" kernels)
%       kind = 'deep' : transform of F_jinf          (deep-water kernels)
%       kind = 'full' : transform of F_j             (= diff + deep, finite depth)
%
%   METHOD (all in internal units scaled so that the deep-water wavenumber kinf = 1):
%     * roots: real roots k0 (finite depth) and kinf=1 (deep), first-quadrant complex roots k1
%       (finite depth) and q1 (deep), imaginary roots i*kappa_n (finite depth); residues.
%     * 'diff', 0 <= r <= R1 (=42):  direct Gauss-Legendre quadrature of the Hankel integral of the
%       smooth, exponentially decaying difference kernel along a contour that dips below the two
%       real poles (depth min(0.3, 2/r), so that |J0(kr)| grows by at most e^2 on the dip); panel
%       widths follow the local scales (1/(2d) near k = 0, the dip depth near the poles, one
%       oscillation period 2pi/r on the real axis).  Values of all 5 derivatives are tabulated on
%       Chebyshev panels in r (widths growing with r inside the wedge of analyticity) and evaluated
%       by Clenshaw recurrence.
%     * 'diff', r > R1: closed form = outgoing waves at k0 and kinf, damped complex-root waves,
%       evanescent K0(kappa_n r) modes, minus the deep-water algebraic tail (asymptotic
%       series in r^-2, optimally truncated).
%     * 'deep', r <= rA (~2/|q|max): convergent power-log series of the deep-water kernel.
%     * 'deep', rA < r <= R1: outgoing + complex-root Hankel terms plus the smooth imaginary-
%       axis integral tau_j(r) tabulated on geometric Chebyshev panels.
%     * 'deep', r > R1: Hankel terms plus asymptotic tail series.
%
%   opts (optional struct): .Ncheb (24), .R1 (42), .verbose (false)

if nargin < 4, opts = struct(); end
Ncheb = getopt(opts, 'Ncheb', 24);
R1    = getopt(opts, 'R1', 42);
rAfac = getopt(opts, 'rAfac', 2);
verbose = getopt(opts, 'verbose', false);
validateattributes_(a, gam, d);

% ---------------------------------------------------------------- scaling
kinf = deep_root(a, gam);                 % a k^5 + gam k = 1
ap = a*kinf^5;  gp = gam*kinf;  dp = d*kinf;
ap = ap/(ap+gp); gp = 1 - ap;             % enforce ap + gp = 1 exactly
if gp < 0
    % gam < 0: the methods below assume (i) a single real root of the deep-water quintic and (ii) one
    % imaginary root per interval ((n-1/2)pi/d, n pi/d), i.e. (a k^4 + gam) k > 0 from k = pi/(2d) on
    km = (-gp/(5*ap))^(1/4);                  % local maximum of ap k^5 + gp k - 1 on k < 0 is at -km
    if ap*(-km)^5 + gp*(-km) - 1 >= 0
        error('fgw_setup:gamma', ['gam = %g is too negative: the deep-water dispersion polynomial has extra ' ...
            'negative real roots (not supported)'], gam);
    end
    if (-gp/ap)^(1/4) >= pi/(2*dp)
        error('fgw_setup:gamma', ['gam = %g is too negative for depth %g: (a k^4 + gam) k is still negative ' ...
            'beyond the first evanescent-root interval (not supported)'], gam, d);
    end
end
S = struct();
S.a = a; S.gam = gam; S.d = d; S.kinf = kinf;
S.ap = ap; S.gp = gp; S.dp = dp;
S.Ncheb = Ncheb; S.R1 = R1;

% ---------------------------------------------------------------- roots (internal units)
Pf   = @(k) (ap*k.^4 + gp).*k;
Pm1  = @(k) (k - 1).*(ap*(k.^4 + k.^3 + k.^2 + k + 1) + gp);   % P(k) - 1 without cancellation near k = 1 (uses ap+gp = 1)
dPf  = @(k) 5*ap*k.^4 + gp;
Ef   = @(k) exp(-2*k*dp);
Dt   = @(k) Pm1(k) - Ef(k).*(Pf(k) + 1);                       % = Delta*(1+E), entire
dDt  = @(k) dPf(k).*(1 - Ef(k)) + 2*dp*Ef(k).*(Pf(k) + 1);
S.k0 = real_root_k0(Dt, dDt);
qs = roots([ap 0 0 0 gp -1]);
for it = 1:5, qs = qs - Pm1(qs)./dPf(qs); end                   % polish
[~, i1] = min(abs(qs - 1)); qs(i1) = 1;                        % the deep-water real root is exactly 1
q1 = qs(real(qs) > 0 & imag(qs) > 0);
if numel(q1) ~= 1, error('fgw_setup:roots', 'unexpected deep-water root structure'); end
S.qs = qs(:); S.q1 = q1;
S.k1 = complex_root_k1(Dt, dDt, q1, ap, gp, dp, S.k0);
S.kappa = imag_roots(ap, gp, dp, R1);                          % kappa_n with kappa_n*R1 <= 50 (>= 3 of them)

% residues of F_j(k) k  (finite depth): p N_j (1+E)/Dt'(p);  N1 = p(1-E)/(1+E), N2 = 1
resF = @(p, j) p .* (j==1) .* (p.*(1-Ef(p))./(1+Ef(p))) .* (1+Ef(p))./dDt(p) + ...
              p .* (j==2) .* (1+Ef(p))./dDt(p);
S.rho0 = [resF(S.k0, 1), resF(S.k0, 2)];
S.rho1 = [resF(S.k1, 1), resF(S.k1, 2)];
S.rhon = [real(resF(1i*S.kappa, 1)), real(resF(1i*S.kappa, 2))];   % evanescent residues (real)
% residues of F_jinf(k) k (deep): q N_jinf / Dinf'(q), N1inf = q, N2inf = 1
resD = @(q, j) q .* ((j==1).*q + (j==2)) ./ dPf(q);
S.Ainf = [resD(1, 1), resD(1, 2)];                             % at kinf = 1
S.Aq1  = [resD(q1, 1), resD(q1, 2)];
S.Aall = [resD(S.qs, 1), resD(S.qs, 2)];                        % all five deep roots (columns j)

% Bessel derivative matrices: d^n/dx^n C0(x) = sum_m DJ(n+1,m+1) C_m(x)
S.DJ = [1 0 0 0 0; 0 -1 0 0 0; -1/2 0 1/2 0 0; 0 3/4 0 -1/4 0; 3/8 0 -1/2 0 1/8];   % J, Y, H
S.DK = [1 0 0 0 0; 0 -1 0 0 0;  1/2 0 1/2 0 0; 0 -3/4 0 -1/4 0; 3/8 0  1/2 0 1/8];  % K

% ---------------------------------------------------------------- deep-water tail series
% tau_j(r) = (1/pi^2) int_0^inf g_j(y) K0(y r) dy,  g1 = y^2/(1+P^2), g2 = P y/(1+P^2),
% 1/(1+P^2) = sum c_m s^m (s = y^2), denominator 1 + gp^2 s + 2 ap gp s^3 + ap^2 s^5
M = 80; c = zeros(1, M+1); c(1) = 1;
for m = 1:M
    v = -gp^2*c(m);
    if m >= 3, v = v - 2*ap*gp*c(m-2); end
    if m >= 5, v = v - ap^2*c(m-4); end
    c(m+1) = v;
end
g1c = c;  g2c = gp*c + ap*[0 0 c(1:end-2)];
mm = 0:M;
fac = 2.^(2*mm+1) .* gamma(mm + 1.5).^2 / pi^2;               % int_0^inf u^{2m+2} K0(u) du / pi^2
S.tailc = [g1c.*fac; g2c.*fac];                                 % tau_j ~ sum_m tailc(j,m+1) r^{-2m-3}


% ---------------------------------------------------------------- deep-water small-r series
S.rA = rAfac/max(abs(S.qs));
S.series = deep_series_coeffs(S, ap, gp);

% ---------------------------------------------------------------- tables
if verbose, fprintf('fgw_setup: kinf=%.6g  k0''=%.12g  k1''=%.6g%+.6gi  q1''=%.6g%+.6gi\n', ...
        kinf, S.k0, real(S.k1), imag(S.k1), real(q1), imag(q1)); end
t0 = tic;
S.difftab = build_diff_table(S, Pf, Pm1, Ef, Dt, Ncheb, R1, verbose);
if verbose, fprintf('fgw_setup: difference-kernel table built in %.1f s (%d panels)\n', toc(t0), numel(S.difftab.edges)-1); end
t0 = tic;
S.tautab = build_tau_table(S, ap, gp, 20, S.rA, R1);
if verbose, fprintf('fgw_setup: deep-water tau table built in %.1f s (%d panels)\n', toc(t0), numel(S.tautab.edges)-1); end
end

% ======================================================================================
function v = getopt(o, name, default)
if isfield(o, name), v = o.(name); else, v = default; end
end

function validateattributes_(a, gam, d)
if ~(isscalar(a) && isscalar(gam) && isscalar(d) && isreal(gam) && a > 0 && d > 0 && isfinite(gam))
    error('fgw_setup:args', 'a and d must be positive scalars and gam a real scalar');
end
end

function k = deep_root(a, gam)
% positive root of a k^5 + gam k - 1 (monotone): bisection + Newton
lo = 0; hi = 1;
while a*hi^5 + gam*hi < 1, hi = 2*hi; end
for it = 1:200
    k = (lo+hi)/2;
    if a*k^5 + gam*k < 1, lo = k; else, hi = k; end
    if hi - lo < 1e-15*hi, break; end
end
for it = 1:4, k = k - (a*k^5 + gam*k - 1)/(5*a*k^4 + gam); end
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

function k1 = complex_root_k1(Dt, dDt, q1, ap, gp, dp, k0)
% first-quadrant complex root of Dt; Newton from several seeds, verified by argument principle
u = roots([ap*dp 0 gp*dp -1]);  u1 = u(imag(u) > 0);           % shallow-water guess: k^2 = u1
seeds = [q1; sqrt(u1(:)); 1.2*exp(1i*pi/3); 2*exp(1i*pi/4); 0.7*exp(1i*1.2)];
best = NaN;
for s = 1:numel(seeds)
    z = seeds(s);
    ok = false;
    for it = 1:60
        dz = Dt(z)/dDt(z); z = z - dz;
        if abs(dz) < 1e-15*abs(z), ok = true; break; end
    end
    if ok && real(z) > 1e-8 && imag(z) > 1e-8
        if isnan(best), best = z; elseif abs(z - best) > 1e-8*abs(z)
            warning('fgw_setup:roots', 'two distinct first-quadrant roots found (%g%+gi, %g%+gi)', real(best), imag(best), real(z), imag(z));
        end
    end
end
if isnan(best), error('fgw_setup:roots', 'no first-quadrant complex root found'); end
k1 = best;
% argument principle on a rectangle in the open first quadrant
Kb = 3*max([abs(k1), abs(q1), k0]) + 2;
nw = winding_number(Dt, 0.05, Kb, 0.05, Kb);
if abs(nw - 1) > 0.1
    warning('fgw_setup:roots', 'winding number in the first quadrant is %g (expected 1): unusual root structure, results may be wrong', nw);
end
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

function kap = imag_roots(ap, gp, dp, R1)
% imaginary roots i*kappa of Delta:  sin(kappa d) + cos(kappa d)/Pt(kappa) = 0, Pt = (ap kappa^4+gp) kappa
Pt = @(x) (ap*x.^4 + gp).*x;
h  = @(x) sin(x*dp) + cos(x*dp)./Pt(x);
dh = @(x) dp*cos(x*dp) - dp*sin(x*dp)./Pt(x) - cos(x*dp).*(5*ap*x.^4+gp)./Pt(x).^2;
kap = [];
n = 0;
while true
    n = n + 1;
    lo = (n-0.5)*pi/dp*(1+1e-14); hi = n*pi/dp*(1-1e-14);
    flo = h(lo);
    for it = 1:200
        x = (lo+hi)/2; fx = h(x);
        if (fx > 0) == (flo > 0), lo = x; flo = fx; else, hi = x; end
        if hi - lo < 4*eps*hi, break; end
    end
    x = (lo+hi)/2;
    for it = 1:3, x = x - h(x)/dh(x); end
    kap(end+1, 1) = x; %#ok<AGROW>
    if x*R1 > 50 && n >= 3, break; end
    if n > 5000, break; end
end
end

% ======================================================================================
function ser = deep_series_coeffs(S, ap, gp)
% G_j^deep(r) = sum_k [alpha_k r^{2k} + beta_k r^{2k} log r] + sum_k sigma_k r^{2k+1}
% built from Phi(q,r) = s_q i pi J0(qr) - (pi/2) H0(qr) - (log(qr/2)+gamma_E) J0(qr) - sum_{k>=1} (-1)^{k+1} H_k (qr/2)^{2k}/(k!)^2
Kmax = 40;
gE = 0.57721566490153286060651209;
Hk = cumsum(1./(1:Kmax));                                        % harmonic numbers
qs = S.qs; A = S.Aall;                                            % 5 x 2
sq = ones(5, 1); sq(imag(qs) < 0) = -1;                           % outgoing real root -> +1
ser = struct('alpha', zeros(Kmax+1, 2), 'beta', zeros(Kmax+1, 2), 'sigma', zeros(Kmax+1, 2), 'Kmax', Kmax);
for j = 1:2
    % exact moment sums S_m = sum_i A_i q_i^m from the large-k expansion of N(k)/Dinf(k)
    Sm = zeros(1, 2*Kmax+3);                                      % index m+1
    for m = 0:2*Kmax+2
        v = 0;
        if m == 2 && j == 1, v = 1; end
        if m == 3 && j == 2, v = 1; end
        if m >= 4, v = v - gp*Sm(m-4+1); end
        if m >= 5, v = v + Sm(m-5+1); end
        Sm(m+1) = v/ap;
    end
    for k = 0:Kmax
        jk = (-1)^k/(4^k*factorial(k)^2);
        Lsum  = sum(A(:,j).*qs.^(2*k).*(log(qs/2) + gE));         % numerical (no exact cancellation)
        Ssum  = sum(sq.*A(:,j).*qs.^(2*k));
        S2k   = Sm(2*k+1);
        alpha = jk*(1i*pi*Ssum - Lsum);
        if k >= 1, alpha = alpha - (-1)^(k+1)*Hk(k)*S2k/(4^k*factorial(k)^2); end
        ser.alpha(k+1, j) = alpha/(2*pi);
        ser.beta(k+1, j)  = -jk*S2k/(2*pi);
        hk = (-1)^k/(2^(2*k+1)*gamma(k+1.5)^2);
        ser.sigma(k+1, j) = -hk*Sm(2*k+2)/4;                      % -(1/2pi)(pi/2) h_k S_{2k+1}
    end
end
% derivative tables, for both forms (f = 1: d^n/dr^n, f = 2: Theta^m = (1/r d/dr)^m):
%   output (j,n) = pe [P_A(s) + log r P_B(s)] + po P_C(s) + singular terms,  s = r^2,
% pe, po = 1 or r depending on the parity of the surviving powers (flags ser.pe/ser.po, 1 = multiply by r);
% row index = power of s + 1.  Theta r^p = p r^{p-2},  Theta (r^p log r) = r^{p-2} (p log r + 1).
ser.PA = zeros(Kmax+1, 5, 2, 2); ser.PB = ser.PA; ser.PC = ser.PA; ser.sing = cell(5, 2, 2);
ser.pe = zeros(5, 2); ser.po = zeros(5, 2);
for f = 1:2
    step = f;                                                          % exponent decrement per derivative
    for j = 1:2
        for n = 0:4
            sg = zeros(0, 3);
            for k = 0:Kmax
                [A, B, p] = dlog_(ser.alpha(k+1, j), ser.beta(k+1, j), 2*k, n, step);   % r^p (A + B log r)
                if A ~= 0 || B ~= 0
                    if p >= 0
                        ip = floor(p/2) + 1;
                        ser.PA(ip, n+1, j, f) = ser.PA(ip, n+1, j, f) + A;
                        ser.PB(ip, n+1, j, f) = ser.PB(ip, n+1, j, f) + B;
                    else
                        sg(end+1, :) = [A, B, p]; %#ok<AGROW>
                    end
                end
                [A, ~, p] = dlog_(ser.sigma(k+1, j), 0, 2*k+1, n, step);
                if A ~= 0
                    if p >= 0
                        ip = floor(p/2) + 1;
                        ser.PC(ip, n+1, j, f) = ser.PC(ip, n+1, j, f) + A;
                    else                                                       % e.g. Theta^2 r^3 = 3/r
                        sg(end+1, :) = [A, 0, p]; %#ok<AGROW>
                    end
                end
            end
            ser.sing{n+1, j, f} = sg;
        end
    end
    if f == 1, ser.pe(:, f) = mod(0:4, 2).'; ser.po(:, f) = 1 - mod(0:4, 2).';   % radial: parity of n
    else,      ser.pe(:, f) = 0;              ser.po(:, f) = 1;                   % theta: even/odd families keep parity
    end
end
end

function [A, B, p] = dlog_(A, B, p, n, step)
% n applications of d/dr (step = 1) or (1/r) d/dr (step = 2) to r^p (A + B log r)
for it = 1:n
    A = p*A + B; B = p*B; p = p - step;
end
end

% ======================================================================================
function tab = build_diff_table(S, Pf, Pm1, Ef, Dt, N, R1, verbose)
dp = S.dp;
Dk = @(k) 2*Ef(k) ./ (Dt(k).*Pm1(k));                             % common factor: D1 = k*Dk, D2 = P*Dk
% panels in r: width grows with r (analyticity wedge |Im r| < 2d + Re r tan(theta))
thw = min(min(angle(S.k1), angle(S.q1)) - 10*pi/180, 40*pi/180);
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
coef = zeros(N, np, 20);                                          % (k, panel, col), complex; col = 5*(j-1)+n+1 (d^n/dr^n), 10+5*(j-1)+m+1 (Theta^m)
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
% contour: [0,ka] real, down to ka-i*delta, across to kb-i*delta, up to kb, then [kb, Kmax];
% delta = min(0.3, 2/rmax) keeps |J0(kr)| <= ~e^2 on the dip for all r <= rmax; it must also stay above the
% lower-half-plane poles conj(k1), conj(q1) and -i kappa_1 (the dispersion function is even in k).
dp = S.dp; k0 = S.k0;
delta = min([0.3, 0.4*imag(S.k1), 0.4*imag(S.q1), 0.5*S.kappa(1), 2/max(rmax, 1e-3)]);
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
% (1/2pi) int_C D_j(k) k^{n+1} J0^{(n)}(k r) dk  (cols 1-10)  and
% Theta^m g = ((-1)^m/2pi) int_C D_j(k) k^{2m+1} [J_m(k r)/(k r)^m] dk  (cols 11-20),  Theta = (1/r) d/dr
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

% ======================================================================================
function tab = build_tau_table(S, ap, gp, N, rA, R1)
% tau_j^{(n)}(r) = (1/pi^2) int_0^inf g_j(y) y^n K0^{(n)}(y r) dy on geometric panels [rA, R1]
ratio = 1.6;
edges = rA;
while edges(end) < R1, edges(end+1) = min(edges(end)*ratio, R1); end %#ok<AGROW>
np = numel(edges) - 1;
xc = cos((2*(1:N)'-1)*pi/(2*N));
Tm = cos((0:N-1)'*acos(xc'));
C2 = (2/N)*Tm; C2(1,:) = C2(1,:)/2;
coef = zeros(N, np, 20);
[xg, wg] = gauss_legendre(16);
for p = 1:np
    rr = (edges(p)+edges(p+1))/2 + (edges(p+1)-edges(p))/2*xc;
    vals = zeros(N, 20);
    for i = 1:N
        vals(i, :) = tau_quad(S, ap, gp, rr(i), xg, wg);
    end
    for q = 1:20, coef(:, p, q) = C2*vals(:, q); end
end
tab = struct('edges', edges, 'coef', coef, 'N', N, 'coefT', permute(coef, [2 1 3]));
end

function vals = tau_quad(S, ap, gp, r, xg, wg)
% substitution y = e^s;  integrand analytic in |Im s| < ~0.3 -> panels of length 0.35
smin = -13; smax = log(42/r);
ns = ceil((smax - smin)/0.35);
bp = smin + (0:ns)/ns*(smax - smin);
s = []; w = [];
for q = 1:ns
    c = (bp(q)+bp(q+1))/2; h = (bp(q+1)-bp(q))/2;
    s = [s; c + h*xg]; w = [w; h*wg]; %#ok<AGROW>
end
y = exp(s);
P = (ap*y.^4 + gp).*y; den = 1 + P.^2;
g1 = y.^2./den; g2 = P.*y./den;
x = y*r;
K = bessel_orders(@besselk, x);                                   % M x 5: K_0..K_4
Kn = K*S.DK.';                                                    % K0^{(n)}(y r)
Kt = K./(x.^(0:4));                                               % K_m(x)/x^m  (Theta^m K0(yr) = (-1)^m y^{2m} K_m/x^m)
vals = zeros(1, 20);
for n = 0:4
    yp = y.^n;
    vals(n+1)    = sum(w.*g1.*yp.*Kn(:, n+1).*y);
    vals(5+n+1)  = sum(w.*g2.*yp.*Kn(:, n+1).*y);
    yq = (-1)^n*y.^(2*n);
    vals(10+n+1) = sum(w.*g1.*yq.*Kt(:, n+1).*y);
    vals(15+n+1) = sum(w.*g2.*yq.*Kt(:, n+1).*y);
end
vals = vals/pi^2;
end

% ======================================================================================
function [x, w] = gauss_legendre(n)
b = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
J = diag(b, 1) + diag(b, -1);
[V, D] = eig(J);
[x, ix] = sort(diag(D));
w = 2*V(1, ix).'.^2;
x = x(:); w = w(:);
end

function B = bessel_orders(fun, x)
% B(:, m+1) = fun(m, x) for m = 0..4, one call per order (MATLAB requires equal-size arguments;
% the old "row orders x column arguments" table syntax is not portable)
x = x(:);
B = zeros(numel(x), 5);
for m = 0:4
    B(:, m+1) = fun(m, x);
end
end
