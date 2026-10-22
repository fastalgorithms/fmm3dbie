function [G1, G2] = fgw_eval(S, r, kind, form)
%FGW_EVAL  Evaluate the point-load kernels prepared by FGW_SETUP and their radial derivatives.
%
%   [G1, G2] = FGW_EVAL(S, r, kind)
%   [G1, G2] = FGW_EVAL(S, r, kind, form)
%
%   Input
%       S    : structure returned by FGW_SETUP
%       r    : radii >= 0 (array of any shape), in the length units of the parameters of FGW_SETUP
%       kind : 'full' (default) - the finite-depth kernels, i.e. the transforms of
%                                     F1(k) = kT/((a k^4+gam)kT-1),  F2(k) = 1/((a k^4+gam)kT-1),  T = tanh(kd)
%              'deep'           - the deep-water kernels (T replaced by 1)
%              'diff'           - the difference of the two (finite depth minus deep water)
%       form : 'radial' (default) - column n+1 holds d^n G / dr^n,          n = 0..4
%              'theta'            - column m+1 holds ((1/r) d/dr)^m G,      m = 0..4
%                                   (the quantities from which Cartesian derivatives are built, see
%                                   FGW_EVAL_XY; they stay finite and accurate as r -> 0)
%
%   Output
%       G1, G2 : complex arrays of size numel(r) x 5; G1 for the deflection kernel F1, G2 for the
%                potential kernel F2.  Row p corresponds to r(p).
%
%   Conventions: G(r) = (1/2pi) int_0^inf F(k) J0(kr) k dk, i.e. the inverse transform
%   (2pi)^-2 int F(|xi|) e^{i xi.x} dxi; real poles are passed with the outgoing prescription
%   (Im k -> 0+, time dependence e^{-i omega t}), so the results are complex.
%   For the point load P e^{-i omega t}:  deflection w(r) = -(P/(rho omega^2)) G1(r),
%   surface potential phi(r,0) = (i P/(rho omega)) G2(r).
%
%   At r = 0 exactly, quantities of the 'deep'/'full' kernels that have no finite limit there
%   (G1 ~ r^2 log r, G2 contains an r^3 term) are returned as NaN.

if nargin < 3 || isempty(kind), kind = 'full'; end
if nargin < 4 || isempty(form), form = 'radial'; end
kind = lower(kind);
switch lower(form)
    case 'radial', f = 1;
    case 'theta',  f = 2;
    otherwise, error('fgw_eval:form', 'form must be ''radial'' or ''theta''');
end
rr = S.kinf*r(:);                                     % internal units (kinf = 1)
if any(rr < 0) || any(~isfinite(rr)), error('fgw_eval:r', 'r must be finite and >= 0'); end
switch kind
    case 'diff'
        V = eval_diff(S, rr, f);
    case 'deep'
        V = eval_deep(S, rr, f);
    case 'full'
        V = eval_diff(S, rr, f) + eval_deep(S, rr, f);
    otherwise
        error('fgw_eval:kind', 'kind must be ''diff'', ''deep'' or ''full''');
end
% back to the user's units: d^n/dr^n G_j = kinf^{s_j+n} (...)',  Theta^m G_j = kinf^{s_j+2m} (...)',  s = [3 2]
sc = S.kinf.^(f*(0:4));
G1 = V(:, 1:5) .* (S.kinf^3*sc);
G2 = V(:, 6:10) .* (S.kinf^2*sc);
end

% ======================================================================================
function V = eval_diff(S, r, f)
V = zeros(numel(r), 10);
m = r <= S.R1;
if any(m),  V(m, :) = cheb_eval_table(S.difftab, r(m), f); end
if any(~m), V(~m, :) = diff_farfield(S, r(~m), f); end
end

function V = diff_farfield(S, r, f)
% closed form for r > R1 (internal units)
M = numel(r);
V = zeros(M, 10);
H0 = hankel_terms(S, S.k0, r, f);          % M x 5 : d^n/dr^n or Theta^n of H0^(1)(k0 r)
Hi = hankel_terms(S, 1, r, f);
H1 = hankel_terms(S, S.k1, r, f);
Hq = hankel_terms(S, S.q1, r, f);
tau = tau_asymptotic(S, r, f);             % M x 10
for j = 1:2
    c = 5*(j-1);
    V(:, c+(1:5)) = 0.5i*(S.rho0(j)*H0 - S.Ainf(j)*Hi) - imag(S.rho1(j)*H1) + imag(S.Aq1(j)*Hq) - tau(:, c+(1:5));
end
% evanescent modes (1/pi) sum_n rho_n K0(kappa_n r)
for q = 1:numel(S.kappa)
    Kq = bessk_terms(S, S.kappa(q), r, f);
    for j = 1:2
        c = 5*(j-1);
        V(:, c+(1:5)) = V(:, c+(1:5)) + (S.rhon(q, j)/pi)*Kq;
    end
end
end

function T = hankel_terms(S, z, r, f)
% f = 1: d^n/dr^n H0^(1)(z r) = z^n H0^(1)(n)(z r);  f = 2: Theta^m H0^(1)(z r) = (-1)^m z^{2m} H_m^(1)(x)/x^m, x = z r
x = z*r;
T = zeros(numel(r), 5);
ok = imag(x) < 700;                                   % zero where the value underflows
if ~any(ok), return; end
H = bessel_orders(@(m, z) besselh(m, 1, z), x(ok));
if f == 1
    T(ok, :) = (H*S.DJ.') .* (z.^(0:4));
else
    T(ok, :) = (H./(x(ok).^(0:4))) .* ((-1).^(0:4) .* z.^(2*(0:4)));
end
end

function T = bessk_terms(S, kap, r, f)
% same for K0(kappa r): f = 1: kappa^n K0^{(n)};  f = 2: (-1)^m kappa^{2m} K_m(x)/x^m
x = kap*r;
T = zeros(numel(r), 5);
ok = x < 700;
if ~any(ok), return; end
K = bessel_orders(@besselk, x(ok));
if f == 1
    T(ok, :) = (K*S.DK.') .* (kap.^(0:4));
else
    T(ok, :) = (K./(x(ok).^(0:4))) .* ((-1).^(0:4) .* kap.^(2*(0:4)));
end
end

function tau = tau_asymptotic(S, r, f)
% tau_j(r) ~ sum_m tailc(j,m+1) r^{-2m-3}, optimally truncated per r, then differentiated term by term:
%   d^n/dr^n r^{-q} = (-1)^n (q)_n r^{-q-n},   Theta^m r^{-q} = (-1)^m q(q+2)...(q+2m-2) r^{-q-2m}.
% The Taylor coefficients of 1/(1+P^2) oscillate (five singularities of comparable modulus), so the
% truncation point is found from a smoothed envelope (running max over 5 consecutive n=0 term sizes);
% the series also stops once the envelope is below 1e-18 of the accumulated n=0 sum.
M = numel(r);
tau = zeros(M, 10);
nterm = size(S.tailc, 2);
W = 5;
rinv = 1./r; rinv2 = rinv.^2;
if f == 1, rstep = rinv; else, rstep = rinv2; end
for j = 1:2
    c = 5*(j-1);
    acc = zeros(M, 5); active = true(M, 1);
    win = zeros(M, W); prevenv = inf(M, 1);
    rp = rinv.^3;                                     % r^{-2m-3}
    for m = 0:nterm-1
        t0 = S.tailc(j, m+1)*rp;                      % n = 0 term
        win(:, mod(m, W)+1) = abs(t0);
        env = max(win, [], 2);
        if m >= W
            active = active & ~(env > prevenv) & (env > 1e-18*abs(acc(:, 1)));
        end
        acc(active, 1) = acc(active, 1) + t0(active);
        tn = t0; fac = 1; q = 2*m + 3;
        for n = 1:4
            fac = fac*(q + f*(n-1)); tn = -tn.*rstep;   % (q)_n or q(q+2)...
            acc(active, n+1) = acc(active, n+1) + fac*tn(active);
        end
        prevenv = env; rp = rp.*rinv2;
        if ~any(active), break; end
    end
    tau(:, c+(1:5)) = acc;
end
end

% ======================================================================================
function V = eval_deep(S, r, f)
V = zeros(numel(r), 10);
m1 = r <= S.rA; m2 = r > S.rA & r <= S.R1; m3 = r > S.R1;
if any(m1), V(m1, :) = deep_series(S, r(m1), f); end
if any(m2), V(m2, :) = deep_hankel(S, r(m2), f) + cheb_eval_table(S.tautab, r(m2), f); end
if any(m3), V(m3, :) = deep_hankel(S, r(m3), f) + tau_asymptotic(S, r(m3), f); end
end

function V = deep_hankel(S, r, f)
M = numel(r);
V = zeros(M, 10);
Hi = hankel_terms(S, 1, r, f);
Hq = hankel_terms(S, S.q1, r, f);
for j = 1:2
    c = 5*(j-1);
    V(:, c+(1:5)) = 0.5i*S.Ainf(j)*Hi - imag(S.Aq1(j)*Hq);
end
end

function V = deep_series(S, r, f)
% output (j,n) = pe [P_A(s) + log(r) P_B(s)] + po P_C(s) + (singular terms), s = r^2,
% pe, po = 1 or r (flags built in fgw_setup for both forms).
ser = S.series; K = ser.Kmax;
M = numel(r);
V = zeros(M, 10);
s = r.^2;
P = ones(M, K+1);
for k = 2:K+1, P(:, k) = P(:, k-1).*s; end           % s^0 .. s^K
lr = log(r); z = r == 0; lr(z) = 0;
for j = 1:2
    c = 5*(j-1);
    VA = P*ser.PA(:, :, j, f); VB = P*ser.PB(:, :, j, f); VC = P*ser.PC(:, :, j, f);   % M x 5 each
    for n = 0:4
        if ser.pe(n+1, f), pe = r; else, pe = 1; end
        if ser.po(n+1, f), po = r; else, po = 1; end
        v = pe.*(VA(:, n+1) + lr.*VB(:, n+1)) + po.*VC(:, n+1);
        sg = ser.sing{n+1, j, f};                      % rows: [coef_A, coef_B, exponent] with exponent < 0
        for q = 1:size(sg, 1)
            v = v + r.^sg(q, 3).*(sg(q, 1) + sg(q, 2)*lr);
        end
        if any(z)
            bad = ~isempty(sg) || (~ser.pe(n+1, f) && ser.PB(1, n+1, j, f) ~= 0);   % log r or r^-m survives at r = 0
            if bad, v(z) = complex(NaN, NaN); end
        end
        V(:, c+n+1) = v;
    end
end
end

% ======================================================================================
function V = cheb_eval_table(tab, r, f)
% f(x) = sum_k c_k T_k(x) on the panel containing r; columns 1-10 (radial) or 11-20 (theta).
% Large batches: Clenshaw per panel (scalar coefficients); small batches: all points at once.
e = tab.edges(:); N = tab.N;
r = r(:);
idx = panel_index(e, r);
x = (2*r - (e(idx) + e(idx+1)))./(e(idx+1) - e(idx));
M = numel(r);
V = zeros(M, 10);
cols = 10*(f-1) + (1:10);
up = unique(idx);
if M > 3000*numel(up)
    for ip = 1:numel(up)
        p = up(ip); m = idx == p; xx = x(m);
        for q = 1:10
            c = tab.coef(:, p, cols(q));
            b1 = zeros(size(xx)); b2 = b1;
            for k = N:-1:2
                b0 = c(k) + 2*xx.*b1 - b2;
                b2 = b1; b1 = b0;
            end
            V(m, q) = c(1) + xx.*b1 - b2;
        end
    end
else
    chunk = 20000;
    for i0 = 1:chunk:M
        ii = i0:min(i0+chunk-1, M); xx = x(ii); id = idx(ii);
        T = ones(numel(ii), N); T(:, 2) = xx;
        for k = 3:N, T(:, k) = 2*xx.*T(:, k-1) - T(:, k-2); end
        for q = 1:10
            V(ii, q) = sum(T.*tab.coefT(id, :, cols(q)), 2);
        end
    end
end
end

function idx = panel_index(edges, r)
edges = edges(:); r = r(:);
n = numel(edges) - 1;
lo = ones(size(r)); hi = n*ones(size(r));
while any(hi > lo)
    mid = floor((lo + hi)/2);
    m = edges(mid + 1) <= r & hi > lo;
    lo(m) = mid(m) + 1;
    hi(~m & hi > lo) = mid(~m & hi > lo);
end
idx = lo;
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
