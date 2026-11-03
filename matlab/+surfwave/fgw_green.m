function [val, grad, hess, third] = fgw_green(S, which, src, targ)
%SURFWAVE.FGW_GREEN   Finite depth correction to the surface
%   wave Green's functions.
%
%   [val,grad,hess,third] = SURFWAVE.FGW_GREEN(S,which,src,targ), with
%   S from SURFWAVE.FGW_SETUP and which one of 'gs', 'gphi' or
%   'bilap_gphi' (value only).
%
% See also SURFWAVE.FGW_SETUP, SURFWAVE.KERNEL

if isstruct(src),  src  = src.r;  end
if isstruct(targ), targ = targ.r; end
ns = size(src, 2);  nt = size(targ, 2);
rx = targ(1,:).' - src(1,:);
ry = targ(2,:).' - src(2,:);

switch lower(which)
    case 'bilap_gphi'
        if nargout > 1
            error('SURFWAVE.FGW_GREEN: ''bilap_gphi'' is value only.');
        end
        % Delta^2 g = r^4 Th^4 g + 8 r^2 Th^3 g + 8 Th^2 g,  Th = (1/r) d/dr
        r = hypot(rx(:), ry(:));
        [~, T2] = fgw_eval(S, r, 'theta');
        val = reshape(r.^4.*T2(:,5) + 8*r.^2.*T2(:,4) + 8*T2(:,3), nt, ns);
        return
    case 'gs'
        [D, ~] = fgw_eval_xy(S, rx, ry);
        D = 2*D;
    case 'gphi'
        [~, D] = fgw_eval_xy(S, rx, ry);
    otherwise
        error('SURFWAVE.FGW_GREEN: unknown kernel ''%s''.', which);
end
Dab = @(a, b) reshape(D(:, a+1, b+1), nt, ns);       % d^(a+b) G / dx^a dy^b
val = Dab(0, 0);
if nargout > 1, grad  = cat(3, Dab(1,0), Dab(0,1)); end
if nargout > 2, hess  = cat(3, Dab(2,0), Dab(1,1), Dab(0,2)); end
if nargout > 3, third = cat(3, Dab(3,0), Dab(2,1), Dab(1,2), Dab(0,3)); end
end


function [D1, D2] = fgw_eval_xy(S, x, y)
% Cartesian derivatives D(:,a+1,b+1) = d^(a+b) G / dx^a dy^b, a+b <= 4, of
% the deflection (D1) and potential (D2) kernels, from FGW_EVAL 'theta'

x = x(:); y = y(:);
if numel(x) ~= numel(y), error('fgw_eval_xy:args', 'x and y must have the same number of elements'); end
r = hypot(x, y);
[T1, T2] = fgw_eval(S, r, 'theta');
M = numel(x);
D1 = nan(M, 5, 5); D2 = nan(M, 5, 5);
c = [1 0 0; 1 0 0; 1 1 0; 1 3 0; 1 6 3];             % c(a+1, i+1) = a!/(i!(a-2i)! 2^i)
xp = x.^(0:4); yp = y.^(0:4);                         % x^0..x^4, y^0..y^4
for a = 0:4
    for b = 0:4-a
        v1 = zeros(M, 1); v2 = zeros(M, 1);
        for i = 0:floor(a/2)
            for j = 0:floor(b/2)
                pf = c(a+1, i+1)*c(b+1, j+1)*xp(:, a-2*i+1).*yp(:, b-2*j+1);
                m = a + b - i - j;
                v1 = v1 + pf.*T1(:, m+1);
                v2 = v2 + pf.*T2(:, m+1);
            end
        end
        D1(:, a+1, b+1) = v1; D2(:, a+1, b+1) = v2;
    end
end
end


function [G1, G2] = fgw_eval(S, r, form)
% finite depth minus deep water kernels G1, G2 (numel(r) x 5): d^n G/dr^n
% (form 'radial') or ((1/r) d/dr)^n G (form 'theta'), n = 0..4

switch lower(form)
    case 'radial', f = 1;
    case 'theta',  f = 2;
    otherwise, error('fgw_eval:form', 'form must be ''radial'' or ''theta''');
end
rr = S.kinf*r(:);                                     % internal units (kinf = 1)
if any(rr < 0) || any(~isfinite(rr)), error('fgw_eval:r', 'r must be finite and >= 0'); end
V = eval_diff(S, rr, f);
% back to the user's units
sc = S.kinf.^(f*(0:4));
G1 = V(:, 1:5) .* (S.kinf^3*sc);
G2 = V(:, 6:10) .* (S.kinf^2*sc);
end

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
tau = tau_asymptotic(S, r, f);             % M x 10
for j = 1:2
    c = 5*(j-1);
    V(:, c+(1:5)) = 0.5i*(S.rho0(j)*H0 - S.Ainf(j)*Hi) - tau(:, c+(1:5));
end
% first-quadrant complex roots, finite depth and deep
for q = 1:numel(S.k1)
    H1 = hankel_terms(S, S.k1(q), r, f);
    for j = 1:2, V(:, 5*(j-1)+(1:5)) = V(:, 5*(j-1)+(1:5)) - imag(S.rho1(q, j)*H1); end
end
for q = 1:numel(S.q1)
    Hq = hankel_terms(S, S.q1(q), r, f);
    for j = 1:2, V(:, 5*(j-1)+(1:5)) = V(:, 5*(j-1)+(1:5)) + imag(S.Aq1(q, j)*Hq); end
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
% derivatives of H0^(1)(z r), radial (f = 1) or ((1/r) d/dr)^n (f = 2)
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
% same for K0(kappa r)
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
% deep water algebraic tail, an optimally truncated series in 1/r
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

function V = cheb_eval_table(tab, r, f)
% Chebyshev table at r, columns 1-10 (radial) or 11-20 (theta)
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
% B(:, m+1) = fun(m, x), m = 0..4
x = x(:);
B = zeros(numel(x), 5);
for m = 0:4
    B(:, m+1) = fun(m, x);
end
end
