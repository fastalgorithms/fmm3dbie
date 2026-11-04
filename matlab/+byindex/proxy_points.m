function proxy = proxy_points(kern_out, width, eps, opts)
%PROXY_POINTS  Proxy surface for byindex.proxyfun, chosen for a kernel.
%
% proxy = byindex.proxy_points(kern_out, width, eps)
% proxy = byindex.proxy_points(kern_out, width, eps, opts)
%
% Adapted from the proxy selection in pcfft's dx_nproxy. Sources fill a
% box with sides width, and the proxy surface is nshell concentric shells,
% spheres stretched with the box, the innermost through crad*sqrt(3)/2
% box sizes along each axis. Points are spread evenly in area over each
% shell:
%
%   1. Each shell has the fewest points that resolve kern_out(.,x),
%      and opts.kern_in(x,.) if given, to eps, judged from their spherical
%      harmonic (Fourier, on rings) tails.
%   2. nshell is the smallest number of shells for which skeletons of
%      random interior sources (and targets, for kern_in) reproduce their
%      exterior interactions to eps.
%      Shells stop being added once 3 more have not halved the error.
%
% Input:
%   kern_out - kernel3d object or array from the box to the proxies, as
%              in byindex.proxyfun (usually the matrix kernel)
%   width    - side lengths of the box, scalar for a cube or 3-vector. A
%              2-vector gives rings in the plane z = 0 instead of shells
%   eps      - tolerance
%   opts     - (optional) options struct
%     opts.kern_in    - incoming kernels, from the proxies to the box, of any
%                       shape; if given, the surface also works for them
%     opts.nshell_max - maximum number of shells (default 20)
%     opts.nproxy_max - maximum number of proxy points per shell (default 6000)
%     opts.crad       - innermost shell at crad*sqrt(3)/2 box sizes
%                       (default sqrt(3)); FLAM needs crad <= sqrt(3)
%     opts.ifaddtrans - (false) proxyfun also returns the transpose block
%                       to rskelf, for non-symmetric matrices
%     opts.rmax       - largest distance from the box center to a point
%                       outside it (default Inf); the shells reach it
%
% Output:
%   proxy - proxy surface for the box centered at 0. For the same shape
%           scaled by s, the points are proxy.r*s + ctr
%     proxy.r, proxy.n - (3, sum(nproxy)) points and outward normals. As
%                        targets, they have zero mean curvature.
%     proxy.du, proxy.dv - (3, sum(nproxy)) d/dphi and -d/dtheta, or
%                        d/dphi and z on rings (proxy.d = du)
%     proxy.w          - (sum(nproxy), 1) area weights, or lengths for rings
%     proxy.radius     - innermost radius in box sizes, i.e. a point x
%                        is inside if norm((x - ctr)./width) < radius
%     proxy.nproxy     - (1, nshell) points on each shell, innermost first
%     proxy.nshell, proxy.crad, proxy.ifaddtrans
%     proxy.width      - (3, 1) or (2, 1) box sides
%     proxy.eps, proxy.opts - the inputs, from which byindex.proxyfun
%                        rebuilds the surface per level
%     proxy.kern_in    - opts.kern_in, if given

if nargin < 4, opts = []; end

nshell_max = 20;
if isfield(opts, 'nshell_max'), nshell_max = opts.nshell_max; end
nproxy_max = 6000;
if isfield(opts, 'nproxy_max'), nproxy_max = opts.nproxy_max; end
crad = sqrt(3);
if isfield(opts, 'crad'), crad = opts.crad; end
nsrc = 117;
if isfield(opts, 'nsrc'), nsrc = opts.nsrc; end
rank_frac = 0.75;
if isfield(opts, 'rank_frac'), rank_frac = opts.rank_frac; end
ifaddtrans = false;
if isfield(opts, 'ifaddtrans'), ifaddtrans = opts.ifaddtrans; end
rmax = Inf;
if isfield(opts, 'rmax'), rmax = opts.rmax; end

if crad > sqrt(3)*(1 + 1e-12)
    error(['BYINDEX.PROXY_POINTS: opts.crad = %g puts non-neighbor boxes ' ...
        'inside the proxy surface, use crad <= sqrt(3)'], crad);
end

% The search compresses interior sources against exterior targets, so the
% incoming kernels are transposed: exterior proxy sources become targets
% and interior targets become sources, with their fields swapped
kin = [];
if isfield(opts, 'kern_in'), kin = opts.kern_in; end
candidates = num2cell(kern_out(:));
for k = 1:numel(kin)
    K = kin(k);
    if ~isempty(K.eval)
        Keval = K.eval;
        K.eval = @(s,t) transpose(Keval(t,s));
    end
    K.opdims = fliplr(K.opdims);
    [K.src_fields, K.targ_fields] = deal(K.targ_fields, K.src_fields);
    candidates{end+1,1} = K;
end

% Kernels to resolve, skipping zero/NaN/empty entries of kernel arrays
dim = 3 - (numel(width) == 2);
carried = {'n', 'du', 'dv', 'mean_curv'};
if dim == 2, carried{end+1} = 'd'; end
kerns = {};
for k = 1:numel(candidates)
    ktmp = candidates{k};
    if isempty(ktmp.eval) || ktmp.iszero || ktmp.isnan, continue; end
    fields = {};
    if ~isempty(ktmp.src_fields),  fields = [fields, ktmp.src_fields(:).']; end
    if ~isempty(ktmp.targ_fields), fields = [fields, ktmp.targ_fields(:).']; end
    fields = setdiff(fields, carried);
    if ~isempty(fields)
        error('BYINDEX.PROXY_POINTS: proxy points only carry r, %s, but kernel %d needs {%s}', ...
            strjoin(carried, ', '), k, strjoin(fields, ', '));
    end
    kerns{end+1} = ktmp;
end
if isempty(kerns)
    error('BYINDEX.PROXY_POINTS: no nonzero kernels given');
end

% Fixed seed, restoring the caller's random state on exit
rng_state = rng;
cleanup = onCleanup(@() rng(rng_state));
rng(0);

% Radii are in box sizes, i.e. a shell of radius rho is rho*(u.*width)
% for u on the unit sphere, or circle for rings. w3 is the box in 3D.
width  = width(:) .* ones(dim, 1);
w3     = [width; zeros(3-dim, 1)];
radius = crad*sqrt(3)/2;
rmax   = max(rmax/min(width), 2*radius);

pm = [];
pm.dim = dim;
if dim == 2
    % Rings are parametrized by the fraction s of their length
    pm.w   = width;
    pm.ph  = linspace(0, 2*pi, 20001);
    pm.s   = cumtrapz(pm.ph, vecnorm(pm.w .* [sin(pm.ph); cos(pm.ph)]));
    pm.len = pm.s(end);
    pm.s   = pm.s/pm.len;
else
    % Shells are parametrized by s = cos(theta) and phi about the side most
    % unlike the other two, cycled to be the last axis. The height t(s) on
    % the shell is chosen so that equal steps in s cover equal areas, as on
    % a sphere, where t = s.
    d = abs(width([2 3 1]) - width([3 1 2]));
    order = [3 1 2];
    [~, i] = min(d(order));
    pm.perm = circshift(1:3, -order(i));
    pm.w    = width(pm.perm);
    pm.t    = linspace(-1, 1, 20001);
    pm.s    = cumtrapz(pm.t, band_area(pm.t, pm.w));
    pm.area = pm.s(end);
    pm.s    = 2*pm.s/pm.area - 1;
end

% Random sources in the box with random frames
src = randinfo((rand(3, nsrc) - 0.5).*w3, width);

%% Step 1: points resolving the kernels on the innermost shell

% Add the box corners before checking resolution
[A, B, C] = ndgrid([-1 1], [-1 1], [-1 1]);
srcres = randinfo([src.r, w3/2 .* [A(:).'; B(:).'; C(:).']], width);

sph = [];
if dim == 2
    % Equispaced grid on the rings, for their Fourier modes
    Lbig  = floor((nproxy_max/1.05 - 1)/2);   % largest mode allowed
    sph.s = (0:2*Lbig+3)/(2*Lbig+4);
    sph.ph = [];
else
    % Tensor grid on the unit sphere, with W{l+1} the weighted Legendre functions
    Lbig = floor(sqrt(nproxy_max/1.05)) - 1;   % largest degree allowed
    Lmax = Lbig + 1;
    sph.ntheta = Lmax + 1;
    sph.nphi   = 2*Lmax + 2;
    [ct, wt] = polytens.lege.pts(sph.ntheta);
    ct = ct(:);
    phi = (0:sph.nphi-1)*(2*pi/sph.nphi);
    sph.s  = reshape(repmat(ct, 1, sph.nphi), 1, []);
    sph.ph = reshape(repmat(phi, sph.ntheta, 1), 1, []);
    sph.W = cell(Lmax+1, 1);
    for l = 0:Lmax
        sph.W{l+1} = wt(:) .* legendre(l, ct.', 'norm').';
    end
end

[nproxy1, tail_big] = nproxy_at(kerns, srcres, radius, pm, sph, eps, Lbig);
if tail_big >= eps
    % Warn if the largest allowed nproxy is within a factor of ten
    if tail_big < 10*eps
        warning('BYINDEX:proxy_points:nproxy', ['BYINDEX.PROXY_POINTS: ' ...
            'nproxy_max = %d only resolves the kernels to %.1e > eps = %.1e'], ...
            nproxy_max, tail_big, eps);
    else
        error('BYINDEX.PROXY_POINTS: nproxy did not converge below nproxy_max = %d', ...
            nproxy_max);
    end
end

%% Step 2: find the number of shells

% Targets spread over the annulus [1.01, 2.02]*radius, with any frame
% and mean curvature
ntarget = 500;
if dim == 2
    th = 2*pi*(0:ntarget-1)/ntarget;
    u  = [cos(th); sin(th); zeros(1, ntarget)];
else
    yv = linspace(1/ntarget - 1, 1 - 1/ntarget, ntarget);
    rv = sqrt(1 - yv.^2);
    th = pi*(sqrt(5) - 1)*(0:ntarget-1);
    u  = [cos(th).*rv; yv; sin(th).*rv];
end
targ = randinfo(1.01*radius*(1 + rand(1, ntarget)) .* u .* w3, width);

nshell   = 0;
err_best = Inf;
err_ref  = Inf;    % error at the last shell count that halved it
ns_ref   = 0;
for ns = 1:nshell_max
    % Chebyshev-like spacing in 1/r from radius to rmax
    sr = 1/rmax + (cos((0:ns-1)/ns*pi) + 1)/2*(1/radius - 1/rmax);
    rads = 1./sr;
    nproxy = nproxy1*ones(1, ns);
    for is = 2:ns
        nproxy(is) = min(nproxy1, nproxy_at(kerns, srcres, rads(is), pm, sph, eps, Lbig));
    end
    pxy = shells(rads, nproxy, pm);

    err = 0;
    for k = 1:numel(kerns)
        % Skeleton of the sources from their proxy interactions, adding
        % sources until the rank is at most rank_frac of the columns
        krank = Inf;
        ncol  = 0;
        nsk   = Inf;
        while krank > rank_frac*ncol
            if ncol > 0
                cps = ncol/nsrc;   % columns per source
                if krank < nsk
                    % The sketch resolved the rank, so take enough sources
                    nsrc = ceil(1.25*krank/(rank_frac*cps));
                else
                    nsrc = 2*nsrc;
                end
                nsrc  = min(nsrc, ceil(nrows/(rank_frac*cps)));
                src = randinfo((rand(3, nsrc) - 0.5).*w3, width);
            end
            m     = kerns{k}.opdims(1);
            nrows = m*size(pxy.r, 2);
            ncol  = kerns{k}.opdims(2)*nsrc;
            % Row sketch; later passes only test rank <= rank_frac*ncol
            if isinf(nsk)
                nsk = min(nrows, ncol);
            else
                nsk = min(nrows, ceil(1.25*rank_frac*ncol) + 10);
            end
            nb = max(1, floor(2e6/(m*ncol)));   % proxy points per block
            Kp = zeros(nsk, ncol);
            for i0 = 1:nb:size(pxy.r, 2)
                idx = i0:min(i0+nb-1, size(pxy.r, 2));
                pb = [];
                pb.r  = pxy.r(:, idx);
                pb.n  = pxy.n(:, idx);
                pb.du = pxy.du(:, idx);
                pb.dv = pxy.dv(:, idx);
                if dim == 2, pb.d = pxy.d(:, idx); end
                pb.mean_curv = zeros(numel(idx), 1);
                Kb = kerns{k}.eval(src, pb);
                check_nan(Kb, k);
                if nsk < nrows
                    Kp = Kp + randn(nsk, size(Kb, 1))*Kb;
                else
                    Kp(m*(idx(1)-1)+1:m*idx(end), :) = Kb;
                end
            end
            % ID well below eps, so the error at the targets judges the
            % shells rather than how much the kernel amplifies the ID
            [~, R, p] = qr(Kp, 0);
            d = abs(diag(R));
            krank = sum(d > eps/100*d(1));
        end

        % Error of the skeleton at the targets
        T  = R(1:krank, 1:krank) \ R(1:krank, krank+1:end);
        Kt = kerns{k}.eval(src, targ);
        check_nan(Kt, k);
        E  = Kt(:, p(krank+1:end)) - Kt(:, p(1:krank))*T;
        err = max(err, max(abs(E(:)))/max(abs(Kt(:))));
    end

    if err < eps
        nshell = ns;
        break
    end
    if err < err_best
        err_best = err;
        ns_best  = ns;
        pxy_best = pxy;
        rads_best   = rads;
        nproxy_best = nproxy;
    end
    if err < err_ref/2
        err_ref = err;
        ns_ref  = ns;
    elseif ns - ns_ref >= 3
        break
    end
end
if nshell == 0
    % Warn if the best number of shells is within a factor of 100
    if err_best < 100*eps
        warning('BYINDEX:proxy_points:nshell', ['BYINDEX.PROXY_POINTS: proxy ' ...
            'error stagnated at %.1e with %d shells, above eps = %.1e; ' ...
            'using it'], err_best, ns_best, eps);
        nshell = ns_best;
        pxy    = pxy_best;
        rads   = rads_best;
        nproxy = nproxy_best;
    else
        error(['BYINDEX.PROXY_POINTS: proxy error stagnated at %.1e with ' ...
            '%d shells, above eps = %g'], err_best, ns_best, eps);
    end
end

proxy = [];
proxy.r      = pxy.r;
proxy.n      = pxy.n;
proxy.du     = pxy.du;
proxy.dv     = pxy.dv;
if dim == 2, proxy.d = pxy.d; end
proxy.w      = pxy.w;
proxy.radius = radius;
proxy.nproxy = nproxy;
proxy.nshell = nshell;
proxy.crad   = crad;
proxy.ifaddtrans = ifaddtrans;
proxy.width      = width;
proxy.eps        = eps;
proxy.opts       = opts;
if isfield(opts, 'kern_in'), proxy.kern_in = kin; end

end


function [n, tail_big] = nproxy_at(kerns, src, radius, pm, sph, eps, Lbig)
% Points resolving the kernels to eps on the shell of the given radius,
% capped at degree Lbig; tail_big is the tail there if the cap is hit.
ngrid = numel(sph.s);
grd = stretch(sph.s, sph.ph, radius, pm);
grd.mean_curv = zeros(ngrid, 1);
n = 0;
tail_big = 0;
for k = 1:numel(kerns)
    F    = kerns{k}.eval(src, grd);
    check_nan(F, k);
    nder = size(F, 1)/ngrid;
    ncol = size(F, 2);

    % Largest coefficient of each degree (Fourier mode on rings), per
    % column and kernel component
    deg_energy = zeros(Lbig+2, ncol);
    for id = 1:nder
        if pm.dim == 2
            G = fft(F(id:nder:end, :))/ngrid;
            e = max(abs(G(1:Lbig+2, :)), abs(G([1, ngrid:-1:ngrid-Lbig], :)));
            deg_energy = max(deg_energy, e);
            continue
        end
        G = fft(reshape(F(id:nder:end, :), sph.ntheta, sph.nphi, ncol), [], 2)/sph.nphi;
        for l = 0:numel(sph.W)-1
            cp = sum(G(:, 1:l+1, :) .* sph.W{l+1}, 1);
            cn = sum(G(:, [1, sph.nphi:-1:sph.nphi-l+1], :) .* sph.W{l+1}, 1);
            e  = max(max(abs(cp), [], 2), max(abs(cn), [], 2));
            deg_energy(l+1, :) = max(deg_energy(l+1, :), reshape(e, 1, ncol));
        end
    end

    scl = max(deg_energy, [], 1);
    scl(scl == 0) = 1;
    tail = max(flipud(cummax(flipud(deg_energy), 1)) ./ scl, [], 2);

    L = find(tail(1:Lbig+1) < eps, 1) - 1;
    if isempty(L)
        L = Lbig;
        tail_big = max(tail_big, tail(Lbig+1));
    end
    if pm.dim == 2
        n = max(n, ceil(1.05*(2*L+1)));
    else
        n = max(n, ceil(1.05*(L+1)^2));
    end
end
end


function pxy = shells(rads, nproxy, pm)
% Golden angle spirals in (s, phi) with nproxy(k) points on shell rads(k),
% or nproxy(k) points evenly spaced in length on ring rads(k), with their
% frames and weights
sv = zeros(1, sum(nproxy));
ph = zeros(1, sum(nproxy));
i0 = 0;
for k = 1:numel(rads)
    n = nproxy(k);
    if pm.dim == 2
        sv(i0+(1:n)) = ((0:n-1) + 0.5)/n;
    else
        sv(i0+(1:n)) = linspace(1/n - 1, 1 - 1/n, n);
        ph(i0+(1:n)) = pi*(sqrt(5) - 1)*(0:n-1);
    end
    i0 = i0 + n;
end
pxy = stretch(sv, ph, repelem(rads(:).', nproxy), pm);
measure = 4*pi;    % of (s, phi), or of s on rings
if pm.dim == 2, measure = 1; end
pxy.w = (pxy.w .* repelem(measure./nproxy(:).', nproxy)).';
end


function pts = stretch(s, ph, rho, pm)
% Points rho*u.*pm.w at height t(s) about the polar axis of pm, with
% outward normals, du = d/dphi, dv = -d/dtheta for s = cos(theta), and w
% the area per unit s and phi. On rings, at the fraction s of the length,
% with du = d = d/dphi, dv = z and w the length per unit s.
if pm.dim == 2
    ph = interp1(pm.s, pm.ph, s, 'spline');
    w  = pm.w(:);
    z  = zeros(size(s));
    pts = [];
    pts.r  = rho .* [w(1)*cos(ph); w(2)*sin(ph); z];
    pts.n  = [w(2)*cos(ph); w(1)*sin(ph); z];
    pts.n  = pts.n ./ vecnorm(pts.n);
    pts.du = rho .* [-w(1)*sin(ph); w(2)*cos(ph); z];
    pts.d  = pts.du;
    pts.dv = [z; z; ones(size(s))];
    pts.w  = rho .* pm.len .* ones(size(s));
    return
end
t  = interp1(pm.s, pm.t, s, 'spline');
st = sqrt(1 - t.^2);
w  = pm.w(:);
u  = [cos(ph).*st; sin(ph).*st; t];
dtds = pm.area/2 ./ band_area(t, w);
x  = rho .* u .* w;
n  = u ./ w;
nn = vecnorm(n);
du = rho .* [-u(2, :); u(1, :); zeros(size(t))] .* w;
dv = rho .* [-cos(ph).*t./st; -sin(ph).*t./st; ones(size(t))] .* w .* dtds .* sqrt(1 - s.^2);
pts = [];
pts.r(pm.perm, :)  = x;
pts.n(pm.perm, :)  = n ./ nn;
pts.du(pm.perm, :) = du;
pts.dv(pm.perm, :) = dv;
pts.w = rho.^2 .* prod(w) .* nn .* dtds;
end


function a = band_area(t, w)
% Area of the shell u.*w per unit t = cos(theta), about the z axis
ph = (0:63).'*(2*pi/64);
st = sqrt(1 - t(:).'.^2);
a  = 2*pi*prod(w)*mean(sqrt((cos(ph).*st/w(1)).^2 + (sin(ph).*st/w(2)).^2 + ...
    (t(:).'/w(3)).^2), 1);
end


function check_nan(F, k)
% Error at once if kernel k is NaN, e.g. beyond the range of a fitted kernel
if any(isnan(F(:)))
    error('BYINDEX:proxy_points:nan', ['BYINDEX.PROXY_POINTS: kernel %d is NaN ' ...
        'at the proxy points; for a kernel fitted up to a distance, such as ' ...
        'KERNEL3D.RADPP, set opts.rmax to that distance less half the box diagonal'], k);
end
end


function pts = randinfo(r, width)
% Points r with random unit normals, orthogonal frames du, dv = n x du and
% mean curvatures on the scale of the box. For rings, frames lie in z=0
% and d is the planar tangent.
n = size(r, 2);
pts = [];
pts.r  = r;
if numel(width) == 2
    nxy = randn(2,n); nxy = nxy./vecnorm(nxy);
    pts.n = [nxy;zeros(1,n)];
    pts.du = [-nxy(2,:);nxy(1,:);zeros(1,n)];
    pts.d = pts.du;
    pts.dv = repmat([0;0;1],1,n);
else
    pts.n  = randn(3, n);
    pts.n  = pts.n ./ vecnorm(pts.n);
    pts.du = randn(3, n);
    pts.du = pts.du - sum(pts.du .* pts.n) .* pts.n;
    pts.du = pts.du ./ vecnorm(pts.du);
    pts.dv = cross(pts.n, pts.du);
end
pts.mean_curv = randn(n, 1)/min(width);
end
