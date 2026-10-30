function proxy = proxy_points(kern, width, eps, nshell_max, opts)
%PROXY_POINTS  Proxy surface for byindex.proxyfun, chosen for a kernel.
%
% proxy = byindex.proxy_points(kern, width, eps)
% proxy = byindex.proxy_points(kern, width, eps, nshell_max)
% proxy = byindex.proxy_points(kern, width, eps, nshell_max, opts)
% proxy = byindex.proxy_points(kern, width, eps, opts)
%
% Adapted from the proxy selection in pcfft's dx_nproxy. Sources fill a
% box of side width, and the proxy surface is nshell concentric spheres of
% nproxy points each, the innermost of radius crad*sqrt(3)*width/2. The
% selection is done in two steps:
%
%   1. Find nproxy that resolves K(.,x) down to eps on a single sphere by
%      looking at the spherical harmonic tails of the kernel. Leaves a
%      margin of at least 3 degrees and at least 25%.
%
%   2. For increasing nshell, check that a skeleton of random sources in
%      the box chosen from their interactions with the proxy points (a
%      column-pivoted QR to eps/10) reproduces their interactions with
%      targets outside the proxy surface to eps.
%
% Proxy points only carry r and n (outward radial normal), so no kernel
% may need fields other than 'n' in src_fields or targ_fields.
%
% Input:
%   kern       - kernel3d object or array of kernel3d objects. The proxy
%                surface is chosen to work for all of them.
%   width      - side length of the box, e.g. the largest side of the
%                bounding box of all points passed to FLAM
%   eps        - tolerance
%   nshell_max - (optional) maximum number of shells (default 20)
%   opts       - (optional) options struct
%     opts.nshell_max - maximum number of shells, if not given as an
%                       argument (default 20)
%     opts.nproxy_max - maximum number of proxy points per shell
%                       (default 6000)
%     opts.crad       - the innermost shell has radius
%                       crad*sqrt(3)*width/2 (default sqrt(3), i.e. 1.5*width).
%                       FLAM requires crad <= sqrt(3) so that every
%                       non-neighbor box is outside the proxy surface.
%     opts.nsrc       - initial number of random sources (default 117)
%     opts.rank_frac  - double the random sources while the rank is more
%                       than rank_frac of them (default 0.5)
%     opts.ifaddtrans - stored in proxy. If true, byindex.proxyfun called
%                       from rskelf also returns the transpose block, as
%                       needed for non-symmetric matrices (default false)
%
% Output:
%   proxy - struct describing the proxy surface in units of the box size,
%           i.e. for a box of size l centered at ctr the proxy points are
%           proxy.r*l + ctr.
%     proxy.r          - (3, nproxy*nshell) proxy points
%     proxy.n          - (3, nproxy*nshell) outward unit normals
%     proxy.w          - (nproxy*nshell, 1) area weights
%     proxy.radius     - radius of the innermost shell; points closer to
%                        the box center are kept as neighbors
%     proxy.nproxy     - number of proxy points per shell
%     proxy.nshell     - number of shells
%     proxy.crad       - see opts.crad
%     proxy.ifaddtrans - see opts.ifaddtrans
%     proxy.width, proxy.eps, proxy.nshell_max, proxy.opts - the inputs,
%                        used by byindex.proxyfun to rebuild the proxy
%                        surface for each level

if nargin < 4, nshell_max = []; end
if nargin < 5, opts = []; end
if isstruct(nshell_max)
    opts = nshell_max;
    nshell_max = [];
end
if isempty(nshell_max)
    nshell_max = 20;
    if isfield(opts, 'nshell_max'), nshell_max = opts.nshell_max; end
end

nproxy_max = 6000;
if isfield(opts, 'nproxy_max'), nproxy_max = opts.nproxy_max; end
crad = sqrt(3);
if isfield(opts, 'crad'), crad = opts.crad; end
nsrc = 117;
if isfield(opts, 'nsrc'), nsrc = opts.nsrc; end
rank_frac = 0.5;
if isfield(opts, 'rank_frac'), rank_frac = opts.rank_frac; end
ifaddtrans = false;
if isfield(opts, 'ifaddtrans'), ifaddtrans = opts.ifaddtrans; end

if crad > sqrt(3)*(1 + 1e-12)
    error(['BYINDEX.PROXY_POINTS: opts.crad = %g puts non-neighbor boxes ' ...
        'inside the proxy surface, use crad <= sqrt(3)'], crad);
end

% Kernels to resolve, skipping zero/NaN/empty entries of kernel arrays
kerns = {};
for k = 1:numel(kern)
    ktmp = kern(k);
    if isempty(ktmp.eval) || ktmp.iszero || ktmp.isnan, continue; end
    fields = {};
    if ~isempty(ktmp.src_fields),  fields = [fields, ktmp.src_fields(:).']; end
    if ~isempty(ktmp.targ_fields), fields = [fields, ktmp.targ_fields(:).']; end
    fields = setdiff(fields, {'n'});
    if ~isempty(fields)
        error(['BYINDEX.PROXY_POINTS: proxy points only carry r and n, ' ...
            'but kernel %d needs {%s}'], k, strjoin(fields, ', '));
    end
    kerns{end+1} = ktmp; %#ok<AGROW>
end
if isempty(kerns)
    error('BYINDEX.PROXY_POINTS: no nonzero kernels given');
end

% Fixed seed, restoring the caller's random state on exit
rng_state = rng;
cleanup = onCleanup(@() rng(rng_state));
rng(0);

halfside = width/2;
radius   = crad*sqrt(3)*halfside;

% Random sources in the box with random unit normals
src = [];
src.r = (rand(3, nsrc) - 0.5)*width;
src.n = randn(3, nsrc);
src.n = src.n ./ vecnorm(src.n);

%% Step 1: find nproxy that resolves the kernels on one sphere

% Add the box corners before checking resolution
[A, B, C] = ndgrid([-1 1], [-1 1], [-1 1]);
srcres = [];
srcres.r = [src.r, halfside*[A(:).'; B(:).'; C(:).']];
srcres.n = randn(3, nsrc + 8);
srcres.n = srcres.n ./ vecnorm(srcres.n);

% Tensor grid resolving every degree up to Lmax
Lmax   = ceil(sqrt(nproxy_max/1.05));
ntheta = Lmax + 1;
nphi   = 2*Lmax + 2;
[ct, wt] = polytens.lege.pts(ntheta);
ct = ct(:);
wt = wt(:);
st  = sqrt(1 - ct.^2);
phi = (0:nphi-1)*(2*pi/nphi);

grd = [];
grd.n = [reshape(st*cos(phi), 1, []); reshape(st*sin(phi), 1, []); ...
         reshape(repmat(ct, 1, nphi), 1, [])];
grd.r = radius*grd.n;
ngrid = ntheta*nphi;

nproxy = 0;
for k = 1:numel(kerns)
    F    = kerns{k}.eval(srcres, grd);
    nder = size(F, 1)/ngrid;
    ncol = size(F, 2);

    % Largest spherical harmonic coefficient of each degree, kept per
    % column so each is normalized by its own peak. Each kernel component
    % is treated as its own function on the sphere.
    deg_energy = zeros(Lmax+1, ncol);
    for id = 1:nder
        G = fft(reshape(F(id:nder:end, :), ntheta, nphi, ncol), [], 2)/nphi;
        for l = 0:Lmax
            W  = wt .* legendre(l, ct.', 'norm').';
            cp = sum(G(:, 1:l+1, :) .* W, 1);
            cn = sum(G(:, [1, nphi:-1:nphi-l+1], :) .* W, 1);
            e  = max(max(abs(cp), [], 2), max(abs(cn), [], 2));
            deg_energy(l+1, :) = max(deg_energy(l+1, :), reshape(e, 1, ncol));
        end
    end

    scl = max(deg_energy, [], 1);
    scl(scl == 0) = 1;
    tail = max(flipud(cummax(flipud(deg_energy), 1)) ./ scl, [], 2);

    nproxy_k = 0;
    for L = 0:Lmax-1
        Lc = L - max(3, ceil(0.25*L));
        if Lc >= 0 && tail(Lc+2) < eps
            nproxy_k = ceil(1.05*(L+1)^2);
            break
        end
    end
    if nproxy_k == 0 || nproxy_k > nproxy_max
        error('BYINDEX.PROXY_POINTS: nproxy did not converge below nproxy_max = %d', ...
            nproxy_max);
    end
    nproxy = max(nproxy, nproxy_k);
end

%% Step 2: find the number of shells

% Unit sphere points, approximately uniform (golden angle spiral)
yv = linspace(1/nproxy - 1, 1 - 1/nproxy, nproxy);
rv = sqrt(1 - yv.^2);
th = pi*(sqrt(5) - 1)*(0:nproxy-1);
p0 = [cos(th).*rv; yv; sin(th).*rv];

% Targets spread over the annulus [1.01, 2.02]*radius
ntarget = 500;
yv = linspace(1/ntarget - 1, 1 - 1/ntarget, ntarget);
rv = sqrt(1 - yv.^2);
th = pi*(sqrt(5) - 1)*(0:ntarget-1);
targ = [];
targ.r = 1.01*radius*(1 + rand(1, ntarget)) .* [cos(th).*rv; yv; sin(th).*rv];
targ.n = randn(3, ntarget);
targ.n = targ.n ./ vecnorm(targ.n);

nshell = 0;
for ns = 1:nshell_max
    % Chebyshev-like spacing in 1/r, starting from the innermost shell
    if ns > 1
        sr   = cos((1:ns-1)/ns*pi);
        sr   = (sr(:) + 1)/2*(1/radius);
        rads = [radius, 1./sr(:).'];
    else
        rads = radius;
    end
    pxy = [];
    pxy.r = reshape(rads .* reshape(p0, 3, 1, []), 3, []);
    pxy.n = pxy.r ./ vecnorm(pxy.r);

    err = 0;
    for k = 1:numel(kerns)
        % Skeleton of the sources from their proxy interactions, doubling
        % the sources while the rank is above rank_frac of them. The rank
        % is at most the number of proxy rows, so this stops by the time
        % there are nrows/rank_frac columns.
        krank = Inf;
        ncol  = 0;
        while krank > rank_frac*ncol
            if ncol > 0
                nsrc  = min(2*nsrc, ceil(size(Kp, 1)/(rank_frac*ncol/nsrc)));
                src.r = (rand(3, nsrc) - 0.5)*width;
                src.n = randn(3, nsrc);
                src.n = src.n ./ vecnorm(src.n);
            end
            Kp = kerns{k}.eval(src, pxy);
            [~, R, p] = qr(Kp, 0);
            d = abs(diag(R));
            krank = sum(d > eps/10*d(1));
            ncol  = size(Kp, 2);
        end

        % Error of the skeleton at the targets
        T  = R(1:krank, 1:krank) \ R(1:krank, krank+1:end);
        Kt = kerns{k}.eval(src, targ);
        E  = Kt(:, p(krank+1:end)) - Kt(:, p(1:krank))*T;
        err = max(err, max(abs(E(:)))/max(abs(Kt(:))));
    end

    if err < eps
        nshell = ns;
        break
    end
end
if nshell == 0
    error('BYINDEX.PROXY_POINTS: no proxy surface with up to %d shells reaches eps = %g', ...
        nshell_max, eps);
end

proxy = [];
proxy.r      = pxy.r/width;
proxy.n      = pxy.n;
proxy.w      = reshape(repmat(4*pi*rads(:).^2/nproxy, 1, nproxy), [], 1)/width^2;
proxy.radius = radius/width;
proxy.nproxy = nproxy;
proxy.nshell = nshell;
proxy.crad   = crad;
proxy.ifaddtrans = ifaddtrans;
proxy.width      = width;
proxy.eps        = eps;
proxy.nshell_max = nshell_max;
proxy.opts       = opts;

end
