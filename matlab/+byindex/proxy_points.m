function proxy = proxy_points(kern, width, eps, opts)
%PROXY_POINTS  Proxy surface for byindex.proxyfun, chosen for a kernel.
%
% proxy = byindex.proxy_points(kern, width, eps)
% proxy = byindex.proxy_points(kern, width, eps, opts)
%
% Adapted from the proxy selection in pcfft's dx_nproxy. Sources fill a
% box of side width, and the proxy surface is nshell concentric spheres,
% the innermost of radius crad*sqrt(3)*width/2:
%
%   1. Each shell has the fewest points that resolve K(.,x) to eps on it,
%      judged from the spherical harmonic tails of the kernel.
%   2. nshell is the smallest number of shells for which a skeleton of
%      random sources reproduces their interactions with targets outside to eps.
%
% Input:
%   kern  - kernel3d object or array; the surface works for all of them
%   width - side length of the box, e.g. the largest side of the
%           bounding box of all points passed to FLAM
%   eps   - tolerance
%   opts  - (optional) options struct
%     opts.nshell_max - maximum number of shells (default 20)
%     opts.nproxy_max - maximum number of proxy points per shell (default 6000)
%     opts.crad       - innermost radius crad*sqrt(3)*width/2 (default
%                       sqrt(3)); FLAM needs crad <= sqrt(3)
%     opts.ifaddtrans - (false) proxyfun also returns the transpose block
%                       to rskelf, for non-symmetric matrices
%     opts.rmax       - largest distance from the box center to a point
%                       outside it (default Inf); the shells stop there
%
% Output:
%   proxy - proxy surface in units of the box size, i.e. for a box of size
%           l centered at ctr the proxy points are proxy.r*l + ctr
%     proxy.r, proxy.n - (3, sum(nproxy)) points and outward normals
%     proxy.w          - (sum(nproxy), 1) area weights
%     proxy.radius     - innermost radius; closer points stay neighbors
%     proxy.nproxy     - (1, nshell) points on each shell, innermost first
%     proxy.nshell, proxy.crad, proxy.ifaddtrans
%     proxy.width, proxy.eps, proxy.opts - the inputs, from which
%                        byindex.proxyfun rebuilds the surface per level

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
    kerns{end+1} = ktmp;
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
rmax     = max(rmax, 2*radius);

% Random sources in the box with random unit normals
src = [];
src.r = (rand(3, nsrc) - 0.5)*width;
src.n = randn(3, nsrc);
src.n = src.n ./ vecnorm(src.n);

%% Step 1: points resolving the kernels on the innermost shell

% Add the box corners before checking resolution
[A, B, C] = ndgrid([-1 1], [-1 1], [-1 1]);
srcres = [];
srcres.r = [src.r, halfside*[A(:).'; B(:).'; C(:).']];
srcres.n = randn(3, nsrc + 8);
srcres.n = srcres.n ./ vecnorm(srcres.n);

% Tensor grid on the unit sphere, with W{l+1} the weighted Legendre functions
Lbig = floor(sqrt(nproxy_max/1.05)) - 1;   % largest degree allowed
Lmax = Lbig + 1;
sph = [];
sph.ntheta = Lmax + 1;
sph.nphi   = 2*Lmax + 2;
[ct, wt] = polytens.lege.pts(sph.ntheta);
ct = ct(:);
st  = sqrt(1 - ct.^2);
phi = (0:sph.nphi-1)*(2*pi/sph.nphi);
sph.n = [reshape(st*cos(phi), 1, []); reshape(st*sin(phi), 1, []); ...
         reshape(repmat(ct, 1, sph.nphi), 1, [])];
sph.W = cell(Lmax+1, 1);
for l = 0:Lmax
    sph.W{l+1} = wt(:) .* legendre(l, ct.', 'norm').';
end

[nproxy1, tail_big] = nproxy_at(kerns, srcres, radius, sph, eps, Lbig);
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

% Targets spread over the annulus [1.01, 2.02]*radius
ntarget = 500;
yv = linspace(1/ntarget - 1, 1 - 1/ntarget, ntarget);
rv = sqrt(1 - yv.^2);
th = pi*(sqrt(5) - 1)*(0:ntarget-1);
targ = [];
targ.r = 1.01*radius*(1 + rand(1, ntarget)) .* [cos(th).*rv; yv; sin(th).*rv];
targ.n = randn(3, ntarget);
targ.n = targ.n ./ vecnorm(targ.n);

nshell   = 0;
err_best = Inf;
for ns = 1:nshell_max
    % Chebyshev-like spacing in 1/r from radius to rmax
    sr = 1/rmax + (cos((0:ns-1)/ns*pi) + 1)/2*(1/radius - 1/rmax);
    rads = 1./sr;
    nproxy = nproxy1*ones(1, ns);
    for is = 2:ns
        nproxy(is) = min(nproxy1, nproxy_at(kerns, srcres, rads(is), sph, eps, Lbig));
    end
    pxy = shells(rads, nproxy);

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
                src.r = (rand(3, nsrc) - 0.5)*width;
                src.n = randn(3, nsrc);
                src.n = src.n ./ vecnorm(src.n);
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
                pb.r = pxy.r(:, idx);
                pb.n = pxy.n(:, idx);
                Kb = kerns{k}.eval(src, pb);
                if nsk < nrows
                    Kp = Kp + randn(nsk, size(Kb, 1))*Kb;
                else
                    Kp(m*(idx(1)-1)+1:m*idx(end), :) = Kb;
                end
            end
            [~, R, p] = qr(Kp, 0);
            d = abs(diag(R));
            krank = sum(d > eps/10*d(1));
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
    if err < err_best
        err_best = err;
        ns_best  = ns;
        pxy_best = pxy;
        rads_best   = rads;
        nproxy_best = nproxy;
    end
end
if nshell == 0
    % Warn if the best number of shells is within a factor of ten
    if err_best < 10*eps
        warning('BYINDEX:proxy_points:nshell', ['BYINDEX.PROXY_POINTS: no ' ...
            'proxy surface with up to %d shells reaches eps = %.1e, using %d ' ...
            'shells with error %.1e'], nshell_max, eps, ns_best, err_best);
        nshell = ns_best;
        pxy    = pxy_best;
        rads   = rads_best;
        nproxy = nproxy_best;
    else
        error('BYINDEX.PROXY_POINTS: no proxy surface with up to %d shells reaches eps = %g', ...
            nshell_max, eps);
    end
end

proxy = [];
proxy.r      = pxy.r/width;
proxy.n      = pxy.n;
proxy.w      = repelem(4*pi*rads(:).^2./nproxy(:), nproxy)/width^2;
proxy.radius = radius/width;
proxy.nproxy = nproxy;
proxy.nshell = nshell;
proxy.crad   = crad;
proxy.ifaddtrans = ifaddtrans;
proxy.width      = width;
proxy.eps        = eps;
proxy.opts       = opts;

end


function [n, tail_big] = nproxy_at(kerns, src, radius, sph, eps, Lbig)
% Points resolving the kernels to eps on a sphere of the given radius,
% capped at degree Lbig; tail_big is the tail there if the cap is hit.
ngrid = sph.ntheta*sph.nphi;
grd = [];
grd.r = radius*sph.n;
grd.n = sph.n;
n = 0;
tail_big = 0;
for k = 1:numel(kerns)
    F    = kerns{k}.eval(src, grd);
    nder = size(F, 1)/ngrid;
    ncol = size(F, 2);

    % Largest coefficient of each degree, per column and kernel component
    deg_energy = zeros(numel(sph.W), ncol);
    for id = 1:nder
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
    n = max(n, ceil(1.05*(L+1)^2));
end
end


function pxy = shells(rads, nproxy)
% Golden angle spirals with nproxy(k) points of radius rads(k)
pxy = [];
pxy.r = zeros(3, sum(nproxy));
i0 = 0;
for k = 1:numel(rads)
    n  = nproxy(k);
    yv = linspace(1/n - 1, 1 - 1/n, n);
    rv = sqrt(1 - yv.^2);
    th = pi*(sqrt(5) - 1)*(0:n-1);
    pxy.r(:, i0+(1:n)) = rads(k)*[cos(th).*rv; yv; sin(th).*rv];
    i0 = i0 + n;
end
pxy.n = pxy.r ./ vecnorm(pxy.r);
end
