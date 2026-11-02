function funs = flamfuns(surfers, kern, eps, targs, opts)
%FLAMFUNS  Inputs for compressing a surfer matrix with FLAM.
%
% funs = byindex.flamfuns(surfers, kern, eps)
% funs = byindex.flamfuns(surfers, kern, eps, targs)
% funs = byindex.flamfuns(surfers, kern, eps, targs, opts)
%
% Builds the quadrature corrections and the proxy surface, and returns the
% matrix entry and proxy functions that FLAM takes. FLAM is not called.
%
% With opts.type = 'ifmm' (default), Afun is the smooth rule on the
% oversampled sources, and the interpolation and corrections are applied
% outside:
%
%   funs = byindex.flamfuns(surfers, kern, eps);          % or (..., targs)
%   F = ifmm(funs.Afun, funs.rx, funs.cx, funs.occ, eps, funs.pxyfun);
%   u = ifmm_mv(F, funs.xinterp*sigma, funs.Afun) + funs.Q*sigma;
%
% which matches surfermat (or surferkernevalmat for targs) with
% opts.unif_nover = true.
%
% With opts.type = 'rskelf', the matrix is square, Afun includes the
% corrections, and nothing is oversampled, i.e. it matches surfermat with
% opts.ifoversamp = 0:
%
%   funs = byindex.flamfuns(surfers, kern, eps, [], struct('type', 'rskelf'));
%   F = rskelf(funs.Afun, funs.x, funs.occ, eps, funs.pxyfun);
%
% The proxy surfaces follow FLAM's own boxes, so do not pass opts.ext.
% The proxy surface is valid on every level. With the fastalgorithms fork
% of FLAM, opts.proxybylevel = true rebuilds it for each level.
%
% Input:
%   surfers - surfer object or array of surfer objects
%   kern    - kernel3d object, scalar or nsurfers x nsurfers, or for
%             targs a vector with one kernel per surfer
%   eps     - quadrature and proxy tolerance
%   targs   - (optional) [] for the on-surface matrix, or targets as a
%             numeric (3 x ntarg) array or a struct with .r (3 x ntarg)
%   opts    - (optional) options struct
%     opts.type       - 'ifmm' (default) or 'rskelf'
%     opts.symm       - FLAM symmetry flag, 'n' (default), 's' or 'h'
%     opts.ifoversamp - oversample the sources for ifmm (default true)
%     opts.quadopts   - options for surfermat or surferkernevalmat
%     opts.proxyopts  - options for byindex.proxy_points
%     opts.occ        - (2048) FLAM leaf occupancy. For rskelf, doubled
%                       until no correction reaches outside a proxy surface
%     opts.occ_max    - (8192) largest occ for rskelf, errors beyond it
%
% Output:
%   funs - struct with fields
%     Afun    - @(i,j) matrix entries
%     pxyfun  - proxy function, see byindex.proxyfun
%     rx, cx  - (3, nrows) and (3, ncols) point of each row and column
%     x       - rx, if rx and cx are the same
%     Q       - sparse quadrature corrections (in Afun for rskelf)
%     xinterp - interpolation from the original to the oversampled nodes
%     proxy   - proxy surfaces for FLAM's levels, see byindex.proxyfun
%     occ     - leaf occupancy to pass to FLAM
%
% See also BYINDEX.KERNBYINDEX, BYINDEX.KERNBYINDEXEVAL, BYINDEX.PROXYFUN,
% BYINDEX.PROXY_POINTS.

if nargin < 4, targs = []; end
if nargin < 5, opts  = []; end

type = 'ifmm';
if isfield(opts, 'type'), type = lower(opts.type); end
symm = 'n';
if isfield(opts, 'symm'), symm = lower(opts.symm); end
ifoversamp = true;
if isfield(opts, 'ifoversamp'), ifoversamp = opts.ifoversamp; end
quadopts = [];
if isfield(opts, 'quadopts'), quadopts = opts.quadopts; end
proxyopts = [];
if isfield(opts, 'proxyopts'), proxyopts = opts.proxyopts; end
occ = 2048;
if isfield(opts, 'occ'), occ = opts.occ; end
occ_max = 8192;
if isfield(opts, 'occ_max'), occ_max = opts.occ_max; end

if ~any(strcmp(type, {'ifmm', 'rskelf'}))
    error('BYINDEX.FLAMFUNS: unknown type ''%s''', type);
end
if strcmp(type, 'rskelf') && ~isempty(targs)
    error('BYINDEX.FLAMFUNS: rskelf needs the on-surface matrix, targs must be []');
end
if isnumeric(targs) && ~isempty(targs), targs = struct('r', targs); end

funs = [];
if strcmp(type, 'rskelf')
    % Square matrix with the corrections in Afun
    quadopts.nonsmoothonly = true;
    quadopts.ifoversamp = 0;
    Q = surfermat(surfers, kern, eps, quadopts);
    funs.Afun    = @(i,j) byindex.kernbyindex(i, j, surfers, kern, Q);
    funs.Q       = Q;
    funs.xinterp = speye(size(Q, 2));
    srcs  = surfers;
    trgs  = [];
    rows  = surfers;
else
    % Smooth rule on the oversampled sources, corrections outside
    quadopts.corrections = true;
    quadopts.ifoversamp  = ifoversamp;
    if isempty(targs)
        % Oversample each source surfer once for all targets
        quadopts.unif_nover = true;
        [funs.Q, objover] = surfermat(surfers, kern, eps, quadopts);
        srcs = [objover{1}{1,:}];
        funs.xinterp = blkdiag(objover{2}{1,:});
        trgs = surfers;
    else
        [funs.Q, objover] = surferkernevalmat(surfers, kern, targs, eps, quadopts);
        srcs = [objover{1}{:}];
        funs.xinterp = blkdiag(objover{2}{:});
        trgs = targs;
    end
    funs.Afun = @(i,j) byindex.kernbyindexeval(i, j, srcs, kern, trgs);
    rows = trgs;
end

% Point of each row and column
nrow = numel(rows);
ncol = numel(srcs);
if numel(kern) == 1
    opdims = repmat(kern.opdims(:), 1, nrow, ncol);
else
    opdims = reshape([kern.opdims], 2, nrow, ncol);
end
rx = cell(1, nrow);
cx = cell(1, ncol);
for k = 1:nrow
    rx{k} = repelem(rows(k).r(:,:), 1, opdims(1, k, 1));
end
for k = 1:ncol
    cx{k} = repelem(srcs(k).r(:,:), 1, opdims(2, 1, k));
end
funs.rx = [rx{:}];
funs.cx = [cx{:}];
if isequal(funs.rx, funs.cx)
    funs.x = funs.rx;
elseif strcmp(type, 'rskelf')
    error(['BYINDEX.FLAMFUNS: rskelf needs the same points for rows and ' ...
        'columns, i.e. opdims(1) = opdims(2) for each surfer']);
end

% FLAM's tree: the root is the bounding box of its points, and each
% level halves the sides at least 1/sqrt(2) of the longest, so every box
% on a level has the same size and the boxes form a grid
if strcmp(type, 'rskelf'), xt = funs.x; else, xt = [funs.rx, funs.cx]; end
n     = size(xt, 2);
lo    = min(xt, [], 2);
lroot = max(xt, [], 2) - lo;
proxylvl = 3;     % FLAM only uses proxies from level 3 down

% box size on each level, and the box of each point and its occupancy,
% down to the leaves for the initial occ
lvls  = lroot;
boxes = {zeros(3, n)};
cnt   = n*ones(1, n);
while (size(lvls, 2) < proxylvl || any(cnt(end, :) > occ)) && size(lvls, 2) < 64
    l = lvls(:, end);
    ldiv = l >= max(l)/sqrt(2);
    l(ldiv) = l(ldiv)/2;
    sp = l < lroot;
    b = zeros(3, n);
    b(sp, :) = min(floor((xt(sp, :) - lo(sp))./l(sp)), round(lroot(sp)./l(sp)) - 1);
    [~, ~, ib] = unique(b.', 'rows');
    c = accumarray(ib, 1);
    lvls(:, end+1)  = l; %#ok<AGROW>
    boxes{end+1}    = b; %#ok<AGROW>
    cnt(end+1, :)   = c(ib).'; %#ok<AGROW>
end
nlvl = size(lvls, 2);

% Proxy shells stretched with the boxes. A side FLAM never split has no
% boxes beyond it, so it is stretched at least as far as the split sides.
leff = lvls;
for k = 1:nlvl
    sp = lvls(:, k) < lroot;
    if any(sp), leff(~sp, k) = max(lvls(~sp, k), max(lvls(sp, k))); end
end

% One surface per box shape, built on the coarsest level with that shape
% and scaled to the finer ones
proxyopts.ifaddtrans = strcmp(type, 'rskelf') && symm == 'n';
proxy = [];
proxy.lvls    = lvls;
proxy.leff    = leff;
proxy.isurf   = zeros(1, nlvl);
proxy.surf    = {};
proxy.surflvl = [];
for k = proxylvl:nlvl
    shape = leff(:, k)/max(leff(:, k));
    if isempty(proxy.surf) || any(abs(shape - leff(:, k-1)/max(leff(:, k-1))) > 1e-8)
        proxy.surf{end+1}    = byindex.proxy_points(kern, leff(:, k), eps, proxyopts);
        proxy.surflvl(end+1) = k;
    end
    proxy.isurf(k) = numel(proxy.surf);
end
% byindex.proxyfun rebuilds them per level, with shells to the root diagonal
proxy.ifaddtrans = proxyopts.ifaddtrans;
proxy.eps        = eps;
proxy.opts       = proxyopts;
proxy.opts.rmax  = norm(lroot);
funs.proxy = proxy;
funs.pxyfun = @(varargin) byindex.proxyfun(srcs, kern, trgs, proxy, varargin{:});

% rskelf keeps a neighbor only inside its box's proxy surface, and would
% drop corrections reaching further. FLAM's boxes are the grid above, so
% this is checked exactly.
if strcmp(type, 'rskelf') && ~isempty(proxy.surf)
    R = proxy.surf{1}.radius;

    % first level on which a point's box has a correction reaching outside
    % its proxy surface
    pbad = inf(1, n);
    for j0 = 1:4096:n
        [iq, jq] = find(Q(:, j0:min(j0+4095, n)));
        jq = jq + j0 - 1;
        for lv = proxylvl:nlvl
            b = boxes{lv};
            for s = 1:2
                if s == 1, p = jq; q = iq; else, p = iq; q = jq; end
                ctr = lo + (b(:, p) + 0.5).*lvls(:, lv);
                out = any(b(:, p) ~= b(:, q), 1) & ...
                    vecnorm((xt(:, q) - ctr)./leff(:, lv)) >= R;
                pbad(p(out)) = min(pbad(p(out)), lv);
            end
        end
    end

    % double occ until every point's leaf is above its first bad level
    [~, leafp] = max(cnt <= occ, [], 1);
    while any(pbad <= leafp)
        if 2*occ > occ_max
            error(['BYINDEX.FLAMFUNS: corrections reach outside the proxy ' ...
                'surfaces of FLAM''s boxes for occ up to occ_max = %d'], occ_max);
        end
        occ = 2*occ;
        [~, leafp] = max(cnt <= occ, [], 1);
    end
end
funs.occ = occ;

end
