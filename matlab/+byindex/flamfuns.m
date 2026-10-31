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
%   F = ifmm(funs.Afun, funs.rx, funs.cx, occ, eps, funs.pxyfun);
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
%   F = rskelf(funs.Afun, funs.x, occ, eps, funs.pxyfun);
%
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
%
% Output:
%   funs - struct with fields
%     Afun    - @(i,j) matrix entries
%     pxyfun  - proxy function, see byindex.proxyfun
%     rx, cx  - (3, nrows) and (3, ncols) point of each row and column
%     x       - rx, if rx and cx are the same
%     Q       - sparse quadrature corrections (in Afun for rskelf)
%     xinterp - interpolation from the original to the oversampled nodes
%     proxy   - proxy surface from byindex.proxy_points
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
    proxylvl = 3;     % rskelf only uses proxies from level 3 down
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
    proxylvl = 1;
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

% Root box, as in FLAM's hypoct
xall  = [funs.rx, funs.cx];
width = max(max(xall, [], 2) - min(xall, [], 2));

proxyopts.ifaddtrans = strcmp(type, 'rskelf') && symm == 'n';
funs.proxy = byindex.proxy_points(kern, width/2^(proxylvl-1), eps, proxyopts);
% byindex.proxyfun rebuilds the proxy surface below proxylvl from the root
funs.proxy.width = width;
funs.proxy.lvl   = proxylvl;

proxy = funs.proxy;
funs.pxyfun = @(varargin) byindex.proxyfun(srcs, kern, trgs, proxy, varargin{:});

end
