function [Kpxy, nbr] = proxyfun(surfers, kern, targs, proxy, varargin)
%PROXYFUN  FLAM proxy function for kernels defined on arrays of surfers.
%
% pxyfun = @(varargin) byindex.proxyfun(surfers, kern, targs, proxy, varargin{:});
%
% The FLAM arguments are passed through varargin, so the same handle works
% for every FLAM routine:
%
%   [Kpxy, nbr] = pxyfun(x, slf, nbr, l, ctr)            % rskelf
%   [Kpxy, nbr] = pxyfun(rc, rx, cx, slf, nbr, l, ctr)    % rskel, ifmm
%   pxyfunlvl   = pxyfun(lvl)                             % opts.proxybylevel
%
% The last form is used by the fastalgorithms fork of FLAM with
% opts.proxybylevel = true, and rebuilds the proxy surface for a box of
% size proxy.width/2^(lvl-1). For a surface built for a smaller box, set
% proxy.lvl to its level (see byindex.flamfuns).
%
% Matrix indices follow byindex.kernbyindex (targs = []) or
% byindex.kernbyindexeval (off-surface targs).
%
% Input:
%   surfers - surfer object or array of surfer objects. Sources, and also
%             targets if targs is empty.
%   kern    - kernel3d object. If targs is empty, scalar or
%             nsurfers x nsurfers; otherwise as in byindex.kernbyindexeval.
%   targs   - [] for on-surface matrices, or off-surface targets as in
%             byindex.kernbyindexeval
%   proxy   - proxy surface from byindex.proxy_points
%
%   FLAM arguments:
%   x, rx, cx - point for each row/column index
%   rc        - 'c' to compress columns (sources), 'r' to compress rows
%               (targets)
%   slf       - indices in the current box
%   nbr       - indices in the neighboring boxes
%   l         - box size, scalar or one size per dimension
%   ctr       - box center
%   lvl       - tree level
%
% Output:
%   Kpxy - interactions of slf with the proxy points, one copy of them per
%          block row (rc = 'c') or column (rc = 'r') of kern, or one copy
%          for a single kernel:
%            rc = 'c': slf sources to proxy targets
%            rc = 'r': proxy sources (weighted by proxy.w) to slf targets
%          For lvl, the proxy function for that level.
%   nbr  - nbr without the indices outside the proxy surface

% Proxy function for level lvl. Levels up to proxy.lvl (default 1) keep
% the given surface, finer levels rebuild it.
if numel(varargin) == 1
    lvl = varargin{1};
    p = proxy;
    lvl0 = 1;
    if isfield(proxy, 'lvl'), lvl0 = proxy.lvl; end
    if lvl > lvl0
        p = byindex.proxy_points(kern, proxy.width/2^(lvl-1), proxy.eps, proxy.opts);
    end
    Kpxy = @(varargin) byindex.proxyfun(surfers, kern, targs, p, varargin{:});
    nbr = [];
    return
end

ifaddtrans = false;
if numel(varargin) == 5
    [rx, slf, nbr, l, ctr] = varargin{:};
    cx = rx;
    rc = 'c';
    ifaddtrans = proxy.ifaddtrans;
elseif numel(varargin) == 7
    [rc, rx, cx, slf, nbr, l, ctr] = varargin{:};
else
    error('BYINDEX.PROXYFUN: expected 1, 5 or 7 FLAM arguments, got %d', ...
        numel(varargin));
end

nsurfers = length(surfers);
if isempty(targs)
    targs = surfers;
elseif isnumeric(targs)
    targs = struct('r', targs);
end
nrow = numel(targs);

% opdims_mat(:, irow, isrc), irow over block rows
if numel(kern) == 1
    opdims_mat = repmat(kern.opdims(:), 1, nrow, nsurfers);
else
    kern = reshape(kern, nrow, nsurfers);
    opdims_mat = reshape([kern.opdims], 2, nrow, nsurfers);
end

% Row and column offsets
irowlocs = ones(nrow+1, 1);
icollocs = ones(nsurfers+1, 1);
for k = 1:nsurfers
    icollocs(k+1) = icollocs(k) + surfers(k).npts * opdims_mat(2, 1, k);
end
for k = 1:nrow
    irowlocs(k+1) = irowlocs(k) + size(targs(k).r(:,:), 2) * opdims_mat(1, k, 1);
end

% Scale proxy points/weights
lmax = max(l);
pxy = [];
pxy.r = proxy.r*lmax + ctr(:);
pxy.n = proxy.n;
pw = proxy.w*lmax^2;
npxy = size(pxy.r, 2);

% Keep only the neighbors inside the proxy surface. For rc = 'c' these
% are rows, for rc = 'r' columns.
if rc == 'c'
    xnbr = rx(:, nbr);
else
    xnbr = cx(:, nbr);
end
nbr = nbr(vecnorm((xnbr - ctr(:))/lmax, 2, 1) < proxy.radius);

% One copy of the proxy points per block row (rc = 'c') or block column
% (rc = 'r'), or a single copy for a single kernel
if numel(kern) == 1
    npr = 1;
    npc = 1;
else
    npr = nrow;
    npc = nsurfers;
end

% Sources in slf -> proxy targets
if rc == 'c'
    prow = ones(npr+1, 1);
    for k = 1:npr
        prow(k+1) = prow(k) + npxy * opdims_mat(1, k, 1);
    end
    Kpxy = zeros(prow(end)-1, length(slf));

    for isrc = 1:nsurfers
        srfj = surfers(isrc);
        f_col = slf >= icollocs(isrc) & slf < icollocs(isrc+1);
        if ~any(f_col), continue; end
        opdims_src = opdims_mat(2, 1, isrc);
        [juni, ijuni2] = decode(double(slf(f_col)) - icollocs(isrc), opdims_src);

        wsrc = repmat(srfj.wts(juni).', opdims_src, 1);
        wsrc = wsrc(:).';

        srcp = ptinfo(srfj, juni, unique([{}, kern(:, min(isrc, end)).src_fields]));

        for irow = 1:npr
            ktmp = kern(min(irow, end), min(isrc, end));
            matuni = ktmp.eval(srcp, pxy) .* wsrc;
            Kpxy(prow(irow):prow(irow+1)-1, f_col) = matuni(:, ijuni2);
        end
    end
end

% Proxy sources -> targets in slf
if rc == 'r' || ifaddtrans
    pcol = ones(npc+1, 1);
    for k = 1:npc
        pcol(k+1) = pcol(k) + npxy * opdims_mat(2, 1, k);
    end
    Ktrg = zeros(length(slf), pcol(end)-1);

    for irow = 1:nrow
        f_row = slf >= irowlocs(irow) & slf < irowlocs(irow+1);
        if ~any(f_row), continue; end
        [iuni, iiuni2] = decode(double(slf(f_row)) - irowlocs(irow), ...
            opdims_mat(1, irow, 1));

        targinfo = ptinfo(targs(irow), iuni, unique([{}, kern(min(irow, end), :).targ_fields]));

        for isrc = 1:npc
            ktmp = kern(min(irow, end), min(isrc, end));

            opdims_src = opdims_mat(2, irow, isrc);
            wpxy = repmat(pw(:).', opdims_src, 1);
            wpxy = wpxy(:).';

            matuni = ktmp.eval(pxy, targinfo) .* wpxy;
            Ktrg(f_row, pcol(isrc):pcol(isrc+1)-1) = matuni(iiuni2, :);
        end
    end

    if rc == 'r'
        Kpxy = Ktrg;
    else
        Kpxy = [Kpxy; Ktrg.'];
    end
end

end


function [uni, ind] = decode(k, opdim)
% Zero-based dof indices k -> unique points uni, and the row/col of each
% dof in the matrix over (uni, components)
k = k(:);
pts = floor(k/opdim) + 1;
[uni, ~, iuni] = unique(pts);
ind = (iuni - 1)*opdim + (k - (pts - 1)*opdim) + 1;
end


function info = ptinfo(obj, pts, fields)
% Point info at pts. Fields may be stored as (:, npts) or (npts, :).
npts = size(obj.r(:,:), 2);
info = [];
info.r = obj.r(:, pts);
for k = 1:numel(fields)
    f = fields{k};
    A = obj.(f);
    if size(A, 2) == npts
        info.(f) = A(:, pts);
    elseif size(A, 1) == npts
        info.(f) = A(pts, :).';
    else
        error('BYINDEX: no dimension of field %s matches npts', f);
    end
end
end
