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
% The last form is used by the fastalgorithms fork of FLAM when
% opts.proxybylevel = true: the proxy surface is rebuilt by
% byindex.proxy_points for a box of size proxy.width/2^(lvl-1). Other
% versions of FLAM ignore opts.proxybylevel and use proxy at every level.
%
% Matrix indices follow byindex.kernbyindex (targs = []) or
% byindex.kernbyindexeval (off-surface targs).
%
% Input:
%   surfers - surfer object or array of surfer objects. Sources, and also
%             targets if targs is empty.
%   kern    - kernel3d object. If targs is empty, scalar or
%             nsurfers x nsurfers; otherwise scalar or vector of length
%             nsurfers.
%   targs   - [] for on-surface matrices, or off-surface targets: struct
%             with .r (3 x ntarg) and any target fields, or numeric
%             (3 x ntarg) array
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
%   Kpxy - interactions between slf and the proxy points, with one copy of
%          the proxy points for each block row (rc = 'c') or block column
%          (rc = 'r') of the kernel array:
%            rc = 'c', rskelf: (npxy*sum opdims(1)) x length(slf), sources
%                              in slf to proxy targets. If
%                              proxy.ifaddtrans, rskelf also gets the
%                              transpose of the rc = 'r' block appended.
%            rc = 'r':         length(slf) x (npxy*sum opdims(2)), proxy
%                              sources (weighted by proxy.w) to targets
%                              in slf
%          For lvl, the proxy function for that level.
%   nbr  - nbr without the indices outside the proxy surface

% Proxy function for one level of the tree
if numel(varargin) == 1
    lvl = varargin{1};
    p = byindex.proxy_points(kern, proxy.width/2^(lvl-1), proxy.eps, ...
        proxy.nshell_max, proxy.opts);
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
ifon = isempty(targs);
if ifon
    nrow = nsurfers;
else
    nrow = 1;
    if isnumeric(targs), targs = struct('r', targs); end
end

% opdims_mat(:, irow, isrc), irow over block rows (one if off-surface)
if numel(kern) == 1
    opdims_mat = repmat(kern.opdims(:), 1, nrow, nsurfers);
else
    opdims_mat = reshape([kern.opdims], 2, nrow, nsurfers);
end

% Row and column offsets
irowlocs = ones(nrow+1, 1);
icollocs = ones(nsurfers+1, 1);
for k = 1:nsurfers
    icollocs(k+1) = icollocs(k) + surfers(k).npts * opdims_mat(2, 1, k);
end
if ifon
    for k = 1:nrow
        irowlocs(k+1) = irowlocs(k) + surfers(k).npts * opdims_mat(1, k, 1);
    end
else
    irowlocs(2) = 1 + size(targs.r, 2) * opdims_mat(1, 1, 1);
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

% Sources in slf -> proxy targets
if rc == 'c'
    prow = ones(nrow+1, 1);
    for k = 1:nrow
        prow(k+1) = prow(k) + npxy * opdims_mat(1, k, 1);
    end
    Kpxy = zeros(prow(end)-1, length(slf));

    for isrc = 1:nsurfers
        srfj = surfers(isrc);
        f_col = slf >= icollocs(isrc) & slf < icollocs(isrc+1);
        if ~any(f_col), continue; end
        mat_col_ind = slf(f_col);

        opdims_src = opdims_mat(2, 1, isrc);
        jpts = idivide(int64(mat_col_ind(:) - icollocs(isrc)), int64(opdims_src)) + 1;
        [juni, ~, ijuni] = unique(jpts);
        ijuni2 = (ijuni - 1)*opdims_src + mod(mat_col_ind(:) - icollocs(isrc), opdims_src) + 1;

        wsrc = repmat(srfj.wts(juni(:)).', opdims_src, 1);
        wsrc = wsrc(:).';

        for irow = 1:nrow
            if numel(kern) == 1
                ktmp = kern;
            elseif ifon
                ktmp = kern(irow, isrc);
            else
                ktmp = kern(isrc);
            end
            srcp = slice_surfer(srfj, juni, ktmp.src_fields);
            matuni = ktmp.eval(srcp, pxy) .* wsrc;
            Kpxy(prow(irow):prow(irow+1)-1, f_col) = matuni(:, ijuni2);
        end
    end
end

% Proxy sources -> targets in slf
if rc == 'r' || ifaddtrans
    pcol = ones(nsurfers+1, 1);
    for k = 1:nsurfers
        pcol(k+1) = pcol(k) + npxy * opdims_mat(2, 1, k);
    end
    Ktrg = zeros(length(slf), pcol(end)-1);

    for irow = 1:nrow
        f_row = slf >= irowlocs(irow) & slf < irowlocs(irow+1);
        if ~any(f_row), continue; end
        mat_row_ind = slf(f_row);

        opdims_row = opdims_mat(1, irow, 1);
        ipts = idivide(int64(mat_row_ind(:) - irowlocs(irow)), int64(opdims_row)) + 1;
        [iuni, ~, iiuni] = unique(ipts);
        iiuni2 = (iiuni - 1)*opdims_row + mod(mat_row_ind(:) - irowlocs(irow), opdims_row) + 1;

        for isrc = 1:nsurfers
            if numel(kern) == 1
                ktmp = kern;
            elseif ifon
                ktmp = kern(irow, isrc);
            else
                ktmp = kern(isrc);
            end
            if ifon
                targinfo = slice_surfer(surfers(irow), iuni, ktmp.targ_fields);
            else
                targinfo = slice_surfer(targs, iuni, ktmp.targ_fields);
            end

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


function info = slice_surfer(obj, pts, fields)
% Point info at pts of a surfer or target struct. Fields may be stored
% as (:, npts) or (npts, :).
    npts = size(obj.r, 2);
    info = [];
    info.r = obj.r(:, pts);
    for k = 1:length(fields)
        f = fields{k};
        A = obj.(f);
        if size(A, 2) == npts
            info.(f) = A(:, pts);
        elseif size(A, 1) == npts
            info.(f) = A(pts, :).';
        else
            error('BYINDEX.PROXYFUN: no dimension of field %s matches npts', f);
        end
    end
end
