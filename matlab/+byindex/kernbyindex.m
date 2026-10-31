function mat = kernbyindex(i, j, surfers, kern, Qsparse, opts)
%KERNBYINDEX  Evaluate smooth-rule matrix entries by index for surfer arrays.
%
% mat = byindex.kernbyindex(i, j, surfers, kern)
% mat = byindex.kernbyindex(i, j, surfers, kern, Qsparse)
% mat = byindex.kernbyindex(i, j, surfers, kern, Qsparse, opts)
%
% Entries use the smooth rule at the nodes of the surfers, without
% oversampling, so with Qsparse this matches surfermat with
% opts.ifoversamp = 0. Self interactions are zero and kern.diag is not
% used. For oversampled sources, see byindex.kernbyindexeval.
%
% For use with FLAM, e.g.
%   Afun   = @(i,j) byindex.kernbyindex(i, j, surfers, kern, Q);
%   pxyfun = @(varargin) byindex.proxyfun(surfers, kern, [], proxy, varargin{:});
%   F = rskelf(Afun, x, occ, rank_or_tol, pxyfun);
% or see byindex.flamfuns.
%
% Input:
%   i, j      - row/col indices to evaluate
%   surfers   - surfer object or array of surfer objects
%   kern      - kernel3d object (scalar or nsurfers x nsurfers)
%   Qsparse   - (optional) sparse matrix of quadrature corrections
%               (same size as the full matrix).  Pass [] to skip.
%               Two modes, selected by opts.replace_quadcorr:
%                 true (default): replacement, i.e. overwrite smooth entries
%                   with Qsparse values where nonzero.  Use with nonsmoothonly=1
%                   output from surfermat.
%                 false: additive, i.e. mat += Qsparse.  Use with corrections=1,
%                   ifoversamp=0 output from surfermat.
%   opts      - options struct
%     opts.replace_quadcorr - logical (default true), see Qsparse above
%
% Output:
%   mat - length(i) x length(j) matrix of smooth-rule entries (with singular
%         corrections applied where Qsparse is provided)

if ~isa(kern, 'kernel3d')
    error('BYINDEX.KERNBYINDEX: kern must be a kernel3d object');
end

if nargin < 5, Qsparse = []; end
if nargin < 6, opts    = []; end

if ~isempty(Qsparse) && ~issparse(Qsparse)
    error(['BYINDEX.KERNBYINDEX: Qsparse must be a sparse matrix (eps and ' ...
        'objover are no longer inputs)']);
end

replace_quadcorr = true;
if isstruct(opts) && isfield(opts, 'replace_quadcorr')
    replace_quadcorr = opts.replace_quadcorr;
end

nsurfers = length(surfers);

% Build opdims_mat (2 x nsurfers x nsurfers) from kern
if numel(kern) == 1
    opdims_mat = repmat(reshape(kern.opdims, 2, 1, 1), 1, nsurfers, nsurfers);
else
    opdims_mat = reshape([kern.opdims], 2, nsurfers, nsurfers);
end

% Build row/col offset tables
irowlocs = ones(nsurfers+1, 1);
icollocs = ones(nsurfers+1, 1);
for k = 1:nsurfers
    irowlocs(k+1) = irowlocs(k) + surfers(k).npts * opdims_mat(1, k, 1);
    icollocs(k+1) = icollocs(k) + surfers(k).npts * opdims_mat(2, 1, k);
end

mat = zeros(length(i), length(j));
if isempty(i) || isempty(j)
    return;
end
i = double(i(:));
j = double(j(:));

for itrg = 1:nsurfers
    srfi = surfers(itrg);

    flag_i = (i >= irowlocs(itrg)) & (i < irowlocs(itrg+1));
    if ~any(flag_i), continue; end
    [iuni, iiuni2] = decode(i(flag_i) - irowlocs(itrg), opdims_mat(1, itrg, 1));

    targinfo = ptinfo(srfi, iuni, unique([{}, kern(min(itrg, end), :).targ_fields]));

    for isrc = 1:nsurfers
        srfj = surfers(isrc);

        flag_j = (j >= icollocs(isrc)) & (j < icollocs(isrc+1));
        if ~any(flag_j), continue; end
        opdims_src = opdims_mat(2, itrg, isrc);
        [juni, ijuni2] = decode(j(flag_j) - icollocs(isrc), opdims_src);

        ktmp = kern(min(itrg, end), min(isrc, end));

        srcp = ptinfo(srfj, juni, ktmp.src_fields);
        wts = repmat(srfj.wts(juni).', opdims_src, 1);
        matuni = ktmp.eval_mask(srcp, targinfo) .* wts(:).';

        mat(flag_i, flag_j) = matuni(iiuni2, ijuni2);
    end
end

% Apply singular quadrature corrections
if ~isempty(Qsparse)
    [isp, jsp, vsp] = find(Qsparse(i, j));
    linsp = isp + (jsp - 1)*length(i);
    if replace_quadcorr
        mat(linsp) = vsp;
    else
        mat(linsp) = mat(linsp) + vsp;
    end
end

end


function [uni, ind] = decode(k, opdim)
% Zero-based dof indices k -> unique points uni, and the row/col of each
% dof in the matrix over (uni, components)
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
