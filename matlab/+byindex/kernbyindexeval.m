function mat = kernbyindexeval(i, j, surfers, kern, targobj, Qsparse, opts)
%KERNBYINDEXEVAL  Evaluate smooth-rule off-surface matrix entries by index.
%
% mat = byindex.kernbyindexeval(i, j, surfers, kern, targobj)
% mat = byindex.kernbyindexeval(i, j, surfers, kern, targobj, Qsparse)
% mat = byindex.kernbyindexeval(i, j, surfers, kern, targobj, Qsparse, opts)
%
% Analogous to byindex.kernbyindex but targets are the points in targobj.
% Entries use the smooth rule at the nodes of the surfers, without
% oversampling. To oversample, pass the oversampled surfers and apply the
% interpolation matrices outside, e.g. with surferkernevalmat:
%
%   [Qc, objover] = surferkernevalmat(surfers, kern, targs, eps, ...
%                                     struct('corrections', 1));
%   Sover   = [objover{1}{:}];
%   xinterp = blkdiag(objover{2}{:});
%   Afun    = @(i,j) byindex.kernbyindexeval(i, j, Sover, kern, targs);
%   pxyfun  = @(varargin) byindex.proxyfun(Sover, kern, targs, proxy, varargin{:});
%   F = ifmm(Afun, rx, cx, occ, rank_or_tol, pxyfun);
%   u = ifmm_mv(F, xinterp*sigma, Afun) + Qc*sigma;
%
% Non-finite entries, from targets that coincide with sources (e.g. the
% original nodes as targets of oversampled sources), are set to zero.
%
% Input:
%   i       - row indices into the flat target DOF array (length ntarg*opdim1)
%   j       - col indices into the source DOF array
%   surfers - surfer object or array of surfer objects (sources)
%   kern    - kernel3d object, scalar or vector of length nsurfers
%   targobj - targets: struct or surfer with .r (3 x ntarg)
%             or numeric (3 x ntarg) array
%   Qsparse - (optional) sparse correction matrix (ntarg*opdim1 x ncols)
%   opts    - options struct
%     opts.replace_quadcorr - logical (default true). If false, Qsparse
%                             is added to the smooth entries, e.g. for
%                             corrections=1 output from surferkernevalmat
%
% Output:
%   mat - length(i) x length(j) submatrix

if ~isa(kern, 'kernel3d')
    error('BYINDEX.KERNBYINDEXEVAL: kern must be a kernel3d object');
end
if nargin < 6, Qsparse = []; end
if nargin < 7, opts    = []; end

if ~isempty(Qsparse) && ~issparse(Qsparse)
    error(['BYINDEX.KERNBYINDEXEVAL: Qsparse must be a sparse matrix (eps ' ...
        'and objover are no longer inputs)']);
end

replace_quadcorr = true;
if isstruct(opts) && isfield(opts, 'replace_quadcorr')
    replace_quadcorr = opts.replace_quadcorr;
end

if isnumeric(targobj)
    targinfo_full.r = targobj;
else
    targinfo_full = targobj;
end

nsurfers = length(surfers);

if numel(kern) == 1
    opdims_mat = repmat(kern.opdims(:), 1, nsurfers);
else
    opdims_mat = reshape([kern.opdims], 2, nsurfers);
end

assert(all(opdims_mat(1,:) == opdims_mat(1,1)), ...
    'BYINDEX.KERNBYINDEXEVAL: opdims(1) must be constant across all source surfers');

opdim1 = opdims_mat(1,1);

icollocs = zeros(nsurfers+1, 1); icollocs(1) = 1;
for k = 1:nsurfers
    icollocs(k+1) = icollocs(k) + surfers(k).npts * opdims_mat(2,k);
end

% Decode row indices -> unique target points
ntarg = size(targinfo_full.r, 2);
ipts  = double(idivide(int64(i(:) - 1), int64(opdim1))) + 1;
[iuni, ~, iiuni] = unique(ipts);
iiuni2 = (iiuni - 1)*opdim1 + mod(i(:) - 1, opdim1) + 1;

% Target fields required by any of the kernels
targ_fields = {};
for k = 1:numel(kern)
    if ~isempty(kern(k).targ_fields)
        targ_fields = union(targ_fields, kern(k).targ_fields);
    end
end

targinfo = [];
targinfo.r = targinfo_full.r(:, iuni);
for k = 1:numel(targ_fields)
    f = targ_fields{k};
    if isfield(targinfo_full, f) || isprop(targinfo_full, f)
        targinfo.(f) = slice_field(targinfo_full.(f), iuni, ntarg);
    end
end

mat = zeros(length(i), length(j));
if isempty(i) || isempty(j), return; end

for isrc = 1:nsurfers
    srfj = surfers(isrc);

    mat_col_start = icollocs(isrc);
    mat_col_end   = icollocs(isrc+1) - 1;
    flag_j = (j >= mat_col_start) & (j <= mat_col_end);
    mat_col_ind = j(flag_j);
    if isempty(mat_col_ind), continue; end

    if numel(kern) == 1, ktmp = kern; else, ktmp = kern(isrc); end
    opdims_src = opdims_mat(2, isrc);

    jpts = idivide(int64(mat_col_ind(:) - mat_col_start), int64(opdims_src)) + 1;
    [juni, ~, ijuni] = unique(jpts);
    ijuni2 = (ijuni - 1)*opdims_src + mod(mat_col_ind(:) - mat_col_start, opdims_src) + 1;

    srcp   = slice_surfer(srfj, juni, ktmp.src_fields);
    wts    = repmat(srfj.wts(juni(:)).', opdims_src, 1);
    wts    = wts(:).';
    matuni = ktmp.eval(srcp, targinfo) .* wts;
    matuni(~isfinite(matuni)) = 0;

    mat(:, flag_j) = matuni(iiuni2, ijuni2);
end

if ~isempty(Qsparse)
    [isp, jsp, vsp] = find(Qsparse(i(:), j(:)));
    linsp = isp + (jsp-1)*length(i(:));
    if replace_quadcorr
        mat(linsp) = vsp;
    else
        mat(linsp) = mat(linsp) + vsp;
    end
end

end


function srcp = slice_surfer(srfj, pts, fields)
    srcp = [];
    srcp.r = srfj.r(:, pts);
    for k = 1:length(fields)
        f = fields{k};
        srcp.(f) = slice_field(srfj.(f), pts, srfj.npts);
    end
end


function v = slice_field(A, idx, npts)
if size(A, 2) == npts
    v = A(:, idx);
elseif size(A, 1) == npts
    v = A(idx, :).';
else
    error('SLICE_FIELD: no dimension of A matches npts');
end
end
