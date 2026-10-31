function mat = kernbyindexeval(i, j, surfers, kern, targobj, Qsparse, opts)
%KERNBYINDEXEVAL  Evaluate smooth-rule off-surface matrix entries by index.
%
% mat = byindex.kernbyindexeval(i, j, surfers, kern, targobj)
% mat = byindex.kernbyindexeval(i, j, surfers, kern, targobj, Qsparse)
% mat = byindex.kernbyindexeval(i, j, surfers, kern, targobj, Qsparse, opts)
%
% Analogous to byindex.kernbyindex but targets are the points in targobj.
% Entries use the smooth rule at the nodes of the surfers, without
% oversampling; to oversample, pass the oversampled surfers and apply the
% interpolation outside, see byindex.flamfuns. Coincident source and
% target points give zero entries.
%
% Input:
%   i       - row indices into the flat target DOF array
%   j       - col indices into the source DOF array
%   surfers - surfer object or array of surfer objects (sources)
%   kern    - kernel3d object, scalar or ntargs x nsurfers
%   targobj - targets: numeric (3 x ntarg) array, or a struct or surfer
%             with .r (3 x ntarg). An array of structs or surfers gives
%             one block row per element.
%   Qsparse - (optional) sparse correction matrix
%   opts    - options struct
%     opts.replace_quadcorr - (default true) as in byindex.kernbyindex
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

if isnumeric(targobj), targobj = struct('r', targobj); end
ntargs   = numel(targobj);
nsurfers = length(surfers);

if numel(kern) == 1
    opdims_mat = repmat(kern.opdims(:), 1, ntargs, nsurfers);
else
    kern = reshape(kern, ntargs, nsurfers);
    opdims_mat = reshape([kern.opdims], 2, ntargs, nsurfers);
end

assert(all(opdims_mat(1,:,:) == opdims_mat(1,:,1), 'all'), ...
    'BYINDEX.KERNBYINDEXEVAL: opdims(1) must be constant across source surfers');

irowlocs = ones(ntargs+1, 1);
for k = 1:ntargs
    irowlocs(k+1) = irowlocs(k) + size(targobj(k).r(:,:), 2) * opdims_mat(1,k,1);
end
icollocs = ones(nsurfers+1, 1);
for k = 1:nsurfers
    icollocs(k+1) = icollocs(k) + surfers(k).npts * opdims_mat(2,1,k);
end

mat = zeros(length(i), length(j));
if isempty(i) || isempty(j), return; end
i = double(i(:));
j = double(j(:));

for itrg = 1:ntargs
    flag_i = (i >= irowlocs(itrg)) & (i < irowlocs(itrg+1));
    if ~any(flag_i), continue; end
    [iuni, iiuni2] = decode(i(flag_i) - irowlocs(itrg), opdims_mat(1, itrg, 1));

    % Only the target fields that targobj has
    tobj = targobj(itrg);
    fields = unique([{}, kern(min(itrg, end), :).targ_fields]);
    fields = fields(cellfun(@(f) isfield(tobj, f) || isprop(tobj, f), fields));
    targinfo = ptinfo(tobj, iuni, fields);

    for isrc = 1:nsurfers
        srfj = surfers(isrc);

        flag_j = (j >= icollocs(isrc)) & (j < icollocs(isrc+1));
        if ~any(flag_j), continue; end
        opdims_src = opdims_mat(2, itrg, isrc);
        [juni, ijuni2] = decode(j(flag_j) - icollocs(isrc), opdims_src);

        ktmp = kern(min(itrg, end), min(isrc, end));

        srcp   = ptinfo(srfj, juni, ktmp.src_fields);
        wts    = repmat(srfj.wts(juni).', opdims_src, 1);
        matuni = ktmp.eval_mask(srcp, targinfo) .* wts(:).';

        mat(flag_i, flag_j) = matuni(iiuni2, ijuni2);
    end
end

if ~isempty(Qsparse)
    [isp, jsp, vsp] = find(Qsparse(i, j));
    linsp = isp + (jsp-1)*length(i);
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
