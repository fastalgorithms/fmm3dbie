function mat = kernbyindex(i, j, surfers, kern, Qsparse, opts)
%KERNBYINDEX  Evaluate smooth-rule matrix entries by index for surfer arrays.
%
% mat = byindex.kernbyindex(i, j, surfers, kern)
% mat = byindex.kernbyindex(i, j, surfers, kern, Qsparse)
% mat = byindex.kernbyindex(i, j, surfers, kern, Qsparse, opts)
%
% Entries use the smooth rule at the nodes of the surfers, without
% oversampling, so together with Qsparse this matches surfermat with
% opts.ifoversamp = 0. Self interactions (same surfer, same node) are
% zero. Matrices that need oversampled sources can be built from
% byindex.kernbyindexeval with the oversampled surfers, applying the
% interpolation matrices outside.
%
% For use with FLAM, e.g.
%   Afun   = @(i,j) byindex.kernbyindex(i, j, surfers, kern, Q);
%   pxyfun = @(varargin) byindex.proxyfun(surfers, kern, [], proxy, varargin{:});
%   F = rskelf(Afun, x, occ, rank_or_tol, pxyfun);
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
%                 false: additive, i.e. mat += Qsparse.  Use with corrections=1
%                   output from surfermat, which stores (near-singular - smooth).
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
irowlocs = zeros(nsurfers+1, 1);
icollocs = zeros(nsurfers+1, 1);
irowlocs(1) = 1;
icollocs(1) = 1;
for k = 1:nsurfers
    irowlocs(k+1) = irowlocs(k) + surfers(k).npts * opdims_mat(1, k, 1);
    icollocs(k+1) = icollocs(k) + surfers(k).npts * opdims_mat(2, 1, k);
end

mat = zeros(length(i), length(j));
if isempty(i) || isempty(j)
    return;
end

for itrg = 1:nsurfers
    srfi = surfers(itrg);

    mat_row_start = irowlocs(itrg);
    mat_row_end   = irowlocs(itrg+1) - 1;
    flag_i = (i >= mat_row_start) & (i <= mat_row_end);
    mat_row_ind = i(flag_i);
    if isempty(mat_row_ind), continue; end

    opdims_row = opdims_mat(1, itrg, 1);

    ipts = idivide(int64(mat_row_ind(:) - mat_row_start), int64(opdims_row)) + 1;
    [iuni, ~, iiuni] = unique(ipts);
    iiuni2 = (iiuni - 1)*opdims_row + mod(mat_row_ind(:) - mat_row_start, opdims_row) + 1;

    for isrc = 1:nsurfers
        srfj = surfers(isrc);

        mat_col_start = icollocs(isrc);
        mat_col_end   = icollocs(isrc+1) - 1;
        flag_j = (j >= mat_col_start) & (j <= mat_col_end);
        mat_col_ind = j(flag_j);
        if isempty(mat_col_ind), continue; end

        opdims_src = opdims_mat(2, itrg, isrc);

        if numel(kern) == 1
            ktmp = kern;
        else
            ktmp = kern(itrg, isrc);
        end

        jpts = idivide(int64(mat_col_ind(:) - mat_col_start), int64(opdims_src)) + 1;
        [juni, ~, ijuni] = unique(jpts);
        ijuni2 = (ijuni - 1)*opdims_src + mod(mat_col_ind(:) - mat_col_start, opdims_src) + 1;

        targinfo = slice_surfer(srfi, iuni, ktmp.targ_fields);
        srcp     = slice_surfer(srfj, juni, ktmp.src_fields);
        wts = repmat(srfj.wts(juni(:)).', opdims_src, 1);
        wts = wts(:).';
        matuni = ktmp.eval(srcp, targinfo) .* wts;

        % Zero self interactions, i.e. source and target at the same node
        if itrg == isrc
            [lia, loc] = ismember(juni, iuni);
            if any(lia)
                q = find(lia);
                p = loc(lia);
                ir = repmat((1:opdims_row).' + opdims_row*(p(:).' - 1), opdims_src, 1);
                ic = repelem((1:opdims_src).' + opdims_src*(q(:).' - 1), opdims_row, 1);
                matuni(ir(:) + size(matuni, 1)*(ic(:) - 1)) = 0;
            end
        end

        mat(flag_i, flag_j) = matuni(iiuni2, ijuni2);
    end
end

% Apply singular quadrature corrections
if ~isempty(Qsparse)
    [isp, jsp, vsp] = find(Qsparse(i(:), j(:)));
    linsp = isp + (jsp - 1)*length(i(:));
    if replace_quadcorr
        % nonsmoothonly mode: Qsparse holds full near-singular values;
        % replace smooth entries wherever Qsparse is structurally nonzero.
        mat(linsp) = vsp;
    else
        % corrections mode: Qsparse holds (near-singular - smooth);
        % add to smooth entries.
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
