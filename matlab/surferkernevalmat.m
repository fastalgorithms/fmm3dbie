function [sysmat,novers] = surferkernevalmat(surferobj,kern, targobj, eps,opts)
%SURFERKERNEVALMAT build evaluation matrix for given kernel, surfer description
% of boundary, and off-surface target points. This is a wrapper for various
% quadrature routines. Optionally, return only those interactions which do
% not use the smooth integration rule in the sparse matrix format.
%
% Syntax: sysmat = surferkernevalmat(S,kern,targinfo,opts)
%
% Input:
%   Surferobj - array of surfer objects describing boundary
%   kern  - kernel function. By default, this should be a function handle
%           accepting input of the form kern(srcinfo,targinfo), where srcinfo
%           and targinfo are in the ptinfo struct format, i.e.
%                ptinfo.r - positions (3,:) array
%                ptinfo.n - unit normals (3,:)
%   targinfo - target point info struct
%                targinfo.r - positions (3,:) array
%                targinfo.n - unit normals (3,:) (if needed by kernel)
%
% Optional input:
%   opts  - options structure. available options (default settings)
%
%           opts.nonsmoothonly = boolean (false), if true, only compute the
%                         entries for which a special quadrature is used
%                         (e.g. self and neighbor interactions) and return
%                         in a sparse array.
%           opts.corrections = boolean (false), if true, only compute the
%                         corrections to the smooth quadrature rule and
%                         return in a sparse array, see opts.nonsmoothonly
%           opts.quad = quad type 'kern' to call kern.get_quad, or 'native'
%               to use the smooth rule
%
% Output:
%   sysmat - the evaluation matrix mapping densities on the boundary
%            defined by S to values at targinfo
%
%
% Examples:
%   sysmat = surferevalmat(S,kern,targinfo); % standard options
%   sysmat = surferevalmat(S,kern,targinfo,opts);
%
% Author: Tristan Goodwill

if iscell(surferobj)
    surfers = surferobj;
else
    surfers = num2cell(surferobj(:).',1);
end

if ~isa(kern,'kernel3d')
        error('SURFEREVALMAT: second input kern not of supported type');
end

if nargin < 5, opts = []; end

% Assign appropriate object to targinfo
if isnumeric(targobj)
    targinfo = targobj.r;
else
    targinfo = targobj;
end

% get opts from struct if available
nonsmoothonly = false;
if isfield(opts,'nonsmoothonly'), nonsmoothonly = opts.nonsmoothonly; end
corrections = false;
if isfield(opts,'corrections'), corrections = opts.corrections; end
if corrections, nonsmoothonly = true; end
adaptive_correction = true;
if isfield(opts,'adaptive_correction'), adaptive_correction = opts.adaptive_correction; end

nsurfers = length(surfers);

opdims_mat = zeros(2,nsurfers);
lsurfer    = zeros(nsurfers,1);


novers = cell(nsurfers, 1);

for j=1:nsurfers
    lsurfer(j) = surfers{j}.npts;

    if (size(kern) == 1)
        ktmp = kern;
    else
        ktmp = kern(j);
    end
    opdims_mat(:,j) = ktmp.opdims;
    if ismethod(surfers{j},'oversample') 
        novers{j} = ktmp.get_overs_orders(surfers{j},targobj,eps);
    else
        novers{j} = NaN;
    end
end

% Assert that opdims(1) is constant across all source blocks (all kernels
% map to the same number of output components per target point), and that
% opdims(2) is consistent for each block-column j.
assert(all(opdims_mat(1,:) == opdims_mat(1,1)), ...
    'SURFEREVALMAT: opdims(1) is not constant across all source blocks');
for j = 1:nsurfers
    assert(opdims_mat(2,j) == opdims_mat(2,j), ...
        'SURFEREVALMAT: opdims(2) is not consistent for block-column %d', j);
end

ntarg   = size(targinfo.r(:,:), 2);
opdim1  = opdims_mat(1,1);

icollocs = zeros(nsurfers+1,1);
icollocs(1) = 1;
for j=1:nsurfers
   icollocs(j+1) = icollocs(j) + lsurfer(j)*opdims_mat(2,j);
end

nrows = ntarg * opdim1;
ncols = icollocs(end)-1;

if (~nonsmoothonly)
    sysmat = zeros(nrows,ncols);
else
    sysmat = sparse(nrows,ncols);
    isysmat = [];
    jsysmat = [];
    vsysmat = [];
end

%% Now build quadratures

for j = 1:nsurfers
    surferj = surfers{j};
    if (surferj.npts < 1 || ntarg < 1)
        sysmat_tmp = [];
        continue
    end

    opdims = opdims_mat(:,j);
    wts = repmat(surferj.wts(:).', opdims(2), 1);
    wts = wts(:).';

    if (size(kern) == 1)
        ktmp = kern;
    else
        ktmp = kern(j);
    end

    noversj = novers{j};
    if ismethod(surferj,'oversample') && ~any(isnan(noversj))
        [surferjover, xinterp] = surferj.oversample(noversj);
        xinterp = kron(xinterp, eye(opdims(2)));
    else
        novers{j} = NaN;
        surferjover = surferj;
        xinterp = eye(numel(wts));
    end

    wtsover = repmat(surferjover.wts(:).', opdims(2), 1);
    wtsover = wtsover(:).';

    if (~nonsmoothonly)
        sysmat_tmp = (ktmp.eval(surferjover,targinfo).*wtsover)*xinterp;
    end

    if adaptive_correction
        Qj = ktmp.getquad(surferj,eps,targinfo);
        if isfield(Qj, 'row_ptr')
            sysmat_quad = conv_rsc_to_spmat(surferj,Qj.row_ptr,Qj.col_ind,Qj.wnear, ktmp.rsc_to_interleave);
        else
            sysmat_quad = Qj;
        end
        if corrections
            if isfield(Qj, 'row_ptr')
                rp = Qj.row_ptr;
                ci = Qj.col_ind;
            else
                [rp, ci] = get_rsc_pattern(surferj, Qj, opdims);
            end
            sysmat_smth = smooth_sparse_quad(ktmp, targinfo, surferj, rp, ci, noversj);
            sysmat_quad = sysmat_quad-sysmat_smth;
        end

        if (~nonsmoothonly)
            inds = abs(sysmat_quad(:))>0;
            sysmat_tmp(inds) = sysmat_quad(inds);
        else
            sysmat_tmp = sysmat_quad;
        end
    end

    if (~nonsmoothonly)
        icolinds = icollocs(j):(icollocs(j+1)-1);
        sysmat(1:nrows,icolinds) = sysmat_tmp;
    else
        if adaptive_correction
            [isys,jsys,vsys] = find(sysmat_tmp);
            isysmat = [isysmat;isys];
            jsysmat = [jsysmat;jsys+icollocs(j)-1];
            vsysmat = [vsysmat;vsys];
        end
    end

end


if (nonsmoothonly)
    sysmat = sparse(isysmat,jsysmat,vsysmat,nrows,ncols);
end

end
