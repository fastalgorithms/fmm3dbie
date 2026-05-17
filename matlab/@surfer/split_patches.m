function [objout,varargout] = split_patches(obj,isplit)
%
%  SPLIT_PATCHES
%
%  Syntax:
%    objout = split_patches(obj,isplit)
%    [objout,val2split] = split_patches(obj,isplit)
%
%  This subroutine refines a subset of patches in obj by splitting each
%  listed patch into four equal sub-patches of the same order, leaving the
%  remaining patches unchanged.
%
%  Input arguments:
%    * obj: surfer object describing the source surface
%    * isplit: (1,m) integer array of patch indices to split
%        (optional, default empty — returns obj unchanged)
%
%  Output arguments:
%    * objout: surfer object with the refined patches. Patches not in
%        isplit appear in their original order, each split patch
%        is replaced by its four sub-patches.
%    * val2split: (optional) sparse prolongation matrix of size
%        (objout.npts, obj.npts) mapping densities on obj to densities
%        on objout
%
%  See also PATCH_SUB, OVERSAMPLE
%

if nargin < 2
    isplit = [];
end

npatches = obj.npatches;

% Convert isplit to a logical mask for O(1) per-patch lookup
lsplit = false(npatches, 1);
lsplit(isplit) = true;

% Build unique (norder, iptype) pairs and precompute patch_sub matrices
ntmp = zeros(npatches, 2);
ntmp(:,1) = obj.norders;
ntmp(:,2) = obj.iptype;

[ntmp_uni, ~, intmp] = unique(ntmp, 'rows');
nuni = size(ntmp_uni, 1);

xinterp_cell = cell(nuni, 1);
dumat_cell   = cell(nuni, 1);
dvmat_cell   = cell(nuni, 1);
npols_uni    = zeros(nuni, 1);

for i = 1:nuni
    norder = ntmp_uni(i,1);
    iptype = ntmp_uni(i,2);

    [~, xinterp, ~, dumat, dvmat] = patch_sub(norder, iptype);
    xinterp_cell{i} = xinterp;
    dumat_cell{i}   = dumat * xinterp;
    dvmat_cell{i}   = dvmat * xinterp;
    npols_uni(i)    = size(xinterp, 2);
end

% Compute total number of output points for preallocation
npts_per_patch = npols_uni(intmp);
npts_per_patch(lsplit) = 4 * npts_per_patch(lsplit);
nptso = sum(npts_per_patch);
ixyzso = cumsum([1; npts_per_patch]);

% Assemble srcmat
srcmat = zeros(12, nptso);

for ii = 1:npatches
    iind  = intmp(ii);
    iout  = ixyzso(ii):ixyzso(ii+1)-1;

    if lsplit(ii)
        rvals = obj.r(:, obj.ixyzs(ii):obj.ixyzs(ii+1)-1);
        rsub  = rvals * xinterp_cell{iind}.';
        dusub = rvals * dumat_cell{iind}.';
        dvsub = rvals * dvmat_cell{iind}.';
        n     = cross(dusub, dvsub);
        n     = n ./ vecnorm(n);
        srcmat(:,iout) = [rsub; dusub; dvsub; n];
    else
        srcmat(:,iout) = obj.srcvals{ii};
    end
end

npatch_out  = npatches + 3*length(isplit);
nsubpatches = npts_per_patch ./ npols_uni(intmp);
norders_out = repelem(obj.norders, nsubpatches);
iptype_out  = repelem(obj.iptype,  nsubpatches);

objout = surfer(npatch_out, norders_out, srcmat, iptype_out);

if nargout > 1
    isysmat_all = cell(nuni, 1);
    jsysmat_all = cell(nuni, 1);
    vsysmat_all = cell(nuni, 1);
    for i = 1:nuni
        [isysmat_all{i},jsysmat_all{i},vsysmat_all{i}] = find(xinterp_cell{i});
    end

    nnz_split   = sum(arrayfun(@(i) length(isysmat_all{intmp(i)}), isplit(:)));
    nnz_unsplit = sum(npts_per_patch(~lsplit));
    isysmat = zeros(nnz_split + nnz_unsplit, 1);
    jsysmat = zeros(nnz_split + nnz_unsplit, 1);
    vsysmat = zeros(nnz_split + nnz_unsplit, 1);

    istart = 1;
    for i = isplit(:).'
        nnz_i = length(isysmat_all{intmp(i)});
        idx = istart:istart+nnz_i-1;
        isysmat(idx) = isysmat_all{intmp(i)} + ixyzso(i) - 1;
        jsysmat(idx) = jsysmat_all{intmp(i)} + obj.ixyzs(i) - 1;
        vsysmat(idx) = vsysmat_all{intmp(i)};
        istart = istart + nnz_i;
    end
    for i = setdiff(1:npatches, isplit)
        np = npts_per_patch(i);
        idx = istart:istart+np-1;
        isysmat(idx) = (ixyzso(i):ixyzso(i+1)-1).';
        jsysmat(idx) = (obj.ixyzs(i):obj.ixyzs(i+1)-1).';
        vsysmat(idx) = 1;
        istart = istart + np;
    end
    varargout{1} = sparse(isysmat, jsysmat, vsysmat, nptso, obj.npts);
end

end
