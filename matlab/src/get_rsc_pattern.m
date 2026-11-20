function [row_ptr, col_ind] = get_rsc_pattern(S, spmat, opdims)
%GET_RSC_PATTERN  Extract the (row_ptr, col_ind) sparsity pattern from a
%  sparse matrix in the RSC (row-sparse-compressed) patch format.
%
%  Syntax:
%    [row_ptr, col_ind] = get_rsc_pattern(S, spmat)
%    [row_ptr, col_ind] = get_rsc_pattern(S, spmat, opdims)
%
%  Finds the set of (target, patch) pairs for which spmat has at least one
%  nonzero entry, and returns them in CSR format.
%
%  Input arguments:
%    * S       : surfer object describing the source surface
%    * spmat   : sparse matrix of size (m*ntarg, n*S.npts)
%    * opdims  : [m, n] operator dimensions (default [1, 1])
%
%  Output arguments:
%    * row_ptr : (ntarg+1, 1) CSR row pointer into col_ind
%    * col_ind : (nnz, 1)     patch indices (1-based), sorted by target

    if nargin < 3 || isempty(opdims)
        opdims = [1, 1];
    end

    m = opdims(1);
    n = opdims(2);
    ntarg = size(spmat, 1) / m;
    npatches = S.npatches;
    patch_of_pt = repelem((1:npatches).', diff(S.ixyzs(:)));

    [itarg, ipatch] = find(spmat);
    itarg = floor((itarg(:) - 1) / m) + 1;
    ipatch = patch_of_pt(floor((ipatch(:) - 1) / n) + 1);

    T = sparse(ipatch, itarg, true, npatches, ntarg);
    clear itarg ipatch
    [col_ind, targ] = find(T);
    col_ind = col_ind(:);
    row_ptr = [1; 1 + cumsum(accumarray(targ(:), 1, [ntarg, 1]))];
end
