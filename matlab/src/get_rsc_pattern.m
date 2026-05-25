function [row_ptr, col_ind] = get_rsc_pattern(S, spmat, opdims)
%GET_RSC_PATTERN  Extract the (row_ptr, col_ind) sparsity pattern from a
%  sparse matrix in the RSC (row-sparse-compressed) patch format.
%
%  Syntax:
%    [row_ptr, col_ind] = get_rsc_pattern(S, spmat)
%    [row_ptr, col_ind] = get_rsc_pattern(S, spmat, opdims)
%
%  Finds the set of (target, patch) pairs for which spmat has at least one
%  nonzero entry, and returns them in CSR format.  For vector/tensor kernels
%  with block size opdims = [m, n], only the (1,1) sub-block is examined
%  (all block components share the same sparsity pattern).
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

    ixyzs    = S.ixyzs(:);
    npatches = S.npatches;
    ntarg    = size(spmat, 1) / m;

    % Work on the (1,1)-block submatrix; sparsity pattern is the same for
    % all block components.
    spmat_11 = spmat(1:m:end, 1:n:end);
    [row_s, col_s] = find(spmat_11);

    if isempty(row_s)
        row_ptr = ones(ntarg+1, 1);
        col_ind = zeros(0, 1);
        return
    end

    % Map point column indices to patch indices.
    patch_of_col = sum(bsxfun(@ge, col_s, ixyzs(1:end-1).'), 2);

    pairs      = unique([row_s, patch_of_col], 'rows');
    targ_inds  = pairs(:,1);
    patch_inds = pairs(:,2);

    counts  = accumarray(targ_inds, 1, [ntarg, 1]);
    row_ptr = [1; 1 + cumsum(counts)];
    col_ind = patch_inds;
end
