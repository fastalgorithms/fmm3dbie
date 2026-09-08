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

    ixyzs    = S.ixyzs(:);
    npatches = S.npatches;
    ntarg    = size(spmat, 1) / m;

    % Union the nonzero pattern across every (block-row, block-col) sub-block.
    row_s = [];
    col_s = [];
    for bi = 1:m
        for bj = 1:n
            spmat_bij = spmat(bi:m:end, bj:n:end);
            [row_b, col_b] = find(spmat_bij);
            row_s = [row_s; row_b]; %#ok<AGROW>
            col_s = [col_s; col_b]; %#ok<AGROW>
        end
    end

    if isempty(row_s)
        row_ptr = ones(ntarg+1, 1);
        col_ind = zeros(0, 1);
        return
    end

    % Map point column indices to patch indices via binary search on the
    % sorted patch-start offsets (avoids materializing an O(nnz x npatches)
    % dense array, which a naive bsxfun(@ge, col_s, ixyzs(1:end-1).')
    % outer-product approach would do -- this blows up memory badly for
    % large nnz and/or large npatches).
    patch_of_col = discretize(col_s, [ixyzs; inf]);

    % Deduplicate (target,patch) pairs by accumulating into a sparse
    % ntarg x npatches indicator matrix.
    T = sparse(row_s, patch_of_col, 1, ntarg, npatches);
    [patch_inds, targ_inds] = find(T.');   % find on transpose -> sorted by target (col-major over T.')
    patch_inds = patch_inds(:);
    targ_inds  = targ_inds(:);

    counts  = accumarray(targ_inds, 1, [ntarg, 1]);
    row_ptr = [1; 1 + cumsum(counts)];
    col_ind = patch_inds;
end
