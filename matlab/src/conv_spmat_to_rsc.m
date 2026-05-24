function rsc = conv_spmat_to_rsc(S, spmat)
%CONV_SPMAT_TO_RSC convert a sparse matrix to patch-based rsc format.
%
%  Syntax:
%    rsc = conv_spmat_to_rsc(S, spmat)
%
%  Given a sparse matrix whose columns correspond to discretisation
%  points on S and whose rows correspond to targets, this routine
%  returns the row-sparse-compressed (rsc) struct used by the
%  layer-potential evaluators.
%
%  The rsc struct stores interactions at the *patch* level: one entry
%  in col_ind per (target, patch) pair that has at least one nonzero
%  column in spmat.  The corresponding quadrature weights (wnear) are
%  stored contiguously in the order prescribed by iquad.
%
%  Input arguments:
%    * S      : surfer object
%    * spmat  : sparse matrix, size (ntarg, S.npts)
%               Columns are ordered consistently with S.ixyzs.
%
%  Output arguments:
%    * rsc.row_ptr : (ntarg+1, 1) double -- standard CSR row pointer
%    * rsc.col_ind : (nnz,    1) double -- patch indices (1-based)
%    * rsc.iquad   : (nnz+1,  1) double -- start of each patch block
%                                          in wnear
%    * rsc.wnear   : (nquad,  1) complex -- quadrature weights
%    * rsc.nquad   : scalar
%    * rsc.nnz     : scalar  (number of target-patch pairs)
%

    ixyzs    = S.ixyzs(:);          % (npatches+1, 1)
    npatches = S.npatches;
    ntarg    = size(spmat, 1);

    npols = ixyzs(2:end) - ixyzs(1:end-1);   % points per patch

    % ----------------------------------------------------------------
    % Build row_ptr and col_ind at the patch level.
    % A (target i, patch p) pair is "active" if spmat(i, ixyzs(p):ixyzs(p+1)-1)
    % contains at least one nonzero.
    % ----------------------------------------------------------------

    % Work column-by-column in the sparse matrix to find active patches
    % per target without expanding to a full dense matrix.

    [row_s, col_s] = find(spmat);   % row = target index, col = point index

    if isempty(row_s)
        rsc.row_ptr = ones(ntarg+1, 1);
        rsc.col_ind = zeros(0,      1);
        rsc.iquad   = ones( 1,      1);
        rsc.wnear   = complex(zeros(0, 1));
        rsc.nquad   = 0;
        rsc.nnz     = 0;
        return
    end

    % Map each nonzero point index to its patch index
    % patch_of_col(k) = patch that contains point col_s(k)
    patch_of_col = zeros(size(col_s));
    for p = 1:npatches
        mask = (col_s >= ixyzs(p)) & (col_s < ixyzs(p+1));
        patch_of_col(mask) = p;
    end

    % Unique (target, patch) pairs
    pairs      = unique([row_s, patch_of_col], 'rows');  % (nnz, 2)
    targ_inds  = pairs(:,1);
    patch_inds = pairs(:,2);
    nnz_pairs  = size(pairs, 1);

    % Build row_ptr
    counts = zeros(ntarg, 1);
    for k = 1:nnz_pairs
        counts(targ_inds(k)) = counts(targ_inds(k)) + 1;
    end
    row_ptr = [1; 1 + cumsum(counts)];

    col_ind = patch_inds;

    % ----------------------------------------------------------------
    % Build iquad: for rsc entry k (the k-th (targ,patch) pair),
    % iquad(k) is where in wnear the weights for that pair start.
    % The block has npols(patch_inds(k)) entries.
    % ----------------------------------------------------------------
    block_sizes = npols(patch_inds);              % (nnz_pairs, 1)
    iquad = [1; 1 + cumsum(block_sizes)];         % (nnz_pairs+1, 1)
    nquad = iquad(end) - 1;

    % ----------------------------------------------------------------
    % Fill wnear: for each (targ, patch) pair extract the row of spmat
    % corresponding to that target and the columns for that patch.
    % ----------------------------------------------------------------
    wnear = complex(zeros(nquad, 1));
    for k = 1:nnz_pairs
        t = targ_inds(k);
        p = patch_inds(k);
        src_cols = ixyzs(p):(ixyzs(p+1)-1);           % point indices for patch p
        wstart   = iquad(k);
        wend     = iquad(k+1) - 1;
        wnear(wstart:wend) = spmat(t, src_cols);
    end

    rsc.row_ptr = row_ptr;
    rsc.col_ind = col_ind;
    rsc.iquad   = iquad;
    rsc.wnear   = wnear;
    rsc.nquad   = nquad;
    rsc.nnz     = nnz_pairs;
end
