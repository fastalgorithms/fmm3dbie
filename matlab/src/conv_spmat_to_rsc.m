function rsc = conv_spmat_to_rsc(S, spmat, ri)
%CONV_SPMAT_TO_RSC convert a sparse matrix to patch-based rsc format.
%
%  Syntax:
%    rsc = conv_spmat_to_rsc(S, spmat)
%    rsc = conv_spmat_to_rsc(S, spmat, ri)
%
%  For scalar kernels (ri omitted): spmat is (ntarg, S.npts), wnear is
%  returned as (nquad,1).  For vector/tensor kernels pass the kernel's
%  rsc_to_interleave struct as ri: spmat is (m*ntarg, n*S.npts) and
%  wnear is returned as (nker, nquad).
%
%  Output rsc fields:
%    .row_ptr  (ntarg+1, 1)  CSR row pointer
%    .col_ind  (nnz, 1)      patch indices (1-based)
%    .iquad    (nnz+1, 1)    start of each patch block in wnear
%    .wnear    (nker, nquad) quadrature weights
%    .nquad    scalar
%    .nnz      scalar

    if nargin < 3 || isempty(ri)
        ri = kernel3d.rsc_interleave_scalar();
    end

    nker   = ri.nker;
    scalar = strcmp(ri.type, 'scalar');
    m      = max(ri.row_ids);
    n      = max(ri.col_ids);

    ixyzs    = S.ixyzs(:);
    npatches = S.npatches;
    npols    = ixyzs(2:end) - ixyzs(1:end-1);
    ntarg    = size(spmat, 1) / m;

    % Find active (target, patch) pairs from the (1,1) block component.
    [row_ptr, col_ind] = get_rsc_pattern(S, spmat, [m, n]);

    if isempty(col_ind)
        rsc.row_ptr = ones(ntarg+1, 1);
        rsc.col_ind = zeros(0, 1);
        rsc.iquad   = ones(1, 1);
        rsc.wnear   = complex(zeros(nker, 0));
        rsc.nquad   = 0;
        rsc.nnz     = 0;
        return
    end

    targ_inds  = repelem((1:ntarg).', diff(row_ptr));
    patch_inds = col_ind;
    nnz_pairs  = numel(col_ind);

    block_sizes = npols(patch_inds);
    iquad = [1; 1 + cumsum(block_sizes)];
    nquad = iquad(end) - 1;

    % Extract wnear: for each (targ, patch) pair and each kernel component
    % read the appropriate block row/col out of spmat.
    kr    = ri.row_ids(:);
    kc    = ri.col_ids(:);
    wnear = complex(zeros(nker, nquad));
    for k = 1:nnz_pairs
        t        = targ_inds(k);
        p        = patch_inds(k);
        src_cols = ixyzs(p):(ixyzs(p+1)-1);
        wstart   = iquad(k);
        wend     = iquad(k+1) - 1;
        for ki = 1:nker
            wnear(ki, wstart:wend) = spmat(m*(t-1)+kr(ki), n*(src_cols-1)+kc(ki));
        end
    end

    if scalar
        wnear = wnear(1,:).';
    end

    rsc.row_ptr = row_ptr;
    rsc.col_ind = col_ind;
    rsc.iquad   = iquad;
    rsc.wnear   = wnear;
    rsc.nquad   = nquad;
    rsc.nnz     = nnz_pairs;
end
