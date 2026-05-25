function spmat = conv_rsc_to_spmat(S, row_ptr, col_ind, wnear, ri)
% CONV_RSC_TO_SPMAT  Make sparse matrix from patch-wise RSC quadrature format.
%
%  Syntax:
%    spmat = conv_rsc_to_spmat(S, row_ptr, col_ind, wnear)
%    spmat = conv_rsc_to_spmat(S, row_ptr, col_ind, wnear, ri)
%
%  For scalar kernels (ri omitted or ri.type == 'scalar'):
%    wnear is (nquad,) or (1,nquad).  Returns sparse(ntarg, S.npts).
%
%  For vector/tensor kernels pass the kernel's rsc_to_interleave struct as ri:
%    wnear is (nker, nquad).
%    Returns sparse(m*ntarg, n*S.npts) where the block layout is read from ri.
%
%  ri fields used:
%    .type     - 'scalar', 'full', or 'symmetric'
%    .nker     - number of rows in wnear
%    .row_ids  - (nker,1) row index within the opdims block for each wnear row
%    .col_ids  - (nker,1) col index within the opdims block for each wnear row
%    For 'symmetric' only:
%    .sym_row_ids - row indices of off-diagonal reflected entries
%    .sym_col_ids - col indices of off-diagonal reflected entries
%    .sym_ker_ids - which wnear row each reflected entry mirrors

    % Default ri: scalar
    if nargin < 5 || isempty(ri)
        ri = kernel3d.rsc_interleave_scalar();
    end

    sz = size(wnear);
    if numel(sz) == 2 && sz(1) > 1 && sz(2) > 1
        nquad = sz(2);
    else
        nquad = numel(wnear);
        wnear = reshape(wnear, 1, nquad);
    end

    nker = ri.nker;

    ixyzs    = S.ixyzs(:);
    npatches = S.npatches;
    npols    = ixyzs(2:end) - ixyzs(1:end-1);   % points per patch
    npts     = S.npts;

    npts_col = npols(col_ind);

    ntarg   = numel(row_ptr) - 1;
    istarts = row_ptr(1:end-1);
    iends   = row_ptr(2:end) - 1;

    % Build point-level (irow_ind_scalar, icol_ind) indices
    isrcinds = cell(npatches, 1);
    for i = 1:npatches
        isrcinds{i} = ixyzs(i):(ixyzs(i+1)-1);
    end

    total_pts = sum(npts_col);
    icol_ind        = zeros(total_pts, 1);
    irow_ind_scalar = zeros(total_pts, 1);

    istart = 1;
    for i = 1:ntarg
        if istarts(i) > iends(i), continue; end
        iinds  = horzcat(isrcinds{col_ind(istarts(i):iends(i))});
        nelem  = numel(iinds);
        icol_ind(istart:istart+nelem-1)       = iinds;
        irow_ind_scalar(istart:istart+nelem-1) = i;
        istart = istart + nelem;
    end

    if strcmp(ri.type, 'scalar')
        % ----- Scalar kernel -----
        spmat = sparse(irow_ind_scalar, icol_ind, wnear(1,:).', ntarg, npts);
        return
    end

    % ----- Vector / tensor kernel -----
    % Infer opdims from the block indices stored in ri.
    m = max(ri.row_ids);
    n = max(ri.col_ids);

    kr = ri.row_ids(:).';   % (1, nker)
    kc = ri.col_ids(:).';   % (1, nker)

    sym = strcmp(ri.type, 'symmetric');
    if sym
        n_entries_per_pair = nker + numel(ri.sym_ker_ids);
    else
        n_entries_per_pair = nker;
    end

    np = numel(icol_ind);   % total scalar (targ, src) pairs
    NZ = np * n_entries_per_pair;
    II = zeros(NZ, 1);
    JJ = zeros(NZ, 1);
    VV = zeros(NZ, 1);

    off = 0;
    for ki = 1:nker
        ir = kr(ki);
        ic = kc(ki);
        II(off+1:off+np) = m*(irow_ind_scalar-1) + ir;
        JJ(off+1:off+np) = n*(icol_ind-1) + ic;
        VV(off+1:off+np) = wnear(ki, :).';
        off = off + np;
    end

    if sym
        sym_kr = ri.sym_row_ids(:).';
        sym_kc = ri.sym_col_ids(:).';
        sym_ki = ri.sym_ker_ids(:).';
        for si = 1:numel(sym_ki)
            ki_src = sym_ki(si);
            ir = sym_kr(si);
            ic = sym_kc(si);
            II(off+1:off+np) = m*(irow_ind_scalar-1) + ir;
            JJ(off+1:off+np) = n*(icol_ind-1) + ic;
            VV(off+1:off+np) = wnear(ki_src, :).';
            off = off + np;
        end
    end

    spmat = sparse(II, JJ, VV, m*ntarg, n*npts);
end
