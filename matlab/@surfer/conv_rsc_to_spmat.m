function spmat = conv_rsc_to_spmat(S, row_ptr, col_ind, wnear, ri)
% CONV_RSC_TO_SPMAT  Make sparse matrix from patch-wise RSC quadrature format.
%
%  Syntax:
%    spmat = conv_rsc_to_spmat(S, row_ptr, col_ind, wnear)
%    spmat = conv_rsc_to_spmat(S, row_ptr, col_ind, wnear, ri)
%
%  For scalar kernels (ri omitted): wnear is (nquad,) or (1,nquad).
%  For vector/tensor kernels pass the kernel's rsc_to_interleave struct as ri:
%    wnear is (nker, nquad).
%  Returns sparse(m*ntarg, n*S.npts) where m,n come from ri.
%
%  ri fields used:
%    .m, .n    - block opdims
%    .nker     - number of rows in wnear
%    .entries  - struct array, each with:
%                  .ker_id  (1-based index into wnear rows)
%                  .row_id  (1-based block row)
%                  .col_id  (1-based block col)
%                  .coef    (scalar coefficient, may be complex)

    % Default ri: scalar
    if nargin < 5 || isempty(ri)
        ri = kernel3d.rsc_interleave_full(1, 1);
    end

    sz = size(wnear);
    if numel(sz) == 2 && sz(1) > 1 && sz(2) > 1
        nquad = sz(2);
    else
        nquad = numel(wnear);
        wnear = reshape(wnear, 1, nquad);
    end

    m    = ri.m;
    n    = ri.n;
    npts = S.npts;

    ixyzs    = S.ixyzs(:);
    npatches = S.npatches;
    npols    = ixyzs(2:end) - ixyzs(1:end-1);   % points per patch

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

    ne = numel(ri.entries);
    np = numel(icol_ind);
    NZ = np * ne;
    II = zeros(NZ, 1);
    JJ = zeros(NZ, 1);
    VV = complex(zeros(NZ, 1));
    off = 0;
    for ei = 1:ne
        e = ri.entries(ei);
        II(off+1:off+np) = m*(irow_ind_scalar-1) + e.row_id;
        JJ(off+1:off+np) = n*(icol_ind-1)        + e.col_id;
        VV(off+1:off+np) = e.coef * wnear(e.ker_id, :).';
        off = off + np;
    end
    spmat = sparse(II, JJ, VV, m*ntarg, n*npts);
end
