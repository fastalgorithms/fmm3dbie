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
    if numel(sz) == 2 && sz(1) == ri.nker
        % wnear is already (nker, nquad) — use as-is even if nquad == 0
        nquad = sz(2);
    else
        nquad = numel(wnear);
        wnear = reshape(wnear, 1, nquad);
    end

    m = ri.m;
    n = ri.n;
    ntarg = numel(row_ptr) - 1;
    if isempty(col_ind)
        spmat = sparse(m*ntarg, n*S.npts);
        return
    end

    ixyzs = S.ixyzs(:);
    npols = ixyzs(2:end) - ixyzs(1:end-1);
    npts_col = npols(col_ind(:));
    istart = [1; 1 + cumsum(npts_col)];

    icol = reshape(repelem(ixyzs(col_ind(:)) - istart(1:end-1), npts_col), [], 1) + (1:nquad).';
    irow = reshape(repelem((1:ntarg).', diff(istart(row_ptr(:)))), [], 1);

    e = ri.entries;
    irow = m*(irow-1) + [e.row_id];
    icol = n*(icol-1) + [e.col_id];
    kid = [e.ker_id];
    coef = [e.coef];
    if isequal(kid, 1:size(wnear, 1)) && all(coef == 1)
        vals = wnear.';
    else
        vals = wnear(kid,:).' .* coef;
    end
    spmat = sparse(irow(:), icol(:), vals(:), m*ntarg, n*S.npts);
end
