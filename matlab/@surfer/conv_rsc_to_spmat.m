function spmat = conv_rsc_to_spmat(S, row_ptr, col_ind, wnear)
% CONV_RSC_TO_SPMAT  Make sparse matrix from patch-wise RSC quadrature format.
%
%  Syntax:
%    spmat = conv_rsc_to_spmat(S, row_ptr, col_ind, wnear)
%
%  For scalar kernels:
%    wnear is (nquad,) or (1,nquad).  Returns sparse(ntarg, S.npts).
%
%  For vector kernels (e.g. Stokes):
%    wnear is (nker, nquad).
%    Returns sparse(m*ntarg, n*S.npts) where the kernel component layout
%    is inferred from nker:
%      nker == m*n  (full tensor, no symmetry)
%      nker == m*(m+1)/2  with m==n  (symmetric tensor, upper-triangle order)
%        Components stored as (1,1),(1,2),...,(1,m),(2,2),(2,3),...,(m,m)
%
%    Currently supported nker values:
%      nker==1  -> scalar (m=n=1)
%      nker==3  -> vector (m=3, n=1)
%      nker==9  -> full 3x3 tensor (m=n=3)
%      nker==6  -> symmetric 3x3 tensor (m=n=3, upper triangle)

    sz = size(wnear);
    if numel(sz) == 2 && sz(1) > 1 && sz(2) >1
        nker  = sz(1);
        nquad = sz(2);
    else
        % scalar: wnear is (nquad,) or (1,nquad) or already flattened
        nker  = 1;
        nquad = numel(wnear);
        wnear = reshape(wnear, 1, nquad);
    end

    % Determine operator dimensions [m, n] from nker
    if nker == 1
        m = 1; n = 1; sym = false;
    elseif nker == 3
        m = 3; n = 1; sym = false;
    elseif nker == 9
        m = 3; n = 3; sym = false;
    elseif nker == 6
        m = 3; n = 3; sym = true;   % symmetric 3x3, upper triangle
    else
        error('conv_rsc_to_spmat: unsupported nker=%d', nker);
    end

    ixyzs    = S.ixyzs(:);
    npatches = S.npatches;
    npols    = ixyzs(2:end) - ixyzs(1:end-1);   % points per patch
    npts     = S.npts;

    npts_col = npols(col_ind);

    ntarg   = numel(row_ptr) - 1;
    istarts = row_ptr(1:end-1);
    iends   = row_ptr(2:end) - 1;

    % Build point-level (irow_ind, icol_ind) for scalar indexing
    isrcinds = cell(npatches, 1);
    for i = 1:npatches
        isrcinds{i} = ixyzs(i):(ixyzs(i+1)-1);
    end

    total_pts = sum(npts_col);
    icol_ind  = zeros(total_pts, 1);
    irow_ind_scalar = zeros(total_pts, 1);  % target index (1..ntarg)

    istart = 1;
    for i = 1:ntarg
        if istarts(i) > iends(i), continue; end
        iinds  = horzcat(isrcinds{col_ind(istarts(i):iends(i))});
        nelem  = numel(iinds);
        icol_ind(istart:istart+nelem-1)       = iinds;
        irow_ind_scalar(istart:istart+nelem-1) = i;
        istart = istart + nelem;
    end

    if m == 1 && n == 1
        % ----- Scalar kernel -----
        spmat = sparse(irow_ind_scalar, icol_ind, wnear(1,:).', ntarg, npts);
        return
    end

    % ----- Vector / tensor kernel -----
    % For each scalar (targ, src) pair we need to fill an [m x n] block.
    % Build triplet lists for sparse().

    if sym
        % Symmetric 3x3: upper-triangle order
        % kernel index -> (row_comp, col_comp) pairs
        % wnear(1,:)=(1,1), (2,:)=(1,2), (3,:)=(1,3),
        %           (4,:)=(2,2), (5,:)=(2,3), (6,:)=(3,3)
        kr = [1 1 1 2 2 3];
        kc = [1 2 3 2 3 3];
        % Symmetry counterparts (off-diagonal reflected entries)
        sym_kr = [2 3 3];
        sym_kc = [1 1 2];
        sym_ki = [2 3 5];   % indices in wnear that give the reflected value
        n_entries_per_pair = nker + numel(sym_ki);  % 6 + 3 = 9
    else
        % Full tensor: col-major order 1..m*n
        [kc_mat, kr_mat] = meshgrid(1:n, 1:m);
        kr = kr_mat(:).';
        kc = kc_mat(:).';
        n_entries_per_pair = nker;
    end

    np = numel(icol_ind);   % total scalar (targ, src) pairs
    NZ = np * n_entries_per_pair;
    II = zeros(NZ, 1);
    JJ = zeros(NZ, 1);
    VV = zeros(NZ, 1);

    % wnear is (nker, nquad); the k-th quadrature weight (in scalar order)
    % for pair index p corresponds to wnear(:, p) (one column per scalar pair).
    % NOTE: wnear is stored with nquad matching the scalar expansion, so
    % column p of wnear corresponds to icol_ind(p), irow_ind_scalar(p).

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