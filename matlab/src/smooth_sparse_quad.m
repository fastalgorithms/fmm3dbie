function Asmth = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,nover)
%
%  SMOOTH_SPARSE_QUAD
%
%  Syntax:
%    Asmth = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,nover)
%
%  This subroutine returns a sparse matrix of smooth quadrature
%  corrections for the interactions between a set of targets and
%  the source patches described by the compressed sparse row
%  structure (row_ptr, col_ind). Each source patch is oversampled
%  to order nover and the kernel is evaluated at every oversampled
%  source point.
%
%  Input arguments:
%    * kern: kernel function handle, or a struct with field kern.eval
%        containing such a handle
%    * targs: target information in ptinfo struct format, or a plain
%        (3,ntarg) array of positions
%    * S: surfer object describing the source surface
%    * row_ptr: (ntarg+1,1) integer array
%        pointer into col_ind; row_ptr(i):row_ptr(i+1)-1 gives the
%        source patches interacting with target i
%    * col_ind: (nnz,1) integer array
%        patch indices of the non-zero target-patch interactions,
%        sorted by target number
%    * nover: integer
%        quadrature order to use when oversampling each source patch
%    * lbat: integer (optional, default 1e3)
%        batch size for kernel evaluations
%
%  Output arguments:
%    * Asmth: sparse matrix of size (opdim(1)*ntarg, opdim(2)*S.npts)
%        smooth quadrature matrix acting on oversampled source
%        densities, projected back to the original S discretization
%
%  See also OVERSAMPLE, GETNEAR
%

    kern = kernel3d(kern);

    ntarg = size(targs.r, 2);
    
    if ismethod(S,'oversample') && ~any(isnan(nover))
        [S_over, val2over] = oversample(S, nover);
    else
        S_over = S; val2over = eye(size(S.r(:,:),2));
    end
    ixyzs = S_over.ixyzs;
    opdims = kern.opdims;

    % Convert rsc -> csc so we can loop over patches (columns) rather
    % than targets (rows), evaluating all oversampled source points on
    % each patch against all interacting targets in one kern.eval call.
    rsc = [];
    rsc.row_ptr = row_ptr;
    rsc.col_ind = col_ind;
    % Build iquad: each rsc entry maps to its own contiguous block of
    % oversampled source points; iquad marks the start of each block.
    npatches = S.npatches;
    nnz_rsc  = length(col_ind);
    npols = ixyzs(2:end) - ixyzs(1:end-1);  % pts per patch (oversampled)
    rsc.iquad = [1; 1 + cumsum(npols(col_ind(:)))];
    csc = conv_rsc_to_csc(npatches, rsc);

    % Accumulate sparse entries patch by patch
    irow_vals = [];
    icol_vals = [];
    vals      = [];

    for p = 1:npatches
        ptstart = csc.col_ptr(p);
        ptend   = csc.col_ptr(p+1) - 1;
        if ptstart > ptend
            continue
        end

        % Targets interacting with patch p
        targ_inds = csc.row_ind(ptstart:ptend);   % (ntarg_p,1)

        % Oversampled source points on patch p
        src_inds  = ixyzs(p):(ixyzs(p+1)-1);      % (npts_p,1)

        % Extract source and target sub-structs
        srcp = [];
        srcp.r = S_over.r(:,src_inds);
        for field = kern.src_fields(:).'
            srcp.(field{1}) = S_over.(field{1})(:,src_inds);
        end


        targp = [];
        targp.r = targs.r(:,targ_inds);
        for field = kern.targ_fields(:).'
            targp.(field{1}) = targs.(field{1})(:,targ_inds);
        end

        % Evaluate kernel. Use eval_mask to zero out self interactions
        Kblock = kern.eval_mask(srcp, targp);

        % Weight columns by quadrature weights of oversampled patch
        wts_p = repmat(S_over.wts(src_inds).', opdims(2), 1);
        wts_p = wts_p(:).';
        Kblock = Kblock .* wts_p;

        % Build global row/col indices (before val2over projection)
        % rows: opdims(1) entries per target
        irows_p = (1:opdims(1)).' + opdims(1)*(targ_inds(:).'-1);
        irows_p = irows_p(:);  % (opdims(1)*ntarg_p, 1)
        % cols: opdims(2) entries per oversampled source point
        icols_p = (1:opdims(2)).' + opdims(2)*(src_inds(:).'-1);
        icols_p = icols_p(:);  % (opdims(2)*npts_p, 1)

        [ii, jj, vv] = find(sparse(Kblock));
        irows_p = irows_p(ii);
        irow_vals = [irow_vals; irows_p(:)];
        icol_vals = [icol_vals; icols_p(jj)];
        vals      = [vals;      vv(:)];
    end

    Asmth = sparse(irow_vals, icol_vals, vals, ...
        opdims(1)*ntarg, opdims(2)*S_over.npts);

    Asmth = Asmth * val2over;
end