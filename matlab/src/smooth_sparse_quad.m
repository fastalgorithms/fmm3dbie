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

    src_field = fieldnames(S)';
    src_field_use = {};
    for i = 1:length(src_field)
        if mod(size(S.(src_field{i})(:,:),2), size(S.r(:,:),2)) == 0 && ~isempty(S.(src_field{i}))
        src_field_use{length(src_field_use)+1} = src_field{i};
        end
    end

    targ_field = fieldnames(targs)';
    targ_field_use = {};
    for i = 1:length(targ_field)
        if mod(size(targs.(targ_field{i})(:,:),2), size(targs.r(:,:),2)) == 0 && ~isempty(targs.(targ_field{i}))
        targ_field_use{length(targ_field_use)+1} = targ_field{i};
        end
    end

    ntarg = size(targs.r, 2);
    
    [S_over, val2over] = oversample(S, nover);

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
        ntarg_p   = length(targ_inds);

        % Oversampled source points on patch p
        src_inds  = ixyzs(p):(ixyzs(p+1)-1);      % (npts_p,1)
        npts_p    = length(src_inds);

        % Extract source and target sub-structs
        srcp = [];
        for field = src_field_use
            srcp.(field{1}) = S_over.(field{1})(:,src_inds);
        end


        targp = [];
        for field = targ_field_use
            targp.(field{1}) = targs.(field{1})(:,targ_inds);
        end

        % Evaluate kernel: result is (opdims(1)*ntarg_p) x (opdims(2)*npts_p)
        Kblock = kern.eval(srcp, targp);

        % Weight columns by quadrature weights of oversampled patch
        wts_p = repmat(S_over.wts(src_inds).', opdims(2), 1);
        wts_p = wts_p(:).';
        Kblock = Kblock .* wts_p;

        % Zero out self-interactions (source == target, relative tol 1e-14)
        % Only relevant when ntarg_p == npts_p (i.e. on-surface self patch).
        % We check each (target, src) pair via norms.
        if ntarg_p == npts_p
            src_norm = max(vecnorm(srcp.r), 1);
            for q = 1:npts_p
                diff_norms = vecnorm(targp.r - srcp.r(:,q));
                self_mask  = diff_norms < 1e-14 * src_norm(q);
                row_off    = (1:opdims(1)) + opdims(1)*(find(self_mask)-1).';
                col_off    = (1:opdims(2)) + opdims(2)*(q-1);
                Kblock(row_off(:), col_off) = 0;
            end
        end

        % Build global row/col indices (before val2over projection)
        % rows: opdims(1) entries per target
        irows_p = (1:opdims(1)).' + opdims(1)*(targ_inds(:).'-1);
        irows_p = irows_p(:);  % (opdims(1)*ntarg_p, 1)
        % cols: opdims(2) entries per oversampled source point
        icols_p = (1:opdims(2)).' + opdims(2)*(src_inds(:).'-1);
        icols_p = icols_p(:);  % (opdims(2)*npts_p, 1)

        [ii, jj, vv] = find(sparse(Kblock));
        irow_vals = [irow_vals; irows_p(ii)];
        icol_vals = [icol_vals; icols_p(jj)];
        vals      = [vals;      vv];
    end

    Asmth = sparse(irow_vals, icol_vals, vals, ...
        opdims(1)*ntarg, opdims(2)*S_over.npts);

    Asmth = Asmth * val2over;
end