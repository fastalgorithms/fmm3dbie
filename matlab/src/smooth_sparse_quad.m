function Asmth = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,nover,S_over,val2over)
%
%  SMOOTH_SPARSE_QUAD
%
%  Syntax:
%    Asmth = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,nover)
%    Asmth = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,nover,S_over,val2over)
%
%  This subroutine returns a sparse matrix of smooth quadrature
%  corrections for the interactions between a set of targets and
%  the source patches described by the compressed sparse row
%  structure (row_ptr, col_ind). Each source patch is oversampled
%  to order nover and the kernel is evaluated at every oversampled
%  source point.
%
%  The kernel is evaluated one patch at a time, in chunks of targets,
%  and each block is mapped back to the nodes of S with that patch's
%  weighted block of the interpolation matrix before it is stored.
%  Memory is therefore proportional to the number of entries of the
%  result, not of the oversampled matrix.
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
%    * nover: integer, or (npatches,1) integers
%        quadrature order to use when oversampling each source patch.
%        If NaN, or if S cannot be oversampled, the nodes of S are used.
%    * S_over, val2over: (optional) S oversampled to nover and the
%        interpolation matrix from S.oversample (not kron'd with the
%        operator dimension), to avoid recomputing them
%
%  Output arguments:
%    * Asmth: sparse matrix of size (opdim(1)*ntarg, opdim(2)*S.npts)
%        smooth quadrature matrix acting on oversampled source
%        densities, projected back to the original S discretization
%
%  See also OVERSAMPLE, GETNEAR
%

    kern = kernel3d(kern);
    m = kern.opdims(1);
    n = kern.opdims(2);

    if isnumeric(targs)
        targs = struct('r', targs);
    end
    ntarg = size(targs.r(:,:), 2);
    npatches = S.npatches;
    ixyzs = S.ixyzs(:);
    npols = ixyzs(2:end) - ixyzs(1:end-1);

    if nargin >= 8
        % precomputed oversampling
    elseif ismethod(S,'oversample') && ~any(isnan(nover))
        [S_over, val2over] = oversample(S, nover);
    else
        S_over = S;
        val2over = [];
    end
    ixyzso = S_over.ixyzs(:);
    wover  = S_over.wts(:);

    % Use the full FMM box size for every patch's self-interaction mask.
    src_r = S_over.r(:,:);
    targ_r = targs.r(:,:);
    targ_r(end+1:size(src_r,1),:) = 0;
    lo = min(min(src_r,[],2), min(targ_r,[],2));
    hi = max(max(src_r,[],2), max(targ_r,[],2));
    L = max(hi-lo);
    clear src_r targ_r

    % targets of each source patch
    rsc = [];
    rsc.row_ptr = row_ptr;
    rsc.col_ind = col_ind;
    rsc.iquad   = [1; 1 + cumsum(npols(col_ind(:)))];
    csc = conv_rsc_to_csc(npatches, rsc);

    % triplets, preallocated to the size of the result
    ntot = m*n*(rsc.iquad(end) - 1);
    irow = zeros(ntot, 1);
    icol = zeros(ntot, 1);
    vals = zeros(ntot, 1);
    ipos = 0;

    % bound the size of a single kernel block (in entries)
    maxblk = 2^22;
    for p = 1:npatches
        cstart = csc.col_ptr(p);
        cend   = csc.col_ptr(p+1) - 1;
        if cstart > cend
            continue
        end

        src_inds = ixyzso(p):(ixyzso(p+1)-1);
        npo = numel(src_inds);
        srcp = [];
        srcp.r = S_over.r(:,src_inds);
        for field = kern.src_fields(:).'
            srcp.(field{1}) = slice_field(S_over.(field{1}), src_inds, S_over.npts);
        end

        % weighted interpolation: oversampled values -> nodes of patch p,
        % expanded to the n components of the density
        icols = ixyzs(p):ixyzs(p+1)-1;
        if isempty(val2over)
            Wp = diag(wover(src_inds));
        else
            Wp = wover(src_inds) .* full(val2over(src_inds, icols));
        end
        Wp = kron(Wp, eye(n));                                   % (n*npo, n*np)
        icols_p = (1:n).' + n*(icols - 1);
        icols_p = icols_p(:);

        nchunk = max(1, floor(maxblk / (m*n*npo)));
        for c0 = cstart:nchunk:cend
            targ_inds = csc.row_ind(c0:min(c0+nchunk-1, cend));
            targ_inds = targ_inds(:).';
            targp = [];
            targp.r = targs.r(:,targ_inds);
            for field = kern.targ_fields(:).'
                targp.(field{1}) = slice_field(targs.(field{1}), targ_inds, ntarg);
            end

            B = kern.eval_mask(srcp, targp, L) * Wp;              % (m*nt, n*np)
            irows_p = (1:m).' + m*(targ_inds - 1);             % (m, nt)

            nb = numel(B);
            idx = ipos + (1:nb);
            [ii, jj] = ndgrid(irows_p(:), icols_p);
            irow(idx) = ii(:);
            icol(idx) = jj(:);
            if ~isreal(B) && isreal(vals)
                vals = complex(vals);
            end
            vals(idx) = B(:);
            ipos = ipos + nb;
        end
    end

    Asmth = sparse(irow(1:ipos), icol(1:ipos), vals(1:ipos), ...
        m*ntarg, n*S.npts);

end


function v = slice_field(A, idx, npts)
if ~ismatrix(A), A = A(:,:); end
if size(A, 2) == npts
    v = A(:, idx);
elseif size(A, 1) == npts
    v = A(idx, :).';
else
    error('SLICE_FIELD: no dimension of A matches npts');
end
end
