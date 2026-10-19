function [wsmth, iquad] = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,nover)
%
%  SMOOTH_SPARSE_QUAD
%
%  Syntax:
%    [wsmth, iquad] = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,nover)
%
%  Oversampled smooth quadrature restricted to the near field
%  (row_ptr, col_ind), interpolated back to the nodes of S patch by patch.
%  Kernels with getquad_sub = true do this in Fortran instead, via
%  getquad opts.corrections (see GETNEARQUAD_MEX).
%
%  Input arguments:
%    * kern: kernel3d object, or a function handle accepted by kernel3d
%    * targs: target information, a struct or surfer with field r and
%        the fields in kern.targ_fields, or a plain (3,ntarg) array
%    * S: surfer object describing the source surface
%    * row_ptr: (ntarg+1,1) integer array
%        row_ptr(i):row_ptr(i+1)-1 index the entries of col_ind
%        for target i
%    * col_ind: (nnz,1) integer array
%        source patches in the near field of each target, sorted by
%        target
%    * nover: oversampling order, scalar or (npatches,1).
%        If NaN, or if S cannot be oversampled, the nodes of S are used.
%
%  Output arguments:
%    * wsmth: (m*n, nquad) array, [m,n] = kern.opdims
%        wnear layout on S, blocks ordered as rsc_interleave_full(m,n);
%        use conv_rsc_to_spmat to get the sparse matrix
%    * iquad: (nnz+1,1) integer array
%
%  See also GETNEAR, GET_RSC_PATTERN, CONV_RSC_TO_SPMAT, GETNEARQUAD_MEX
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

    if ismethod(S,'oversample') && ~any(isnan(nover))
        [S_over, val2over] = oversample(S, nover);
    else
        S_over = S;
        val2over = [];
    end
    ixyzso = S_over.ixyzs(:);
    wover  = S_over.wts(:);

    % one interpolation matrix per (order, oversampled order, patch type)
    keys = [S.norders(:), S_over.norders(:), S.iptype(:)];
    [~, irep, ikey] = unique(keys, 'rows');
    pmats = cell(numel(irep), 1);
    for k = 1:numel(irep)
        p = irep(k);
        if isempty(val2over)
            pmats{k} = eye(npols(p));
        else
            pmats{k} = full(val2over(ixyzso(p):ixyzso(p+1)-1, ...
                ixyzs(p):ixyzs(p+1)-1));
        end
    end

    % Layout of the output, indexed by the nodes of S
    iquad = [1; 1 + cumsum(npols(col_ind(:)))];
    nquad = iquad(end) - 1;

    rsc = [];
    rsc.row_ptr = row_ptr;
    rsc.col_ind = col_ind;
    rsc.iquad   = iquad;
    csc = conv_rsc_to_csc(npatches, rsc);

    wsmth = zeros(m*n, nquad);

    % Bound the size of a single kernel block (in entries)
    maxblk = 2^24;

    for p = 1:npatches
        cstart = csc.col_ptr(p);
        cend   = csc.col_ptr(p+1) - 1;
        if cstart > cend
            continue
        end

        % Oversampled source nodes on patch p
        src_inds = ixyzso(p):(ixyzso(p+1)-1);
        npo = numel(src_inds);
        np  = npols(p);
        srcp = [];
        srcp.r = S_over.r(:,src_inds);
        for field = kern.src_fields(:).'
            srcp.(field{1}) = slice_field(S_over.(field{1}), src_inds, S_over.npts);
        end

        % Weighted interpolation: oversampled values -> nodes of patch p
        Wp = kron(wover(src_inds) .* pmats{ikey(p)}, eye(n));  % (n*npo, n*np)

        nchunk = max(1, floor(maxblk / (m*n*npo)));
        for c0 = cstart:nchunk:cend
            cinds = c0:min(c0+nchunk-1, cend);
            targ_inds = csc.row_ind(cinds);
            nt = numel(targ_inds);

            targp = [];
            targp.r = targs.r(:,targ_inds);
            for field = kern.targ_fields(:).'
                targp.(field{1}) = slice_field(targs.(field{1}), targ_inds, ntarg);
            end

            % (m*nt, n*npo) -> (m*nt, n*np)
            B = kern.eval_mask(srcp, targp) * Wp;
            if ~isreal(B) && isreal(wsmth)
                wsmth = complex(wsmth);
            end

            % (a,t) x (b,k) -> (b,a) x (k,t)
            B = reshape(permute(reshape(B, m, nt, n, np), [3 1 4 2]), m*n, np*nt);
            cols = (0:np-1).' + csc.iquad_r(cinds).';
            wsmth(:, cols(:)) = B;
        end
    end
end

function v = slice_field(A, idx, npts)
if size(A, 2) == npts
    v = A(:, idx);
elseif size(A, 1) == npts
    v = A(idx, :).';
else
    error('SLICE_FIELD: no dimension of A matches npts');
end
end
