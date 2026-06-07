function Asmth = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,nover,lbat)
%
%  SMOOTH_SPARSE_QUAD
%
%  Syntax:
%    Asmth = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,nover)
%    Asmth = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,nover,lbat)
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

    if nargin < 7
        lbat = 1e3;
    end
    try
        kernuse = kern.eval;
    catch
        kernuse = kern;
    end
    try
        targs.r; 
    catch 
        targs = struct('r',targs(1:3,:)); 
    end

    [S_over,val2over] = oversample(S,nover);

    ixyzs = S_over.ixyzs;
    npols = ixyzs(2:end)-ixyzs(1:end-1);
    npts_col = npols(col_ind);

    [nt,~] = size(row_ptr);
    ntarg = nt-1;
    nrep = zeros(ntarg,1);
    istarts = row_ptr(1:end-1);
    iends = row_ptr(2:end)-1;
    icol_ind = zeros(sum(npts_col),1);
    isrcinds = cell(S_over.npatches,1);
    for i=1:S_over.npatches
        isrcinds{i} = ixyzs(i):ixyzs(i+1)-1;
    end

    istart = 1;
    for i=1:ntarg
        nrep(i) = sum(npts_col(istarts(i):iends(i))); 
        iinds = horzcat(isrcinds{col_ind(istarts(i):iends(i))});
        nelem = length(iinds);
        icol_ind(istart:istart+nelem-1) = iinds;
        istart = istart+nelem;
    end
    irow_ind = repelem((1:ntarg)',nrep);

    nnz = length(irow_ind);
    
    nbat = ceil(nnz/lbat);
    
    loc_field = fieldnames(targs)';
    loc_field_use = {};
    for i = 1:length(loc_field)
        if mod(size(targs.(loc_field{i})(:,:),2), size(targs.r(:,:),2)) == 0 && ~isempty(targs.(loc_field{i}))
        loc_field_use{length(loc_field_use)+1} = loc_field{i};
        end
    end

    targdim = [];
    for field = loc_field_use
        targdim.(field{1}) = targs.(field{1})(:,1);
    end
    tmp = kernuse(struct('r',[0;0;0]),targdim);
    opdims = size(tmp);


    irow_indopdim = (1:opdims(1)).' + opdims(1)*(irow_ind(:).'-1);
    icol_indopdim = 0*(1:opdims(1)).' + (icol_ind.');
    irow_indopdim = irow_indopdim(:);
    icol_indopdim = icol_indopdim(:);
    irow_indopdim = 0*(1:opdims(2)).' + (irow_indopdim.');
    icol_indopdim = (1:opdims(2)).' + opdims(2)*(icol_indopdim.'-1);
    irow_indopdim = irow_indopdim(:);
    icol_indopdim = icol_indopdim(:);

    vals = zeros(1, prod(opdims)*nnz);

    for k = 1:nbat
        ks = (lbat*(k-1)+1):min(lbat*k,nnz);
        ksloc = (prod(opdims)*lbat*(k-1)+1):(prod(opdims)*min(lbat*k,nnz));
        rsrc = S_over.r(:,icol_ind(ks));
        targk = [];
        for field = loc_field_use
            targk.(field{1}) = targs.(field{1})(:,irow_ind(ks));
        end
        targk.r = [targk.r;zeros(3-size(targk.r,1),size(targk.r,2))];
        targk.r = targk.r - rsrc;
        kernvals = reshape(kernuse(struct('r',[0;0;0]),targk), opdims(1), length(ks), opdims(2)) ...
                   .* reshape(S_over.wts(icol_ind(ks)), 1, length(ks), 1);
        % Zero out entries where source and target coincide (relative norm < 1e-14)
        src_norm = max(vecnorm(rsrc), 1);
        self_mask = vecnorm(targk.r) < 1e-14 * src_norm;
        kernvals(:, self_mask, :) = 0;
        % Permute to (opdims(1), opdims(2), length(ks)) and flatten
        vals(ksloc) = reshape(permute(kernvals, [1 3 2]), 1, []);
    end

   Asmth = sparse(irow_indopdim,icol_indopdim, vals, opdims(1)*ntarg, opdims(2)*S_over.npts);

   Asmth = Asmth*kron(val2over, speye(opdims(2)));
end
