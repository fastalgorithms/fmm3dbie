function [csc] = conv_rsc_to_csc(npatches, rsc)
%
%  CONV_RSC_TO_CSC(npatches, rsc)
%
%
%  Syntax:  
%    [csc] = conv_rsc_to_csc(npatches, rsc)
%
%  This subroutine reorders the sparse indices of 
%  a row sparse compressed (rsc) representation to the 
%  column sparse compressed (csc) representation. 
%
%  Input arguments:
%    * npatches: number of columns in the sparse representation
%    * rsc: row sparse compressed struct
%       rsc.row_ptr: (nt+1,1)
%         pointer to col_ind array where list of relevant source
%         patches for target i starts
%       rsc.col_ind: (nnz,1)
%         list of source patches relevant for all targets,
%         sorted by target number
%       rsc.iquad: (nnz+1,1)
%         The rsc format stores relevant
%         indices going from patches on the surface to the target
%         while the quadrature corrections are sparse corrections
%         from discretization points on the surface to the target. 
%         iquad(i) is the location in quadrature correction array 
%         where quadrature for interaction corresponding to
%         col_ind(i) starts.
%
%  Output arguments:
%    * csc: column sparse compressed struct
%       csc.col_ptr: (npatches+1,1)
%         pointer to row_ind array where list of relevant targets
%         for patch i start
%       csc.row_ind: (nnz,1)
%         list of targets relevant for all patches, sorted by patch
%         number
%       csc.iquad_r: (nnz+1,1)
%         rsc.iquad(iper) where iper is the permutation that goes
%         from rsc.row_ind to csc.col_ind
%

    row_ptr = rsc.row_ptr;
    col_ind = rsc.col_ind;
    ntarg = size(row_ptr,1)-1;
    nnz = length(col_ind);
    [jsort, iper] = sort(col_ind);

    row_ind_exp = zeros(nnz,1);
    for i = 1:ntarg
        row_ind_exp(row_ptr(i):(row_ptr(i+1)-1)) = i;
    end
    row_ind = row_ind_exp(iper);

    col_ptr = zeros(npatches+1,1);
    
    icur = 1;
    col_ptr(1) = 1;
    for i = 1:nnz
        if (jsort(i) > icur) 
            icur0 = icur;
            icur = jsort(i);
            col_ptr((icur0+1):icur) = i;
        end
    end
    col_ptr((icur+1):end) = nnz+1;

    csc = [];
    csc.col_ptr = col_ptr;
    csc.row_ind = row_ind;
    csc.iquad_r = rsc.iquad(iper);
end
