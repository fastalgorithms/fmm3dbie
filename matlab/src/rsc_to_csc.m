function [col_ptr, row_ind, iper] = rsc_to_csc(npatches,ntarg,row_ptr,col_ind,nnz)
% reorder rsc format to csc format
% row_ind(col_ptr(i)):(row_ind(col_ptr(i+1))-1) are the indices of targets near patch i

    col_ptr = zeros(npatches+1,1);
    row_ind = zeros(nnz,1);
    iper = zeros(nnz,1);
    npatp1 = npatches+1;
    mex_id_ = 'rsc_to_csc(i int[x], i int[x], i int[x], i int[x], i int[x], io int[x], io int[x], io int[x])';
[col_ptr, row_ind, iper] = fmm3dbie_routs(mex_id_, npatches, ntarg, nnz, row_ptr, col_ind, col_ptr, row_ind, iper, 1, 1, 1, ntp1, nnz, npatp1, nnz, nnz);
        
end