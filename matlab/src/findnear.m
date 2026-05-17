function [row_ptr, col_ind] = findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,nnz)
    row_ptr = zeros(ntarg+1,1);
    col_ind = zeros(nnz,1);
    n3 = 3;

    ntp1 = ntarg+1;
    nnzp1 = nnz+1;
    mex_id_ = 'findnear(i double[xx], i int64_t[x], i double[x], i int64_t[x], i double[xx], i int64_t[x], io int64_t[x], io int64_t[x])';
[row_ptr, col_ind] = kern_routs(mex_id_, cms, npatches, rad_near, ndtarg, targs, ntarg, row_ptr, col_ind, n3, npatches, 1, npatches, 1, ndtarg, ntarg, 1, ntp1, nnz);
end

