function iquad = get_iquad_rsc(npatches,ixyzs,npts,ntarg,nnz,row_ptr,col_ind)
    iquad = zeros(nnz+1,1);
    npp1 = npatches+1;
    ntp1 = ntarg+1;
    nnzp1 = nnz+1;
    mex_id_ = 'get_iquad_rsc(i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], io int64_t[x])';
[iquad] = kern_routs(mex_id_, npatches, ixyzs, npts, nnz, row_ptr, col_ind, iquad, 1, npp1, 1, 1, ntp1, nnz, nnzp1);
end


%
%
%
