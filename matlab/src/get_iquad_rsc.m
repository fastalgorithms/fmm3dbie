function iquad = get_iquad_rsc(npatches,ixyzs,npts,ntarg,nnz,row_ptr,col_ind)
% get location in wnear_ij array where quadrature for col_ind(i) starts for a single kernel. 
    iquad = zeros(nnz+1,1);
    npp1 = npatches+1;
    ntp1 = ntarg+1;
    nnzp1 = nnz+1;
    mex_id_ = 'get_iquad_rsc(i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], io int[x])';
[iquad] = fmm3dbie_routs(mex_id_, npatches, ixyzs, npts, nnz, row_ptr, col_ind, iquad, 1, npp1, 1, 1, ntp1, nnz, nnzp1);
end


