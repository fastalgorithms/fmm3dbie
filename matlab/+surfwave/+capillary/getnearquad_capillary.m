function wnear = getnearquad_capillary(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,eps,iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker)
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;
    ntarg = npts;
    ntp1 = ntarg+1;
    nnzp1 = nnz+1;

    wnear = zeros(nquad,1,'like',1i);

    mex_id_ = 'getnearquad_capillary_all(i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i double[xx], i double[x], i dcomplex[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[x], i int64_t[x], i int64_t[x], io dcomplex[x])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 6, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, nquad);
end
%
%
%
%
%
%
