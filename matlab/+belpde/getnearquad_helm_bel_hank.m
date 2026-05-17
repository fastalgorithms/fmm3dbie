function wnear = getnearquad_helm_bel_hank(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,eps,iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,mean_curv,zk,nquad,iktype)
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;
    ntarg = npts;
    ntp1 = ntarg+1;
    nnzp1 = nnz+1;

    wnear = zeros(nquad,1,'like',1i);

    mex_id_ = 'getnearquad_helm_bel_hank(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c io dcomplex[x])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, mean_curv, nquad, zk, iktype, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, 1, ntp1, nnz, nnzp1, 1, npts, 1, 1, 1, nquad);

end
%
%
%
%----------------------------------
%  SURFWAVE routs
%----------------------------------
%
%
