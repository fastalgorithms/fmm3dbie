function wnear = getnearquad_radpp_d(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id,uvs_targ,eps, ...
          ipvs,nnz,row_ptr,col_ind,iquad,rfac0,nquad,dpars,ipars)
%
%  Near quadrature for the length(ipvs) rows of a real radpp kernel.
%  Returns wnear (nrow,nquad).
%
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;
    nnzp1 = nnz+1;
    nrow = length(ipvs);
    ndd = length(dpars);
    ndi = length(ipars);

    wnear = zeros(nrow,nquad);
    mex_id_ = 'getnearquad_radpp_d(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io double[xx])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, nrow, ipvs, nnz, row_ptr, col_ind, iquad, rfac0, nquad, ndd, dpars, ndi, ipars, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, nrow, 1, ntargp1, nnz, nnzp1, 1, 1, 1, ndd, 1, ndi, nrow, nquad);
end

%
%

