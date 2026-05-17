function wnear = getnearquad(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,type)
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;

    targs = [targs.r;zeros(6,size(targs.r,2));targs.n];
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;
    nnzp1 = nnz+1;
    ndz = length(zpars);

    wnear = zeros(1,nquad,'like',1i);

    if strcmp(type, 'v2v') || strcmp(type, 'v2bd')
    %# FORTRAN getnearquad_helm2d_gv2v(int[1] npatches, int[npatches] norders, int[npp1] ixyzs, int[npatches] iptype, int[1] npts, double[n9,npts] srccoefs, double[n12,npts] srcvals, int[1] ndtarg, int[1] ntarg, double[ndtarg,ntarg] targs, int[ntarg] ipatch_id, double[2,ntarg] uvs_targ, double[1] eps, dcomplex[ndz] zpars, int[1] iquadtype, int[1] nnz, int[ntargp1] row_ptr, int[nnz] col_ind, int[nnzp1] iquad, double[1] rfac0, int[1] nquad, inout dcomplex[1,nquad] wnear);
    mex_id_ = 'getnearquad_helm2d_dir(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
    elseif strcmp(type, 'v2bn')
mex_id_ = 'getnearquad_helm2d_gdn(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[xx])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, ndz, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 1, nquad);
    end

end
%
%
%
