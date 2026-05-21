function [xmat,novers] = get_quad_cor_sub(S, zk, eps)

    helm2d_kern = kernel('h','s',zk);

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npatches = S.npatches;
    ndtarg = 3;
    npp1 = npatches+1;

    iptype_avg = floor(sum(iptype)/(npatches+0.0d0));
    norder_avg = floor(sum(norders)/(npatches+0.0d0));

    [rfac, rfac0] = get_rfacs(norder_avg,iptype_avg);

    [cms, rads] = get_centroid_rads(npatches,norders,ixyzs,iptype,npts, ...
       srccoefs);

    rad_near = rads*rfac;

    targs = srcvals(1:3,:);
    ntarg = npts;
    ntp1 = ntarg+1;

    nnz = findnearmem(cms, npatches, rad_near, ndtarg, targs, npts);
    nnzp1 = nnz+1;

    [row_ptr,col_ind] = findnear(cms,npatches,rad_near,ndtarg,targs,npts,nnz);

    iquad = get_iquad_rsc(npatches,ixyzs,npts,ntarg,nnz,row_ptr,col_ind);

    nquad = iquad(nnz+1)-1;

    iquadtype = 1;

    A = zeros(nquad,1,'like',1);

    zpars = complex(double(real(zk)), double(imag(zk)));

    mex_id_ = 'getnearquad_helm2d_gv2v(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);

    xmat = conv_rsc_to_spmat(S, row_ptr, col_ind, A);

    Q = []; Q.targinfo = S; Q.rfac = rfac; Q.wavenumber = zk;
    Q.row_ptr = row_ptr; Q.col_ind = col_ind; Q.kernel_order = -1;
    novers = get_oversampling_parameters(S,Q,eps);
    
    Asmth_over = smooth_sparse_quad(helm2d_kern,targs,S,row_ptr,col_ind,novers);

    xmat = xmat - Asmth_over;

end

