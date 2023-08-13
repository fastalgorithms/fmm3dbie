function [spmat] = helm_trans_near_corr(S,zpars,eps)
    npts = S.npts;
    npatches = S.npatches;
    ipatch_id = zeros(npts,1);
    uvs_targ = zeros(2,npts);
    norders = S.norders;
    ixyzs = S.ixyzs;
    iptype = S.iptype;
    srcvals = S.srcvals;
    srccoefs = S.srccoefs;
    npp1 = npatches+1;

    n3 = 3;
    n9 = 9;
    n12 = 12;
    mex_id_ = 'get_patch_id_uvs(i int[x], i int[x], i int[x], i int[x], i int[x], io int[x], io double[xx])';
[ipatch_id, uvs_targ] = fmm3dbierouts(mex_id_, npatches, norders, ixyzs, iptype, npts, ipatch_id, uvs_targ, 1, npatches, npp1, npatches, 1, npts, 2, npts);
    
    iptype0 = iptype(1);
    norder0 = norders(1);
    rfac = 0.0;
    rfac0 = 0.0;
    mex_id_ = 'get_rfacs(i int[x], i int[x], io double[x], io double[x])';
[rfac, rfac0] = fmm3dbierouts(mex_id_, norder0, iptype0, rfac, rfac0, 1, 1, 1, 1);
    

    cms = zeros(3,npatches);
    rads = zeros(npatches,1);
    mex_id_ = 'get_centroid_rads(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], io double[xx], io double[x])';
[cms, rads] = fmm3dbierouts(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, cms, rads, 1, npatches, npp1, npatches, 1, n9, npts, n3, npatches, npatches);

    rad_near = rads*rfac;
    ndtarg = 12;
    nnz = 0;
    mex_id_ = 'findnearmem(i double[xx], i int[x], i double[x], i int[x], i double[xx], i int[x], io int[x])';
[nnz] = fmm3dbierouts(mex_id_, cms, npatches, rad_near, ndtarg, srcvals, npts, nnz, n3, npatches, 1, npatches, 1, ndtarg, npts, 1, 1);

    row_ptr = zeros(npts+1,1);
    col_ind = zeros(nnz,1);
    nptsp1 = npts+1;
    nnzp1 = nnz+1;
    mex_id_ = 'findnear(i double[xx], i int[x], i double[x], i int[x], i double[xx], i int[x], io int[x], io int[x])';
[row_ptr, col_ind] = fmm3dbierouts(mex_id_, cms, npatches, rad_near, ndtarg, srcvals, npts, row_ptr, col_ind, n3, npatches, 1, npatches, 1, ndtarg, npts, 1, nptsp1, nnz);

    iquad = zeros(nnz+1,1);
    mex_id_ = 'get_iquad_rsc(i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], io int[x])';
[iquad] = fmm3dbierouts(mex_id_, npatches, ixyzs, npts, nnz, row_ptr, col_ind, iquad, 1, npp1, 1, 1, nptsp1, nnz, nnzp1);

    nquad = iquad(nnz+1)-1;
    wnear = complex(zeros(nquad,4));
    irowind = zeros(nquad,1);
    icolind = zeros(nquad,1);
    iquadtype = 1;
    
    mex_id_ = 'getnearquad_helm_comb_trans_spmat(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[x], i int[x], io dcomplex[xx], io int[x], io int[x])';
[wnear, irowind, icolind] = fmm3dbierouts(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, irowind, icolind, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 5, 1, 1, nptsp1, nnz, nnzp1, 1, 1, nquad, 4, nquad, nquad);
    spmat = cell(1,4);
    spmat{1} = sparse(irowind,icolind,wnear(:,1),npts,npts);
    spmat{2} = sparse(irowind,icolind,wnear(:,2),npts,npts);
    spmat{3} = sparse(irowind,icolind,wnear(:,3),npts,npts);
    spmat{4} = sparse(irowind,icolind,wnear(:,4),npts,npts);
end

