function [pot] = lpcomp_helm_comb_dir_boundary(S,zpars,sigma,targs,eps)
    [ndtarg,ntarg] = size(targs);
    npts = S.npts;
    npatches = S.npatches;
    ipatch_id = zeros(ntarg,1);
    uvs_targ = zeros(2,ntarg);
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
    
    pot = zeros(ntarg,1);
    mex_id_ = 'lpcomp_helm_comb_dir(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i int[x], i double[xx], i double[x], i dcomplex[x], i dcomplex[x], io dcomplex[x])';
[pot] = fmm3dbierouts(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, sigma, pot, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 3, npts, ntarg);
    

end
