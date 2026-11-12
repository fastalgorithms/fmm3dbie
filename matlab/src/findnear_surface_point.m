function [nearest_pnts,nearest_uvs,dists,flags] = ...
    findnear_surface_point(S,targs,idx_init,tol,maxiter)

    [srcvals,srccoefs,norders,ixyzs,iptype,~] = extract_arrays(S);

    ipatch_id = S.patch_id;
    uvs_src   = S.uvs_targ;
    ipts_init = idx_init(:);

    npts     = S.npts;
    npatches = S.npatches;
    npatp1   = npatches + 1;

    [ndtarg,ntarg] = size(targs);

    nearest_pnts = zeros(3,ntarg);
    nearest_uvs  = zeros(2,ntarg);
    dists        = zeros(ntarg,1);
    flags        = zeros(ntarg,1);

    mex_id_ = 107;
[nearest_pnts, nearest_uvs, dists, flags] = fmm3dbie_routs(mex_id_, ntarg, ndtarg, targs, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ipatch_id, uvs_src, ipts_init, maxiter, tol, nearest_pnts, nearest_uvs, dists, flags, 1, 1, ndtarg, ntarg, 1, npatches, npatp1, npatches, 1, 9, npts, 12, npts, npts, 2, npts, ntarg, 1, 1, 3, ntarg, 2, ntarg, ntarg, ntarg);

return

end