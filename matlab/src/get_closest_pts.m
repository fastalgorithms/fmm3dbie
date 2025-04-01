function [sxyz, uvsloc, dists, flags] = get_closest_pts(S, targinfo)
%
%  get_closest_pts
%    This subroutine takes in a surface (S) and collection of target points
%    (targinfo) and returns the closest point on the surface (sxyz), the
%    local uv coordinates of that point (uvsloc), the distance of that
%    point (dists), and flags which indicate whether Newton method was
%    successful.
%
%  Syntax
%   [sxyz, uvsloc, dists, flags] = get_closest_pts(S,targinfo)
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * targinfo: target info
%
%  Output arguments:
%    * sxyz: closest point
%    * uvsloc: local (u,v) coordinates of closest point on each patch 
%    * dists: distance between target point and closest point
%    * flags: whether or not Newton method succeeded

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;
    npp1 = npatches+1;
    n3 = 3;

    cms = zeros(3,npatches);
    rads = zeros(npatches,1);
    mex_id_ = 'get_centroid_rads(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], io double[xx], io double[x])';
[cms, rads] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, cms, rads, 1, npatches, npp1, npatches, 1, n9, npts, n3, npatches, npatches);

    targs = extract_targ_array(targinfo); 
    [ndtarg,ntarg] = size(targs);

    nnz = 0;
    mex_id_ = 'findnearmem(i double[xx], i int[x], i double[x], i int[x], i double[xx], i int[x], io int[x])';
[nnz] = fmm3dbie_routs(mex_id_, cms, npatches, rads, ndtarg, targs, ntarg, nnz, n3, npatches, 1, npatches, 1, ndtarg, ntarg, 1, 1);

    row_ptr = zeros(ntarg+1,1);
    col_ind = zeros(nnz,1);
    ntp1 = ntarg+1;
    mex_id_ = 'findnear(i double[xx], i int[x], i double[x], i int[x], i double[xx], i int[x], io int[x], io int[x])';
[row_ptr, col_ind] = fmm3dbie_routs(mex_id_, cms, npatches, rads, ndtarg, targs, ntarg, row_ptr, col_ind, n3, npatches, 1, npatches, 1, ndtarg, ntarg, 1, ntp1, nnz);

    istarts = row_ptr(1:end-1);
    iends = row_ptr(2:end)-1;
    isrcinds = cell(npatches,1);
    for i=1:npatches
        isrcinds{i} = ixyzs(i):ixyzs(i+1)-1;
    end

    patch_inds = zeros(ntarg,1);
    local_inds = zeros(ntarg,1);

    for i=1:ntarg
        iinds = horzcat(isrcinds{col_ind(istarts(i):iends(i))});
        dr2 = (targinfo.r(1,i) - S.r(1,iinds)).^2 + (targinfo.r(2,i) - S.r(2,iinds)).^2 + (targinfo.r(3,i) - S.r(3,iinds)).^2;
        [~,mind] = min(dr2);
        glob_ind = iinds(mind);
        patch_inds(i) = S.patch_id(glob_ind);
        local_inds(i) = glob_ind - ixyzs(S.patch_id(glob_ind)) + 1;
    end

    porder = S.norders(1);
    snpols = ixyzs(2) - ixyzs(1);
    maxiter = 10;
    tol = 1e-8;
    lpatchidxvec = local_inds;
    itvec = ixyzs(patch_inds);
    sxyz = zeros(ndtarg,ntarg);
    uvsloc = zeros(2,ntarg);
    dists = zeros(ntarg,1);
    flags = zeros(ntarg,1);

mex_id_ = 'findnearsrfcpnt(i int[x], i double[xx], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[x], i int[x], i int[x], io double[xx], io double[xx], io double[x], io int[x])';
[sxyz, uvsloc, dists, flags] = fmm3dbie_routs(mex_id_, ntarg, targs, npts, porder, snpols, maxiter, srcvals, srccoefs, tol, lpatchidxvec, itvec, sxyz, uvsloc, dists, flags, 1, ndtarg, ntarg, 1, 1, 1, 1, n12, npts, n9, npts, 1, ntarg, ntarg, ndtarg, ntarg, 2, ntarg, ntarg, ntarg);
    

end
%
%

