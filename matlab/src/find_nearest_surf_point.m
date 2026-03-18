function [nearest_pnts, nearest_uvs,dists,flags] = find_nearest_surf_point(S,targs,idx_init,tol,maxiter)

    [srcvals,srccoefs,norders,ixyzs,~,~] = extract_arrays(S);
        
    npts = S.npts;

    
    ndtarg = size(targs,1); 
    ntarg = size(targs,2);

    n2 = 2;
    n9 = size(srccoefs,1);
    n12 = size(srcvals,1);

%   assume the same order for all patches
    norder = S.norders(1);
%   assume triangular patch with VR-nodes    
    npols = (norder+1)*(norder+2)/2;

% redo with ixyzs
%    idx_init_patch_pnt = mod(idx_init + (npols-1),npols) + 1
    %idx_init_patch_idx = idx_init - idx_init_patch_pnt + 1;
    idx_init_patch_pnt = mod(idx_init-1,npols) + 1;
    idx_init_patch_idx = ixyzs(S.patch_id(idx_init));

% itvec = ixyzs, lpatchidxvec = index on that patch

% output
  nearest_pnts = zeros(ndtarg,ntarg);
  nearest_uvs = zeros(2,ntarg);
  dists = zeros(ntarg,1);
  flags = zeros(ntarg,1);

  mex_id_ = 104;
[nearest_pnts, nearest_uvs, dists, flags] = fmm3dbie_routs(mex_id_, ntarg, targs, npts, norder, npols, maxiter, srcvals, srccoefs, tol, idx_init_patch_pnt, idx_init_patch_idx, nearest_pnts, nearest_uvs, dists, flags, 1, ndtarg, ntarg, 1, 1, 1, 1, n12, npts, n9, npts, 1, ntarg, ntarg, ndtarg, ntarg, n2, ntarg, ntarg, ntarg);

return

end