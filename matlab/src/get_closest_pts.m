function [sxyz, patch_inds, uvsloc, dists, flags] = get_closest_pts(S, targinfo, opts)
%
%  get_closest_pts
%    This subroutine takes in a surface (S) and collection of target points
%    (targinfo) and returns the closest point on the surface (sxyz), the
%    local uv coordinates of that point (uvsloc), the distance of that
%    point (dists), and flags which indicate whether Newton method was
%    successful.
%
%  For efficiency reasons, this routine only returns the closest 
%  points for targets that lie in a tubular neighborhood of the surface. 
%  In particular, suppose c_{j} and r_{j} are the centroid and radius of 
%  patch j, then the tubular neighborhood is the union of the spheres 
%  centered at c_{j} with radius opts.rfac*r_{j}.
%
%  Syntax
%   [sxyz, patch_inds, uvsloc, dists, flags] = get_closest_pts(S,targinfo)
%   [sxyz, patch_inds, uvsloc, dists, flags] = get_closest_pts(S,targinfo,opts)
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * targinfo: target info
%    * opts: options struct
%        opts.tol - tolerance to be used for Newton (1e-10)
%        opts.maxiter - maximum number of iterations for Newton (10)
%        opts.rfac - rfac*r_{j} is the radius of the bounding sphere
%         centered at the patch centroid. Closest points will
%         be determined for only targets that lie in the 
%         union of these spheres (1.25)

%
%  Output arguments:
%    * sxyz: closest point (only relevant if flags(i) = 0)
%    * patchinds: patch number of closest point (only relevant if flags(i) = 0)
%    * uvsloc: local (u,v) coordinates of closest point (only relevant if flags(i) = 0)
%    * dists: distance between target point and closest point (only relevant if flags(i) = 0)
%    * flags: whether or not Newton method succeeded
%        flags(i) = -1, if target not in tubular neighborhood
%        flags(i) = 0, if target in tubular neighborhood and newton for closest
%                      patch succeeded
%        flags(i) = 1, if target in tubular neighborhood and newton for closest
%                      patch failed

    if (nargin < 3)
        opts = [];
    end
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    patch_id = S.patch_id;
    uvs_src = S.uvs_targ;
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;
    npp1 = npatches+1;
    n3 = 3;

    if numel(targinfo) == 1
        targs = extract_targ_array(targinfo);
    else
        targs = targinfo;
    end
    [ndtarg,ntarg] = size(targs);

    sxyz = zeros(ndtarg,ntarg);
    patch_inds = zeros(ntarg,1);
    uvsloc = zeros(2,ntarg);
    dists = inf*ones(ntarg,1);
    flags = zeros(ntarg,1);

    tol = 1e-10;    if isfield(opts, 'tol'), tol = opts.tol; end
    maxiter = 20;   if isfield(opts, 'maxiter'), maxiter = opts.maxiter; end
    rfac = 1.25;    if isfield(opts, 'rfac'), rfac = opts.rfac; end

    mex_id_ = 'get_closest_points(c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[x], c i double[x], c io double[xx], c io int64_t[x], c io double[xx], c io double[x], c io int64_t[x])';
[sxyz, patch_inds, uvsloc, dists, flags] = fmm3dbie_routs(mex_id_, ndtarg, ntarg, targs, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, patch_id, uvs_src, maxiter, tol, rfac, sxyz, patch_inds, uvsloc, dists, flags, 1, 1, ndtarg, ntarg, 1, npatches, npatp1, npatches, 1, 9, npts, 12, npts, npts, 2, npts, 1, 1, 1, 3, ntarg, ntarg, 2, ntarg, ntarg, ntarg);

end
%
%
%
%
