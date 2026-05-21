function wnear = getnearquad_gravity(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker,ivpp,S3d_scal)
%SURFWAVE.GRAVITY.GETNEARQUAD_GRAVITY near-field quadrature for gravity wave kernels
%
% Syntax:
%   wnear = surfwave.gravity.getnearquad_gravity(npatches,norders,ixyzs, ...
%             iptype,npts,srccoefs,srcvals,targs,ipatch_id,uvs_targ,eps, ...
%             iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker,ivpp,S3d_scal)
%
% Computes the near-field quadrature correction wnear for a single
% gravity wave kernel selected by iker, combined with a scaled Laplace
% single-layer (S3d) contribution.  Calls the Fortran MEX routines
% getnearquad_gravity_all or getnearquad_gravity_all_vpp depending on ivpp.
%
% Input:
%   npatches  - number of patches
%   norders   - (npatches,1) polynomial order on each patch
%   ixyzs     - (npatches+1,1) start indices into srccoefs/srcvals
%   iptype    - (npatches,1) patch type (1=RV triangles, 11/12=quads)
%   npts      - total number of discretization points
%   srccoefs  - (9,npts) basis expansion coefficients (xyz, dxyz/du, dxyz/dv)
%   srcvals   - (12,npts) xyz, derivatives, and normals at nodes
%   targs     - (ndtarg,ntarg) target point coordinates, or a struct with field .r
%   ipatch_id - (ntarg,1) patch index for each target (-1 if off-surface)
%   uvs_targ  - (2,ntarg) local uv coordinates for on-surface targets
%   eps       - requested precision
%   iquadtype - quadrature type (1 = ggq + adaptive)
%   nnz       - number of source-patch/target interactions in the near field
%   row_ptr   - (ntarg+1,1) row pointers into col_ind
%   col_ind   - (nnz,1) source patch indices for all near-field interactions
%   iquad     - (nnz+1,1) start indices into wnear for each interaction
%   rfac0     - radius parameter for switching to adaptive quadrature (rec. 1.25)
%   zpars     - complex array of kernel parameters (dispersion root and residue)
%   nquad     - total number of near-field quadrature entries
%   iker      - kernel selector: 0 = G_s, 1 = G_phi
%   ivpp      - flag: 1 = use variable-periodization (vpp) Fortran routine,
%               0 = standard routine
%   S3d_scal  - scalar weight for the Laplace S3d correction added to wnear
%
% Output:
%   wnear     - (nquad,1) complex near-field quadrature weights
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;
    if isfield(targs,'r') || isprop(targs,'r')
        targs = targs.r;
    end
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;
    nnzp1 = nnz+1;
    ndz = length(zpars);

    wnear = zeros(1,nquad,'like',1i);

    if ivpp
        max_x = max([max(srcvals(1,:)),max(targs(1,:))]);
        min_x = min([min(srcvals(1,:)),min(targs(1,:))]);

        max_y = max([max(srcvals(2,:)),max(targs(2,:))]);
        min_y = min([min(srcvals(2,:)),min(targs(2,:))]);

        maxdist = 1.5*sqrt((max_x - min_x).^2 + (max_y - min_y).^2);
        mex_id_ = 'getnearquad_gravity_all_vpp(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i double[x], c io dcomplex[xx])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, maxdist, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, ndz, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 1, 1, 1, nquad);
    else
        mex_id_ = 'getnearquad_gravity_all(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io dcomplex[xx])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, ndz, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 1, 1, nquad);
    end
    wnear_lap = zeros(1,nquad);

    mex_id_ = 'getnearquad_lap_s_neu_eval(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[wnear_lap] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear_lap, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
    wnear = wnear(:) + S3d_scal*wnear_lap(:);
end
%
%
