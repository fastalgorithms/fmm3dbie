function wnear = getnearquad_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker,ivpp)
%SURFWAVE.FLEX.GETNEARQUAD_FLEX near-field quadrature for flexural-gravity
%  wave kernels at arbitrary (possibly off-surface) targets.
%
% Syntax:
%   wnear = surfwave.flex.getnearquad_flex(npatches,norders,ixyzs, ...
%             iptype,npts,srccoefs,srcvals,targs,ipatch_id,uvs_targ,eps, ...
%             iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker,ivpp)
%
% For iker < 6 calls getnearquad_flex_all (or its vpp variant).
% For iker == 6 calls getnearquad_lap_s_neu_eval (Laplace S_{3d}).
% For iker == 7 recurses as iker=5 plus iker=6.
%
% Input:
%   npatches  - number of patches
%   norders   - (npatches,1) polynomial order on each patch
%   ixyzs     - (npatches+1,1) start indices into srccoefs/srcvals
%   iptype    - (npatches,1) patch type (1=RV triangles, 11/12=quads)
%   npts      - total number of discretization points
%   srccoefs  - (9,npts)  basis expansion coefficients
%   srcvals   - (12,npts) xyz, derivatives, and normals at nodes
%   targs     - (ndtarg,ntarg) target coordinates, or struct with field .r
%   ipatch_id - (ntarg,1) patch index for on-surface targets (-1 otherwise)
%   uvs_targ  - (2,ntarg) local uv coordinates for on-surface targets
%   eps       - requested precision
%   iquadtype - quadrature type (1 = ggq + adaptive)
%   nnz       - number of source-patch/target near-field interactions
%   row_ptr   - (ntarg+1,1) row pointers into col_ind
%   col_ind   - (nnz,1) source patch indices
%   iquad     - (nnz+1,1) start indices into wnear for each interaction
%   rfac0     - radius parameter for adaptive quadrature (rec. 1.25)
%   zpars     - complex array of kernel parameters:
%                 zpars(1:5)  = dispersion roots
%                 zpars(6:10) = partial-fraction residues
%   nquad     - total number of near-field quadrature entries
%   iker      - kernel selector:
%                 1 = G_s,  2 = G_phi,  3 = Delta^2 G_s,
%                 4 = Delta^2 G_phi,  5 = S_{3d} G_phi,
%                 6 = Laplace S_{3d},  7 = S_{3d} G_phi + Laplace S_{3d}
%   ivpp      - flag: 1 = use variable-periodization MEX routine, 0 = standard
%
% Output:
%   wnear - (nquad,1) complex near-field quadrature weights
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;
    %if isfield(targs,'r') || isprop(targs,'r')
    if ~isfloat(targs)
        targs = targs.r;
    end
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;
    ntp1 = ntarg+1;
    nnzp1 = nnz+1;
    ndz = length(zpars);

    wnear = zeros(1,nquad,'like',1i);

    if iker < 7

    if iker < 6
    if ivpp
        max_x = max([max(srcvals(1,:)),max(targs(1,:))]);
        min_x = min([min(srcvals(1,:)),min(targs(1,:))]);

        max_y = max([max(srcvals(2,:)),max(targs(2,:))]);
        min_y = min([min(srcvals(2,:)),min(targs(2,:))]);

        maxdist = 1.5*sqrt((max_x - min_x).^2 + (max_y - min_y).^2);
        mex_id_ = 'getnearquad_flex_all_vpp(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i double[x], c io dcomplex[x])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, maxdist, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, 1, nquad);
    else
        mex_id_ = 'getnearquad_flex_all(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io dcomplex[x])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, nquad);
    end
    else
        wnear_lap = zeros(1,nquad);    
        mex_id_ = 'getnearquad_lap_s_neu_eval(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[wnear_lap] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear_lap, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        wnear = wnear_lap;
    end


    else
    wnear1 = surfwave.flex.getnearquad_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,5,ivpp);
    wnear2 = surfwave.flex.getnearquad_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,6,ivpp);    
    wnear = wnear1 + wnear2;
    end
    wnear = wnear(:);
end

