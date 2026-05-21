function wnear = getnearquad_v2b_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targ,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad)
%SURFWAVE.FLEX.GETNEARQUAD_V2B_FLEX near-field quadrature for the flexural-gravity
%  wave volume-to-boundary (v2b) boundary condition kernels.
%
% Syntax:
%   wnear = surfwave.flex.getnearquad_v2b_flex(npatches,norders,ixyzs, ...
%             iptype,npts,srccoefs,srcvals,targ,ipatch_id,uvs_targ,eps, ...
%             iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,zpars,nquad)
%
% Calls getnearquad_flex_bcs, which computes quadrature weights for both
% v2b boundary condition kernels simultaneously. 
%
% Input:
%   npatches  - number of patches
%   norders   - (npatches,1) polynomial order on each patch
%   ixyzs     - (npatches+1,1) start indices into srccoefs/srcvals
%   iptype    - (npatches,1) patch type (1=RV triangles, 11/12=quads)
%   npts      - total number of discretization points
%   srccoefs  - (9,npts)  basis expansion coefficients
%   srcvals   - (12,npts) xyz, derivatives, and normals at nodes
%   targ      - target geometry struct with fields:
%                 targ.r     - (2,ntarg) positions
%                 targ.n     - (2,ntarg) unit outward normals
%                 targ.d     - (2,ntarg) parameterization tangent derivative
%                 targ.kappa - (ntarg,1) signed curvature
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
%
% Output:
%   wnear - (2,nquad) complex near-field quadrature weights:
%             wnear(1,:) = weights for v2b bc kernel 1
%             wnear(2,:) = weights for v2b bc kernel 2

    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;
    nnzp1 = nnz+1;
    ndz = length(zpars);

    
    rx = targ.r(1,:); 
    ry = targ.r(2,:); 
    
    nx = targ.n(1,:); 
    ny = targ.n(2,:);
    
    dx = targ.d(1,:);
    dy = targ.d(2,:);
    
    ds = sqrt(dx.*dx+dy.*dy);
    taux = (dx./ds); % normalization
    tauy = (dy./ds);
    
    targs = zeros(13,length(targ.r(1,:)));
    targs(1:2,:) = [rx; ry];
    targs(4:5,:) = [taux; tauy];
    targs(10:11,:) = [nx; ny];
    targs(13,:) = targ.kappa(:);

    [ndtarg,ntarg] = size(targs);
    ntp1 = ntarg+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    wnear = zeros(2,nquad,'like',1i);

    mex_id_ = 'getnearquad_flex_bcs(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[xx])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 2, nquad);


end


%
%
