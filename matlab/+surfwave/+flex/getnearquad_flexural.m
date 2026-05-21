function wnear = getnearquad_flexural(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,eps,iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker)
%SURFWAVE.FLEX.GETNEARQUAD_FLEXURAL near-field quadrature for flexural-gravity
%  wave kernels, on-surface targets only.
%
% Syntax:
%   wnear = surfwave.flex.getnearquad_flexural(npatches,norders,ixyzs, ...
%             iptype,npts,srccoefs,srcvals,eps,iquadtype,nnz, ...
%             row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker)
%
% Calls getnearquad_flex_all (iker < 7) or recurses as iker=5 plus iker=6
% (iker == 7) to compute the near-field quadrature weights for an
% on-surface evaluation (targets = source nodes).
%
% Input:
%   npatches  - number of patches
%   norders   - (npatches,1) polynomial order on each patch
%   ixyzs     - (npatches+1,1) start indices into srccoefs/srcvals
%   iptype    - (npatches,1) patch type (1=RV triangles, 11/12=quads)
%   npts      - total number of discretization points
%   srccoefs  - (9,npts)  basis expansion coefficients
%   srcvals   - (12,npts) xyz, derivatives, and normals at nodes
%   eps       - requested precision
%   iquadtype - quadrature type (1 = ggq + adaptive)
%   nnz       - number of source-patch/target interactions in near field
%   row_ptr   - (npts+1,1) row pointers into col_ind
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
%
% Output:
%   wnear - (nquad,1) complex near-field quadrature weights
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;
    ntarg = npts;
    ntp1 = ntarg+1;
    nnzp1 = nnz+1;

    if iker < 7
    wnear = zeros(nquad,1,'like',1i);

    mex_id_ = 'getnearquad_flex_all(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io dcomplex[x])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, nquad);
    else
    wnear1 = surfwave.flex.getnearquad_flexural(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,eps,iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,5);
    wnear2 = surfwave.flex.getnearquad_flexural(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,eps,iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,6);    
    wnear = wnear1 + wnear2;
    end
end
%
%
%
%
%
