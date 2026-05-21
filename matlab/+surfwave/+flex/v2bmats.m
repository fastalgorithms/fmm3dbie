function wnear = v2bmats(S,targinfo,zpars,eps)
%SURFWAVE.FLEX.V2BMATS dense near-field quadrature weights for the
%  flexural-gravity wave volume-to-boundary (v2b) boundary condition kernels.
%
% Syntax:
%   wnear = surfwave.flex.v2bmats(S, targinfo, zpars, eps)
%
% Builds the full (all target-source pairs) quadrature weight array for the
% two v2b boundary condition kernels by calling the Fortran MEX routine
% getnearquad_flex_bcs with a fully-dense row_ptr/col_ind/iquad structure.
% The target struct must carry position, normal, tangent-derivative, and
% curvature information; these are packed into a 13-row array internally
% (rows 1:2 = position, 4:5 = unit tangent, 10:11 = unit outward normal,
% 13 = signed curvature).
%
% Input:
%   S        - surfer object describing the surface discretization
%   targinfo - target geometry struct with fields:
%                targinfo.r     - (2,ntarg) positions
%                targinfo.n     - (2,ntarg) unit outward normals
%                targinfo.d     - (2,ntarg) parameterization tangent derivative
%                targinfo.kappa - (ntarg,1) signed curvature
%   zpars    - complex array of kernel parameters:
%                zpars(1:5)  = dispersion roots
%                zpars(6:10) = partial-fraction residues
%   eps      - requested quadrature precision
%
% Output:
%   wnear - (2,nquad) complex near-field quadrature weights:
%             wnear(1,:) = weights for v2b bc kernel 1
%             wnear(2,:) = weights for v2b bc kernel 2

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npp1 = npatches+1;

    [targs] = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;

    row_ptr = 1:npatches:(npatches*ntarg+1);
    col_ind = repmat(1:npatches,[1,ntarg]);
    row_ptr = row_ptr(:);
    col_ind = col_ind(:);
    nnz = npatches*ntarg;
    nnzp1 = nnz + 1;
    iquadtype = 1;
    ntp1 = ntarg + 1;
    npols = ixyzs(2) - ixyzs(1);

    rfac0 = 1.25;

    ipatch_id = zeros(ntarg,1);
    uvs_targ = zeros(2,ntarg);

    iquad = 1:npols:(npts*ntarg+1);
    nquad = iquad(nnzp1)-1;
    wnear = zeros(2,nquad);

    mex_id_ = 'getnearquad_flex_bcs(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[xx])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 2, nquad);

    
end
%
%
%
%
