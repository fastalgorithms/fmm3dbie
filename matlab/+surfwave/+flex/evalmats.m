function wnear = evalmats(S,targinfo,zpars,eps,ipatch_id,uvs_targ)
%SURFWAVE.FLEX.EVALMATS dense near-field quadrature weights for the
%  flexural-gravity wave post-processing kernels at off-surface targets.
%
% Syntax:
%   wnear = surfwave.flex.evalmats(S, targinfo, zpars, eps)
%   wnear = surfwave.flex.evalmats(S, targinfo, zpars, eps, ipatch_id, uvs_targ)
%
% Builds the full (all target-source pairs) quadrature weight array for the
% three evaluation kernels G_s, G_phi, and S_{3d}G_phi
%
% Input:
%   S         - surfer object describing the surface discretization
%   targinfo  - target geometry struct or array:
%                 targinfo.r  - (3,ntarg) target positions (required)
%   zpars     - complex array of kernel parameters:
%                 zpars(1:5)  = dispersion roots
%                 zpars(6:10) = partial-fraction residues
%   eps       - requested quadrature precision
%   ipatch_id - (optional) (ntarg,1) patch index for on-surface targets
%               (-1 for off-surface); default zeros(ntarg,1)
%   uvs_targ  - (optional) (2,ntarg) local uv coordinates for on-surface
%               targets; default zeros(2,ntarg)
%
% Output:
%   wnear - (3,nquad) complex near-field quadrature weights:
%             wnear(1,:) = weights for G_s
%             wnear(2,:) = weights for G_phi
%             wnear(3,:) = weights for S_{3d}G_phi
%

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
    ntp1 = npts + 1;
    npols = ixyzs(2) - ixyzs(1);

    rfac0 = 1.25;

    if nargin < 5
    ipatch_id = zeros(ntarg,1);
    uvs_targ = zeros(2,ntarg);
    end

    iquad = 1:npols:(npts*ntarg+1);
    nquad = iquad(nnzp1)-1;
    wnear = zeros(3,nquad);

    mex_id_ = 'getnearquad_flex_eval(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[xx])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 13, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 3, nquad);


end
%
%
