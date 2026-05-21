function [xmat1,xmat2,wnear] = get_quad_cor_v2b(S, targ, kern, eps, zpars, uv_bndry, patch_id)
%SURFWAVE.FLEX.GET_QUAD_COR_V2B near-field quadrature correction matrices
%  for the flexural-gravity wave volume-to-boundary (v2b) kernels.
%
% Syntax:
%   [xmat1,xmat2,wnear] = surfwave.flex.get_quad_cor_v2b( ...
%       S, targ, kern, eps, zpars)
%   [xmat1,xmat2,wnear] = surfwave.flex.get_quad_cor_v2b( ...
%       S, targ, kern, eps, zpars, uv_bndry, patch_id)
%
% Computes the near-field quadrature correction for the two v2b boundary
% condition kernels (v2b bc kernel 1 and v2b bc kernel 2) by calling the
% Fortran MEX routine getnearquad_flex_bcs. The smooth oversampled
% quadrature is subtracted from each correction matrix via smooth_sparse_quad.
%
% Input:
%   S        - surfer object describing the surface discretization
%   targ     - target geometry struct with fields:
%                targ.r     - (2,ntarg) positions
%                targ.n     - (2,ntarg) unit outward normals
%                targ.d     - (2,ntarg) parameterization tangent derivative
%                targ.kappa - (ntarg,1) signed curvature
%   kern     - (2,1) cell array of kernel function handles for
%              smooth_sparse_quad subtraction:
%                kern{1} = v2b bc kernel 1 handle
%                kern{2} = v2b bc kernel 2 handle
%   eps      - requested quadrature precision
%   zpars    - complex array of kernel parameters:
%                zpars(1:5)  = dispersion roots
%                zpars(6:10) = partial-fraction residues
%   uv_bndry - (optional) (2,ntarg) local uv coordinates for on-surface
%              targets; default zeros(2,ntarg)
%   patch_id - (optional) (ntarg,1) patch index for on-surface targets
%              (-1 for off-surface); default zeros(ntarg,1)
%
% Output:
%   xmat1 - (ntarg, S.npts) sparse complex correction matrix for v2b bc kernel 1
%   xmat2 - (ntarg, S.npts) sparse complex correction matrix for v2b bc kernel 2
%   wnear - (2,nquad) raw near-field quadrature weights before matrix assembly


    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npatches = S.npatches;
    npp1 = npatches+1;
    
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

    % [patch_id, uvs_targ] = get_patch_id_uvs(S);

    %this might need fixing

    iptype_avg = floor(sum(iptype)/(npatches+0.0d0));
    norder_avg = floor(sum(norders)/(npatches+0.0d0));

    [rfac, rfac0] = get_rfacs(norder_avg,iptype_avg);

    [cms, rads] = get_centroid_rads(npatches,norders,ixyzs,iptype,npts, ... 
       srccoefs);

    rad_near = rads*rfac;

    % 
    % find near quadrature correction interactions
    % 

    nnz = findnearmem(cms, npatches, rad_near, ndtarg, targs, ntarg);
    nnzp1 = nnz+1;

    [row_ptr,col_ind] = findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,nnz);

    iquad = get_iquad_rsc(npatches,ixyzs,npts,ntarg,nnz,row_ptr,col_ind);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wts = get_qwts_f(npatches,norders,ixyzs,iptype,ntarg,srcvals);

    targinfo = [];
    targinfo.r = targs(1:3,:);

    if nargin < 6
        ipatch_id = zeros(ntarg,1);
        uvs_targ = zeros(2,ntarg);
    else
        ipatch_id = patch_id; 
        uvs_targ = uv_bndry;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    nquad = iquad(nnz+1)-1;

    iquadtype = 1;

    wnear = zeros(2,nquad,'like',1i);

    mex_id_ = 'getnearquad_flex_bcs(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[xx])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 2, nquad);

    xmat1 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(1,:).');
    xmat2 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(2,:).');

    nover = S.norders(1) + 4;

    xmat1 = xmat1 - smooth_sparse_quad(kern(1), targinfo, S, row_ptr, col_ind, nover);
    xmat2 = xmat2 - smooth_sparse_quad(kern(2), targinfo, S, row_ptr, col_ind, nover);

end
%
%
%
%
%
%
%
%
%
