function [xmat1,xmat2,xmat3,xmat4] = get_quad_cor_targ(S, targinfo, gkerns, eps, zpars,ipatch_id,uvs_targ)
%SURFWAVE.FLEX.GET_QUAD_COR_TARG near-field quadrature correction matrices
%  for flexural-gravity wave post-processing kernels at off-surface targets.
%
% Syntax:
%   [xmat1,xmat2,xmat3,xmat4] = surfwave.flex.get_quad_cor_targ( ...
%       S, targinfo, gkerns, eps, zpars)
%   [xmat1,xmat2,xmat3,xmat4] = surfwave.flex.get_quad_cor_targ( ...
%       S, targinfo, gkerns, eps, zpars, ipatch_id, uvs_targ)
%
% Computes four sparse near-field quadrature correction matrices for the
% evaluation kernels G_s, G_phi, S_{3d}G_phi, and Laplace S_{3d} at
% arbitrary (typically off-surface) targets.  Calls the Fortran MEX routine
% getnearquad_flex_eval for the first three kernels simultaneously and
% lap3d.neumann.get_quadrature_correction for the fourth, then subtracts
% the smooth oversampled quadrature via smooth_sparse_quad.
%
% Input:
%   S         - surfer object describing the surface discretization
%   targinfo  - target geometry, struct with field .r (3,ntarg)
%   gkerns    - (4,1) cell array of kernel function handles for
%               smooth_sparse_quad subtraction:
%                 gkerns{1} = G_s kernel handle
%                 gkerns{2} = G_phi kernel handle
%                 gkerns{3} = S_{3d} G_phi kernel handle
%                 gkerns{4} = Laplace S_{3d} kernel handle
%   eps       - requested quadrature precision
%   zpars     - complex array of kernel parameters:
%                 zpars(1:5)  = dispersion roots
%                 zpars(6:10) = partial-fraction residues
%   ipatch_id - (optional) (ntarg,1) patch index for on-surface targets
%               (-1 for off-surface); default zeros(ntarg,1)
%   uvs_targ  - (optional) (2,ntarg) local uv coordinates for on-surface
%               targets; default zeros(2,ntarg)
%
% Output:
%   xmat1 - (ntarg, S.npts) sparse complex correction matrix for G_s
%   xmat2 - (ntarg, S.npts) sparse complex correction matrix for G_phi
%   xmat3 - (ntarg, S.npts) sparse complex correction matrix for S_{3d}G_phi
%   xmat4 - (ntarg, S.npts) sparse complex correction matrix for Laplace S_{3d}

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npatches = S.npatches;
    npp1 = npatches+1;

    [targs] = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;


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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    nquad = iquad(nnz+1)-1;

    iquadtype = 1;

    nover = S.norders(1);

    if nargin < 6
        ipatch_id = zeros(ntarg,1);
        uvs_targ = zeros(2,ntarg);
        nover = S.norders(1) + 4; 
    else
        targinfo.uvs_targ = uvs_targ;
        targinfo.patch_id = ipatch_id;
    end

    wnear = zeros(3,nquad,'like',1i);

    mex_id_ = 'getnearquad_flex_eval(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[xx])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 13, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 3, nquad);

    opts.rep = 'eval';
    opts.format = 'sparse';
    wstruc = lap3d.neumann.get_quadrature_correction(S,eps,targinfo,opts);
    xmat4 = wstruc.spmat;

    xmat1 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(1,:).');
    xmat2 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(2,:).');
    xmat3 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(3,:).');

    xmat1 = xmat1 - smooth_sparse_quad(gkerns(1), targinfo, S, row_ptr, col_ind, nover);
    xmat2 = xmat2 - smooth_sparse_quad(gkerns(2), targinfo, S, row_ptr, col_ind, nover);
    xmat3 = xmat3 - smooth_sparse_quad(gkerns(3), targinfo, S, row_ptr, col_ind, nover);
    xmat4 = xmat4 - smooth_sparse_quad(gkerns(4), targinfo, S, row_ptr, col_ind, nover);

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
%
