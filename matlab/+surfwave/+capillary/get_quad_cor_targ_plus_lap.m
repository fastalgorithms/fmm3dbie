function [xmat1,xmat2,xmat3] = get_quad_cor_targ_plus_lap(S, targinfo, kern, eps, zpars,ipatch_id,uvs_targ)
% returns quadrature corrections for G_S, G_phi, AND lap G_phi evaluation / postproc

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
    % prin2('rfac = *',rfac,1)
    % prin2('rfac0 = *',rfac0,1)

        % allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

    [cms, rads] = get_centroid_rads(npatches,norders,ixyzs,iptype,npts, ... 
       srccoefs);

    rad_near = rads*rfac;

    % 
    % find near quadrature correction interactions
    % 

    % prinf('entering find near mem',0,0)
    nnz = findnearmem(cms, npatches, rad_near, ndtarg, targs, ntarg);
    nnzp1 = nnz+1;
    % prinf('nnz = *',nnz,1);

    [row_ptr,col_ind] = findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,nnz);

    iquad = get_iquad_rsc(npatches,ixyzs,npts,ntarg,nnz,row_ptr,col_ind);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wts = get_qwts_f(npatches,norders,ixyzs,iptype,ntarg,srcvals);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    nquad = iquad(nnz+1)-1;

    iquadtype = 1;

    if nargin < 6
        ipatch_id = zeros(ntarg,1);
        uvs_targ = zeros(2,ntarg);
    end
    wnear = zeros(2,nquad);

    mex_id_ = 'getnearquad_capillary_eval_all(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i int[x], i double[xx], i double[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[x], i int[x], io dcomplex[xx])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 6, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 2, nquad);
    
    xmat1 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(1,:).');
    xmat2 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(2,:).');

    wnear = zeros(1,nquad);

    mex_id_ = 'getnearquad_capillary_lapgphi_eval(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i int[x], i double[xx], i double[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[x], i int[x], io dcomplex[xx])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 6, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 1, nquad);
    
    xmat3 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(1,:).');

    nover = S.norders(1) + 4;
    targ_struct = struct('r', targs(1:2,:));

    xmat1 = xmat1 - smooth_sparse_quad(kern(1), targ_struct, S, row_ptr, col_ind, nover);
    xmat2 = xmat2 - smooth_sparse_quad(kern(2), targ_struct, S, row_ptr, col_ind, nover);
    xmat3 = xmat3 - smooth_sparse_quad(kern(3), targ_struct, S, row_ptr, col_ind, nover);

end