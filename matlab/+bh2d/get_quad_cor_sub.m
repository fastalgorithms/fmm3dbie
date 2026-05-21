function [xmat,novers] = get_quad_cor_sub(S, zk, eps)

    if abs(zk)>1e-8
        error('No support for flexural problem')
    end 

    kern = @(s,t) chnk.flex2d.kern(0, s, t, 's');

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    npatches = S.npatches;
    npts = S.npts;
    ndtarg = 3;
    npp1 = npatches+1;
    ntarg = npts;
    ntp1 = ntarg+1;
    n9 = 9;
    n12 = 12;

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

    targs = srcvals(1:3,:);
    ntarg = npts;

    nnz = findnearmem(cms, npatches, rad_near, ndtarg, targs, npts);
    nnzp1 = nnz+1;

    [row_ptr,col_ind] = findnear(cms,npatches,rad_near,ndtarg,targs,npts,nnz);

    iquad = get_iquad_rsc(npatches,ixyzs,npts,ntarg,nnz,row_ptr,col_ind);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % wts = get_qwts_f(npatches,norders,ixyzs,iptype,ntarg,srcvals);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    nquad = iquad(nnz+1)-1;

    iquadtype = 1;

    wnear = zeros(nquad,1,'like',1);

    if zk <= 1e-8
    mex_id_ = 'getnearquad_bh2d_gv2v(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
    end

    xmat = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear);

    Q = []; Q.targinfo = S; Q.rfac = rfac; Q.wavenumber = 0;
    Q.row_ptr = row_ptr; Q.col_ind = col_ind; Q.kernel_order = -1;
    novers = get_oversampling_parameters(S,Q,eps);

    Asmth_over = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,novers);

    xmat = xmat - Asmth_over;

end
%
%
%
%
%
