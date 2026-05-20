function [xmat,novers] = get_quad_cor_v2b_dir(S, zk, targinfo, eps, uv_bndry, patch_id)

    h2d_s = kernel('h','s',zk);

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
    % wts = get_qwts_f(npatches,norders,ixyzs,iptype,ntarg,srcvals);

    if nargin < 5
        ipatch_id = zeros(ntarg,1);
        uvs_targ = zeros(2,ntarg);
    else
        if ~isempty(patch_id)
            ipatch_id = patch_id;
            uvs_targ = uv_bndry;
        else
            ipatch_id = zeros(ntarg,1);
            uvs_targ = zeros(2,ntarg);
        end
    end

    targinfo = [];
    targinfo.r = targs(1:3,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    nquad = iquad(nnz+1)-1;

    iquadtype = 1;

    A = zeros(1,nquad,'like',1);
    zpars = complex(double(real(zk)), double(imag(zk)));

    mex_id_ = 'getnearquad_helm2d_dir(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);

    xmat = conv_rsc_to_spmat(S,row_ptr,col_ind,A);

    Q = []; Q.targinfo = targinfo; Q.rfac = rfac; Q.wavenumber = zk;
    Q.row_ptr = row_ptr; Q.col_ind = col_ind; Q.kernel_order = -1;
    novers = get_oversampling_parameters(S,Q,1e2*eps);

    Asmth_over = smooth_sparse_quad(h2d_s,targs,S,row_ptr,col_ind,novers);

    xmat = xmat - Asmth_over;
end
%
%
%
%
%
