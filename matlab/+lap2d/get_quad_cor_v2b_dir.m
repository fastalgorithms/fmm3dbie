function [xmat, novers] = get_quad_cor_v2b_dir(S, targinfo, eps, uv_bndry, patch_id)
%
%  lap2d.get_quad_cor_v2b_dir
%
%  Near-quadrature correction for evaluating the Laplace SLP volume
%  potential at boundary targets:
%
%    val(r) = \int G(r,r') mu(r') dA(r'),   r on boundary
%
%  Inputs:
%    S        - surfer object (source domain)
%    targinfo - target struct with field r (2 or 3 x ntarg), and
%               optionally patch_id, uvs_targ for on-surface targets
%    eps      - precision
%    uv_bndry - (optional) (2,ntarg) local uv coords of on-surface targets
%    patch_id - (optional) (ntarg,1) patch ids of on-surface targets
%
%  Outputs:
%    xmat     - sparse correction matrix (ntarg x S.npts)
%    novers   - oversampling orders per patch

    lap2d_kern = kernel('l', 's');

    [srcvals, srccoefs, norders, ixyzs, iptype, wts] = extract_arrays(S);
    [n12, npts] = size(srcvals);
    [n9, ~]     = size(srccoefs);
    npatches    = S.npatches;
    npp1        = npatches + 1;

    [targs]       = extract_targ_array(targinfo);
    [ndtarg, ntarg] = size(targs);
    ntargp1       = ntarg + 1;

    iptype_avg = floor(sum(iptype) / (npatches + 0.0d0));
    norder_avg = floor(sum(norders) / (npatches + 0.0d0));

    [rfac, rfac0] = get_rfacs(norder_avg, iptype_avg);

    [cms, rads] = get_centroid_rads(npatches, norders, ixyzs, iptype, npts, srccoefs);
    rad_near = rads * rfac;

    nnz   = findnearmem(cms, npatches, rad_near, ndtarg, targs, ntarg);
    nnzp1 = nnz + 1;

    [row_ptr, col_ind] = findnear(cms, npatches, rad_near, ndtarg, targs, ntarg, nnz);

    iquad = get_iquad_rsc(npatches, ixyzs, npts, ntarg, nnz, row_ptr, col_ind);

    if nargin < 4
        ipatch_id = zeros(ntarg, 1);
        uvs_targ  = zeros(2, ntarg);
    else
        if ~isempty(patch_id)
            ipatch_id = patch_id;
            uvs_targ  = uv_bndry;
        else
            ipatch_id = zeros(ntarg, 1);
            uvs_targ  = zeros(2, ntarg);
        end
    end

    nquad     = iquad(nnz+1) - 1;
    iquadtype = 1;

    A    = zeros(1, nquad, 'like', 1);
    zpars = 0;

    mex_id_ = 'getnearquad_lap2d_dir(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
    [A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);

    xmat = conv_rsc_to_spmat(S, row_ptr, col_ind, real(A));

    Q = []; Q.targinfo = targinfo; Q.rfac = rfac; Q.wavenumber = 0;
    Q.row_ptr = row_ptr; Q.col_ind = col_ind; Q.kernel_order = -1;
    novers = get_oversampling_parameters(S, Q, eps);

    Asmth_over = smooth_sparse_quad(lap2d_kern, targs, S, row_ptr, col_ind, novers);

    xmat = xmat - Asmth_over;
end
