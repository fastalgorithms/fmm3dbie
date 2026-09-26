function wnear = getnearquad_capillary(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id,uvs_targ,eps, ...
          iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker)
%
%  surfwave.capillary.getnearquad_capillary
%    Near-field quadrature correction for the capillary surface wave
%    kernels, at an arbitrary target array.
%
%    iker selects the kernel:
%      0 = G_S            1 = G_phi        3 = lap G_phi
%      5 = S3d G_phi      6 = Laplace S3d  7 = 5 + 6
%      8 = S'_S           9 = S'_phi
%
%    Targets on the source surface carry their patch id in ipatch_id and
%    local coordinates in uvs_targ; off-surface targets pass -1.
%
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;
    if ~isnumeric(targs)   % may already be a flattened array (recursive calls)
        targs = extract_targ_array(targs);
    end
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;
    nnzp1 = nnz+1;

    wnear = zeros(1,nquad,'like',1i);

    if iker < 6

        mex_id_ = 'getnearquad_capillary_all(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io dcomplex[x])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 6, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 1, nquad);

    elseif iker == 6

        wnear_lap = zeros(1,nquad);
        mex_id_ = 'getnearquad_lap_s_neu_eval(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[wnear_lap] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear_lap, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        wnear = wnear_lap;

    elseif iker == 7

        wnear1 = surfwave.capillary.getnearquad_capillary(npatches,norders, ...
            ixyzs,iptype,npts,srccoefs,srcvals,targs,ipatch_id,uvs_targ, ...
            eps,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,zpars,nquad,5);
        wnear2 = surfwave.capillary.getnearquad_capillary(npatches,norders, ...
            ixyzs,iptype,npts,srccoefs,srcvals,targs,ipatch_id,uvs_targ, ...
            eps,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,zpars,nquad,6);
        wnear = wnear1 + wnear2;

    else
        % iker 8, 9: S' kernels.  S' depends on the target normal, so it
        % is not radial and is handled by the direct routine only.
        mex_id_ = 'getnearquad_capillary_all(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io dcomplex[x])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 6, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 1, nquad);

    end

end

