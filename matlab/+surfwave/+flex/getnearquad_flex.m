function wnear = getnearquad_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker)
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;
    %if isfield(targs,'r') || isprop(targs,'r')
    if ~isfloat(targs)
        targs = targs.r;
    end
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;
    ntp1 = ntarg+1;
    nnzp1 = nnz+1;
    ndz = length(zpars);

    wnear = zeros(1,nquad,'like',1i);

    if iker < 7

    if iker < 6
        mex_id_ = 'getnearquad_flex_all_targ(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io dcomplex[x])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, nquad);
    else
        wnear_lap = zeros(1,nquad);    
        mex_id_ = 'getnearquad_lap_s_neu_eval(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[wnear_lap] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear_lap, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        wnear = wnear_lap;
    end


    else
    wnear1 = surfwave.flex.getnearquad_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,5);
    wnear2 = surfwave.flex.getnearquad_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,6);    
    wnear = wnear1 + wnear2;
    end
    wnear = wnear(:);
end

