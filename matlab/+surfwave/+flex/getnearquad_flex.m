function wnear = getnearquad_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker,ivpp)
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
    if ivpp
        max_x = max([max(srcvals(1,:)),max(targs(1,:))]);
        min_x = min([min(srcvals(1,:)),min(targs(1,:))]);

        max_y = max([max(srcvals(2,:)),max(targs(2,:))]);
        min_y = min([min(srcvals(2,:)),min(targs(2,:))]);

        maxdist = 1.5*sqrt((max_x - min_x).^2 + (max_y - min_y).^2);
        mex_id_ = 'getnearquad_flex_all_vpp(i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i double[xx], i double[x], i dcomplex[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[x], i int64_t[x], i int64_t[x], i double[x], io dcomplex[x])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, maxdist, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, 1, nquad);
    else
        mex_id_ = 'getnearquad_flex_all(i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i double[xx], i double[x], i dcomplex[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[x], i int64_t[x], i int64_t[x], io dcomplex[x])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, nquad);
    end
    else
        wnear_lap = zeros(1,nquad);    
        mex_id_ = 'getnearquad_lap_s_neu_eval(i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i double[xx], i int64_t[x], i int64_t[x], i double[xx], i int64_t[x], i double[xx], i double[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[x], i int64_t[x], io double[x])';
[wnear_lap] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear_lap, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        wnear = wnear_lap;
    end


    else
    wnear1 = surfwave.flex.getnearquad_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,5,ivpp);
    wnear2 = surfwave.flex.getnearquad_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,6,ivpp);    
    wnear = wnear1 + wnear2;
    end
    wnear = wnear(:);
end

