function nnz = findnearmem(cms, npatches, rad_near, ndtarg, targs, ntarg)
% find the number of near field points
    nnz = 0;
    n3 = 3;

    mex_id_ = 'findnearmem(i double[xx], i int[x], i double[x], i int[x], i double[xx], i int[x], io int[x])';
[nnz] = fmm3dbie_routs(mex_id_, cms, npatches, rad_near, ndtarg, targs, ntarg, nnz, n3, npatches, 1, npatches, 1, ndtarg, ntarg, 1, 1);
end

