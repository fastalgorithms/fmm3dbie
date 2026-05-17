function nnz = findnearmem(cms, npatches, rad_near, ndtarg, targs, ntarg)
    nnz = 0;
    n3 = 3;

    mex_id_ = 'findnearmem(c i double[xx], c i int64_t[x], c i double[x], c i int64_t[x], c i double[xx], c i int64_t[x], c io int64_t[x])';
[nnz] = kern_routs(mex_id_, cms, npatches, rad_near, ndtarg, targs, ntarg, nnz, n3, npatches, 1, npatches, 1, ndtarg, ntarg, 1, 1);
end

