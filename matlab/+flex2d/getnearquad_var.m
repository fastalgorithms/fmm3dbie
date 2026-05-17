function A = getnearquad_var(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targinfo,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,dpars,zpars,nquad,type)
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;

    [targs] = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;
    nnzp1 = nnz+1;
    ndz = length(zpars);

    A = zeros(1,nquad,'like',1i);
    wmat = zeros(1,nquad,'like',1i);

    if strcmp(type, 'v2v') 
        for kernel_id = 1:7
        mex_id_ = 'getnearquad_flex2d_var(i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i double[xx], i int64_t[x], i int64_t[x], i double[xx], i int64_t[x], i double[xx], i double[x], i dcomplex[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[x], i int64_t[x], i int64_t[x], io dcomplex[x])';
[wmat] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, kernel_id, wmat, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 2, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 1, nquad);

            for i = 1:length(row_ptr)-1
                for k = row_ptr(i):row_ptr(i+1)-1
                    A(:,iquad(k):iquad(k+1)-1) = A(:,iquad(k):iquad(k+1)-1) + ...
                        targinfo.icoefs(kernel_id,i).*wmat(iquad(k):iquad(k+1)-1).';
                end
            end
        end
    end
end
%
%
