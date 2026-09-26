function wnear = getnearquad_gravity_grad(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker,S3d_scal)
%
%  surfwave.gravity.getnearquad_gravity_grad
%
%  Near field quadrature for the TARGET Cartesian gradient (d/dx, d/dy)
%  of gsgrav (iker=0) or gphigrav (iker=1), with the same calling
%  convention as getnearquad_gravity.m, returning a 2-row wnear
%  (x-component, y-component). getnearquad_gravity_grad_all returns the
%  bulk part only; the gradient of the algebraic 1/r term,
%  S3d_scal*d/dx[1/(4*pi*r)], is added here via
%  getnearquad_lap_grad_s_neu_eval.
%
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;
    if ~isnumeric(targs)   % may already be a flattened array (recursive calls)
        targs = targs.r;
    end
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;
    nnzp1 = nnz+1;
    ndz = length(zpars);

    wnear = zeros(2,nquad,'like',1i);

    mex_id_ = 'getnearquad_gravity_grad_all(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io dcomplex[xx])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, ndz, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 1, 2, nquad);

    wnear_lap_grad = zeros(2,nquad);

    mex_id_ = 'getnearquad_lap_grad_s_neu_eval(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[xx])';
[wnear_lap_grad] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear_lap_grad, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 2, nquad);

    wnear = wnear + S3d_scal*wnear_lap_grad;
end
%
%


%
%%  Flexural 2D near-field quadrature
%%
%

