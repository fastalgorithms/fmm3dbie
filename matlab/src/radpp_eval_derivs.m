function comps = radpp_eval_derivs(dx,dy,dz,ipars,dpars)
%RADPP_EVAL_DERIVS radpp fit and derivatives at displacements (dx,dy,dz)
%
%  comps = radpp_eval_derivs(dx,dy,dz,ipars,dpars), with ipars (including
%  ipars(6:9) = nderiv, ncomp, dim, isrc) and dpars from KERNEL3D.RADPP.
%  dx, dy, dz: target - source, same size. comps: (numel(dx),nc), ordered
%  as RADPP_ROTATE. Outside the fit, the end panels are extrapolated.
%
%  See also RADPP_ROTATE, RADPP_EVAL
%
    dx = dx(:);  dy = dy(:);  dz = dz(:);
    n = length(dx);
    ndd = length(dpars);
    ndi = length(ipars);
    nfit = ipars(6) + ipars(9);
    switch ipars(8)
      case 2,  nc = [1 3 6 10];   nc = nc(nfit+1);
      case 3,  nc = [1 4 10 20];  nc = nc(nfit+1);
      case 13, nc = 6 + 18*(ipars(6) > 0);
    end
    comps = complex(zeros(nc,n));
    mex_id_ = 'radpp_eval_derivs_vec(c i int64_t[x], c i double[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io dcomplex[xx])';
[comps] = fmm3dbie_routs(mex_id_, n, dx, dy, dz, ndd, dpars, ndi, ipars, nc, comps, 1, n, n, n, 1, ndd, 1, ndi, 1, nc, n);
    comps = comps.';
end

%
%

