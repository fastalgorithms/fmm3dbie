function val = radpp_eval_fortran(r,ipars,dpars)
%RADPP_EVAL_FORTRAN evaluate a radpp_fit expansion with the Fortran
%  evaluator vpp_eval_vec
%
%  val = radpp_eval_fortran(r,ipars,dpars), with ipars and dpars from
%  RADPP_FIT. r is any array; val is (numel(r),m), as RADPP_EVAL. Radii
%  outside the fit give NaN.
%
%  See also RADPP_EVAL, RADPP_FIT
%
    r = r(:);
    nr = length(r);
    ndd = length(dpars);
    ndi = length(ipars);
    nv = ipars(5);
    vals = zeros(nv,nr);
    mex_id_ = 'vpp_eval_vec(c i int64_t[x], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io double[xx])';
[vals] = fmm3dbie_routs(mex_id_, nr, r, ndd, dpars, ndi, ipars, vals, 1, nr, 1, ndd, 1, ndi, nv, nr);
    val = complex(vals(1:2:end,:),vals(2:2:end,:)).';
    val(r < dpars(1) | r > dpars(ipars(1)+1),:) = NaN;
end

%
%

