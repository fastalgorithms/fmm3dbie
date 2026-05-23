function Q = cprime_getquad(S, eps, zk, rep_pars, varargin)
%CPRIME_GETQUAD  Helper: call get_quadrature_correction with iprime=1.
%
%   Accepts the same varargin as the standard getquad lambda,
%   i.e. varargin = {targinfo} or {targinfo, opts}.
%   Ensures opts.iprime = 1 before forwarding.

if isempty(varargin)
    targinfo = [];
    opts = struct('iprime', 1);
elseif length(varargin) == 1
    targinfo = varargin{1};
    opts = struct('iprime', 1);
else
    targinfo = varargin{1};
    opts = varargin{2};
    opts.iprime = 1;
end

Q = helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars, targinfo, opts);

end
