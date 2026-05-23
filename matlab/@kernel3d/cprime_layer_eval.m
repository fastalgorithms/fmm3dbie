function p = cprime_layer_eval(S, sigma, targ, eps, zk, rep_pars, varargin)
%CPRIME_LAYER_EVAL  Helper: call helm3d.dirichlet.eval with iprime=1.
%
%   Accepts optional trailing opts struct and ensures opts.iprime = 1.

if isempty(varargin)
    opts = struct('iprime', 1);
else
    opts = varargin{1};
    opts.iprime = 1;
end

p = helm3d.dirichlet.eval(S, sigma, targ, eps, zk, rep_pars, opts);

end
