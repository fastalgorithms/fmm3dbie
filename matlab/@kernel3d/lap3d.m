function obj = lap3d(type, coefs)
%KERNEL3D.LAP3D   Construct a Laplace kernel in 3D.
%
%   KERNEL3D.LAP3D(type) or KERNEL3D.LAP3D(type, coefs), where type is:
%      's'  - single layer,   S(x,y) = 1/(4*pi*|x-y|)
%      'd'  - double layer,   D(x,y) = d/dn_y G(x,y)
%      'sp' - sprime,         S'(x,y) = d/dn_x G(x,y)
%      'c'  - combined layer, coefs(1)*S + coefs(2)*D  (default coefs=[1 1])
%
% Kernel evaluation:
%   obj.eval(srcinfo, targinfo)
%      srcinfo / targinfo are structs with fields .r (3,:) and .n (3,:).
%      Returns an (nt x ns) matrix.
%
% Raw FMM call:
%   obj.fmm(eps, srcinfo, targinfo, sigma)
%      Calls lap3d.fmm directly (no quadrature corrections).
%
% FMM + quadrature layer-potential eval:
%   obj.layer_eval(S, sigma, targinfo, eps)
%   obj.layer_eval(S, sigma, targinfo, eps, opts)
%      Calls lap3d.dirichlet.eval or lap3d.neumann.eval (FMM + near
%      corrections) with dpars baked in from the kernel object.
%
% Quadrature corrections:
%   obj.getquad(S, eps, varargin)
%      Calls lap3d.dirichlet.get_quadrature_correction or
%      lap3d.neumann.get_quadrature_correction.
%
% Kernel orders (obj.sing):  s -> -1,  d -> 0,  sp -> 0,  c -> -1
%
% Oversampling orders:
%   obj.get_overs_orders(S, t, eps)
%      Returns novers, the per-patch oversampling orders needed for
%      accurate near-field quadrature.  Wraps get_oversampling_parameters,
%      using getnear(S,t) to build the required rsc struct.
%
% See also LAP3D.KERN, LAP3D.DIRICHLET.EVAL, LAP3D.NEUMANN.EVAL

if ( nargin < 1 )
    error('KERNEL3D.LAP3D: missing Laplace kernel type.');
end

obj           = kernel3d();
obj.name      = 'laplace';
obj.opdims    = [1 1];
obj.zk        = 0;
obj.ifcomplex = 0;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.sing = -1;

        obj.eval = @(s,t) lap3d.kern(s, t, 's');

        % dpars for dirichlet eval: [alpha_S, beta_D] = [1, 0]
        dpars_s = [1.0; 0.0];
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 's', sigma);
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            lap3d.dirichlet.eval(S, sigma, targ, eps, dpars_s, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          lap3d.dirichlet.get_quadrature_correction(S, eps, dpars_s, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.sing);

    case {'d', 'double'}
        obj.type = 'd';
        obj.sing = -1;

        obj.eval = @(s,t) lap3d.kern(s, t, 'd');

        % dpars for dirichlet eval: [alpha_S, beta_D] = [0, 1]
        dpars_d = [0.0; 1.0];
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'd', sigma);
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            lap3d.dirichlet.eval(S, sigma, targ, eps, dpars_d, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          lap3d.dirichlet.get_quadrature_correction(S, eps,dpars_d, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.sing);

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.sing = -1;

        obj.eval = @(s,t) lap3d.kern(s, t, 'sprime');

        % sprime is the normal-derivative of S at the target: use
        % neumann eval with dpars = [1, 0] (S only, no D)
        dpars_sp = [1.0; 0.0];
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'sp', sigma);
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            lap3d.neumann.eval(S, sigma, targ, eps, dpars_sp, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          lap3d.neumann.get_quadrature_correction(S, eps, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.sing);

    case {'c', 'combined'}
        if ( nargin < 2 )
            warning('KERNEL3D.LAP3D: missing coefs for combined layer. Defaulting to [1 1].');
            coefs = [1; 1];
        end
        coefs = coefs(:);

        obj.type        = 'c';
        obj.sing         = -1;
        obj.params.coefs = coefs;

        obj.eval = @(s,t) lap3d.kern(s, t, 'c', coefs(1), coefs(2));

        % dpars for dirichlet eval: [alpha_S, beta_D] = coefs
        dpars_c = coefs;
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'c', sigma, coefs);
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            lap3d.dirichlet.eval(S, sigma, targ, eps, dpars_c, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          lap3d.dirichlet.get_quadrature_correction(S, eps, dpars_c, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.sing);

    otherwise
        error('KERNEL3D.LAP3D: unknown Laplace kernel type ''%s''.', type);

end

end
