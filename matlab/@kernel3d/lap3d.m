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
% Kernel orders (obj.kernel_order):  s -> -1,  d -> 0,  sp -> 0,  dp -> 1,  c -> 0,  cp -> 1
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
        obj.kernel_order = -1;

        obj.eval = @(s,t) lap3d.kern(s, t, 's');

        % dpars for dirichlet eval: [alpha_S, beta_D] = [1, 0]
        dpars_s = [1.0; 0.0];
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 's', sigma);
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            lap3d.dirichlet.eval(S, sigma, targ, eps, dpars_s, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          lap3d.dirichlet.get_quadrature_correction(S, eps, dpars_s, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

    case {'d', 'double'}
        obj.type = 'd';
        obj.kernel_order = 0;

        obj.eval = @(s,t) lap3d.kern(s, t, 'd');

        % dpars for dirichlet eval: [alpha_S, beta_D] = [0, 1]
        dpars_d = [0.0; 1.0];
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'd', sigma);
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            lap3d.dirichlet.eval(S, sigma, targ, eps, dpars_d, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          lap3d.dirichlet.get_quadrature_correction(S, eps,dpars_d, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.kernel_order = 0;

        obj.eval = @(s,t) lap3d.kern(s, t, 'sprime');

        % sprime is the normal-derivative of S at the target: use
        % neumann eval with dpars = [1, 0] (S only, no D)
        dpars_sp = [1.0; 0.0];
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'sp', sigma);
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            lap_combprime_eval(S, sigma, targ, eps, dpars_sp, varargin);
        obj.getquad  = @(S,eps,varargin) ...
                          lap_combprime_getquad(S, eps, dpars_sp, varargin);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.targ_fields = {'n'};

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.kernel_order = 1;

        obj.eval = @(s,t) lap3d.kern(s, t, 'dprime');

        % dpars for dirichlet eval: [alpha_S, beta_D] = [0, 1]
        dpars_dp = [0.0; 1.0];
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'dp', sigma);
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            lap_combprime_eval(S, sigma, targ, eps, dpars_dp, varargin);
        obj.getquad  = @(S,eps,varargin) ...
                          lap_combprime_getquad(S, eps, dpars_dp, varargin);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};
        obj.targ_fields = {'n'};

    case {'c', 'combined'}
        if ( nargin < 2 )
            warning('KERNEL3D.LAP3D: missing coefs for combined layer. Defaulting to [1 1].');
            coefs = [1; 1];
        end
        coefs = coefs(:);
        assert(isreal(coefs), 'KERNEL3D.LAP3D: coefs must be real for combined field.');

        obj.type         = 'c';
        obj.kernel_order = 0;
        obj.params.coefs = coefs;

        obj.eval = @(s,t) lap3d.kern(s, t, 'c', coefs(1), coefs(2));

        % dpars for dirichlet eval: [alpha_S, beta_D] = coefs
        dpars_c = coefs;
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'c', sigma, coefs);
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            lap3d.dirichlet.eval(S, sigma, targ, eps, dpars_c, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          lap3d.dirichlet.get_quadrature_correction(S, eps, dpars_c, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};

    case {'cp', 'cprime'}
        if ( nargin < 2 )
            warning('KERNEL3D.LAP3D: missing coefs for cprime. Defaulting to [1 1].');
            coefs = [1; 1];
        end
        coefs = coefs(:);
        assert(isreal(coefs), 'KERNEL3D.LAP3D: coefs must be real for combined field.');
        

        obj.type         = 'cprime';
        obj.kernel_order = 1;
        obj.params.coefs = coefs;

        obj.eval = @(s,t) lap3d.kern(s, t, 'cprime', coefs);
        obj.fmm  = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'cp', sigma, coefs);

        dpars_c = coefs;
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            lap_combprime_eval(S, sigma, targ, eps, dpars_c, varargin);
        obj.getquad  = @(S,eps,varargin) ...
                          lap_combprime_getquad(S, eps, dpars_c, varargin);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};
        obj.targ_fields = {'n'};
        
    otherwise
        error('KERNEL3D.LAP3D: unknown Laplace kernel type ''%s''.', type);

end

obj.rsc_to_interleave = kernel3d.rsc_interleave_scalar();

end

function p = lap_combprime_eval(S, sigma, targ, eps, rep_pars, args)
%LAP_COMBPRIME_EVAL  Call lap3d.dirichlet.eval with iprime=1.
if nargin < 6 || isempty(args), args = {}; end
% Ensure opts struct is present and has iprime=1
if isempty(args) || ~isstruct(args{end})
    args{end+1} = struct();
end
args{end}.iprime = 1;
p = lap3d.dirichlet.eval(S, sigma, targ, eps, rep_pars, args{:});
end

function Q = lap_combprime_getquad(S, eps, rep_pars, args)
%LAP_COMBPRIME_GETQUAD  Call get_quadrature_correction with iprime=1.
if nargin < 4 || isempty(args), args = {}; end
% Separate targinfo and opts from args
if length(args) >= 2
    targinfo = args{1};
    opts = args{2};
elseif length(args) == 1
    targinfo = args{1};
    opts = struct();
else
    % No targets provided: use S as on-surface self-quadrature
    targinfo = S;
    opts = struct();
end
% If targinfo is empty use S
if isempty(targinfo)
    targinfo = S;
end
% Set iprime flag
opts.iprime = 1;
Q = lap3d.dirichlet.get_quadrature_correction(S, eps, rep_pars, targinfo, opts);
end