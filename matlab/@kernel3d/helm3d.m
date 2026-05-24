function obj = helm3d(type,zk, coefs)
%KERNEL3D.HELM3D   Construct a Helmholtz kernel in 3D.
%
%   KERNEL3D.HELM3D(type,zk) or KERNEL3D.HELM3D(zk, type, coefs), where
%   zk is the (complex) wavenumber and type is one of:
%      's'     - single layer,   S(x,y) = exp(ik|x-y|)/(4*pi*|x-y|)
%      'd'     - double layer,   D(x,y) = d/dn_y G(x,y)
%      'sp'    - sprime,         S'(x,y) = d/dn_x G(x,y)
%      'dprime'- dprime,         D'(x,y) = d^2/dn_x dn_y G(x,y)
%      'c'     - combined layer, coefs(1)*S + coefs(2)*D  (default coefs=[1 1])
%      'cprime'- combined prime, coefs(1)*S' + coefs(2)*D'
%      'trans' - transmission,   two-density transmission representation
%                                (requires rep_params = [alpha0;beta0;alpha1;beta1])
%
% Kernel evaluation:
%   obj.eval(srcinfo, targinfo)
%      srcinfo / targinfo are structs with fields .r (3,:) and .n (3,:).
%      Returns an (nt x ns) matrix, or (2*nt x 2*ns) for transmission.
%
% Raw FMM call (s, d, sp, c only):
%   obj.fmm(eps, srcinfo, targinfo, sigma)
%      Calls helm3d.fmm directly (no quadrature corrections).
%
% FMM + quadrature layer-potential eval:
%   obj.layer_eval(S, sigma, targinfo, eps)
%   obj.layer_eval(S, sigma, targinfo, eps, opts)
%
% Quadrature corrections:
%   obj.getquad(S, eps, varargin)
%
% Kernel orders (obj.kernel_order):  s -> -1,  d -> 0,  sp -> 0,  dp -> 1,  c -> 0,  cp -> 1,  trans -> 1
%
% See also HELM3D.KERN, HELM3D.DIRICHLET.EVAL, HELM3D.NEUMANN.EVAL,
%          HELM3D.TRANSMISSION.EVAL

if ( nargin < 2 )
    error('KERNEL3D.HELM3D: requires at least zk and type arguments.');
end

obj           = kernel3d();
obj.name      = 'helmholtz';
obj.opdims    = [1 1];
obj.zk        = zk;
obj.ifcomplex = 1;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.kernel_order = -1;

        obj.eval = @(s,t) helm3d.kern(zk, s, t, 's');

        obj.fmm = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 's', sigma);

        % rep_pars for dirichlet eval: [alpha_S, beta_D] = [1, 0]
        rep_pars_s = [1.0; 0.0];
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            helm3d.dirichlet.eval(S, sigma, targ, eps, zk, rep_pars_s, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars_s, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

    case {'d', 'double'}
        obj.type = 'd';
        obj.kernel_order = 0;

        obj.eval = @(s,t) helm3d.kern(zk, s, t, 'd');

        obj.fmm = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'd', sigma);

        % rep_pars for dirichlet eval: [alpha_S, beta_D] = [0, 1]
        rep_pars_d = [0.0; 1.0];
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            helm3d.dirichlet.eval(S, sigma, targ, eps, zk, rep_pars_d, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars_d, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.kernel_order = 0;

        obj.eval = @(s,t) helm3d.kern(zk, s, t, 'sprime');

        obj.fmm = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'sp', sigma);

        % sprime = cprime with coefs = [1, 0]
        rep_pars_sp = [1.0; 0.0];
        obj.layer_eval = @(S,sigma,targ,eps) ...
                            helm_combprime_eval(S, sigma, targ, eps, zk, rep_pars_sp);
        obj.getquad  = @(S,eps,varargin) ...
                          helm_combprime_getquad(S, eps, zk, rep_pars_sp, varargin);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.targ_fields = {'n'};

    case {'dp', 'dprime'}
        obj.type = 'dprime';
        obj.kernel_order = 1;

        obj.eval = @(s,t) helm3d.kern(zk, s, t, 'dprime');

        % dprime = cprime with coefs = [0, 1], iprime = 1
        rep_pars_dp = [0.0; 1.0];
        obj.layer_eval = @(S,sigma,targ,eps) ...
                            helm_combprime_eval(S, sigma, targ, eps, zk, rep_pars_dp);
        obj.getquad  = @(S,eps,varargin) ...
                          helm_combprime_getquad(S, eps, zk, rep_pars_dp, varargin);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};
        obj.targ_fields = {'n'};

    case {'c', 'combined'}
        if ( nargin < 3 )
            warning('KERNEL3D.HELM3D: missing coefs for combined layer. Defaulting to [1 1].');
            coefs = [1; 1];
        end
        coefs = coefs(:);

        obj.type         = 'c';
        obj.kernel_order = 0;
        obj.params.coefs = coefs;

        obj.eval = @(s,t) helm3d.kern(zk, s, t, 'c', coefs);

        obj.fmm = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'c', sigma, coefs);

        % rep_pars for dirichlet eval: [alpha_S, beta_D] = coefs
        rep_pars_c = coefs;
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            helm3d.dirichlet.eval(S, sigma, targ, eps, zk, rep_pars_c, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars_c, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};

    case {'cp', 'cprime'}
        if ( nargin < 3 )
            warning('KERNEL3D.HELM3D: missing coefs for cprime. Defaulting to [1 1].');
            coefs = [1; 1];
        end
        coefs = coefs(:);

        obj.type         = 'cprime';
        obj.kernel_order = 1;
        obj.params.coefs = coefs;

        obj.eval = @(s,t) helm3d.kern(zk, s, t, 'cprime', coefs);

        rep_pars_cp = coefs;
        obj.layer_eval = @(S,sigma,targ,eps) ...
                            helm_combprime_eval(S, sigma, targ, eps, zk, rep_pars_cp);
        obj.getquad  = @(S,eps,varargin) ...
                          helm_combprime_getquad(S, eps, zk, rep_pars_cp, varargin);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};
        obj.targ_fields = {'n'};

    case {'trans', 'transmission'}
        if ( nargin < 3 )
            error('KERNEL3D.HELM3D: transmission kernel requires rep_params = [alpha0;beta0;alpha1;beta1].');
        end
        rep_params = coefs(:);   % reuse coefs slot for rep_params
        if ( length(rep_params) ~= 4 )
            error('KERNEL3D.HELM3D: rep_params must have 4 entries [alpha0;beta0;alpha1;beta1].');
        end

        obj.type           = 'trans';
        obj.kernel_order   = 1;
        obj.opdims         = [2 2];
        obj.params.rep_params = rep_params;

        % eval: 2ns columns interleaved [rho, sigma], 2nt rows interleaved
        % Transmission kernel evaluated via helm3d.kern 'eval' type
        % coef = i*zk (the S weight in the transmission representation)
        coef_trans = 1i * zk;
        obj.eval = @(s,t) helm3d.kern(zk, s, t, 'eval', coef_trans);

        % layer_eval and getquad delegate to helm3d.transmission
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            helm3d.transmission.eval(S, sigma, targ, eps, [zk; zk], rep_params, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          helm3d.transmission.get_quadrature_correction(S, eps, [zk; zk], rep_params, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};
        obj.targ_fields = {'n'};

    otherwise
        error('KERNEL3D.HELM3D: unknown Helmholtz kernel type ''%s''.', type);

end

end

function p = helm_combprime_eval(S, sigma, targ, eps, zk, rep_pars,args)
%HELM_COMBPRIME_EVAL  Call helm3d.dirichlet.eval with iprime=1.
if isempty(args) || ~isstruct(args{end})
    args{end+1} = [];
end
args{end}.iprime = 1;
p = helm3d.dirichlet.eval(S, sigma, targ, eps, zk, rep_pars, args{:});
end

function Q = helm_combprime_getquad(S, eps, zk, rep_pars, args)
%HELM_COMBPRIME_GETQUAD  Call get_quadrature_correction with iprime=1.
%   args is the varargin cell from the lambda: {} or {targinfo} or {targinfo, opts}
if isempty(args) || ~isstruct(args{end})
    args{end+1} = [];
end
args{end}.iprime = 1;
Q = helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars, args{:});
end
