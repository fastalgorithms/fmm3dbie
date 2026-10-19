function obj = lap3d(type, coefs)
%KERNEL3D.LAP3D   Construct a Laplace kernel in 3D.
%
%   KERNEL3D.LAP3D(type) or KERNEL3D.LAP3D(type, coefs), where type is:
%      's'  - single layer,   S(x,y) = 1/(4*pi*|x-y|)
%      'd'  - double layer,   D(x,y) = d/dn_y G(x,y)
%      'sp' - sprime,         S'(x,y) = d/dn_x G(x,y)
%      'c'  - combined layer, coefs(1)*S + coefs(2)*D  (default coefs=[1 1])
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
obj.getquad_sub = true;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.kernel_order = -1;

        obj.eval = @(s,t) lap3d.kern(s, t, 's');

        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 's', sigma);
        obj.getquad  = @(S,eps,varargin) surf_quad(S, eps, 'lap3d_s', [], varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

    case {'d', 'double'}
        obj.type = 'd';
        obj.kernel_order = 0;

        obj.eval = @(s,t) lap3d.kern(s, t, 'd');

        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'd', sigma);
        obj.getquad  = @(S,eps,varargin) surf_quad(S, eps, 'lap3d_d', [], varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.kernel_order = 0;

        obj.eval = @(s,t) lap3d.kern(s, t, 'sprime');

        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'sp', sigma);
        obj.getquad  = @(S,eps,varargin) surf_quad(S, eps, 'lap3d_sp', [], varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.targ_fields = {'n'};

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.kernel_order = 1;

        obj.eval = @(s,t) lap3d.kern(s, t, 'dprime');

        % dprime = 0*S' + 1*D'
        pars_dp      = struct('dpars', [0.0; 1.0]);
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'dp', sigma);
        obj.getquad  = @(S,eps,varargin) surf_quad(S, eps, 'lap3d_cp', pars_dp, varargin{:});
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

        % dpars = [alpha_S; beta_D] = coefs
        pars_c       = struct('dpars', coefs);
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'c', sigma, coefs);
        obj.getquad  = @(S,eps,varargin) surf_quad(S, eps, 'lap3d_c', pars_c, varargin{:});
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

        % dpars = [alpha_S'; beta_D'] = coefs; without the hypersingular
        % part a principal value self-quadrature is enough
        pars_cp      = struct('dpars', coefs);
        name_cp      = 'lap3d_cp';
        if coefs(2) == 0, name_cp = 'lap3d_cp_pv'; end
        obj.getquad  = @(S,eps,varargin) surf_quad(S, eps, name_cp, pars_cp, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};
        obj.targ_fields = {'n'};
        
    otherwise
        error('KERNEL3D.LAP3D: unknown Laplace kernel type ''%s''.', type);

end

end

function spmat = surf_quad(S, eps, name, pars, targinfo, opts)
%SURF_QUAD  getquad via GETNEARQUAD_MEX.
if nargin < 5, targinfo = []; end
if nargin < 6 || isempty(opts), opts = struct(); end
opts.format = 'sparse';
Q = getnearquad_mex(S, name, eps, pars, targinfo, opts);
spmat = Q.spmat;
end
