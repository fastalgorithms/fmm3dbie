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

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.kernel_order = -1;

        obj.eval = @(s,t) lap3d.kern(s, t, 's');

        % dpars for dirichlet eval: [alpha_S, beta_D] = [1, 0]
        dpars_s = [1.0; 0.0];
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 's', sigma);
        obj.getquad  = @(S,eps,varargin) rsc_to_sparse( ...
                          lap3d.dirichlet.get_quadrature_correction(S, eps, dpars_s, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

    case {'d', 'double'}
        obj.type = 'd';
        obj.kernel_order = 0;

        obj.eval = @(s,t) lap3d.kern(s, t, 'd');

        % dpars for dirichlet eval: [alpha_S, beta_D] = [0, 1]
        dpars_d = [0.0; 1.0];
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'd', sigma);
        obj.getquad  = @(S,eps,varargin) rsc_to_sparse( ...
                          lap3d.dirichlet.get_quadrature_correction(S, eps,dpars_d, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.kernel_order = 0;

        obj.eval = @(s,t) lap3d.kern(s, t, 'sprime');

        % sprime is the normal-derivative of S at the target: use
        % neumann eval with dpars = [1, 0]
        dpars_sp = [1.0; 0.0];
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'sp', sigma);
        obj.getquad  = @(S,eps,varargin) rsc_to_sparse( ...
                          lap_combprime_getquad(S, eps, dpars_sp, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.targ_fields = {'n'};

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.kernel_order = 1;

        obj.eval = @(s,t) lap3d.kern(s, t, 'dprime');

        % dpars for dirichlet eval: [alpha_S, beta_D] = [0, 1]
        dpars_dp = [0.0; 1.0];
        obj.fmm      = @(eps,src,targ,sigma) lap3d.fmm(eps, src, targ, 'dp', sigma);
        obj.getquad  = @(S,eps,varargin) rsc_to_sparse( ...
                          lap_combprime_getquad(S, eps, dpars_dp, varargin{:}), S);
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
        obj.getquad  = @(S,eps,varargin) rsc_to_sparse( ...
                          lap3d.dirichlet.get_quadrature_correction(S, eps, dpars_c, varargin{:}), S);
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
        obj.getquad  = @(S,eps,varargin) rsc_to_sparse( ...
                          lap_combprime_getquad(S, eps, dpars_c, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};
        obj.targ_fields = {'n'};
        
    otherwise
        error('KERNEL3D.LAP3D: unknown Laplace kernel type ''%s''.', type);

end

end

function spmat = rsc_to_sparse(Q, S)
%RSC_TO_SPARSE  Convert Laplace getquad RSC output to a sparse matrix.
% Laplace kernels are scalar (opdims=[1,1]); no ri argument needed.
spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear);
end

function p = lap_combprime_eval(S, sigma, targ, eps, rep_pars, opts)
%LAP_COMBPRIME_EVAL  Call lap3d.dirichlet.eval with iprime=1.
if nargin < 6 || isempty(opts), opts = struct(); end
opts.iprime = 1;
p = lap3d.dirichlet.eval(S, sigma, targ, eps, rep_pars, opts);
end

function Q = lap_combprime_getquad(S, eps, rep_pars, targinfo, opts)
%LAP_COMBPRIME_GETQUAD  Call get_quadrature_correction with iprime=1.
if nargin < 4 || isempty(targinfo), targinfo = S; end
if nargin < 5 || isempty(opts),     opts     = struct(); end
opts.iprime = 1;
Q = lap3d.dirichlet.get_quadrature_correction(S, eps, rep_pars, targinfo, opts);
end