function obj = helm3d(type,zk, coefs)
%KERNEL3D.HELM3D   Construct a Helmholtz kernel in 3D.
%
%   KERNEL3D.HELM3D(type,zk) or KERNEL3D.HELM3D(zk, type, coefs), where
%   zk is the (complex) wavenumber and type is one of:
%
%   Single-wavenumber kernels  (zk is scalar):
%      's'           - single layer,   S(x,y) = exp(ik|x-y|)/(4*pi*|x-y|)
%      'd'           - double layer,   D(x,y) = d/dn_y G(x,y)
%      'sp'          - sprime,         S'(x,y) = d/dn_x G(x,y)
%      'dprime'      - dprime,         D'(x,y) = d^2/dn_x dn_y G(x,y)
%      'c'           - combined layer, coefs(1)*S + coefs(2)*D  (default coefs=[1 1])
%      'cprime'      - combined prime, coefs(1)*S' + coefs(2)*D'
%
%   Transmission kernels  (single wavenumber, zk scalar):
%      'trans_sys'   - [2x2] BIE system kernel using the transmission
%                      representation with the same wavenumber on both sides.
%                      coefs = [alpha0;beta0;alpha1;beta1]  (default [1;1;1;1])
%                      eval  = helm3d.kern 'eval' (representation matrix)
%                      layer_eval/getquad via helm3d.transmission
%      'trans_rep'   - [1x2] representation kernel: evaluates
%                        u(x) = D_{zk}[rho](x) + i*zk*S_{zk}[sigma](x)
%                      built from interleave of d and (i*zk)*s kernels
%      's2trans'     - [2x1]  [S_{zk}; S'_{zk}]  (single → transmission)
%      'd2trans'     - [2x1]  [D_{zk}; D'_{zk}]  (double → transmission)
%      'c2trans'     - [2x1]  [C_{zk}; C'_{zk}]  coefs=[alpha,beta] default [1;1]
%
%   Two-wavenumber difference kernels  (zk = [zk0; zk1]):
%      'trans_sys_diff' - [2x2] BIE system kernel built from difference kernels
%                         S_{k0}-S_{k1}, D_{k0}-D_{k1}, etc.
%                         coefs = [alpha0;beta0;alpha1;beta1]  (default [1;1;1;1])
%                         eval subtracts two single-zk representation matrices;
%                         layer_eval/getquad via helm3d.transmission with [zk0;zk1]
%      's2trans_diff'   - [2x1]  [a0*S_{k0}-a1*S_{k1}; a0*S'_{k0}-a1*S'_{k1}]
%                         coefs = [a0;a1]  (default [1;1])
%      'd2trans_diff'   - [2x1]  [a0*D_{k0}-a1*D_{k1}; a0*D'_{k0}-a1*D'_{k1}]
%                         coefs = [a0;a1]  (default [1;1])
%      'c2trans_diff'   - [2x1]  [a0*C_{k0}-a1*C_{k1}; a0*C'_{k0}-a1*C'_{k1}]
%                         coefs = [alpha;beta;a0;a1]  (default [1;1;1;1])
%
% Kernel evaluation:
%   obj.eval(srcinfo, targinfo)
%      srcinfo / targinfo are structs with fields .r (3,:) and .n (3,:).
%      Returns an (opdims(1)*nt) x (opdims(2)*ns) matrix.
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
% Kernel orders (obj.kernel_order):
%   s -> -1,  d/sp -> 0,  dp/c/cp -> 1,
%   trans_sys/trans_sys_diff -> 1,  trans_rep -> 0
%   s2trans/d2trans/c2trans/diff variants -> 1
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
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            helm_combprime_eval(S, sigma, targ, eps, zk, rep_pars_sp, varargin);
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
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            helm_combprime_eval(S, sigma, targ, eps, zk, rep_pars_dp, varargin);
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
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            helm_combprime_eval(S, sigma, targ, eps, zk, rep_pars_cp, varargin);
        obj.getquad  = @(S,eps,varargin) ...
                          helm_combprime_getquad(S, eps, zk, rep_pars_cp, varargin);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};
        obj.targ_fields = {'n'};

    % ------------------------------------------------------------------
    %  trans_sys  -- [2x2] BIE system kernel, single wavenumber
    %
    %  Representation (both sides use the same zk):
    %    u = 1/beta * (D_{zk}[rho] + i*zk*S_{zk}[sigma])
    %
    %  BIE system (from transmission BC):
    %    Row 1: alpha0/beta0*D_{zk} - alpha1/beta1*D_{zk}  (Dirichlet)
    %    Row 2: D'_{zk} - D'_{zk}  ... (Neumann)
    %  (all differences vanish for identical zk; the meaningful quadrature
    %   is provided by helm3d.transmission with zks=[zk;zk])
    %
    %  eval returns the [2x2] representation matrix (helm3d.kern 'eval').
    %  layer_eval/getquad use helm3d.transmission.{eval,get_quadrature_correction}.
    % ------------------------------------------------------------------
    case {'trans_sys', 'transmission_sys', 'trans', 'transmission'}
        if ( nargin < 3 )
            warning('KERNEL3D.HELM3D: missing rep_params for trans_sys. Defaulting to [1;1;1;1].');
            coefs = [1; 1; 1; 1];
        end
        rep_params = coefs(:);
        if ( length(rep_params) ~= 4 )
            error('KERNEL3D.HELM3D: trans_sys coefs must have 4 entries [alpha0;beta0;alpha1;beta1].');
        end

        obj.type           = lower(type);
        obj.kernel_order   = 1;
        obj.opdims         = [2 2];
        obj.params.rep_params = rep_params;

        % eval: (2*nt) x (2*ns) representation matrix
        % Columns interleaved [rho, sigma]: odd = coef*D_{zk}, even = S_{zk}
        coef_trans = 1i * zk;
        obj.eval = @(s,t) helm3d.kern(zk, s, t, 'eval', coef_trans);

        % layer_eval and getquad via helm3d.transmission (zks=[zk;zk])
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            helm3d.transmission.eval(S, sigma, targ, eps, [zk; zk], rep_params, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          helm3d.transmission.get_quadrature_correction(S, eps, [zk; zk], rep_params, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields  = {'n'};
        obj.targ_fields = {'n'};

        % rsc_to_interleave for the 4-kernel transmission quadrature output.
        % From lpcomp_helm_comb_trans_addsub (Fortran lines 667-673):
        %   wnear(1,:) acts on sigma(2,:)=lambda  -> block (1,2)
        %   wnear(2,:) acts on sigma(1,:)=rho     -> block (1,1)
        %   wnear(3,:) acts on sigma(2,:)=lambda  -> block (2,2)
        %   wnear(4,:) acts on sigma(1,:)=rho     -> block (2,1)
        e = struct('ker_id',{},'row_id',{},'col_id',{},'coef',{});
        e(end+1) = struct('ker_id',1,'row_id',1,'col_id',2,'coef',1);
        e(end+1) = struct('ker_id',2,'row_id',1,'col_id',1,'coef',1);
        e(end+1) = struct('ker_id',3,'row_id',2,'col_id',2,'coef',1);
        e(end+1) = struct('ker_id',4,'row_id',2,'col_id',1,'coef',1);
        obj.rsc_to_interleave = kernel3d.rsc_interleave_basis(2, 2, 4, e);

    % ------------------------------------------------------------------
    %  trans_sys_diff  -- [2x2] BIE system kernel, two wavenumbers
    %
    %  zk = [zk0; zk1]  (exterior and interior wavenumbers)
    %  coefs = [alpha0; beta0; alpha1; beta1]  (default [1;1;1;1])
    %
    %  eval subtracts the two single-wavenumber representation matrices:
    %    M_diff = M_{zk0}(coef0) - M_{zk1}(coef1)
    %  where coef_j = i*zk_j.
    %
    %  layer_eval/getquad use helm3d.transmission with the two distinct zks.
    % ------------------------------------------------------------------
    case {'trans_sys_diff', 'transmission_sys_diff'}
        zk = zk(:);
        if ( length(zk) ~= 2 )
            error('KERNEL3D.HELM3D: trans_sys_diff requires zk = [zk0; zk1].');
        end
        zk0 = zk(1);
        zk1 = zk(2);
        if ( nargin < 3 )
            warning('KERNEL3D.HELM3D: missing rep_params for trans_sys_diff. Defaulting to [1;1;1;1].');
            coefs = [1; 1; 1; 1];
        end
        rep_params = coefs(:);
        if ( length(rep_params) ~= 4 )
            error('KERNEL3D.HELM3D: trans_sys_diff coefs must have 4 entries [alpha0;beta0;alpha1;beta1].');
        end

        obj.type           = 'trans_sys_diff';
        obj.kernel_order   = 1;
        obj.opdims         = [2 2];
        obj.zk             = zk0;   % store exterior wavenumber as primary
        obj.params.zks        = [zk0; zk1];
        obj.params.rep_params = rep_params;

        % eval: difference of two representation matrices
        coef0 = 1i * zk0;
        coef1 = 1i * zk1;
        obj.eval = @(s,t) helm3d.kern(zk0, s, t, 'eval', coef0) - ...
                          helm3d.kern(zk1, s, t, 'eval', coef1);

        % layer_eval/getquad via helm3d.transmission with [zk0; zk1]
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            helm3d.transmission.eval(S, sigma, targ, eps, [zk0; zk1], rep_params, varargin{:});
        obj.getquad  = @(S,eps,varargin) ...
                          helm3d.transmission.get_quadrature_correction(S, eps, [zk0; zk1], rep_params, varargin{:});
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, zk0, obj.kernel_order);
        obj.src_fields  = {'n'};
        obj.targ_fields = {'n'};

        % Same rsc_to_interleave as trans_sys (same Fortran routine)
        e = struct('ker_id',{},'row_id',{},'col_id',{},'coef',{});
        e(end+1) = struct('ker_id',1,'row_id',1,'col_id',2,'coef',1);
        e(end+1) = struct('ker_id',2,'row_id',1,'col_id',1,'coef',1);
        e(end+1) = struct('ker_id',3,'row_id',2,'col_id',2,'coef',1);
        e(end+1) = struct('ker_id',4,'row_id',2,'col_id',1,'coef',1);
        obj.rsc_to_interleave = kernel3d.rsc_interleave_basis(2, 2, 4, e);

    % ------------------------------------------------------------------
    %  trans_rep  -- [1x2] representation kernel (single wavenumber)
    %
    %  u(x) = D_{zk}[rho](x) + i*zk * S_{zk}[sigma](x)
    %  opdims = [1, 2]: row = scalar field value, cols = [rho, sigma]
    %  Built from interleave of D and (i*zk)*S scalar kernels.
    % ------------------------------------------------------------------
    case {'trans_rep', 'transmission_rep'}
        coef_s = 1i * zk;
        Ks_rep = kernel3d.helm3d('s', zk);
        Kd_rep = kernel3d.helm3d('d', zk);
        KikS   = coef_s .* Ks_rep;
        obj = kernel3d.interleave([Kd_rep, KikS]);
        obj.type = 'trans_rep';

    % ------------------------------------------------------------------
    %  s2trans  -- [2x1] single-layer to transmission: [S_{zk}; S'_{zk}]
    % ------------------------------------------------------------------
    case {'s2trans'}
        Ks  = kernel3d.helm3d('s',  zk);
        Ksp = kernel3d.helm3d('sp', zk);
        obj = kernel3d.interleave([Ks; Ksp]);
        obj.type = 's2trans';

    % ------------------------------------------------------------------
    %  d2trans  -- [2x1] double-layer to transmission: [D_{zk}; D'_{zk}]
    % ------------------------------------------------------------------
    case {'d2trans'}
        Kd  = kernel3d.helm3d('d',  zk);
        Kdp = kernel3d.helm3d('dp', zk);
        obj = kernel3d.interleave([Kd; Kdp]);
        obj.type = 'd2trans';

    % ------------------------------------------------------------------
    %  c2trans  -- [2x1] combined-layer to transmission
    %  coefs = [alpha; beta]  (default [1;1])
    %  Output: [alpha*S+beta*D; alpha*S'+beta*D']
    % ------------------------------------------------------------------
    case {'c2trans'}
        if ( nargin < 3 )
            warning('KERNEL3D.HELM3D: missing coefs for c2trans. Defaulting to [1 1].');
            coefs = [1; 1];
        end
        coefs = coefs(:);
        Kc  = kernel3d.helm3d('c',  zk, coefs);
        Kcp = kernel3d.helm3d('cp', zk, coefs);
        obj = kernel3d.interleave([Kc; Kcp]);
        obj.type = 'c2trans';

    % ------------------------------------------------------------------
    %  s2trans_diff  -- [2x1] difference version, two wavenumbers
    %  zk = [zk0; zk1]
    %  coefs = [a0; a1]  (default [1;1])
    %  Output: [a0*S_{k0} - a1*S_{k1}; a0*S'_{k0} - a1*S'_{k1}]
    %  eval only (no layer_eval/getquad: no dedicated Fortran quadrature)
    % ------------------------------------------------------------------
    case {'s2trans_diff'}
        zk = zk(:);
        if ( length(zk) ~= 2 )
            error('KERNEL3D.HELM3D: s2trans_diff requires zk = [zk0; zk1].');
        end
        zk0 = zk(1);  zk1 = zk(2);
        if ( nargin < 3 )
            warning('KERNEL3D.HELM3D: missing coefs for s2trans_diff. Defaulting to [1;1].');
            coefs = [1; 1];
        end
        coefs = coefs(:);
        if ( length(coefs) ~= 2 )
            error('KERNEL3D.HELM3D: s2trans_diff coefs must have 2 entries [a0;a1].');
        end
        a0 = coefs(1);  a1 = coefs(2);

        Ks0  = a0 .* kernel3d.helm3d('s',  zk0);
        Ksp0 = a0 .* kernel3d.helm3d('sp', zk0);
        Ks1  = a1 .* kernel3d.helm3d('s',  zk1);
        Ksp1 = a1 .* kernel3d.helm3d('sp', zk1);

        Kdiff_s  = Ks0  - Ks1;
        Kdiff_sp = Ksp0 - Ksp1;
        obj = kernel3d.interleave([Kdiff_s; Kdiff_sp]);
        obj.type             = 's2trans_diff';
        obj.zk               = zk0;
        obj.params.zks       = [zk0; zk1];
        obj.params.coefs     = coefs;

    % ------------------------------------------------------------------
    %  d2trans_diff  -- [2x1] difference version, two wavenumbers
    %  zk = [zk0; zk1]
    %  coefs = [a0; a1]  (default [1;1])
    %  Output: [a0*D_{k0} - a1*D_{k1}; a0*D'_{k0} - a1*D'_{k1}]
    %  eval only (no layer_eval/getquad)
    % ------------------------------------------------------------------
    case {'d2trans_diff'}
        zk = zk(:);
        if ( length(zk) ~= 2 )
            error('KERNEL3D.HELM3D: d2trans_diff requires zk = [zk0; zk1].');
        end
        zk0 = zk(1);  zk1 = zk(2);
        if ( nargin < 3 )
            warning('KERNEL3D.HELM3D: missing coefs for d2trans_diff. Defaulting to [1;1].');
            coefs = [1; 1];
        end
        coefs = coefs(:);
        if ( length(coefs) ~= 2 )
            error('KERNEL3D.HELM3D: d2trans_diff coefs must have 2 entries [a0;a1].');
        end
        a0 = coefs(1);  a1 = coefs(2);

        Kd0  = a0 .* kernel3d.helm3d('d',  zk0);
        Kdp0 = a0 .* kernel3d.helm3d('dp', zk0);
        Kd1  = a1 .* kernel3d.helm3d('d',  zk1);
        Kdp1 = a1 .* kernel3d.helm3d('dp', zk1);

        Kdiff_d  = Kd0  - Kd1;
        Kdiff_dp = Kdp0 - Kdp1;
        obj = kernel3d.interleave([Kdiff_d; Kdiff_dp]);
        obj.type             = 'd2trans_diff';
        obj.zk               = zk0;
        obj.params.zks       = [zk0; zk1];
        obj.params.coefs     = coefs;

    % ------------------------------------------------------------------
    %  c2trans_diff  -- [2x1] difference version, two wavenumbers
    %  zk = [zk0; zk1]
    %  coefs = [alpha; beta; a0; a1]  (default [1;1;1;1])
    %    alpha, beta: combined-layer coefficients (same for both sides)
    %    a0, a1:      subtraction weights
    %  Output: [a0*C_{k0} - a1*C_{k1}; a0*C'_{k0} - a1*C'_{k1}]
    %  where C = alpha*S + beta*D
    %  eval only (no layer_eval/getquad)
    % ------------------------------------------------------------------
    case {'c2trans_diff'}
        zk = zk(:);
        if ( length(zk) ~= 2 )
            error('KERNEL3D.HELM3D: c2trans_diff requires zk = [zk0; zk1].');
        end
        zk0 = zk(1);  zk1 = zk(2);
        if ( nargin < 3 )
            warning('KERNEL3D.HELM3D: missing coefs for c2trans_diff. Defaulting to [1;1;1;1].');
            coefs = [1; 1; 1; 1];
        end
        coefs = coefs(:);
        if ( length(coefs) ~= 4 )
            error('KERNEL3D.HELM3D: c2trans_diff coefs must have 4 entries [alpha;beta;a0;a1].');
        end
        comb_coefs = coefs(1:2);
        a0 = coefs(3);  a1 = coefs(4);

        Kc0  = a0 .* kernel3d.helm3d('c',  zk0, comb_coefs);
        Kcp0 = a0 .* kernel3d.helm3d('cp', zk0, comb_coefs);
        Kc1  = a1 .* kernel3d.helm3d('c',  zk1, comb_coefs);
        Kcp1 = a1 .* kernel3d.helm3d('cp', zk1, comb_coefs);

        Kdiff_c  = Kc0  - Kc1;
        Kdiff_cp = Kcp0 - Kcp1;
        obj = kernel3d.interleave([Kdiff_c; Kdiff_cp]);
        obj.type             = 'c2trans_diff';
        obj.zk               = zk0;
        obj.params.zks       = [zk0; zk1];
        obj.params.coefs     = coefs;

    otherwise
        error('KERNEL3D.HELM3D: unknown Helmholtz kernel type ''%s''.', type);

end

if isempty(obj.rsc_to_interleave)
    obj.rsc_to_interleave = kernel3d.rsc_interleave_scalar();
end

end

function p = helm_combprime_eval(S, sigma, targ, eps, zk, rep_pars, args)
%HELM_COMBPRIME_EVAL  Call helm3d.dirichlet.eval with iprime=1.
if nargin < 7 || isempty(args), args = {}; end
% Ensure opts struct is present and has iprime=1
if isempty(args) || ~isstruct(args{end})
    args{end+1} = struct();
end
args{end}.iprime = 1;
p = helm3d.dirichlet.eval(S, sigma, targ, eps, zk, rep_pars, args{:});
end

function Q = helm_combprime_getquad(S, eps, zk, rep_pars, args)
%HELM_COMBPRIME_GETQUAD  Call get_quadrature_correction with iprime=1.
if nargin < 5 || isempty(args), args = {}; end
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
Q = helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars, targinfo, opts);
end
