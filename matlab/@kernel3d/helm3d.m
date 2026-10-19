function obj = helm3d(type,zk, coefs)
%KERNEL3D.HELM3D   Construct a Helmholtz kernel in 3D.
%
%   KERNEL3D.HELM3D(type,zk) or KERNEL3D.HELM3D(type,zk,coefs), where
%   zk is the (complex) wavenumber and type is one of:
%
%   Single-wavenumber scalar kernels  (zk scalar):
%      's'       - single layer,   G(x,y) = exp(ik|x-y|)/(4*pi*|x-y|)
%      'd'       - double layer,   D(x,y) = d/dn_y G(x,y)
%      'sp'      - sprime,         S'(x,y) = d/dn_x G(x,y)
%      'dprime'  - dprime,         D'(x,y) = d^2/dn_x dn_y G(x,y)
%      'c'       - combined,       coefs(1)*S + coefs(2)*D  (default [1;1])
%      'cprime'  - combined prime, coefs(1)*S' + coefs(2)*D'
%
%   Two-wavenumber difference scalar kernels  (zk = [zk0; zk1]):
%      's_diff'     - a0*S_{k0} - a1*S_{k1},    coefs=[a0;a1] default [1;1]
%      'd_diff'     - a0*D_{k0} - a1*D_{k1}
%      'sp_diff'    - a0*S'_{k0} - a1*S'_{k1}
%      'dp_diff'    - a0*D'_{k0} - a1*D'_{k1}
%      'c_diff'     - a0*C_{k0} - a1*C_{k1},    coefs=[comb_alpha;comb_beta;a0;a1]
%      'cp_diff'    - a0*C'_{k0} - a1*C'_{k1},  same coefs as c_diff
%
%   Transmission kernels (single wavenumber, zk scalar):
%      'trans_sys'  - [2x2] BIE kernel built from interleave of D,S,D',S'.
%                    coefs is 2x2, default ones(2,2).
%                    Block (i,j) = coefs(i,j) * K_ij where the kernel
%                    assignments are:
%                       (1,1) -> D    (1,2) -> S
%                       (2,1) -> D'   (2,2) -> S'
%                    Inherits layer_eval/getquad from scalar sub-kernels.
%      'trans_rep'  - [1x2] representation: u = D[rho] + S[sigma].
%                    Built from interleave([coefs(1)*D, coefs(2)*S]),
%                    coefs default [1;1].
%      's2trans'    - [2x1]  [S; S']     (col 5 of the BIE system)
%      'd2trans'    - [2x1]  [D; D']
%      'c2trans'    - [2x1]  [C; C'],    coefs=[alpha;beta] default [1;1]
%
%   Transmission diff kernels  (zk = [zk0; zk1]):
%      'trans_sys_diff' ('trans_diff') - [2x2] difference of trans_sys.
%                    coefs is 2x2x2, default ones(2,2,2).
%                    Block (i,j) = coefs(i,j,1)*K_ij_{k0} - coefs(i,j,2)*K_ij_{k1}
%                    where K assignments are the same as trans_sys.
%      's2trans_diff'   - [2x1]  [s_diff; sp_diff],  coefs=[a0;a1] default [1;1],
%                         with an fmm
%      'd2trans_diff'   - [2x1]  [d_diff; dp_diff],  coefs=[a0;a1] default [1;1]
%      'c2trans_diff'   - [2x1]  [c_diff; cp_diff],
%                         coefs=[comb_alpha;comb_beta;a0;a1] default [1;1;1;1]
%
% See also HELM3D.KERN, HELM3D.DIRICHLET.EVAL, HELM3D.TRANSMISSION.EVAL

if ( nargin < 2 )
    error('KERNEL3D.HELM3D: requires at least type and zk arguments.');
end

obj           = kernel3d();
obj.name      = 'helmholtz';
obj.opdims    = [1 1];
obj.zk        = zk;
obj.ifcomplex = 1;

if (nargin < 3)
    coefs = [];
end

switch lower(type)

    % ----------------------------------------------------------------
    %  Single-wavenumber scalar kernels
    % ----------------------------------------------------------------

    case {'s', 'single', 'slp'}
        obj.type         = 's';
        obj.kernel_order = -1;
        obj.eval         = @(s,t) helm3d.kern(zk, s, t, 's');
        obj.fmm          = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 's', sigma);
        rep_pars         = [1.0; 0.0];
        obj.getquad      = @(S,eps,varargin) rsc_to_sparse( ...
                               helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

    case {'d', 'double', 'dlp'}
        obj.type         = 'd';
        obj.kernel_order = 0;
        obj.eval         = @(s,t) helm3d.kern(zk, s, t, 'd');
        obj.fmm          = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'd', sigma);
        rep_pars         = [0.0; 1.0];
        obj.getquad      = @(S,eps,varargin) rsc_to_sparse( ...
                               helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields   = {'n'};

    case {'sp', 'sprime'}
        obj.type         = 'sp';
        obj.kernel_order = 0;
        obj.eval         = @(s,t) helm3d.kern(zk, s, t, 'sprime');
        obj.fmm          = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'sp', sigma);
        rep_pars         = [1.0; 0.0];
        obj.getquad      = @(S,eps,varargin) rsc_to_sparse( ...
                               helm_combprime_getquad(S, eps, zk, rep_pars, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.targ_fields  = {'n'};

    case {'dp', 'dprime'}
        obj.type         = 'dprime';
        obj.kernel_order = 1;
        obj.eval         = @(s,t) helm3d.kern(zk, s, t, 'dprime');
        obj.fmm          = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'dp', sigma);
        rep_pars         = [0.0; 1.0];
        obj.getquad      = @(S,eps,varargin) rsc_to_sparse( ...
                               helm_combprime_getquad(S, eps, zk, rep_pars, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields   = {'n'};
        obj.targ_fields  = {'n'};

    case {'c', 'combined', 'comb'}
        if isempty(coefs), coefs = [1; 1]; end
        coefs = coefs(:);
        obj.type         = 'c';
        obj.kernel_order = 0;
        obj.params.coefs = coefs;
        obj.eval         = @(s,t) helm3d.kern(zk, s, t, 'c', coefs);
        obj.fmm          = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'c', sigma, coefs);
        obj.getquad      = @(S,eps,varargin) rsc_to_sparse( ...
                               helm3d.dirichlet.get_quadrature_correction(S, eps, zk, coefs, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields   = {'n'};

    case {'cp', 'cprime', 'combprime'}
        if isempty(coefs), coefs = [1; 1]; end
        coefs = coefs(:);
        obj.type         = 'cprime';
        obj.kernel_order = 1;
        obj.params.coefs = coefs;
        obj.eval         = @(s,t) helm3d.kern(zk, s, t, 'cprime', coefs);
        obj.fmm          = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'cp', sigma, coefs);
        obj.getquad      = @(S,eps,varargin) rsc_to_sparse( ...
                               helm_combprime_getquad(S, eps, zk, coefs, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields   = {'n'};
        obj.targ_fields  = {'n'};

    % ----------------------------------------------------------------
    %  Two-wavenumber difference scalar kernels
    %  zk = [zk0; zk1],  coefs = [a0; a1]  (default [1;1])
    %  Result: a0*K_{k0} - a1*K_{k1}
    %  Built with kernel3d(g, type) on g = helm3d.helmdiff_green.
    % ----------------------------------------------------------------

    case {'s_diff'}
        [zk0, zk1, a0, a1] = parse_diff2(zk, coefs, 's_diff');
        obj = diff_kernel('s', [zk0; zk1], [a0; a1], []);
        obj.type = 's_diff';  obj.params.zks = [zk0;zk1];  obj.params.coefs = [a0;a1];

    case {'d_diff'}
        [zk0, zk1, a0, a1] = parse_diff2(zk, coefs, 'd_diff');
        obj = diff_kernel('d', [zk0; zk1], [a0; a1], []);
        obj.type = 'd_diff';  obj.params.zks = [zk0;zk1];  obj.params.coefs = [a0;a1];

    case {'sp_diff'}
        [zk0, zk1, a0, a1] = parse_diff2(zk, coefs, 'sp_diff');
        obj = diff_kernel('sp', [zk0; zk1], [a0; a1], []);
        obj.type = 'sp_diff';  obj.params.zks = [zk0;zk1];  obj.params.coefs = [a0;a1];

    case {'dp_diff'}
        [zk0, zk1, a0, a1] = parse_diff2(zk, coefs, 'dp_diff');
        obj = diff_kernel('dp', [zk0; zk1], [a0; a1], []);
        obj.type = 'dp_diff';  obj.params.zks = [zk0;zk1];  obj.params.coefs = [a0;a1];

    case {'c_diff'}
        % coefs = [comb_alpha; comb_beta; a0; a1]
        [zk0, zk1, comb, a0, a1] = parse_diff4(zk, coefs, 'c_diff');
        obj = diff_kernel('c', [zk0; zk1], [a0; a1], comb);
        obj.type = 'c_diff';  obj.params.zks = [zk0;zk1];  obj.params.coefs = [comb;a0;a1];

    case {'cp_diff'}
        [zk0, zk1, comb, a0, a1] = parse_diff4(zk, coefs, 'cp_diff');
        obj = diff_kernel('cp', [zk0; zk1], [a0; a1], comb);
        obj.type = 'cp_diff';  obj.params.zks = [zk0;zk1];  obj.params.coefs = [comb;a0;a1];

    % ----------------------------------------------------------------
    %  trans_sys  -- [2x2] BIE system kernel, single wavenumber
    %
    %  coefs is 2x2, default ones(2,2).
    %  Block layout:
    %     (1,1) = coefs(1,1) * D
    %     (1,2) = coefs(1,2) * S
    %     (2,1) = coefs(2,1) * D'
    %     (2,2) = coefs(2,2) * S'
    %
    %  Built via interleave.
    % ----------------------------------------------------------------
    case {'trans_sys', 'transmission_sys', 'trans', 'transmission', 'tsys'}
        if nargin < 3 || isempty(coefs), coefs = ones(2,2); end
        if ~isequal(size(coefs), [2 2])
            error('KERNEL3D.HELM3D: trans_sys coefs must be 2x2.');
        end
        Kblocks = [coefs(1,1) .* kernel3d.helm3d('d',  zk), ...
                   coefs(1,2) .* kernel3d.helm3d('s',  zk); ...
                   coefs(2,1) .* kernel3d.helm3d('dp', zk), ...
                   coefs(2,2) .* kernel3d.helm3d('sp', zk)];
        obj              = kernel3d.interleave(Kblocks);
        obj.type         = lower(type);
        obj.zk           = zk;
        obj.params.coefs = coefs;
        obj.fmm          = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'trans_sys', sigma, coefs);

    % ----------------------------------------------------------------
    %  trans_rep  -- [1x2] representation kernel
    %
    %  u(x) = coefs(1)*D[rho](x) + coefs(2)*S[sigma](x)
    %  coefs default [1;1].
    % ----------------------------------------------------------------
    case {'trans_rep', 'transmission_rep', 'trep'}
        if nargin < 3 || isempty(coefs), coefs = [1; 1]; end
        coefs = coefs(:);
        obj              = kernel3d.interleave([coefs(1) .* kernel3d.helm3d('d', zk), ...
                                                coefs(2) .* kernel3d.helm3d('s', zk)]);
        obj.type         = 'trans_rep';
        obj.zk           = zk;
        obj.params.coefs = coefs;
        obj.fmm          = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'trans_rep', sigma, coefs);

    % ----------------------------------------------------------------
    %  s2trans / d2trans / c2trans  -- [2x1] single-zk stacking kernels
    % ----------------------------------------------------------------
    case {'s2trans', 's2t'}
        obj      = kernel3d.interleave([kernel3d.helm3d('s', zk); kernel3d.helm3d('sp', zk)]);
        obj.type = 's2trans';
        obj.zk   = zk;
        obj.fmm  = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 's2trans', sigma);

    case {'d2trans', 'd2t'}
        obj      = kernel3d.interleave([kernel3d.helm3d('d', zk); kernel3d.helm3d('dp', zk)]);
        obj.type = 'd2trans';
        obj.zk   = zk;
        obj.fmm  = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'd2trans', sigma);

    case {'c2trans', 'c2t', 'c2tr'}
        if nargin < 3 || isempty(coefs), coefs = [1; 1]; end
        coefs = coefs(:);
        obj      = kernel3d.interleave([kernel3d.helm3d('c',  zk, coefs); ...
                                        kernel3d.helm3d('cp', zk, coefs)]);
        obj.type         = 'c2trans';
        obj.zk           = zk;
        obj.params.coefs = coefs;
        obj.fmm          = @(eps,src,targ,sigma) helm3d.fmm(eps, zk, src, targ, 'c2trans', sigma, coefs);

    % ----------------------------------------------------------------
    %  trans_sys_diff  -- [2x2], two wavenumbers
    %
    %  zk = [zk0; zk1]
    %  coefs is 2x2x2, default ones(2,2,2).
    %  Block (i,j) = coefs(i,j,1)*K_ij_{k0} - coefs(i,j,2)*K_ij_{k1}
    %  where K assignments are the same as trans_sys.
    % ----------------------------------------------------------------
    case {'trans_sys_diff', 'transmission_sys_diff', 'trans_diff', 'transmission_diff'}
        zk = zk(:);
        if length(zk) ~= 2
            error('KERNEL3D.HELM3D: trans_sys_diff requires zk = [zk0; zk1].');
        end
        zk0 = zk(1);  zk1 = zk(2);
        if nargin < 3 || isempty(coefs), coefs = ones(2,2,2); end
        if ~isequal(size(coefs), [2 2 2])
            error('KERNEL3D.HELM3D: trans_sys_diff coefs must be 2x2x2.');
        end
        c0 = coefs(:,:,1);   % weights on k0
        c1 = coefs(:,:,2);   % weights on k1
        zks = [zk0; zk1];
        rho = c1./c0;
        if all(c0(:) ~= 0) && all(abs(rho(:) - rho(1)) <= 1e-14*abs(rho(1)))
            % one difference green for all blocks: one kernel, one fit
            obj = diff_kernel('trans_sys', zks, [1; rho(1)], c0);
        else
            K11 = diff_kernel('d',  zks, [c0(1,1); c1(1,1)], []);
            K12 = diff_kernel('s',  zks, [c0(1,2); c1(1,2)], []);
            K21 = diff_kernel('dp', zks, [c0(2,1); c1(2,1)], []);
            K22 = diff_kernel('sp', zks, [c0(2,2); c1(2,2)], []);
            obj = kernel3d.interleave([K11, K12; K21, K22]);
        end
        obj.type         = 'trans_sys_diff';
        obj.zk           = max(abs(zk));
        obj.params.zks   = [zk0; zk1];
        obj.params.coefs = coefs;

    % ----------------------------------------------------------------
    %  s2trans_diff / d2trans_diff / c2trans_diff  -- [2x1], two wavenumbers
    %  zk = [zk0; zk1],  coefs = [a0; a1]  (or [alpha;beta;a0;a1] for c)
    % ----------------------------------------------------------------
    case {'s2trans_diff'}
        [zk0, zk1, a0, a1] = parse_diff2(zk, coefs, 's2trans_diff');
        obj = diff_kernel('s2trans', [zk0; zk1], [a0; a1], []);
        obj.type = 's2trans_diff';  obj.params.zks = [zk0;zk1];  obj.params.coefs = [a0;a1];

    case {'d2trans_diff'}
        [zk0, zk1, a0, a1] = parse_diff2(zk, coefs, 'd2trans_diff');
        obj = diff_kernel('d2trans', [zk0; zk1], [a0; a1], []);
        obj.type = 'd2trans_diff';  obj.params.zks = [zk0;zk1];  obj.params.coefs = [a0;a1];

    case {'c2trans_diff'}
        [zk0, zk1, comb, a0, a1] = parse_diff4(zk, coefs, 'c2trans_diff');
        obj = diff_kernel('c2trans', [zk0; zk1], [a0; a1], comb);
        obj.type = 'c2trans_diff';  obj.params.zks = [zk0;zk1];  obj.params.coefs = [comb;a0;a1];

    otherwise
        error('KERNEL3D.HELM3D: unknown Helmholtz kernel type ''%s''.', type);

end

end % function helm3d

function spmat = rsc_to_sparse(Q, S)
%RSC_TO_SPARSE  Convert Helmholtz getquad RSC output to a sparse matrix.
% Helmholtz kernels are scalar (opdims=[1,1]); no ri argument needed.
spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear);
end

% =========================================================================
%  Local helpers
% =========================================================================

function K = diff_kernel(type, zks, a, comb, g)
%DIFF_KERNEL  a(1)*K_{zks(1)} - a(2)*K_{zks(2)} for a helmholtz type, as
%  kernel3d(g, type) on g = helm3d.helmdiff_green, with the fmm
%  a(1)*helm3d.fmm(zks(1)) - a(2)*helm3d.fmm(zks(2)) where helm3d.fmm has
%  the type(s). comb is the coefs of 'c', 'cp', 'c2trans', 'trans_sys'.
    if nargin < 5 || isempty(g)
        g = @(s,t) helm3d.helmdiff_green(zks, a, s, t);
    end
    args = {'zk', max(abs(zks))};
    if ~isempty(comb)
        if isvector(comb), comb = comb(:).'; end
        args = [args, {'coefs', comb}];
    end
    K = kernel3d(g, type, args{:});

    % the difference of helmholtz fmms, one call per wavenumber
    fa = {};
    if ~isempty(comb), fa = {comb}; end
    if any(strcmp(type, {'c', 'cp', 'c2trans'})), fa = {comb(:)}; end
    K.fmm = @(eps, s, t, sigma) a(1)*helm3d.fmm(eps, zks(1), s, t, type, sigma, fa{:}) ...
                              - a(2)*helm3d.fmm(eps, zks(2), s, t, type, sigma, fa{:});
end

function [zk0, zk1, a0, a1] = parse_diff2(zk, coefs, name)
%PARSE_DIFF2  Parse zk=[zk0;zk1] and coefs=[a0;a1] for simple diff kernels.
    zk = zk(:);
    if length(zk) ~= 2
        error('KERNEL3D.HELM3D: %s requires zk = [zk0; zk1].', name);
    end
    zk0 = zk(1);  zk1 = zk(2);
    if isempty(coefs), coefs = [1; 1]; end
    coefs = coefs(:);
    if length(coefs) ~= 2
        error('KERNEL3D.HELM3D: %s coefs must have 2 entries [a0; a1].', name);
    end
    a0 = coefs(1);  a1 = coefs(2);
end

function [zk0, zk1, comb, a0, a1] = parse_diff4(zk, coefs, name)
%PARSE_DIFF4  Parse zk=[zk0;zk1] and coefs=[alpha;beta;a0;a1].
    zk = zk(:);
    if length(zk) ~= 2
        error('KERNEL3D.HELM3D: %s requires zk = [zk0; zk1].', name);
    end
    zk0 = zk(1);  zk1 = zk(2);
    if isempty(coefs), coefs = [1; 1; 1; 1]; end
    coefs = coefs(:);
    if length(coefs) ~= 4
        error('KERNEL3D.HELM3D: %s coefs must have 4 entries [alpha;beta;a0;a1].', name);
    end
    comb = coefs(1:2);
    a0   = coefs(3);
    a1   = coefs(4);
end

function Q = helm_combprime_getquad(S, eps, zk, rep_pars, targinfo, opts)
%HELM_COMBPRIME_GETQUAD  Call get_quadrature_correction with iprime=1.
    if nargin < 5 || isempty(targinfo), targinfo = S; end
    if nargin < 6 || isempty(opts),     opts     = struct(); end
    opts.iprime = 1;
    Q = helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars, targinfo, opts);
end
