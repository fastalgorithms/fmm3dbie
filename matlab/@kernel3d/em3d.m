function obj = em3d(type, zk, rep_params)
%KERNEL3D.EM3D   Construct a Maxwell (EM3D) kernel in 3D.
%
%   KERNEL3D.EM3D(type, zk) or KERNEL3D.EM3D(type, zk, rep_params)
%
%   Currently supported PEC (perfect electric conductor) kernels:
%
%     'nrccie-bc'    Non-resonant charge-current integral equation,
%                    boundary condition (system) operator.
%                    Representation:
%                      H = \nabla \times S_{k}[J]
%                      E = ik S_{k}[J] - \nabla S_{k}[\rho]
%                    Density: [j_ru, j_rv, rho]  (3 components per point)
%                      j_ru = coeff of ru = du/|du| (orthonormal source frame)
%                      j_rv = coeff of rv = n x ru
%                    Output:  [pot_ru, pot_rv, pot_rho]  (3 equations per point)
%                      pot_ru/rv projected onto orthonormal target frame
%                    opdims = [3, 3]
%                    rep_params = alpha  (scalar regularisation parameter)
%
%     'nrccie-eval'  NRCCIE field evaluation kernel.
%                    Computes E and H from the Cartesian density
%                    [J_x, J_y, J_z, rho]  (4 components per point).
%                    Output: [E_x, E_y, E_z, H_x, H_y, H_z] (6 components).
%                    opdims = [6, 4]
%                    rep_params unused (pass [] or omit).
%
%   Parameters
%   ----------
%     zk         - complex wavenumber
%     rep_params - for 'nrccie-bc': scalar alpha
%                  for 'nrccie-eval': unused
%
%   Kernel evaluation  (K.eval(srcinfo, targinfo))
%   ------------------------------------------------
%     srcinfo / targinfo are structs with fields .r (3,ns/nt),
%     .n (3,:), .du (3,:), .dv (3,:).
%     Returns an (opdims(1)*nt) x (opdims(2)*ns) matrix.
%
%   FMM + quadrature layer-potential eval  (K.layer_eval)
%   -------------------------------------------------------
%     K.layer_eval(S, sigma, targinfo, eps)
%     Delegates to em3d.pec.eval.
%
%   Quadrature corrections  (K.getquad)
%   --------------------------------------
%     K.getquad(S, eps, varargin)
%     Delegates to em3d.pec.get_quadrature_correction.
%
% See also EM3D.PEC.EVAL, EM3D.PEC.GET_QUADRATURE_CORRECTION

if nargin < 2
    error('KERNEL3D.EM3D: requires at least type and zk arguments.');
end
if nargin < 3
    rep_params = [];
end

obj           = kernel3d();
obj.name      = 'maxwell';
obj.zk        = real(zk);   % use real part for oversampling decisions
obj.ifcomplex = 1;

switch lower(type)

    case 'nrccie-bc'
        % ----------------------------------------------------------------
        % System (boundary condition) operator for NRCCIE.
        %   Density : [j_u, j_v, rho]  (3 components / point)
        %   Output  : [M_u, M_v, E_n]  (3 equations / point)
        %   Quadrature: 9 kernels stored in wnear(9, nquad), column-major
        %               3x3: wnear(k) = block entry (mod(k-1,3)+1, ceil(k/3))
        % ----------------------------------------------------------------
        if isempty(rep_params)
            error('KERNEL3D.EM3D: nrccie-bc requires rep_params = alpha.');
        end
        alpha = rep_params(1);

        obj.type         = 'nrccie-bc';
        obj.kernel_order = 0;
        obj.opdims       = [3, 3];
        obj.params.alpha = alpha;
        obj.params.zk    = zk;

        obj.eval = @(s,t) em3d.kern(zk, s, t, 'nrccie-bc', alpha);

        obj.fmm = @(eps,src,targ,sigma) em3d.fmm(eps, zk, src, targ, 'nrccie-bc', sigma, alpha);

        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
            em3d_nrccie_bc_layer_eval(S, sigma, targ, eps, zk, alpha, varargin);

        obj.getquad = @(S,eps,varargin) ...
            em3d_nrccie_bc_getquad(S, eps, zk, alpha, varargin);

        obj.get_overs_orders = @(S,t,eps) ...
            kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

        % wnear(9,nquad): full 3x3 in column-major order within block
        obj.rsc_to_interleave = kernel3d.rsc_interleave_full(3, 3);

        obj.src_fields = {'n', 'du', 'dv'};
        obj.targ_fields = {'n', 'du', 'dv'};

    case 'nrccie-eval'
        % ----------------------------------------------------------------
        % Field evaluation operator for NRCCIE.
        %   Density : [J_x, J_y, J_z, rho]  (4 Cartesian components / point)
        %   Output  : [E_x, E_y, E_z, H_x, H_y, H_z]  (6 components)
        %   Quadrature: 4 kernels (S_k and its 3 Cartesian gradients)
        % ----------------------------------------------------------------
        obj.type         = 'nrccie-eval';
        obj.kernel_order = -1;
        obj.opdims       = [6, 4];
        obj.params.zk    = zk;

        obj.eval = @(s,t) em3d.kern(zk, s, t, 'nrccie-eval');

        obj.fmm = @(eps,src,targ,sigma) em3d.fmm(eps, zk, src, targ, 'nrccie-eval', sigma);

        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
            em3d_nrccie_eval_layer_eval(S, sigma, targ, eps, zk, varargin);

        obj.getquad = @(S,eps,varargin) ...
            em3d_nrccie_eval_getquad(S, eps, zk, varargin);

        obj.get_overs_orders = @(S,t,eps) ...
            kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

        obj.src_fields = {'n', 'du', 'dv'};
        obj.targ_fields = {};
        % wnear(4,nquad) stores scalar basis kernels [S_k, dx S_k, dy S_k, dz S_k]
        % combined with frequency-dependent coefficients to fill the 6x4 block.
        obj.rsc_to_interleave = kernel3d.rsc_interleave_nrccie_eval(zk);

    otherwise
        error('KERNEL3D.EM3D: unknown kernel type ''%s''.', type);
end

end

% =========================================================================
% NRCCIE-BC helpers
% =========================================================================

function Q = em3d_nrccie_bc_getquad(S, eps, zk, alpha, args)
%EM3D_NRCCIE_BC_GETQUAD  Near-quadrature for nrccie-bc.
%  Wraps em3d.pec.get_quadrature_correction with rep = 'nrccie-bc'.
if nargin < 5 || isempty(args), args = {}; end
opts_use = struct();
if length(args) >= 2
    opts_use = args{2};
elseif length(args) >= 1 && isstruct(args{1}) && ~isfield(args{1},'r')
    % single arg and it looks like an opts struct, not targinfo
    opts_use = args{1};
end
opts_use.rep = 'nrccie-bc';
% nrccie-bc uses surface itself as targets (self-quadrature)
Q = em3d.pec.get_quadrature_correction(S, eps, zk, alpha, S, opts_use);
end

function p = em3d_nrccie_bc_layer_eval(S, sigma, targ, eps, zk, alpha, args)
%EM3D_NRCCIE_BC_LAYER_EVAL  Layer-potential for nrccie-bc via em3d.pec.eval.
%
%  Uses the orthonormal-frame convention matching the Fortran routines and
%  em3d.kern 'nrccie-bc'.
%
%  sigma : (3, npts) or (3*npts, 1) density [j_ru; j_rv; rho]
%    j_ru = coefficient of ru_s = du_s / |du_s|
%    j_rv = coefficient of rv_s = n_s x ru_s
%
%  Output: (3*nt, 1) interleaved [pot_ru(1); pot_rv(1); pot_rho(1); ...]
%    pot_ru = zvec3 . ru_t,  pot_rv = zvec3 . rv_t
%    where ru_t = du_t/|du_t|,  rv_t = n_t x ru_t
%
%  em3d.pec.eval with rep='nrccie-bc' calls lpcomp_em_nrccie_pec_addsub_targ
%  which expects and returns exactly this orthonormal-frame convention, so
%  no frame conversion is needed here.
if nargin < 7 || isempty(args), args = {}; end

npts      = S.npts;
densities = reshape(sigma, 3, npts);   % (3, npts) = [j_ru; j_rv; rho]

opts_pe   = struct('rep', 'nrccie-bc');
[E_pe, ~] = em3d.pec.eval(S, densities, targ, eps, zk, alpha, opts_pe);
% E_pe = (3, ntarg) = [pot_ru; pot_rv; pot_rho] in orthonormal target frame

% Output interleaved: [pot_ru(1); pot_rv(1); pot_rho(1); pot_ru(2); ...]
p = E_pe(:);
end

% =========================================================================
% NRCCIE-EVAL helpers
% =========================================================================

function Q = em3d_nrccie_eval_getquad(S, eps, zk, args)
%EM3D_NRCCIE_EVAL_GETQUAD  Near-quadrature for nrccie-eval.
if nargin < 4 || isempty(args), args = {}; end
% Separate targinfo and opts
if length(args) >= 2
    targinfo = args{1};
    opts_use = args{2};
elseif length(args) == 1
    if isstruct(args{1}) && isfield(args{1}, 'r')
        targinfo = args{1};
        opts_use = struct();
    elseif isa(args{1}, 'surfer')
        targinfo = args{1};
        opts_use = struct();
    else
        targinfo = S;
        opts_use = args{1};
    end
else
    targinfo = S;
    opts_use = struct();
end
if isempty(targinfo), targinfo = S; end
opts_use.rep = 'nrccie-eval';
Q = em3d.pec.get_quadrature_correction(S, eps, zk, [], targinfo, opts_use);
end

function p = em3d_nrccie_eval_layer_eval(S, sigma, targ, eps, zk, args)
%EM3D_NRCCIE_EVAL_LAYER_EVAL  FMM layer-potential for nrccie-eval.
%
%  sigma : (4, npts) or (4*npts, 1) density [Jx; Jy; Jz; rho].
%  Output: (6*nt, 1) interleaved [Ex(1);Ey(1);Ez(1);Hx(1);Hy(1);Hz(1); Ex(2);...].
%
%  Delegates to em3d.pec.eval with rep='nrccie'.
%  densities expected as (4, npts): rows = [Jx; Jy; Jz; rho].
if nargin < 6 || isempty(args), args = {}; end

npts      = S.npts;
densities = reshape(sigma, 4, npts);   % (4, npts)

opts     = struct('rep', 'nrccie');
[E, H]   = em3d.pec.eval(S, densities, targ, eps, zk, 0, opts);
% E, H are each (3, nt)

EH = [E; H];   % (6, nt): rows = Ex,Ey,Ez,Hx,Hy,Hz
p  = EH(:);    % (6*nt, 1) interleaved
end

