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
%                    Density: [j_u, j_v, rho]  (3 components per point)
%                    Output:  [M_u, M_v, E_n]  (3 equations per point)
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
obj.src_fields = {'n', 'du', 'dv'};
obj.targ_fields = {'n', 'du', 'dv'};

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

        % nrccie-eval is a pure evaluation kernel (not a BIE system operator);
        % its wnear(4,nquad) stores 4 scalar basis kernels [S_k, dx, dy, dz] that
        % are combined with a frequency-dependent coefficient pattern to fill the
        % 6x4 block.  rsc_to_interleave is left empty until a 'basis' assembler
        % is implemented in conv_rsc_to_spmat.
        % obj.rsc_to_interleave = [];   % (intentionally unset)

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
%EM3D_NRCCIE_BC_LAYER_EVAL  FMM + quadrature for nrccie-bc.
if nargin < 7 || isempty(args), args = {}; end
opts_use = struct();
if ~isempty(args) && isstruct(args{end})
    opts_use = args{end};
end
opts_use.rep = 'nrccie';
% sigma is (3, npts): [ju; jv; rho] per point
% em3d.pec.eval expects (4, npts) in Cartesian representation
% For BC evaluation on the surface we call the solver path via pec.eval
[p_cell] = em3d.pec.eval(S, sigma, targ, eps, zk, alpha, opts_use);
% p_cell contains [E; H] (6 x ntarg), return it
p = p_cell;
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
%EM3D_NRCCIE_EVAL_LAYER_EVAL  FMM + quadrature for nrccie-eval.
if nargin < 6 || isempty(args), args = {}; end
opts_use = struct();
if ~isempty(args) && isstruct(args{end})
    opts_use = args{end};
end
opts_use.rep = 'nrccie';
[E, H] = em3d.pec.eval(S, sigma, targ, eps, zk, [], opts_use);
p = [E; H];
end

