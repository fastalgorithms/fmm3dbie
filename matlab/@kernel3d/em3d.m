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

        obj.getquad = @(S,eps,varargin) rsc_to_sparse_bc( ...
            em3d_nrccie_bc_getquad(S, eps, zk, alpha, varargin{:}), S);

        obj.get_overs_orders = @(S,t,eps) ...
            kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

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

        obj.getquad = @(S,eps,varargin) rsc_to_sparse_eval( ...
            em3d_nrccie_eval_getquad(S, eps, zk, varargin{:}), S, zk);

        obj.get_overs_orders = @(S,t,eps) ...
            kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

        obj.src_fields = {'n', 'du', 'dv'};
        obj.targ_fields = {};

    otherwise
        error('KERNEL3D.EM3D: unknown kernel type ''%s''.', type);
end

end

% =========================================================================
% NRCCIE-BC helper
% =========================================================================

function Q = em3d_nrccie_bc_getquad(S, eps, zk, alpha, targinfo, opts)
%EM3D_NRCCIE_BC_GETQUAD  Near-quadrature for nrccie-bc.
if nargin < 5 || isempty(targinfo), targinfo = S; end
if nargin < 6 || isempty(opts),     opts     = struct(); end
opts.rep = 'nrccie-bc';
Q = em3d.pec.get_quadrature_correction(S, eps, zk, alpha, targinfo, opts);
end

% =========================================================================
% NRCCIE-EVAL helper
% =========================================================================

function Q = em3d_nrccie_eval_getquad(S, eps, zk, targinfo, opts)
%EM3D_NRCCIE_EVAL_GETQUAD  Near-quadrature for nrccie-eval.
if nargin < 4 || isempty(targinfo), targinfo = S; end
if nargin < 5 || isempty(opts),     opts     = struct(); end
opts.rep = 'nrccie-eval';
Q = em3d.pec.get_quadrature_correction(S, eps, zk, [], targinfo, opts);
end

function spmat = rsc_to_sparse_bc(Q, S)
%RSC_TO_SPARSE_BC  Convert nrccie-bc getquad RSC output to sparse.
% wnear(9,nquad) -> 3x3 block spmat.
spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear, ...
    kernel3d.rsc_interleave_full(3, 3));
end

function spmat = rsc_to_sparse_eval(Q, S, zk)
%RSC_TO_SPARSE_EVAL  Convert nrccie-eval getquad RSC output to sparse.
% wnear(4,nquad) -> 6x4 block spmat.
spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear, ...
    kernel3d.rsc_interleave_nrccie_eval(zk));
end



