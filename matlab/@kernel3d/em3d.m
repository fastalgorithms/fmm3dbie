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

        obj.eval = @(s,t) em3d_nrccie_bc_eval(zk, alpha, s, t);

        obj.fmm = [];   % no standalone FMM for this kernel

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

        obj.eval = @(s,t) em3d_nrccie_eval_kern(zk, s, t);

        obj.fmm = [];   % no standalone FMM for this kernel

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

function mat = em3d_nrccie_bc_eval(zk, alpha, srcinfo, targinfo)
%EM3D_NRCCIE_BC_EVAL  Evaluate 3x3 block NRCCIE system kernel pointwise.
%
%  Calls the Fortran kernel em3d_nrccie_kers (nd=9) for each (targ,src)
%  pair and assembles the (3*nt) x (3*ns) matrix.
%
%  Density layout: col index = 3*(j-1) + {1=ju, 2=jv, 3=rho}
%  Output  layout: row index = 3*(i-1) + {1=Mu, 2=Mv, 3=En}

src  = srcinfo.r;
targ = targinfo.r;

[~, ns] = size(src);
[~, nt] = size(targ);

% Pack srcinfo and targinfo into the 12-component arrays expected by
% em3d_nrccie_kers.  (Only the first 12 entries are used by the kernel.)
srcinfo12  = pack_ptinfo12(srcinfo);
targinfo12 = pack_ptinfo12(targinfo);

zpars = complex([zk; alpha]);
ndz = int64(2);
ndd = int64(0);
dpars = zeros(0,1);
ndi = int64(0);
ipars = zeros(0,1,'int64');
nd = int64(9);
ndt = int64(12);

mat = complex(zeros(3*nt, 3*ns));

for isrc = 1:ns
    for itarg = 1:nt
        val = complex(zeros(9,1));
        mex_id_ = 'em3d_nrccie_kers(c i int64_t[x], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c io dcomplex[x])';
        [val] = fmm3dbie_routs(mex_id_, nd, srcinfo12(:,isrc), ndt, targinfo12(:,itarg), ndd, dpars, ndz, zpars, ndi, ipars, val, 1, 12, 1, 12, 1, 1, 1, 2, 1, 1, 9);
        % val(1:9) -> [A11 A21 A31 A12 A22 A32 A13 A23 A33] (col-major 3x3)
        % Map: val = [row1_col1, row2_col1, row3_col1, row1_col2, ...]
        roff = 3*(itarg-1);
        coff = 3*(isrc-1);
        mat(roff+1, coff+1) = val(1);  % Mu from ju
        mat(roff+2, coff+1) = val(4);  % Mv from ju
        mat(roff+3, coff+1) = val(7);  % En from ju
        mat(roff+1, coff+2) = val(2);  % Mu from jv
        mat(roff+2, coff+2) = val(5);  % Mv from jv
        mat(roff+3, coff+2) = val(8);  % En from jv
        mat(roff+1, coff+3) = val(3);  % Mu from rho
        mat(roff+2, coff+3) = val(6);  % Mv from rho
        mat(roff+3, coff+3) = val(9);  % En from rho
    end
end
end

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

function mat = em3d_nrccie_eval_kern(zk, srcinfo, targinfo)
%EM3D_NRCCIE_EVAL_KERN  Evaluate (6*nt) x (4*ns) nrccie-eval kernel matrix.
%
%  Representation:
%    E = ik*S_k[J] - grad(S_k)[rho]
%    H = curl(S_k)[J]  =  (grad S_k) x J
%
%  Density column order: [J_x, J_y, J_z, rho]  (4 cols per source point)
%  Output row order:     [E_x, E_y, E_z, H_x, H_y, H_z]  (6 rows per target)
%
%  For a single (targ t, src s) pair the 6x4 block is:
%
%    [ik*G   0      0      -dGdx ]
%    [0      ik*G   0      -dGdy ]
%    [0      0      ik*G   -dGdz ]
%    [0      -dGdz  dGdy   0     ]   <- H_x from {J_x, J_y, J_z}
%    [dGdz   0      -dGdx  0     ]   <- H_y
%    [-dGdy  dGdx   0      0     ]   <- H_z
%
%  where G = S_k(t,s), dGdx = d/dx G(t,s)  (derivative w.r.t. target).

[~, ns] = size(srcinfo.r);
[~, nt] = size(targinfo.r);

[G, gradG] = helm3d.green(zk, srcinfo.r, targinfo.r);
% G     : (nt, ns)
% gradG : (nt, ns, 3)  -- gradient w.r.t. target position

% Build (6, nt, 4, ns) tensor then reshape to (6*nt, 4*ns)
T = complex(zeros(6, nt, 4, ns));

for isrc = 1:ns
    gx = gradG(:, isrc, 1);   % (nt,1) col vectors
    gy = gradG(:, isrc, 2);
    gz = gradG(:, isrc, 3);
    sk = G(:, isrc);

    % col 1: J_x component of density
    T(1,:,1,isrc) = 1i*zk*sk;   % E_x
    T(5,:,1,isrc) = gz;          % H_y = +dGdz * J_x
    T(6,:,1,isrc) = -gy;         % H_z = -dGdy * J_x

    % col 2: J_y
    T(2,:,2,isrc) = 1i*zk*sk;   % E_y
    T(4,:,2,isrc) = -gz;         % H_x = -dGdz * J_y
    T(6,:,2,isrc) = gx;          % H_z = +dGdx * J_y

    % col 3: J_z
    T(3,:,3,isrc) = 1i*zk*sk;   % E_z
    T(4,:,3,isrc) = gy;          % H_x = +dGdy * J_z
    T(5,:,3,isrc) = -gx;         % H_y = -dGdx * J_z

    % col 4: rho  -> E = -grad(S_k)[rho], H = 0
    T(1,:,4,isrc) = -gx;
    T(2,:,4,isrc) = -gy;
    T(3,:,4,isrc) = -gz;
end

% Permute to (nt, 6, ns, 4) then reshape to (6*nt, 4*ns)
% Row index: 6*(itarg-1) + {1..6}   Col index: 4*(isrc-1) + {1..4}
mat = reshape(permute(T, [2,1,4,3]), [6*nt, 4*ns]);
end

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

% =========================================================================
% Utility
% =========================================================================

function arr12 = pack_ptinfo12(ptinfo)
%PACK_PTINFO12  Pack a ptinfo struct into a (12, n) real array.
n   = size(ptinfo.r, 2);
r   = ptinfo.r;
du  = zeros(3, n);
dv  = zeros(3, n);
nrm = zeros(3, n);
if isfield(ptinfo, 'du'),  du  = ptinfo.du;  end
if isfield(ptinfo, 'dv'),  dv  = ptinfo.dv;  end
if isfield(ptinfo, 'n'),   nrm = ptinfo.n;   end
arr12 = [r; du; dv; nrm];
end
