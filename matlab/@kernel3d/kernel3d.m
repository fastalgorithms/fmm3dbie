classdef kernel3d
%KERNEL3D   Class describing an integral kernel in 3D, usually
% related to the solution of a PDE on a surface.
%
%   K = KERNEL3D(NAME, TYPE) constructs a kernel of the specified name
%   and type. The currently supported kernel names and types are:
%
%      NAME                              TYPE
%      ----                              ----
%      'laplace'   ('lap', 'l')          's', 'd', 'sp', 'dp', 'c', 'cp'
%      'helmholtz' ('helm', 'h')         's', 'd', 'sp', 'dprime', 'c', 'cprime',
%                                        's_diff', 'd_diff', 'sp_diff', 'dp_diff',
%                                        'c_diff', 'cp_diff',
%                                        'trans_sys' ([2x2], single zk),
%                                        'trans_sys_diff' ([2x2], two zks),
%                                        'trans_rep' ([1x2]),
%                                        's2trans', 'd2trans', 'c2trans',
%                                        's2trans_diff', 'd2trans_diff', 'c2trans_diff'
%      'maxwell'   ('em3d', 'em')        'nrccie-bc', 'nrccie-eval'
%      'stokes'    ('stok3d', 'stok')    's', 'd', 'c'
%      'zero'/'zeros'                    (no type needed)
%
%   Some kernels accept extra parameters as trailing arguments, e.g.
%   the combined-layer Laplace kernel:
%
%      K = KERNEL3D('laplace', 'c', [alpha, beta])
%
%   The kernel class stores the following data:
%
%      K.eval(srcinfo, targinfo)
%         Function handle.  Evaluates the kernel matrix between sources
%         and targets. Returns an (opdims(1)*nt) x (opdims(2)*ns) matrix.
%
%      K.getquads(S, eps, dpars, targinfo, opts)
%         Function handle.  Returns the near-quadrature struct Q.
%
%      K.layer_eval(S, sigma, targinfo, eps, dpars)
%         Function handle.  Evaluates the layer potential via the
%         FMM-accelerated Fortran evaluator.
%
%      K.fmm(eps, s, t, sigma)
%         Function handle.  K.fmm(eps, s, t, sigma) evaluates the FMM for
%         the corresponding kernel with density sigma from sources s to
%         targets t with accuracy eps.
%
%      K.kernel_order - kernel order, used to get oversampling orders. 
%         Options are  -1 (single layer), 0 (double layer), 1 (derivative
%         of double layer).
%      K.zk        - wavenumber
%      K.ifcomplex - 0 if kernel is real-valued, 1 if complex-valued
%      K.opdims    - [m n], operator dimensions (scalar kernels: [1 1])
%      K.params    - struct of kernel parameters (e.g. coefs)
%
%      K.src_fields  - source fields expected by K.eval. ('r' is assumed)
%      K.targ_fields - target fields expected by K.eval. ('r' is assumed)
%
    properties

        name           % Name of the kernel
        type           % Type of the kernel
        eval           % Function handle for kernel evaluation
        getquad        % Function handle to get quadrature corrections
        fmm            % Function handle for FMM call
        kernel_order   % Kernel order: -1 (single layer), 0 (double layer), 1 (derivative of double layer)
        opdims = [0 0] % Dimension of the operator [m n]
        src_fields     % Ptinfo fields required at sources, i.e. {'n'} for double layer
                       % the location 'r' is implicitly added
        targ_fields    % Ptinfo fields required at targets
                       % the location 'r' is implicitly added

        rsc_to_interleave  % Struct describing how to map wnear(nker,nquad) to the
                           % (opdims(1)*nt, opdims(2)*ns) sparse matrix.
                           %
                           % All kernels use the same canonical form (see rsc_interleave_full):
                           %   .m       - block row dimension (opdims(1))
                           %   .n       - block col dimension (opdims(2))
                           %   .nker    - number of rows in wnear
                           %   .entries - struct array, one element per (wnear row, block position)
                           %              pair.  Each element has:
                           %                .ker_id  1-based index into wnear rows
                           %                .row_id  1-based block row
                           %                .col_id  1-based block col
                           %                .coef    scalar coefficient (may be complex)
                           %
                           %   Sparse assembly:  B(row_id, col_id) += coef * wnear(ker_id, :)
                           %   Inversion:        wnear(ker_id, :)   = B(row_id, col_id) / coef
                           %                    (using the first nonzero-coef entry per ker_id)

    end

    properties (Hidden = true)
        zk      = 0    % Wavenumber for oversampling purposes (largest wavenumber in the probelm)
        ifcomplex = 0  % 0 = real-valued kernel, 1 = complex-valued kernel
        layer_eval       % Function handle for FMM + quadrature layer-potential eval
        get_overs_orders % Function handle to get oversampling orders
        isnan  = false % Boolean, true for NaN kernels
        iszero = false % Boolean, true for zero kernels
        params         % Structure of kernel parameters
    end
    methods

        function obj = kernel3d(kern, varargin)

            if ( nargin < 1 )
                return
            end

            if ( isa(kern, 'string') || isa(kern, 'char') )
                switch lower(kern)
                    case {'laplace', 'lap', 'l'}
                        obj = kernel3d.lap3d(varargin{:});
                    case {'helmholtz', 'helm', 'h'}
                        obj = kernel3d.helm3d(varargin{:});
                    case {'maxwell', 'em3d', 'em'}
                        obj = kernel3d.em3d(varargin{:});
                    case {'stokes', 'stok3d', 'stok'}
                        obj = kernel3d.stok3d(varargin{:});
                    case {'z', 'zero', 'zeros'}
                        if ~isempty(varargin)
                            obj = kernel3d.zeros(varargin{1});
                        else
                            obj = kernel3d.zeros();
                        end
                    otherwise
                        error('KERNEL3D: kernel ''%s'' not found.', kern);
                end
            elseif ( isa(kern, 'function_handle') )
                obj.eval = kern;
                try
                    s = []; s.r = randn(3,1); s.n = randn(3,1);
                    s.du = randn(3,1); s.dv = randn(3,1);
                    s.d = randn(3,1); s.d2 = randn(3,1);
                    t = []; t.r = randn(3,1); t.n = randn(3,1);
                    t.du = randn(3,1); t.dv = randn(3,1);
                    t.d = randn(3,1); t.d2 = randn(3,1);
                    obj.opdims = size(kern(s,t));
                catch
                    % unable to determine opdims automatically
                end
            elseif ( isa(kern, 'kernel3d') )
                if ( numel(kern) == 1 )
                    obj = kern;
                else
                    obj = kernel3d.interleave(kern);
                end
            else
                error('KERNEL3D: first input not of a supported type.');
            end

        end

        val    = eval_mask(obj,src,targ);

    end

    methods ( Static )

        obj    = lap3d(varargin);
        obj    = helm3d(varargin);
        obj    = em3d(varargin);
        obj    = stok3d(varargin);
        obj    = zeros(opdims);
        K      = interleave(kerns);
        novers = kernel3d_getnear_overs(S,t,eps,zk,sing);
        Q      = addquad(Qf,Qg,S,sign,ri);
        Q      = scalequad(Qf,S,c,ri);
        ri     = rsc_interleave_full(m, n, nker, entries);
        ri     = rsc_interleave_symmetric3();
        ri     = rsc_interleave_nrccie_eval(zk);
        obj    = tangent_kern(kernin, ids, jds);

    end

end
