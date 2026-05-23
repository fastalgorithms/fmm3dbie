classdef kernel3d
%KERNEL3D   Class describing an integral kernel in 3D, usually
% related to the solution of a PDE on a surface.
%
%   K = KERNEL3D(NAME, TYPE) constructs a kernel of the specified name
%   and type. The currently supported kernel names and types are:
%
%      NAME                              TYPE
%      ----                              ----
%      'laplace'   ('lap', 'l')          's', 'd', 'sp', 'c'
%      'helmholtz' ('helm', 'h')         's', 'd', 'sp', 'dprime',
%                                        'c', 'cprime', 'trans'
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
%         and targets.  srcinfo / targinfo are structs with fields
%            .r   - (3,:) positions
%            .n   - (3,:) unit normals
%            .du  - (3,:) first u-derivative (if needed)
%            .dv  - (3,:) first v-derivative (if needed)
%         Returns an (opdims(1)*nt) x (opdims(2)*ns) matrix.
%
%      K.getquads(S, eps, dpars, targinfo, opts)
%         Function handle.  Returns the near-quadrature struct
%         Q used by the layer-potential evaluator.
%
%      K.layer_eval(S, sigma, targinfo, eps, dpars)
%         Function handle.  Evaluates the layer potential via the
%         FMM-accelerated evaluator.
%
%      K.fmm(eps, s, t, sigma)
%         Function handle.  Evaluates the FMM for the corresponding
%         kernel with density sigma from sources s to targets t with
%         accuracy eps.
%
%      K.sing      - kernel order: -1 (weakly singular), 0 (pv), 1 (hypersingular)
%      K.zk        - wavenumber
%      K.ifcomplex - 0 if kernel is real-valued, 1 if complex-valued
%      K.opdims    - [m n], operator dimensions (scalar kernels: [1 1])
%      K.params    - struct of kernel parameters (e.g. coefs)
%
    properties

        name           % Name of the kernel
        type           % Type of the kernel
        params         % Structure of kernel parameters
        eval           % Function handle for kernel evaluation
        getquad        % Function handle to get quadrature corrections
        fmm              % Function handle for raw FMM call
        layer_eval       % Function handle for FMM + quadrature layer-potential eval
        get_overs_orders % Function handle to get oversampling orders
        sing           % Kernel order: -1 (weakly singular), 0 (pv), 1 (hypersingular)
        zk      = 0    % Wavenumber (0 for Laplace)
        ifcomplex = 0  % 0 = real-valued kernel, 1 = complex-valued kernel
        opdims = [0 0] % Dimension of the operator [m n]
        isnan  = false % Boolean, true for NaN kernels
        iszero = false % Boolean, true for zero kernels

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
                    otherwise
                        error('KERNEL3D: kernel ''%s'' not found.', kern);
                end
            elseif ( isa(kern, 'function_handle') )
                obj.eval = kern;
                try
                    s = []; s.r = randn(3,1); s.n = randn(3,1);
                    s.du = randn(3,1); s.dv = randn(3,1);
                    t = []; t.r = randn(3,1); t.n = randn(3,1);
                    t.du = randn(3,1); t.dv = randn(3,1);
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

    end

    methods ( Static )

        obj    = lap3d(varargin);
        obj    = helm3d(varargin);
        K      = interleave(kerns);
        novers = kernel3d_getnear_overs(S,t,eps,zk,sing);
        Q      = addquad(Qf,Qg,S,sign);
        Q      = scalequad(Qf,S,c);

    end

end
