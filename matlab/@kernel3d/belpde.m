function obj = belpde(type, zk)
%KERNEL3D.BELPDE   Construct a Beltrami parametrix/remainder kernel in 3D.
%
%   KERNEL3D.BELPDE(type) or KERNEL3D.BELPDE(type, zk), where type is:
%      'klb' - Laplace-Beltrami parametrix,   K(x,y) = log|x-y|^2/(4*pi)
%      'rlb' - Laplace-Beltrami remainder,    R(x,y) = \Delta_\Gamma K(x,y)
%      'khb' - Helmholtz-Beltrami parametrix, K(x,y) = H_0(zk|x-y|)/(4i)
%      'rhb' - Helmholtz-Beltrami remainder,
%                                    R(x,y) = (\Delta_\Gamma + zk^2) K(x,y)
%
%   zk is the Helmholtz wavenumber, and is ignored for 'klb' and 'rlb'.
%
%   The remainder kernels depend on the mean curvature at the target, which
%   is read from targinfo.mean_curv.
%
% See also BELPDE.KERN, BELPDE.GET_QUADRATURE_CORRECTION, KERNEL3D

if ( nargin < 1 )
    error('KERNEL3D.BELPDE: missing Beltrami kernel type.');
end

if ( nargin < 2 )
    zk = 0;
end

islap = strcmpi(type, 'klb') || strcmpi(type, 'rlb');
if ( islap )
    zk = 0;
end

obj           = kernel3d();
obj.name      = 'belpde';
obj.opdims    = [1 1];
obj.zk        = zk;
obj.ifcomplex = double(~islap);

switch lower(type)

    case {'klb'}
        obj.type = 'klb';
        obj.kernel_order = -1;

    case {'rlb'}
        obj.type = 'rlb';
        obj.kernel_order = -1;
        obj.targ_fields = {'n', 'mean_curv'};

    case {'khb'}
        obj.type = 'khb';
        obj.kernel_order = -1;

    case {'rhb'}
        obj.type = 'rhb';
        obj.kernel_order = -1;
        obj.targ_fields = {'n', 'mean_curv'};

    otherwise
        error('KERNEL3D.BELPDE: unknown Beltrami kernel type ''%s''.', type);

end

% belpde.kern takes the wavenumber squared
zk2 = zk^2;
obj.eval = @(s,t) belpde.kern(zk2, s, t, obj.type);

obj.getquad = @(S,eps,varargin) rsc_to_sparse( ...
    belpde.get_quadrature_correction(S, obj.type, zk, eps, varargin{:}), S);
obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

end

function spmat = rsc_to_sparse(Q, S)
%RSC_TO_SPARSE  Convert belpde getquad RSC output to a sparse matrix.
% Beltrami kernels are scalar (opdims=[1,1]); no ri argument needed.
spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear);
end
