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
%   The remainder kernels carry the identity term of the second kind
%   equation in obj.diag, so surfermat(S, kernel3d('bel','rlb')) builds
%   I + R directly.
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
        obj.diag = @(t) ones(size(t.r(:,:),2), 1);

    case {'khb'}
        obj.type = 'khb';
        obj.kernel_order = -1;

    case {'rhb'}
        obj.type = 'rhb';
        obj.kernel_order = -1;
        obj.targ_fields = {'n', 'mean_curv'};
        obj.diag = @(t) ones(size(t.r(:,:),2), 1);

    otherwise
        error('KERNEL3D.BELPDE: unknown Beltrami kernel type ''%s''.', type);

end

% belpde.kern takes the wavenumber squared
zk2 = zk^2;
obj.eval = @(s,t) belpde.kern(zk2, s, t, obj.type);

% registry names belpde_klb, belpde_rlb, belpde_khb, belpde_rhb
qname = ['belpde_', lower(obj.type)];
obj.getquad = @(S,eps,varargin) surf_quad(S, eps, qname, zk, varargin{:});
obj.getquad_sub = true;
obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

end

function spmat = surf_quad(S, eps, name, zk, targinfo, opts)
%SURF_QUAD  getquad via GETNEARQUAD_MEX; targets need the mean curvature.
if nargin < 5 || isempty(targinfo), targinfo = S; end
if nargin < 6 || isempty(opts), opts = struct(); end

if isa(targinfo, 'surfer')
    targinfo = struct('r', targinfo.r, 'du', targinfo.du, ...
        'dv', targinfo.dv, 'n', targinfo.n, ...
        'patch_id', targinfo.patch_id, 'uvs_targ', targinfo.uvs_targ, ...
        'kappa', targinfo.mean_curv);
end
if ~isfield(targinfo, 'kappa')
    if isfield(targinfo, 'mean_curv')
        targinfo.kappa = targinfo.mean_curv;
        targinfo = rmfield(targinfo, 'mean_curv');
    else
        error(['KERNEL3D.BELPDE: mean curvature at the targets must be ' ...
            'supplied in targinfo.kappa or targinfo.mean_curv.']);
    end
end
targinfo.kappa = reshape(targinfo.kappa, 1, []);

opts.format = 'sparse';
Q = getnearquad_mex(S, name, eps, struct('zpars', zk), targinfo, opts);
spmat = Q.spmat;
end
