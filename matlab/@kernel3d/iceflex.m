function obj = iceflex(type, alpha, gamma, green, varargin)
%KERNEL3D.ICEFLEX   Construct a flexural (ice plate) surface wave kernel
% arising from the free surface condition on the velocity potential
%      (alpha Delta^2 + gamma) d phi/dz - phi = 0   on z = 0.
%
% KERNEL3D.ICEFLEX(type,alpha,gamma,green,'rmax',rmax,...) 
%   green is a weight vector on [G_S, G_phi, S3d G_phi,
%   Lap^2 G_phi] or one of 'gs' (default), 'gphi', 's3d', 'gphi_bilap'.
%
%   Can provide 'depth', for a fluid of depth d. (default inf)
%
%   KERNEL3D.ICEFLEX('flex_varalpha',alpha,gamma,'nu',nu,'pfun',pfun,
%   'rmax',rmax) is the kernel K of the plate with variable rigidity a,
%      Delta(a Delta w) - (1-nu)(a_xx w_yy - 2 a_xy w_xy + a_yy w_xx)
%         + gamma w - phi = 0,   w = d phi/dz,
%   for w = G_S[sigma], phi = 2 G_phi[sigma]. The kernel is scaled by
%   alpha/2 (see KERNEL3D.FLEX2D).
%
% See also SURFWAVE.KERNEL, KERNEL3D.CAPILLARY, KERNEL3D.GRAVITY

if ( nargin < 3 )
    error('KERNEL3D.ICEFLEX: requires type, alpha and gamma arguments.');
end

if ( any(strcmpi(type, {'flex_varalpha', 'varalpha', 'varcoef'})) )
    if ( nargin > 3 ), varargin = [{green}, varargin]; end
    [nu,   varargin] = getarg(varargin, 'nu', []);
    [pfun, varargin] = getarg(varargin, 'pfun', []);
    if ( isempty(nu) || isempty(pfun) )
        error('KERNEL3D.ICEFLEX: ''flex_varalpha'' requires nu and pfun.');
    end
    [depth, varargin] = getarg(varargin, 'depth', []);
    if ( ~isempty(depth) )
        varargin = [varargin, {'fgw', surfwave.fgw_setup(alpha, gamma, depth, 5)}];
    end
    [depth, varargin] = getarg(varargin, 'depth', []);
    if ( ~isempty(depth) )
        varargin = [varargin, {'fgw', surfwave.fgw_setup(alpha, gamma, depth, 5)}];
    end
    K1 = surfwave.kernel(alpha, gamma, 5, [1 0], ...
             @(t) wvar(t, alpha, nu, pfun), 'nderiv', 3, varargin{:});
    K2 = surfwave.kernel(alpha, gamma, 5, [-gamma/2 1], ...
             @(t) reshape(pfun(t).alpha, [], 1) - alpha, 'nderiv', 0, varargin{:});
    obj = K1 + K2;
    obj.name = 'iceflex';
    obj.type = 'flex_varalpha';
    obj.zk   = K1.zk;
    obj.diag = @(t) reshape(pfun(t).alpha, [], 1);
    return
end

if ( nargin < 4 || isempty(green) ), green = 'gs'; end
if ( ischar(green) )
    green = double(strcmpi(green, {'gs', 'gphi', 's3d', 'gphi_bilap'}));
end
obj = surfwave.kernel(alpha, gamma, 5, green, type, varargin{:});

end


function W = wvar(t, al, nu, pfun)
% alpha/2 times the weights on [G | G_x G_y | G_xx G_xy G_yy |
% G_xxx G_xxy G_xyy G_yyy]
p = pfun(t);
nt = numel(p.alpha);
da = reshape(p.dalpha, 2, nt).';  d2a = reshape(p.d2alpha, 3, nt).';
la = d2a(:,1) + d2a(:,3);
z = zeros(nt, 1);
W = al/2*[z, z, z, la - (1-nu)*d2a(:,3), 2*(1-nu)*d2a(:,2), ...
          la - (1-nu)*d2a(:,1), 2*da(:,1), 2*da(:,2), 2*da(:,1), 2*da(:,2)];
end


function [v, args] = getarg(args, name, v)
k = find(strcmpi(args(1:2:end), name), 1);
if ( ~isempty(k) )
    v = args{2*k};
    args(2*k-1:2*k) = [];
end
end
