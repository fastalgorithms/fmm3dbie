function obj = capillary(type, beta, gamma, green, varargin)
%KERNEL3D.CAPILLARY   Construct a capillary surface wave kernel
% arising from the free surface condition on the velocity potential
%      (-beta Delta + gamma) d phi/dz - phi = 0   on z = 0.
%
% KERNEL3D.CAPILLARY(type,beta,gamma,green,'rmax',rmax,...)
%   green is a weight vector on [G_S, G_phi, S3d G_phi] or one of
%   'gs' (default), 'gphi', 's3d'.
%
%   Can provide 'depth', for a fluid of depth d. (default inf)
%
%   KERNEL3D.CAPILLARY('cap_varbeta',beta,gamma,'pfun',pfun,'rmax',rmax)
%   is the kernel K of the surface with variable tension b,
%      -div(b grad w) + gamma w - phi = 0,   w = d phi/dz,
%   for w = G_S[sigma], phi = 2 G_phi[sigma]. The kernel is scaled by
%   beta/2. pfun(t) returns beta (1,nt) and dbeta (2,nt).
%
% See also SURFWAVE.KERNEL, KERNEL3D.ICEFLEX, KERNEL3D.GRAVITY

if ( nargin < 3 )
    error('KERNEL3D.CAPILLARY: requires type, beta and gamma arguments.');
end

if ( any(strcmpi(type, {'cap_varbeta', 'varbeta', 'varcoef'})) )
    if ( nargin > 3 ), varargin = [{green}, varargin]; end
    [pfun, varargin] = getarg(varargin, 'pfun', []);
    if ( isempty(pfun) )
        error('KERNEL3D.CAPILLARY: ''cap_varbeta'' requires pfun.');
    end
    [depth, varargin] = getarg(varargin, 'depth', []);
    if ( ~isempty(depth) )
        varargin = [varargin, {'fgw', surfwave.fgw_setup(beta, gamma, depth, 3)}];
    end
    K1 = surfwave.kernel(beta, gamma, 3, [1 0], ...
             @(t) wvar(t, beta, pfun), 'nderiv', 1, varargin{:});
    K2 = surfwave.kernel(beta, gamma, 3, [-gamma/2 1], ...
             @(t) reshape(pfun(t).beta, [], 1) - beta, 'nderiv', 0, varargin{:});
    obj = K1 + K2;
    obj.name = 'capillary';
    obj.type = 'cap_varbeta';
    obj.zk   = K1.zk;
    obj.diag = @(t) reshape(pfun(t).beta, [], 1);
    return
end

if ( nargin < 4 || isempty(green) ), green = 'gs'; end
if ( ischar(green) )
    green = double(strcmpi(green, {'gs', 'gphi', 's3d'}));
end
obj = surfwave.kernel(beta, gamma, 3, green, type, varargin{:});

end


function W = wvar(t, be, pfun)
% beta/2 times the weights on [G | G_x G_y]
p = pfun(t);
db = reshape(p.dbeta, 2, []).';
W = -be/2*[zeros(size(db, 1), 1), db];
end


function [v, args] = getarg(args, name, v)
k = find(strcmpi(args(1:2:end), name), 1);
if ( ~isempty(k) )
    v = args{2*k};
    args(2*k-1:2*k) = [];
end
end
