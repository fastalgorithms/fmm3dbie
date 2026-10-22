function obj = gravity(type, g, green, varargin)
%KERNEL3D.GRAVITY   Gravity surface wave kernels, for a flat surfer in
%   z = 0.
%
%   K = KERNEL3D.GRAVITY(type, g)
%   K = KERNEL3D.GRAVITY(type, g, green)
%   K = KERNEL3D.GRAVITY(type, g, green, ...)
%
%   G is the gravity parameter (= 2*rho): the free surface Green's
%   function has one dispersion root rho = g/2 with residue 1, and
%   G_phi = G_S/g. GREEN is a weight vector w on [G_S, G_phi] (see
%   SURFWAVE.GRAVITY.GSGRAV), or 'gs' (default) or 'gphi'.
%
%   TYPE is any type of KERNEL3D.RADPP with 'dim' = 2 ('s', 'sp', 'grad',
%   'free_bc', ...), or 'lap', see SURFWAVE.KERNEL, or
%
%      'vol'   (a Lap^2 + b) G_S - g G_phi, with [a b] given as
%              'coefs', [a b]. This needs four derivatives, so it is
%              eval only (no near quadrature) and GREEN is ignored.
%
%   By default this is a KERNEL3D.RADGREEN kernel (exact eval); with
%   'rmax', rmax it is a KERNEL3D.RADPP kernel (eval from a fit on
%   [1e-12 rmax], fast), see SURFWAVE.KERNEL.
%
%   Further arguments go to KERNEL3D.RADGREEN. The plate types take the
%   Poisson ratio as 'nu', nu.
%
% See also SURFWAVE.GRAVITY.GREEN, SURFWAVE.KERNEL, KERNEL3D.CAPILLARY

if ( nargin < 2 || isempty(g) )
    error('KERNEL3D.GRAVITY: need type and the gravity parameter g.');
end
if ( nargin < 3 || isempty(green) ), green = 'gs'; end

if ( any(strcmpi(type, {'vol', 'volume'})) )
    obj = volkern(g, varargin{:});
    return
end

obj = surfwave.kernel(surfwave.gravity.green(g, green), type, abs(g/2), ...
                     varargin{:});

end


function obj = volkern(g, varargin)
% (a Lap^2 + b) G_S - g G_phi, eval only
k = find(strcmpi(varargin(1:2:end), 'coefs'), 1);
if ( isempty(k) )
    error('KERNEL3D.GRAVITY: ''vol'' needs ''coefs'', [a b].');
end
coefs = varargin{2*k};
obj              = kernel3d();
obj.name         = 'gravity';
obj.type         = 'vol';
obj.opdims       = [1 1];
obj.kernel_order = 1;
obj.zk           = abs(g/2);
obj.ifcomplex    = 1;
obj.params       = struct('g', g, 'coefs', coefs);
obj.eval         = @(s,t) surfwave.gravity.kern(g, s, t, 'vol', coefs);
obj.getquad      = [];
obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, ...
                           eps, obj.zk, obj.kernel_order);
end
