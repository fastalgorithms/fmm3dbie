function obj = gravity(type, g, green, varargin)
%KERNEL3D.GRAVITY   Construct a gravity surface wave kernel
% arising from the free surface condition on the velocity potential
%      (2/g) d phi/dz - phi = 0   on z = 0.
%
% KERNEL3D.GRAVITY(type,g,green,'rmax',rmax,...)
%   green is a weight vector on [G_S, G_phi] or one of 'gs' (default),
%   'gphi'. For type 's', obj.diag is the delta part of G_S.
%
%   Can provide 'depth', for a fluid of depth d. (default inf)
%
%   KERNEL3D.GRAVITY('floe',g,'coefs',[a b],'rmax',rmax) is the plate
%   condition (a Lap^2 + b) w - (g/2) phi applied to w = G_S, phi =
%   2 G_phi. obj.diag is g b.
%
% See also SURFWAVE.KERNEL, KERNEL3D.CAPILLARY, KERNEL3D.ICEFLEX

if ( nargin < 2 || isempty(g) )
    error('KERNEL3D.GRAVITY: requires type and g arguments.');
end

if ( any(strcmpi(type, {'floe', 'plate'})) )
    if ( nargin > 2 ), varargin = [{green}, varargin]; end
    [coefs, varargin] = getarg(varargin, 'coefs', []);
    if ( isempty(coefs) )
        error('KERNEL3D.GRAVITY: ''floe'' requires ''coefs'', [a b].');
    end
    obj = surfwave.kernel(2/g, 0, 1, g*[0, coefs(2) - 1, 0, coefs(1)], ...
                          's', varargin{:});
    obj.name = 'gravity';
    obj.type = 'floe';
    obj.diag = @(t) g*coefs(2)*ones(size(t.r(:,:), 2), 1);
    return
end

% G_S = g (delta + G_phi), G_phi of SURFWAVE.KERNEL with a = 2/g, n = 1
if ( nargin < 3 || isempty(green) ), green = 'gs'; end
if ( ischar(green) ), green = double(strcmpi(green, {'gs', 'gphi'})); end
green(end+1:2) = 0;
obj = surfwave.kernel(2/g, 0, 1, [0, g*green(1) + green(2)], type, ...
                      varargin{:});
if ( green(1) ~= 0 && ischar(type) && any(strcmpi(type, {'s', 'single', 'slp'})) )
    obj.diag = @(t) g*green(1)*ones(size(t.r(:,:), 2), 1);
end

end


function [v, args] = getarg(args, name, v)
k = find(strcmpi(args(1:2:end), name), 1);
if ( ~isempty(k) )
    v = args{2*k};
    args(2*k-1:2*k) = [];
end
end
