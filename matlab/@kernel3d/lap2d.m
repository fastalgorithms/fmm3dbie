function obj = lap2d(type, coefs)
%KERNEL3D.LAP2D   Construct a Laplace volume potential kernel in 2D.
%
%   KERNEL3D.LAP2D(type) or KERNEL3D.LAP2D(type, coefs), for densities on
%   a flat surfer in the z = 0 plane, where type is:
%      's'  - single layer,   S(x,y) = -log|x-y|/(2*pi)
%      'sp' - sprime,         S'(x,y) = d/dn_x S(x,y)
%      'sg' - gradient,       [d/dx S; d/dy S], [2 1], rows interleaved
%      's2trans' - [2x1]      [coefs(1)*S; coefs(2)*S'], default [1;1]
%
%   The kernel is a KERNEL3D.RADGREEN kernel on LAP2D.GREEN, with an fmm
%   if rfmm2d is on the path.
%
% See also LAP2D.GREEN, KERNEL3D.RADGREEN

if ( nargin < 1 )
    error('KERNEL3D.LAP2D: missing Laplace kernel type.');
end
if ( nargin < 2 || isempty(coefs) ), coefs = [1; 1]; end

extra = {};
switch lower(type)
    case {'s', 'single', 'slp'}
        rtype = 's';
    case {'sp', 'sprime'}
        rtype = 'sp';
    case {'sg', 'sgrad'}
        rtype = 'grad';
    case {'s2trans', 's2t'}
        rtype = 's2trans';
        extra = {'coefs', coefs(:).'};
    otherwise
        error('KERNEL3D.LAP2D: unknown Laplace kernel type ''%s''.', type);
end

args = {'dim', 2};
if ( exist('rfmm2d', 'file') )
    args = [args, {'fmm', @(varargin) fmm2dto3d(@lap2d.fmm, varargin{:})}];
end
obj = kernel3d(@(s,t) lap2d.green(s, t), rtype, args{:}, extra{:});

end
