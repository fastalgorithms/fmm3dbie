function obj = helm2d(type, zk, coefs)
%KERNEL3D.HELM2D   Construct a Helmholtz volume potential kernel in 2D.
%
%   KERNEL3D.HELM2D(type, zk) or KERNEL3D.HELM2D(type, zk, coefs), for
%   densities on a flat surfer in the z = 0 plane, where zk is the
%   (complex) wavenumber and type is:
%      's'  - single layer,   S(x,y) = i/4 H_0^(1)(zk|x-y|)
%      'sp' - sprime,         S'(x,y) = d/dn_x S(x,y)
%      'sg' - gradient,       [d/dx S; d/dy S], [2 1], rows interleaved
%      's2trans' - [2x1]      [coefs(1)*S; coefs(2)*S'], default [1;1]
%
%   The kernel is a KERNEL3D.RADGREEN kernel on HELM2D.GREEN, with an fmm
%   if hfmm2d is on the path.
%
% See also HELM2D.GREEN, KERNEL3D.RADGREEN

if ( nargin < 2 )
    error('KERNEL3D.HELM2D: requires type and zk arguments.');
end
if ( nargin < 3 || isempty(coefs) ), coefs = [1; 1]; end

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
        error('KERNEL3D.HELM2D: unknown Helmholtz kernel type ''%s''.', type);
end

args = {'dim', 2, 'zk', zk};
if ( exist('hfmm2d', 'file') )
    fmm = @(eps,s,pg,t,pgt) hfmm2d(eps, zk, s, pg, t, pgt);
    args = [args, {'fmm', @(varargin) fmm2dto3d(fmm, varargin{:})}];
end
obj = kernel3d(@(s,t) helm2d.green(zk, s, t), rtype, args{:}, extra{:});

end
