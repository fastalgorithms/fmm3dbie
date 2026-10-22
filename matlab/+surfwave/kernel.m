function obj = kernel(green, type, zk, varargin)
%SURFWAVE.KERNEL   Radial kernel on a surface wave Green's function, for a
%   flat surfer in z = 0.
%
%   K = SURFWAVE.KERNEL(green, type, zk, ...)
%   K = SURFWAVE.KERNEL(green, type, zk, 'rmax', rmax, ...)
%
%   GREEN is a radial Green's function in the 2D format of
%   SURFWAVE.HELMGREEN, i.e. [val, grad, hess, third] = green(src, targ)
%   with (2,n) points, e.g. from SURFWAVE.CAPILLARY.GREEN. ZK is its
%   largest wavenumber (for oversampling). TYPE is a KERNEL3D.RADPP type, or
%
%      'lap'   - surface Laplacian, G_xx + G_yy
%
%   Without 'rmax' the kernel is a KERNEL3D.RADGREEN kernel: exact eval,
%   and a radpp fit for the near quadrature, built at getquad time. With
%   'rmax', rmax (or [rmin rmax]) it is a KERNEL3D.RADPP kernel: eval and
%   near quadrature both from one fit on [rmin rmax], which is much faster
%   to evaluate (e.g. for pcfft) but NaN outside [rmin rmax].
%
%   Further arguments go to KERNEL3D.RADGREEN or KERNEL3D.RADPP (as
%   name-value pairs, e.g. 'nu', nu or 'opts', struct('eps', 1e-13)); the
%   fit tolerance defaults to 1e-12.
%   Kernels of one type add under +, - and scalar *.
%
% See also KERNEL3D.CAPILLARY, KERNEL3D.ICEFLEX, KERNEL3D.GRAVITY

if ( isempty(type) ), type = 's'; end

rmax = [];
k = find(strcmpi(varargin(1:2:end), 'rmax'), 1);
if ( ~isempty(k) )
    rmax = varargin{2*k};
    varargin(2*k-1:2*k) = [];
end

% SURFWAVE.STRUVER switches to its asymptotic form at |z| = 95, which
% leaves ~1e-13 noise in the profiles, so the fits default to eps = 1e-12
if ( ~any(strcmpi(varargin(1:2:end), 'opts')) )
    varargin = [varargin, {'opts', struct('eps', 1e-12)}];
end
args = [{'dim', 2, 'zk', zk}, varargin];
if ( ischar(type) && any(strcmpi(type, {'lap', 'lap_s', 'laplacian'})) )
    type = @lapweights;
    args = [{'nderiv', 2, 'name', 'lap'}, args];
end

if ( isempty(rmax) )
    obj = kernel3d(green, type, args{:});
else
    obj = kernel3d.radpp(green, rmax, type, args{:});
end

end


function W = lapweights(t)
% weights on [G | G_x G_y | G_xx G_xy G_yy]
nt = size(t.r(:,:), 2);
W = zeros(nt, 6);
W(:, [4 6]) = 1;
end
