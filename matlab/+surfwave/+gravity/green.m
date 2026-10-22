function g = green(gpar, w)
%SURFWAVE.GRAVITY.GREEN   gravity Green's function handle.
%
%   g = SURFWAVE.GRAVITY.GREEN(gpar, w)
%
%   returns [val, grad, hess, third, fourth] = g(src, targ) for
%   w(1) G_S + w(2) G_phi, with gravity parameter gpar (= 2*rho) and
%   G_phi = G_S/gpar. With ej = w(1) + w(2)/gpar,
%
%      G = ej rho/(2 pi r) + ej rho^2/4 K(rho r),
%
%   K the radial profile of SURFWAVE.ROOTSUM. W may also be 'gs'
%   (default) or 'gphi'.
%
% See also KERNEL3D.GRAVITY, SURFWAVE.ROOTSUM

if ( nargin < 2 || isempty(w) ), w = 'gs'; end
names = {'gs', 'gphi'};
if ( ischar(w) ), w = double(strcmpi(w, names)); end
w(end+1:numel(names)) = 0;
rho = gpar/2;
ej = w(1) + w(2)/gpar;
g = @(s,t) surfwave.gravity.gsgrav_der(rho, ej, s, t);
end
