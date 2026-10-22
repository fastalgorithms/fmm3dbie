function g = green(rts, ejs, w)
%SURFWAVE.CAPILLARY.GREEN   capillary Green's function handle.
%
%   g = SURFWAVE.CAPILLARY.GREEN(rts, ejs, w)
%
%   returns [val, grad, hess] = g(src, targ) for
%
%      w(1) G_S + w(2) G_phi + w(3) S3d
%
%   with G_S, G_phi the capillary Green's functions and S3d = S3d G_phi +
%   the Laplace single layer 1/(4 pi r). All three are sums over the
%   dispersion roots of one radial profile with different weights, so the
%   combination is one sum (SURFWAVE.ROOTSUM), with weights
%
%      ejs.*(w(1) rts.^2/2 + w(2) rts/4 + w(3)/8).
%
%   Trailing zeros of w may be left out. W may also be 'gs' (default),
%   'gphi' or 's3d'.
%
% See also KERNEL3D.CAPILLARY, SURFWAVE.ROOTSUM

if ( nargin < 3 || isempty(w) ), w = 'gs'; end
names = {'gs', 'gphi', 's3d'};
if ( ischar(w) ), w = double(strcmpi(w, names)); end
w(end+1:numel(names)) = 0;
rts = rts(:);  ejs = ejs(:);
cjs = ejs.*(w(1)*rts.^2/2 + w(2)*rts/4 + w(3)/8);
g = @(s,t) combo(rts, cjs, w(3), s, t);
end


function varargout = combo(rts, cjs, w3, s, t)
nout = max(nargout, 1);
[varargout{1:nout}] = surfwave.rootsum(rts, cjs, s, t);
if ( w3 ~= 0 )
    lap = cell(1, nout);
    [lap{:}] = surfwave.lap3dgreen(s, t);
    for k = 1:nout, varargout{k} = varargout{k} + w3*lap{k}; end
end
end
