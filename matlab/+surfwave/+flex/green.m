function g = green(rts, ejs, w, add)
%SURFWAVE.FLEX.GREEN   flexural (ice plate) Green's function handle.
%
%   g = SURFWAVE.FLEX.GREEN(rts, ejs, w)
%   g = SURFWAVE.FLEX.GREEN(rts, ejs, w, add)
%
%   returns [val, grad, hess, third] = g(src, targ) for
%
%      w(1) G_S + w(2) G_phi + w(3) S3d + w(4) Lap^2 G_phi
%
%   with G_S, G_phi the flexural Green's functions, S3d = S3d G_phi + the
%   Laplace single layer 1/(4 pi r) and the bilaplacian of G_phi (values
%   only). The first three are sums over the dispersion roots of one
%   radial profile with different weights, so the combination is one sum
%   (SURFWAVE.ROOTSUM), with weights
%
%      ejs.*(w(1) rts.^2/2 + w(2) rts/4 + w(3)/8).
%
%   Trailing zeros of w may be left out. W may also be 'gs' (default),
%   'gphi', 's3d' or 'gphi_bilap'.
%
%   ADD = {G_S, G_phi, S3d, Lap^2 G_phi} are smooth corrections added to
%   the four Green's functions with the same weights w, e.g. finite depth
%   minus deep water. Entries are [] or handles in the format of G.
%
% See also KERNEL3D.ICEFLEX, SURFWAVE.ROOTSUM

if ( nargin < 3 || isempty(w) ), w = 'gs'; end
if ( nargin < 4 ), add = {}; end
names = {'gs', 'gphi', 's3d', 'gphi_bilap'};
if ( ischar(w) ), w = double(strcmpi(w, names)); end
w(end+1:numel(names)) = 0;
add(end+1:numel(names)) = {[]};
rts = rts(:);  ejs = ejs(:);
cjs = ejs.*(w(1)*rts.^2/2 + w(2)*rts/4 + w(3)/8);
g = @(s,t) combo(rts, ejs, cjs, w, add, s, t);
end


function varargout = combo(rts, ejs, cjs, w, add, s, t)
nout = max(nargout, 1);
if ( w(4) ~= 0 && nout > 1 )
    error('SURFWAVE.FLEX.GREEN: the bilaplacian of G_phi has values only.');
end
varargout = cell(1, nout);
if ( any(w(1:3)) )
    [varargout{:}] = surfwave.rootsum(rts, cjs, s, t, 1e-8);
else
    [varargout{:}] = deal(0);
end
if ( w(3) ~= 0 )
    lap = cell(1, nout);
    [lap{:}] = surfwave.lap3dgreen(s, t);
    for k = 1:nout, varargout{k} = varargout{k} + w(3)*lap{k}; end
end
if ( w(4) ~= 0 )
    [~, ~, ~, ~, bl] = surfwave.flex.gphiflex(rts, ejs, s, t);
    varargout{1} = varargout{1} + w(4)*bl;
end
for j = find(w ~= 0 & ~cellfun(@isempty, add))
    d = cell(1, nout);
    [d{:}] = add{j}(s, t);
    for k = 1:nout, varargout{k} = varargout{k} + w(j)*d{k}; end
end
end
