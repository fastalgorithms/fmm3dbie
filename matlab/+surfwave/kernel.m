function obj = kernel(a, gam, n, w, type, varargin)
%SURFWAVE.KERNEL   Construct a surface wave kernel arising from the
% free surface condition on the velocity potential
%      (a (-Delta)^((n-1)/2) + gam) d phi/dz - phi = 0   on z = 0,
% with n = 1, 3 or 5.
%
% SURFWAVE.KERNEL(a,gam,n,w,type,'rmax',rmax,...)
%   w is a weight vector on [G_S, G_phi, S3d G_phi, Lap^2 G_phi] and
%   type is a KERNEL3D.RADPP type or 'lap'. The kernels are compressed
%   out to rmax. obj.zk is the propagating wavenumber.
%
%   Can provide 'depth', for a fluid of depth d. (default inf)
%
% See also KERNEL3D.CAPILLARY, KERNEL3D.ICEFLEX, KERNEL3D.GRAVITY

if ( isempty(type) ), type = 's'; end

[rmax, varargin] = getarg(varargin, 'rmax', []);
if ( isempty(rmax) )
    error('SURFWAVE.KERNEL: requires ''rmax''.');
end
[depth, varargin] = getarg(varargin, 'depth', []);
[fgw,   varargin] = getarg(varargin, 'fgw', []);
if ( ~isempty(depth) )
    fgw = surfwave.fgw_setup(a, gam, depth, n);
end

[g, zk] = green(a, gam, n, w, fgw);

% the Struve asymptotics leave ~1e-13 noise
if ( ~any(strcmpi(varargin(1:2:end), 'opts')) )
    varargin = [varargin, {'opts', struct('eps', 1e-12)}];
end
args = [{'dim', 2, 'zk', zk}, varargin];
if ( ischar(type) && any(strcmpi(type, {'lap', 'lap_s', 'laplacian'})) )
    type = @lapweights;
    args = [{'nderiv', 2, 'name', 'lap'}, args];
end

obj = kernel3d.radpp(g, rmax, type, args{:});

end


function [g, zk] = green(a, gam, n, w, fgw)
% one root sum, plus c3d/(4 pi r), powers c r^-m and the finite depth
% corrections; away from the source Lap^2 G_phi = (2/a)(S3d - gam/2 G_phi)
% for n = 5 and rho^4 G_phi + (9 rho/r^5 - rho^3/r^3)/(2 pi) for n = 1
[rts, ejs] = surfwave.find_roots(a, gam, n);
rts = rts(:);  ejs = ejs(:);
w(end+1:4) = 0;
if ( (n == 1 && w(1) ~= 0) || (n == 3 && w(4) ~= 0) )
    error('SURFWAVE.KERNEL: weights not available for n = %d.', n);
end
if ( ~isempty(fgw) && w(3) ~= 0 )
    error('SURFWAVE.KERNEL: S3d has no finite depth form.');
end
w2 = w(2);  w3 = w(3);  pw = zeros(0, 2);
if ( n == 5 )
    w2 = w2 - gam/a*w(4);
    w3 = w3 + 2/a*w(4);
elseif ( n == 1 && w(4) ~= 0 )
    rho = rts(1);
    w2 = w2 + rho^4*w(4);
    pw = [w(4)/(2*pi)*[9*rho; -rho^3], [5; 3]];
end
cjs = ejs.*(w(1)*rts.^2/2 + w2*rts/4 + w3/8);
c3d = w3;
if ( n == 1 ), c3d = c3d + 2*w2*sum(ejs); end     % sum(ejs) = 0 for n > 1

zk = real(rts(1));
add = {};
if ( ~isempty(fgw) )
    zk = fgw.k0*fgw.kinf;
    fw = [w(1), w(2), w(4)];
    names = {'gs', 'gphi', 'bilap_gphi'};
    add = [num2cell(fw(fw ~= 0)); names(fw ~= 0)];
end
g = @(s,t) combo(rts, cjs, c3d, pw, fgw, add, s, t);
end


function varargout = combo(rts, cjs, c3d, pw, fgw, add, s, t)
nout = max(nargout, 1);
varargout = cell(1, nout);
[varargout{:}] = surfwave.rootsum(rts, cjs, s, t, 1e-8, c3d);
for j = 1:size(pw, 1)
    d = cell(1, nout);
    [d{:}] = rpow(pw(j,2), s, t);
    for k = 1:nout, varargout{k} = varargout{k} + pw(j,1)*d{k}; end
end
for j = 1:size(add, 2)
    d = cell(1, nout);
    [d{:}] = surfwave.fgw_green(fgw, add{2,j}, s, t);
    for k = 1:nout, varargout{k} = varargout{k} + add{1,j}*d{k}; end
end
end


function [v, grad, hess] = rpow(m, s, t)
% r^-m and its derivatives
dx = t(1,:).' - s(1,:);  dy = t(2,:).' - s(2,:);
r2 = dx.^2 + dy.^2;
v = r2.^(-m/2);
if ( nargout > 1 )
    f1 = -m*v./r2;
    grad = cat(3, f1.*dx, f1.*dy);
end
if ( nargout > 2 )
    f2 = m*(m+2)*v./r2.^2;
    hess = cat(3, f1 + f2.*dx.^2, f2.*dx.*dy, f1 + f2.*dy.^2);
end
end


function W = lapweights(t)
% weights on [G | G_x G_y | G_xx G_xy G_yy]
nt = size(t.r(:,:), 2);
W = zeros(nt, 6);
W(:, [4 6]) = 1;
end


function [v, args] = getarg(args, name, v)
k = find(strcmpi(args(1:2:end), name), 1);
if ( ~isempty(k) )
    v = args{2*k};
    args(2*k-1:2*k) = [];
end
end
