function varargout = lap3dgreen(src, targ)
%SURFWAVE.LAP3DGREEN   the Laplace single layer 1/(4 pi r) between points
%   in the plane, in the 2D format of SURFWAVE.HELMGREEN:
%   [val, grad, hess, third, fourth], with target derivatives in the
%   plane packed as G_x, G_y | G_xx, G_xy, G_yy | G_xxx, ... | G_xxxx, ...
if ( isstruct(src) ),  src  = src.r;  end
if ( isstruct(targ) ), targ = targ.r; end
x = targ(1,:).' - src(1,:);
y = targ(2,:).' - src(2,:);
r = sqrt(x.^2 + y.^2);
c = 1./(4*pi*r);
nout = max(nargout, 1);
varargout = cell(1, nout);
varargout{1} = c;
if ( nout > 1 ), varargout{2} = -c./r.^2.*cat(3, x, y); end
if ( nout > 2 )
    varargout{3} = c./r.^4.*cat(3, 2*x.^2 - y.^2, 3*x.*y, 2*y.^2 - x.^2);
end
if ( nout > 3 )
    varargout{4} = c./r.^6.*cat(3, -3*x.*(2*x.^2 - 3*y.^2), ...
        -3*y.*(4*x.^2 - y.^2), 3*x.*(x.^2 - 4*y.^2), 3*y.*(3*x.^2 - 2*y.^2));
end
if ( nout > 4 )
    varargout{5} = c./r.^8.*cat(3, 3*(8*x.^4 - 24*x.^2.*y.^2 + 3*y.^4), ...
        15*x.*y.*(4*x.^2 - 3*y.^2), -3*(4*x.^4 - 27*x.^2.*y.^2 + 4*y.^4), ...
        -15*x.*y.*(3*x.^2 - 4*y.^2), 3*(3*x.^4 - 24*x.^2.*y.^2 + 8*y.^4));
end
end
