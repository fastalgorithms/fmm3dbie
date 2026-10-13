function [val,grad,hess] = green(src,targ,nolog)
%LAP2D.GREEN  2D Laplace Green's function G(x,y) = -log|x-y|/(2 pi).
%
%  [val,grad,hess] = lap2d.green(src,targ)
%  [val,grad,hess] = lap2d.green(src,targ,nolog)
%
%  src (2,ns) and targ (2,nt) are points in the plane (further rows are
%  ignored). val is (nt,ns); grad(:,:,1:2) and hess(:,:,1:3) hold the
%  target derivatives d/dx, d/dy and d2/dx2, d2/dxdy, d2/dy2. If nolog
%  is true, val is returned empty.

if nargin < 3, nolog = false; end

switch nargout
    case {0,1}
        val = green2d.lap(src,targ,nolog);
    case 2
        [val,grad] = green2d.lap(src,targ,nolog);
    otherwise
        [val,grad,hess] = green2d.lap(src,targ,nolog);
end

end
