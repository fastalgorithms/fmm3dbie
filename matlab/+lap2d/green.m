function [val,grad,hess] = green(src,targ)
%LAP2D.GREEN  2D Laplace Green's function G(x,y) = -log|x-y|/(2 pi).
%
%  [val,grad,hess] = lap2d.green(src,targ)
%
%  src (2,ns) and targ (2,nt) are points in the plane (further rows are
%  ignored). val is (nt,ns); grad(:,:,1:2) and hess(:,:,1:3) hold the
%  target derivatives d/dx, d/dy and d2/dx2, d2/dxdy, d2/dy2.
%
%  Adapted from chunkie (BSD 3-clause).

rx = targ(1,:).' - src(1,:);
ry = targ(2,:).' - src(2,:);
r2 = rx.^2 + ry.^2;

val = -log(r2)/(4*pi);

if nargout > 1
    grad = zeros([size(r2),2]);
    grad(:,:,1) = -rx./(2*pi*r2);
    grad(:,:,2) = -ry./(2*pi*r2);
end

if nargout > 2
    r4 = r2.*r2;
    hess = zeros([size(r2),3]);
    hess(:,:,1) = rx.^2./(pi*r4) - 1./(2*pi*r2);
    hess(:,:,2) = rx.*ry./(pi*r4);
    hess(:,:,3) = ry.^2./(pi*r4) - 1./(2*pi*r2);
end

end
