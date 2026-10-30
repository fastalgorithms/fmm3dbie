function [val,grad,hess] = green(zk,src,targ)
%HELM2D.GREEN  2D Helmholtz Green's function G(x,y) = (i/4) H_0^(1)(zk|x-y|).
%
%  [val,grad,hess] = helm2d.green(zk,src,targ)
%
%  src (2,ns) and targ (2,nt) are points in the plane (further rows are
%  ignored). val is (nt,ns); grad(:,:,1:2) and hess(:,:,1:3) hold the
%  target derivatives d/dx, d/dy and d2/dx2, d2/dxdy, d2/dy2.
%
%  Adapted from chunkie (BSD 3-clause).

rx = targ(1,:).' - src(1,:);
ry = targ(2,:).' - src(2,:);
rx2 = rx.*rx;
ry2 = ry.*ry;
r2 = rx2 + ry2;
r  = sqrt(r2);

h0  = besselh(0,1,zk*r);
val = 0.25*1i*h0;

if nargout > 1
    h1 = besselh(1,1,zk*r);
    grad = zeros([size(r),2],'like',val);
    grad(:,:,1) = -1i*zk*0.25*h1.*rx./r;
    grad(:,:,2) = -1i*zk*0.25*h1.*ry./r;
end

if nargout > 2
    r3 = r.^3;
    h2 = 2*h1./(zk*r) - h0;
    hess = zeros([size(r),3],'like',val);
    hess(:,:,1) = 0.25*1i*zk*((rx-ry).*(rx+ry).*h1./r3 - zk*rx2.*h0./r2);
    hess(:,:,2) = 0.25*1i*zk*zk*rx.*ry.*h2./r2;
    hess(:,:,3) = 0.25*1i*zk*((ry-rx).*(rx+ry).*h1./r3 - zk*ry2.*h0./r2);
end

end
