function [val,grad,hess] = green(zk,src,targ)
%HELM2D.GREEN  2D Helmholtz Green's function G(x,y) = (i/4) H_0^(1)(zk|x-y|).
%
%  [val,grad,hess] = helm2d.green(zk,src,targ)
%
%  src (2,ns) and targ (2,nt) are points in the plane (further rows are
%  ignored). val is (nt,ns); grad(:,:,1:2) and hess(:,:,1:3) hold the
%  target derivatives d/dx, d/dy and d2/dx2, d2/dxdy, d2/dy2.

switch nargout
    case {0,1}
        val = green2d.helm(zk,src,targ);
    case 2
        [val,grad] = green2d.helm(zk,src,targ);
    otherwise
        [val,grad,hess] = green2d.helm(zk,src,targ);
end

end
