function [val,grad,hess] = lap(src,targ,nolog)
%GREEN2D.LAP  2D Laplace Green's function -log|x-y|/(2 pi).
%
%  [val,grad,hess] = green2d.lap(src,targ)
%  [val,grad,hess] = green2d.lap(src,targ,nolog)
%
%  Same layout as GREEN2D.HELM. If nolog is true, val is returned empty.
%
%  Adapted from chunkie (BSD 3-clause).
%
%  See also GREEN2D.HELM

if nargin < 3, nolog = false; end

rx = targ(1,:).' - src(1,:);
ry = targ(2,:).' - src(2,:);
r2 = rx.^2 + ry.^2;

if nolog
    val = [];
else
    val = -log(r2)/(4*pi);
end

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
