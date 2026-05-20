function [val,grad,hess] = struveK(rhoj,src,targ)
% evaluator for struveK(0,rhoj*x) and its derivatives with singularity
% subtraction
%
% K(x,y) = -i*R_0(rhoj|x-y|) + i*H^{(1)}_0(rhoj|x-y|) + 2/pi log(|x-y|)
%
% where the singular parts of the derivatives of R_0 are subtracted off. 
% 
% output :
% - val has the value of K(x,y) 
% - grad(:,:,1) has K_{x}, grad(:,:,2) has K_{y}
% - hess(:,:,1) has K_{xx}, hess(:,:,2) has K_{xy}, 
%   hess(:,:,3) has K_{yy}
%
% input: 
%
% rhoj - complex number
% src - [2,ns] array of source locations
% targ - [2,nt] array of target locations
%

rslf = 1e-14;
eulergamma = 0.57721566490153286060651209008240243;

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

r = sqrt((xt-xs).^2 + (yt-ys).^2);

ilow = (imag(rhoj) < 0);

if ilow
    rhoj = conj(rhoj); % flip rhoj up into upper half plane
end

zt = r*rhoj;
[cr0,cr1] = surfwave.capillary.struveR(zt);
h0 = chnk.helm2d.green(rhoj,src,targ);

h0(r < rslf) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
h0 = -4i*h0;

val = -1i*cr0+1i*h0;

if ilow
    val = conj(val);
end

if nargout > 1

dx = xt-xs;
dy = yt-ys;

[~,gradh0] = chnk.helm2d.green(rhoj,src,targ);

h0x = gradh0(:,:,1);
h0y = gradh0(:,:,2);

h0x = -4i*h0x;
h0y = -4i*h0y;

cr0x = 2i*rhoj*dx./(pi*r) - rhoj*dx./r.*cr1;
cr0y = 2i*rhoj*dy./(pi*r) - rhoj*dy./r.*cr1;

gradx = -1i*cr0x+1i*h0x;
grady = -1i*cr0y+1i*h0y;

gradx(r < rslf) = 0;
grady(r < rslf) = 0;

if ilow
    gradx = conj(gradx);
    grady = conj(grady);
end

grad(:,:,1) = gradx;
grad(:,:,2) = grady;
end

if nargout > 2

dx2 = dx.*dx;
dy2 = dy.*dy;

r2 = dx2 + dy2;
r3 = r.^3;

[~,~,hessh0] = chnk.helm2d.green(rhoj,src,targ);

h0xx = hessh0(:,:,1);
h0xy = hessh0(:,:,2);
h0yy = hessh0(:,:,3);

h0xx = -4i*h0xx;
h0xy = -4i*h0xy;
h0yy = -4i*h0yy;

cr1x = rhoj*dx./r.*cr0 - dx./r2.*cr1;
cr1y = rhoj*dy./r.*cr0 - dy./r2.*cr1;
cr0xx = 2i*rhoj*dy2./(pi*r3) - rhoj*dy2./r3.*cr1 - rhoj*dx./r.*cr1x;
cr0xy = -2i*rhoj*dx.*dy./(pi*r3) + rhoj*dx.*dy./r3.*cr1 - rhoj*dx./r.*cr1y;
cr0yy = 2i*rhoj*dx2./(pi*r3) - rhoj*dx2./r3.*cr1 - rhoj*dy./r.*cr1y;

hessxx = -1i*cr0xx+1i*h0xx;
hessxy = -1i*cr0xy+1i*h0xy;
hessyy = -1i*cr0yy+1i*h0yy;

hessxx(r < rslf) = 1i*rhoj^2-1i*rhoj^2/2-4i*1i*rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));
hessxy(r < rslf) = 0;
hessyy(r < rslf) = 1i*rhoj^2-1i*rhoj^2/2-4i*1i*rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));

if ilow
    hessxx = conj(hessxx);
    hessxy = conj(hessxy);
    hessyy = conj(hessyy);
end

hess(:,:,1) = hessxx;
hess(:,:,2) = hessxy;
hess(:,:,3) = hessyy;
end



end