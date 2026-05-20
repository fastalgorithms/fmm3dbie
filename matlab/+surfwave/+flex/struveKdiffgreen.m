function [val,grad,hess,third,fourth] = struveKdiffgreen(rhoj,src,targ,ifr2logr)
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

if nargin < 4
    ifr2logr = false;
end

[~,ns] = size(src);
[~,nt] = size(targ);

% grad = zeros([nt,ns,2]);
% hess = zeros([nt,ns,3]);
% third = zeros([nt,ns,4]);
% fourth = zeros([nt,ns,5]);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

dx = xt-xs;
dy = yt-ys;

r2 = dx.*dx + dy.*dy;
r = sqrt(r2);

eulergamma = 0.5772156649015328606065120900824;

ilow = (imag(rhoj) < 0);

if ilow
    rhoj = conj(rhoj); % flip rhoj up into upper half plane
end

zt = r*rhoj;
[cr0,cr1] = surfwave.flex.struveR(zt);

if nargout > 4
    [h0,gradh0,hessh0,thirdh0,fourthh0] = surfwave.flex.helmdiffgreen(rhoj,src,targ,ifr2logr);
elseif nargout > 3
    [h0,gradh0,hessh0,thirdh0] = surfwave.flex.helmdiffgreen(rhoj,src,targ,ifr2logr);
elseif nargout > 2
    [h0,gradh0,hessh0] = surfwave.flex.helmdiffgreen(rhoj,src,targ,ifr2logr);
elseif nargout > 1
    [h0,gradh0] = surfwave.flex.helmdiffgreen(rhoj,src,targ,ifr2logr);
else
    [h0] = surfwave.flex.helmdiffgreen(rhoj,src,targ,ifr2logr);
end

h0(r < rslf) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
h0 = -4i*h0;

val = -1i*cr0+1i*h0;

if ilow
    val = conj(val);
end

if nargout > 1

h0x = gradh0(:,:,1);
h0y = gradh0(:,:,2);

h0x = -4i*h0x;
h0y = -4i*h0y;

cr0x = - rhoj*dx./r.*cr1;
cr0y = - rhoj*dy./r.*cr1;

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

r3 = r.^3;

dx2 = dx.*dx;
dy2 = dy.*dy;
    
h0xx = hessh0(:,:,1);
h0xy = hessh0(:,:,2);
h0yy = hessh0(:,:,3);

h0xx = -4i*h0xx;
h0xy = -4i*h0xy;
h0yy = -4i*h0yy;

cr1x = rhoj*dx./r.*cr0 - dx./r2.*cr1;
cr1y = rhoj*dy./r.*cr0 - dy./r2.*cr1;
cr0xx = - rhoj*dy2./r3.*cr1 - rhoj*dx./r.*cr1x;
cr0xy = rhoj*dx.*dy./r3.*cr1 - rhoj*dx./r.*cr1y;
cr0yy = - rhoj*dx2./r3.*cr1 - rhoj*dy./r.*cr1y;

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


if nargout > 3

r4 = r.^4;
r5 = r.^5;

h0xxx = thirdh0(:,:,1);
h0xxy = thirdh0(:,:,2);
h0xyy = thirdh0(:,:,3);
h0yyy = thirdh0(:,:,4);

h0xxx = -4i*h0xxx;
h0xxy = -4i*h0xxy;
h0xyy = -4i*h0xyy;
h0yyy = -4i*h0yyy;

cr1xx = rhoj*dy2./r3.*cr0 + rhoj*dx./r.*cr0x + (dx2-dy2)./r4.*cr1 - dx./r2.*cr1x;
cr1xy = - rhoj*dx.*dy./r3.*cr0 + 2*dx.*dy./r4.*cr1 + rhoj*dx./r.*cr0y - dx./r2.*cr1y;
cr1yy = rhoj*dx2./r3.*cr0 + rhoj*dy./r.*cr0y + (dy2-dx2)./r4.*cr1 - dy./r2.*cr1y;
cr0xxx = 3*rhoj*dy2.*dx./r5.*cr1 - 2*rhoj*dy2./r3.*cr1x - rhoj*dx./r.*cr1xx;
cr0yyy = 3*rhoj*dx2.*dy./r5.*cr1 - 2*rhoj*dx2./r3.*cr1y - rhoj*dy./r.*cr1yy;
cr0xxy = rhoj*dy.*(dy2 - 2*dx2)./r5.*cr1 + rhoj*dx.*dy./r3.*cr1x - rhoj*dy2./r3.*cr1y - rhoj*dx./r.*cr1xy;
cr0xyy = rhoj*dx.*(dx2 - 2*dy2)./r5.*cr1 + rhoj*dy.*dx./r3.*cr1y - rhoj*dx2./r3.*cr1x - rhoj*dy./r.*cr1xy;


thirdxxx = -1i*cr0xxx+1i*h0xxx;
thirdxxy = -1i*cr0xxy+1i*h0xxy;
thirdxyy = -1i*cr0xyy+1i*h0xyy;
thirdyyy = -1i*cr0yyy+1i*h0yyy;

thirdxxx(r < rslf) = 0;
thirdxxy(r < rslf) = 0;
thirdxyy(r < rslf) = 0;
thirdyyy(r < rslf) = 0;

if ilow
    thirdxxx = conj(thirdxxx);
    thirdxxy = conj(thirdxxy);
    thirdxyy = conj(thirdxyy);
    thirdyyy = conj(thirdyyy);
end

third(:,:,1) = thirdxxx;
third(:,:,2) = thirdxxy;
third(:,:,3) = thirdxyy;
third(:,:,4) = thirdyyy;
 
% gradlap(:,:,1) = thirdxxx + thirdxyy;
% gradlap(:,:,2) = thirdxxy + thirdyyy; 

end

if nargout > 4

r5 = r.^5;
r6 = r.^6;
r7 = r.^7;

dx3 = dx.*dx.*dx;
dy3 = dy.*dy.*dy;

dx4 = dx.*dx.*dx.*dx;
dy4 = dy.*dy.*dy.*dy;

h0xxxx = fourthh0(:,:,1);
h0xxxy = fourthh0(:,:,2);
h0xxyy = fourthh0(:,:,3);
h0xyyy = fourthh0(:,:,4);
h0yyyy = fourthh0(:,:,5);

h0xxxx = -4i*h0xxxx;
h0xxxy = -4i*h0xxxy;
h0xxyy = -4i*h0xxyy;
h0xyyy = -4i*h0xyyy;
h0yyyy = -4i*h0yyyy;

cr1xxx = -3*rhoj*dy2.*dx./r5.*cr0 + rhoj*dy2./r3.*cr0x - 2*dx.*(dx2-3*dy2)./r6.*cr1 + (dx2-dy2)./r4.*cr1x + ...
    rhoj*dy2./r3.*cr0x + rhoj*dx./r.*cr0xx + (dx2-dy2)./r4.*cr1x - dx./r2.*cr1xx;
cr1yyy = -3*rhoj*dx2.*dy./r5.*cr0 + rhoj*dx2./r3.*cr0y - 2*dy.*(dy2-3*dx2)./r6.*cr1 + (dy2-dx2)./r4.*cr1y + ...
    rhoj*dx2./r3.*cr0y + rhoj*dy./r.*cr0yy + (dy2-dx2)./r4.*cr1y - dy./r2.*cr1yy;
cr1xxy = rhoj*(2*dx2.*dy-dy3)./r5.*cr0 - rhoj*dx.*dy./r3.*cr0x + 2*dy.*(dy2-3*dx2)./r6.*cr1 + 2*dx.*dy./r4.*cr1x + ...
    rhoj*dy2./r3.*cr0y + rhoj*dx./r.*cr0xy + (dx2-dy2)./r4.*cr1y - dx./r2.*cr1xy;
cr1xyy = rhoj*(2*dy2.*dx-dx3)./r5.*cr0 - rhoj*dy.*dx./r3.*cr0y + 2*dx.*(dx2-3*dy2)./r6.*cr1 + 2*dy.*dx./r4.*cr1y + ...
    rhoj*dx2./r3.*cr0x + rhoj*dy./r.*cr0xy + (dy2-dx2)./r4.*cr1x - dy./r2.*cr1xy;

cr0xxxx = 3*rhoj*(dy4-4*dx2.*dy2)./r7.*cr1 + 6*rhoj*dy2.*dx./r5.*cr1x - rhoj*dy2./r3.*cr1xx + ...
    3*rhoj*dy2.*dx./r5.*cr1x - 2*rhoj*dy2./r3.*cr1xx - rhoj*dx./r.*cr1xxx;
cr0yyyy = 3*rhoj*(dx4-4*dy2.*dx2)./r7.*cr1 + 6*rhoj*dx2.*dy./r5.*cr1y - rhoj*dx2./r3.*cr1yy + ...
    3*rhoj*dx2.*dy./r5.*cr1y - 2*rhoj*dx2./r3.*cr1yy - rhoj*dy./r.*cr1yyy;
cr0xxxy = 3*rhoj*(2*dx3.*dy-3*dx.*dy3)./r7.*cr1 - 2*rhoj*(2*dx2.*dy-dy3)./r5.*cr1x + rhoj*dx.*dy./r3.*cr1xx+...
    3*rhoj*dy2.*dx./r5.*cr1y - 2*rhoj*dy2./r3.*cr1xy - rhoj*dx./r.*cr1xxy;
cr0xyyy = 3*rhoj*(2*dy3.*dx-3*dy.*dx3)./r7.*cr1 - 2*rhoj*(2*dy2.*dx-dx3)./r5.*cr1y + rhoj*dy.*dx./r3.*cr1yy+...
    3*rhoj*dx2.*dy./r5.*cr1x - 2*rhoj*dx2./r3.*cr1xy - rhoj*dy./r.*cr1xyy;

cr0xxyy = rhoj*(-2*dx4+11*dx2.*dy2-2*dy4)./r7.*cr1 + rhoj*(dy3-2*dy.*dx2)./r5.*cr1y + rhoj*(dx3-2*dy2.*dx)./r5.*cr1x + rhoj*dx.*dy./r3.*cr1xy+...
    rhoj*dx.*(dx2 - 2*dy2)./r5.*cr1x + rhoj*dy.*dx./r3.*cr1xy - rhoj*dx2./r3.*cr1xx - rhoj*dy./r.*cr1xxy;

fourthxxxx = -1i*cr0xxxx+1i*h0xxxx;
fourthxxxy = -1i*cr0xxxy+1i*h0xxxy;
fourthxxyy = -1i*cr0xxyy+1i*h0xxyy;
fourthxyyy = -1i*cr0xyyy+1i*h0xyyy;
fourthyyyy = -1i*cr0yyyy+1i*h0yyyy;

if ilow
    fourthxxxx = conj(fourthxxxx);
    fourthxxxy = conj(fourthxxxy);
    fourthxxyy = conj(fourthxxyy);
    fourthxyyy = conj(fourthxyyy);
    fourthyyyy = conj(fourthyyyy);
end

fourth(:,:,1) = fourthxxxx;
fourth(:,:,2) = fourthxxxy;
fourth(:,:,3) = fourthxxyy;
fourth(:,:,4) = fourthxyyy;
fourth(:,:,5) = fourthyyyy;

% bilap = fourthxxxx + 2*fourthxxyy + fourthyyyy;

end

end