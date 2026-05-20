function [val,grad,hess,third,bilap] = gphiflex(rts,ejs,src,targ)
%
% computes the green's function centered at (x,y) = 0 for the 
% integro-differential equation determined by the roots of the polynomial:
%             z^5 - beta*z + gamma = 0
%
% output is a cell array with all the kernels needed to solve the adjointed
% Lippman Schwinger equation: {val,hessxx,hessxy,hessyy,gradlapx,gradlapy}
% where:
% - val is the value of the Green's function centered at zero and
%   evaluated at (x,y)
% - hessxx is G_{xx}, hessxy is G_{xy}, 
%   hessyy is G_{yy}
% - gradlap is the gradient of the Laplacian, namely 
%   gradlapx is G_{xxx} + G_{xyy}, gradlapy is G_{yxx} + G_{yyy}
%
% input:
%
% x - x-coordinates array
% y - y-coordinates array
% beta - coefficient beta in the equation
% gamma - coefficient gamma in the equation
%
% optional input:
%
% opt - bool, default: false. 
%         Possible options are:
%         opt = false => Green's function (and derivatives)
%         opt = true => kernel used to evaluate phi on surface
%

rslf = 1e-14;

if isstruct(src)
    src = src.r;
    targ = targ.r;
end

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

dx = xt-xs;
dy = yt-ys;

r = sqrt(dx.*dx + dy.*dy);

val = 0;
gradx = 0;
grady = 0;
hessxx = 0;
hessxy = 0;
hessyy = 0;
thirdxxx = 0;
thirdxxy = 0;
thirdxyy = 0;
thirdyyy = 0;
bilap = 0;

oneoveralpha = sum(ejs.*rts.^4);

eulergamma = 0.5772156649015328606065120900824;


for i = 1:5
    
    rhoj = rts(i);
    ej = ejs(i);

    if (abs(angle(rhoj)) < 1e-8) && (abs(rhoj) > rslf)

       [sk0,gradsk0,hesssk0,thirdsk0,fourthsk0] = surfwave.flex.struveKdiffgreen(rhoj,src,targ);
       [h0,gradh0,hessh0,thirdh0,fourthh0] = surfwave.flex.helmdiffgreen(rhoj,src,targ);

       h0(r < rslf) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));

       h0 = -4i*h0;
       gradh0 = -4i*gradh0;
       
       h0x = gradh0(:,:,1);
       h0y = gradh0(:,:,2);

       h0x(r < rslf) = 0;
       h0y(r < rslf) = 0;
       
       h0xx = hessh0(:,:,1);
       h0xy = hessh0(:,:,2);
       h0yy = hessh0(:,:,3);

       h0xx(r < rslf) = rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));
       h0xy(r < rslf) = 0;
       h0yy(r < rslf) = rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));
       
       h0xx = -4i*h0xx;
       h0xy = -4i*h0xy;
       h0yy = -4i*h0yy;

       h0xxx = thirdh0(:,:,1);
       h0xxy = thirdh0(:,:,2);
       h0xyy = thirdh0(:,:,3);
       h0yyy = thirdh0(:,:,4);

       h0xxx(r < rslf) = 0;
       h0xxy(r < rslf) = 0;
       h0xyy(r < rslf) = 0;
       h0yyy(r < rslf) = 0;

       h0xxx = -4i*h0xxx;
       h0xxy = -4i*h0xxy;
       h0xyy = -4i*h0xyy;
       h0yyy = -4i*h0yyy;

       h0xxxx = fourthh0(:,:,1);
       h0xxxy = fourthh0(:,:,2);
       h0xxyy = fourthh0(:,:,3);
       h0xyyy = fourthh0(:,:,4);
       h0yyyy = fourthh0(:,:,5);

       h0xxxx(r < rslf) = 0;
       h0xxxy(r < rslf) = 0;
       h0xxyy(r < rslf) = 0;
       h0xyyy(r < rslf) = 0;
       h0yyyy(r < rslf) = 0;

       h0xxxx = -4i*h0xxxx;
       h0xxxy = -4i*h0xxxy;
       h0xxyy = -4i*h0xxyy;
       h0xyyy = -4i*h0xyyy;
       h0yyyy = -4i*h0yyyy;

       sk0x = gradsk0(:,:,1);
       sk0y = gradsk0(:,:,2);

       sk0xx = hesssk0(:,:,1);
       sk0xy = hesssk0(:,:,2);
       sk0yy = hesssk0(:,:,3);

       sk0xxx = thirdsk0(:,:,1);
       sk0xxy = thirdsk0(:,:,2);
       sk0xyy = thirdsk0(:,:,3);
       sk0yyy = thirdsk0(:,:,4);

       sk0xxxx = fourthsk0(:,:,1);
       sk0xxxy = fourthsk0(:,:,2);
       sk0xxyy = fourthsk0(:,:,3);
       sk0xyyy = fourthsk0(:,:,4);
       sk0yyyy = fourthsk0(:,:,5);

       val = val + ej*rhoj*(-sk0 + 2i*h0);
       % phi = phi + ej*rhoj*(-sk0 + 2i*h0);

       gradx = gradx + ej*rhoj*(-sk0x + 2i*h0x);
       grady = grady + ej*rhoj*(-sk0y + 2i*h0y);
       
       hessxx = hessxx + ej*rhoj*(-sk0xx + 2i*h0xx);
       hessxy = hessxy + ej*rhoj*(-sk0xy + 2i*h0xy);
       hessyy = hessyy + ej*rhoj*(-sk0yy + 2i*h0yy);

       thirdxxx = thirdxxx + ej*rhoj*(-sk0xxx + 2i*h0xxx);
       thirdxxy = thirdxxy + ej*rhoj*(-sk0xxy + 2i*h0xxy);
       thirdxyy = thirdxyy + ej*rhoj*(-sk0xyy + 2i*h0xyy);
       thirdyyy = thirdyyy + ej*rhoj*(-sk0yyy + 2i*h0yyy);

       bilap = bilap + ej*rhoj*(-sk0xxxx - 2*sk0xxyy - sk0yyyy ...
           + 2i*h0xxxx +  4i*h0xxyy +  2i*h0yyyy);


    elseif abs(rhoj) > rslf

       [sk0,gradsk0,hesssk0,thirdsk0,fourthsk0] = surfwave.flex.struveKdiffgreen(-rhoj,src,targ);

       sk0x = gradsk0(:,:,1);
       sk0y = gradsk0(:,:,2);

       sk0xx = hesssk0(:,:,1);
       sk0xy = hesssk0(:,:,2);
       sk0yy = hesssk0(:,:,3);

       sk0xxx = thirdsk0(:,:,1);
       sk0xxy = thirdsk0(:,:,2);
       sk0xyy = thirdsk0(:,:,3);
       sk0yyy = thirdsk0(:,:,4);

       sk0xxxx = fourthsk0(:,:,1);
       sk0xxxy = fourthsk0(:,:,2);
       sk0xxyy = fourthsk0(:,:,3);
       sk0xyyy = fourthsk0(:,:,4);
       sk0yyyy = fourthsk0(:,:,5);

       val = val + ej*rhoj*sk0;

       gradx = gradx + ej*rhoj*sk0x;
       grady = grady + ej*rhoj*sk0y;

       hessxx = hessxx + ej*rhoj*sk0xx;
       hessxy = hessxy + ej*rhoj*sk0xy;
       hessyy = hessyy + ej*rhoj*sk0yy;

       thirdxxx = thirdxxx + ej*rhoj*(sk0xxx);
       thirdxxy = thirdxxy + ej*rhoj*(sk0xxy);
       thirdxyy = thirdxyy + ej*rhoj*(sk0xyy);
       thirdyyy = thirdyyy + ej*rhoj*(sk0yyy);

       bilap = bilap + ej*rhoj*(sk0xxxx + 2*sk0xxyy + sk0yyyy);

    end

end

if nargout > 3
    r3 = r.^3;
    dx2 = dx.*dx;
    dx3 = dx.*dx2;

    dy2 = dy.*dy;
    dy3 = dy.*dy2;

    thirdxxx = oneoveralpha/(2*pi)*dx3./r3 + 1/4*thirdxxx;
    thirdxxy = oneoveralpha/(2*pi)*dx2.*dy./r3 + 1/4*thirdxxy;
    thirdxyy = oneoveralpha/(2*pi)*dx.*dy2./r3 + 1/4*thirdxyy;
    thirdyyy = oneoveralpha/(2*pi)*dy3./r3 + 1/4*thirdyyy;
    bilap = oneoveralpha/(2*pi)./r + 1/4*bilap;
    bilap(r < rslf) = 0;
end

val = 1/4*val;
gradx = 1/4*gradx;
grady = 1/4*grady;
hessxx = 1/4*hessxx;
hessxy = 1/4*hessxy;
hessyy = 1/4*hessyy;


grad = cat(3,gradx,grady);
hess = cat(3,hessxx,hessxy,hessyy);
third = cat(3,thirdxxx,thirdxxy,thirdxyy,thirdyyy);

end