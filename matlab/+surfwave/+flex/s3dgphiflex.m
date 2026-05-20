function val = s3dgphiflex(rts,ejs,src,targ)
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

       sk0 = surfwave.flex.struveKdiffgreen(rhoj,src,targ);
       h0 = surfwave.flex.helmdiffgreen(rhoj,src,targ);

       h0(r < rslf) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));

       h0 = -4i*h0;
       
       val = val + ej*(-sk0 + 2i*h0);


    elseif abs(rhoj) > rslf

       sk0 = surfwave.flex.struveKdiffgreen(-rhoj,src,targ);

       val = val + ej*sk0;

    end

end

val = 1/8*val;

end