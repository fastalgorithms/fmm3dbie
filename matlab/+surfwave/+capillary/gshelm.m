function [val,grad,hess] = gshelm(rts,ejs,src,targ)
%
% computes the green's function for the integro-differential equation 
% determined by the polynomial:
%             -alpha z^2 + beta + gamma/|z| = -2
%
% outputs are:
% - val is the value of the Green's function centered at zero and
%   evaluated at (x,y)
% - grad(:,:,1) is G_{x}, grad(:,:,2) is G_{y} 
% - hess(:,:,1) is G_{xx}, hess(:,:,2) is G_{xy}, 
%   hess(:,:,3) is G_{yy}
%
% input:
%
% x - x-coordinates array
% y - y-coordinates array
% rts - roots of cubic polynomial
% ejs - residues (see notes)
%

% src = src.r;
% targ = targ.r;

rslf = 1e-14;

if isstruct(src)
    src = src.r;
end

if isstruct(targ)
    targ = targ.r;
end

eulergamma = 0.57721566490153286060651209008240243;

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

r = sqrt((xt-xs).^2 + (yt-ys).^2);

% clear xs ys xt yt

val = 0;
gradx = 0;
grady = 0;
hessxx = 0;
hessxy = 0;
hessyy = 0;

for i = 1:3
    
    rhoj = rts(i);
    ej = ejs(i);

    if (abs(angle(rhoj)) < rslf) && (abs(rhoj) > rslf)

       zt = r*rhoj;
       [cr0] = surfwave.capillary.struveR(zt);

       h0 = chnk.helm2d.green(rhoj,src,targ);
       h0(r < rslf) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
       h0 = -4i*h0;
        
       sk0 = -1i*cr0+1i*h0;
       
       val = val + ej*rhoj^2*(-sk0 + 2i*h0);

    elseif abs(rhoj) > rslf

       ilow = (imag(rhoj) > 0);

        if ilow
            rhoj = conj(rhoj); % flip -rhoj up into upper half plane
        end

       zt = -r*rhoj;
       [cr0] = surfwave.capillary.struveR(zt);

       h0 = chnk.helm2d.green(-rhoj,src,targ);

       h0(r < rslf) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
       h0 = -4i*h0;

       sk0 = -1i*cr0+1i*h0;

       if ilow
           sk0 = conj(sk0);
           rhoj = conj(rhoj); % flip -rhoj up into upper half plane
       end

       val = val + ej*rhoj^2*sk0;

    end

end

val = 1/2*val;

if nargout > 1

for i = 1:3
    
    rhoj = rts(i);
    ej = ejs(i);

    if (abs(angle(rhoj)) < rslf) && (abs(rhoj) > rslf)

       [~,gradsk0] = surfwave.capillary.struveK(rhoj,src,targ);
       [~,gradh0] = chnk.helm2d.green(rhoj,src,targ);

       gradh0 = -4i*gradh0;
       
       h0x = gradh0(:,:,1);
       h0y = gradh0(:,:,2);

       h0x(r < rslf) = 0;
       h0y(r < rslf) = 0;
       
       sk0x = gradsk0(:,:,1);
       sk0y = gradsk0(:,:,2);

       gradx = gradx + ej*rhoj^2*(-sk0x + 2i*h0x);
       grady = grady + ej*rhoj^2*(-sk0y + 2i*h0y);

    elseif abs(rhoj) > rslf

       [~,gradsk0] = surfwave.capillary.struveK(-rhoj,src,targ);

       sk0x = gradsk0(:,:,1);
       sk0y = gradsk0(:,:,2);

       gradx = gradx + ej*rhoj^2*sk0x;
       grady = grady + ej*rhoj^2*sk0y;

    end

end

gradx = 1/2*gradx;
grady = 1/2*grady;

grad = cat(3,gradx,grady);

end


if nargout > 2 

for i = 1:3
    
    rhoj = rts(i);
    ej = ejs(i);

    if (abs(angle(rhoj)) < rslf) && (abs(rhoj) > rslf)

       [~,~,hesssk0] = surfwave.capillary.struveK(rhoj,src,targ);
       [~,~,hessh0] = chnk.helm2d.green(rhoj,src,targ);
       
       h0xx = hessh0(:,:,1);
       h0xy = hessh0(:,:,2);
       h0yy = hessh0(:,:,3);

       h0xx(r < rslf) = rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));
       h0xy(r < rslf) = 0;
       h0yy(r < rslf) = rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));
       
       h0xx = -4i*h0xx;
       h0xy = -4i*h0xy;
       h0yy = -4i*h0yy;

       sk0xx = hesssk0(:,:,1);
       sk0xy = hesssk0(:,:,2);
       sk0yy = hesssk0(:,:,3);
       
       hessxx = hessxx + ej*rhoj^2*(-sk0xx + 2i*h0xx);
       hessxy = hessxy + ej*rhoj^2*(-sk0xy + 2i*h0xy);
       hessyy = hessyy + ej*rhoj^2*(-sk0yy + 2i*h0yy);


    elseif abs(rhoj) > rslf

       [~,~,hesssk0] = surfwave.capillary.struveK(-rhoj,src,targ);

       sk0xx = hesssk0(:,:,1);
       sk0xy = hesssk0(:,:,2);
       sk0yy = hesssk0(:,:,3);

       hessxx = hessxx + ej*rhoj^2*sk0xx;
       hessxy = hessxy + ej*rhoj^2*sk0xy;
       hessyy = hessyy + ej*rhoj^2*sk0yy;

    end

end

hessxx = 1/2*hessxx;
hessxy = 1/2*hessxy;
hessyy = 1/2*hessyy;

hess = cat(3,hessxx,hessxy,hessyy);

end




end