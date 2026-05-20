function [val,grad,hess,der3] = gsflexstab(rts,ejs,src,targ)
%GSFLEXSTAB evaluates the derivatives of the flexural-gravity wave green's  
% function in a stable manner using power series
%
% - grad(:,:,1) has G_{x1}, grad(:,:,2) has G_{x2}
% - hess(:,:,1) has G_{x1x1}, hess(:,:,2) has G_{x1x2}, 
% hess(:,:,3) has G_{x2x2}
% - der3 has the third derivatives in the order G_{x1x1x1}, G_{x1x1x2}, 
% G_{x1x2x2}, G_{x2x2x2}
% - der4 has the fourth derivatives in the order G_{x1x1x1x1}, 
% G_{x1x1x1x2}, G_{x1x1x2x2}, G_{x1x2x2x2}, G_{x2x2x2x2}
%
% derivatives are on the *target* variables
%
% input:
%
% src - (2,ns) array of source locations
% targ - (2,nt) array of target locations
% k - wave number, as above
%
% optional input:
%
% ifr2logr - boolean, default: false. If true, also subtract off the 
%             k^2/(8pi) r^2 log r kernel


[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

dx = xt-xs;
dy = yt-ys;

dx2 = dx.*dx;
dy2 = dy.*dy;

r = sqrt(dx2 + dy2);

% get value and r derivatives
      
[g0,g1,g21,g3,g4,g5] = gs_and_rders(rts,ejs,r);

%     evaluate potential and derivatives

if nargout > 0
    val = g0;  
end
if nargout > 1
    rm1 = 1./r;

    grad(:,:,1) = dx.*g1.*rm1;
    grad(:,:,2) = dy.*g1.*rm1;
end
if nargout > 2
    rm2 = rm1.*rm1;

    hess(:,:,1) = dx2.*g21.*rm2+g1.*rm1;
    hess(:,:,2) = dx.*dy.*g21.*rm2;
    hess(:,:,3) = dy2.*g21.*rm2+g1.*rm1;
end
if nargout > 3
    rm3 = rm1.*rm2;
    rm4 = rm1.*rm3;
    dx3 = dx2.*dx;
    dy3 = dy2.*dy;

    der3(:,:,1) = (dx3.*g3+3*dy2.*dx.*g21.*rm1).*rm3;
    der3(:,:,2) = dx2.*dy.*(g3.*rm3-3*g21.*rm4) + ...
             dy.*g21.*rm2;
    der3(:,:,3) = dx.*dy2.*(g3.*rm3-3*g21.*rm4) + ...
             dx.*g21.*rm2;
    der3(:,:,4) = (dy3.*g3+3*dx2.*dy.*g21.*rm1).*rm3;
end
if nargout > 4
    dx4 = dx3.*dx;
    dy4 = dy3.*dy;

    der4(:,:,1) = (dx4.*(g4-6*g3.*rm1+15*g21.*rm2)).*rm4 + ...
             (6*dx2.*(g3-3*g21.*rm1)).*rm3 + ...
             3*g21.*rm2;
    der4(:,:,2) = (dx3.*dy.*(g4-6*g3.*rm1+15*g21.*rm2)).*rm4 + ...
             (3*dx.*dy.*(g3-3*g21.*rm1)).*rm3;
    der4(:,:,3) = dx2.*dy2.*(g4-6*g3.*rm1+15*g21.*rm2).*rm4 + ...
             g3.*rm1 - 2*g21.*rm2;
    der4(:,:,4) = dx.*dy3.*(g4-6*g3.*rm1+15*g21.*rm2).*rm4 + ...
             3*dx.*dy.*(g3-3*g21.*rm1).*rm3;
    der4(:,:,5) = dy4.*(g4-6*g3.*rm1+15*g21.*rm2).*rm4 + ...
             6*dy2.*(g3-3*g21.*rm1).*rm3 + ...
             3*g21.*rm2;
end
if nargout > 5
    rm5 = rm1.*rm4;
    dx5 = dx4.*dx;
    dy5 = dy4.*dy;

    der5(:,:,1) = (dx5.*g5+10*dy2.*dx3.*g4.*rm1 + ...
          (15*dy4.*dx-30*dy2.*dx3).*g3.*rm2 + ...
         (60*dy2.*dx3-45*dy4.*dx).*g21.*rm3).*rm5;
    der5(:,:,2) = (dy.*dx4.*g5+(6*dy3.*dx2-4*dy.*dx4).*g4.*rm1 + ...
      (3*dy5+12*dy.*dx4-30*dy3.*dx2).*g3.*rm2 + ...
     (72*dy3.*dx2-9*dy5-24*dy.*dx4).*g21.*rm3).*rm5;
    der5(:,:,3) = (dy2.*dx3.*g5+(3*dy4.*dx-6*dy2.*dx3+dx5).*g4.*rm1 + ...
      (27*dy2.*dx3-15*dy4.*dx-3*dx5).*g3.*rm2 + ...
     (36*dy4.*dx-63*dy2.*dx3+6*dx5).*g21.*rm3).*rm5;
    der5(:,:,4) = (dx2.*dy3.*g5+(3*dx4.*dy-6*dx2.*dy3+dy5).*g4.*rm1 + ...
      (27*dx2.*dy3-15*dx4.*dy-3*dy5).*g3.*rm2 + ...
     (36*dx4.*dy-63*dx2.*dy3+6*dy5).*g21.*rm3).*rm5;
    der5(:,:,5) = (dx.*dy4.*g5+(6*dx3.*dy2-4*dx.*dy4).*g4.*rm1 + ...
      (3*dx5+12*dx.*dy4-30*dx3.*dy2).*g3.*rm2 + ...
     (72*dx3.*dy2-9*dx5-24*dx.*dy4).*g21.*rm3).*rm5;
    der5(:,:,6) = (dy5.*g5+10*dx2.*dy3.*g4.*rm1 + ...
      (15*dx4.*dy-30*dx2.*dy3).*g3.*rm2 + ...
     (60*dx2.*dy3-45*dx4.*dy).*g21.*rm3).*rm5;
end

end

function [g0,g1,g21,g3,g4,g5] = gs_and_rders(rts,ejs,r)
% g0 = g
% g1 = g'
% g21 = g'' - g'/r
%
% maybe later:
% g321 = g''' - 3*g''/r + 3g'/r^2
% g4321 = g'''' - 6*g'''/r + 15*g''/r^2 - 15*g'/r^3

k = rts(1); % assuming real root is the first one

g0 = zeros(size(r));
g1 = zeros(size(r));
g2 = zeros(size(r));
g3 = zeros(size(r));
g4 = zeros(size(r));
g5 = zeros(size(r));

io4 = 1i*0.25;
o2p = 1/(2*pi);

isus = abs(k)*r < 1;
%isus = false(size(r));
%isus = true(size(r));

% straightforward formulas for sufficiently large
% isus = false;
rnot = r(~isus);

[g00,g01,g02,g03,g04] = surfwave.flex.gsflex(rts,ejs,[0;0],[rnot;0*rnot]) ;

g0(~isus) = g00(:,:,1);
g1(~isus) = g01(:,:,1);
g2(~isus) = g02(:,:,1);
g3(~isus) = g03(:,:,1);
g4(~isus) = g04(:,:,1);

% manually cancel when small

rsus = r(isus);
rm1 = 1./rsus;
rm2 = rm1.*rm1;
rm3 = rm1.*rm2;
rm4 = rm1.*rm3;
rm5 = rm1.*rm4;

gam = 0.57721566490153286060651209;
nterms = 14;
const1 = (io4+(log(2)-gam-log(k))*o2p);

% relevant parts of hankel function represented as power series
[cj0m1,cy0] = surfwave.flex.besseldiff_etc_pscoefs(nterms);
ch0 = surfwave.flex.struve_pscoefs(nterms);

ejrj2 = ejs.*rts.^(2*(1:nterms)+2);
ejrj3 = ejs.*rts.^(2*(1:nterms)+1);

ct1 = 1/pi*(log(2./rts(1)) - gam + 1i*pi)*ejrj2(1,:).*cj0m1.';
ct2 = 1/pi*sum((log(-2./rts(2:end)) - gam).*(ejrj2(2:end,:))).*cj0m1.';
ct3 = 1/pi*sum(ejrj2) .* cy0.';
ct4 = -1/pi*sum(ejrj2) .* cj0m1.'; 
ct5 = -1/2*sum(ejrj3) .* ch0.' ;

logr = log(rsus);

t1 = surfwave.flex.even_pseval(ct1,rsus);
t2 = surfwave.flex.even_pseval(ct2,rsus);
t3 = surfwave.flex.even_pseval(ct3,rsus);
t4 = surfwave.flex.even_pseval(ct4,rsus);
t5 = surfwave.flex.odd_pseval(ct5,rsus);

t1 = t1 + 1/pi*ejs(1).*rts(1).^2*(log(2./rts(1)) + 1i*pi);
t2 = t2 + 1/pi*sum(ejs(2:end).*rts(2:end).^2.*(log(-2./rts(2:end))));

g0(isus) = t1+t2+t3+t4.*logr+t5


% differentiate power series to get derivatives
fac = 2*(1:nterms);
d21 = fac(:).*(fac(:)-1)-fac(:);
fd21 = surfwave.flex.even_pseval(cf2(:).*d21,rsus).*rm2;
cf1 = cf1.*fac(:); cf2 = cf2.*fac(:);
j0m1d1 = surfwave.flex.even_pseval(cf1,rsus).*rm1;
fd1 = surfwave.flex.even_pseval(cf2,rsus).*rm1;
cf1 = cf1.*(fac(:)-1); cf2 = cf2.*(fac(:)-1);
j0m1d2 = surfwave.flex.even_pseval(cf1,rsus).*rm2;
fd2 = surfwave.flex.even_pseval(cf2,rsus).*rm2;

cf1 = cf1(:).*(fac(:)-2); cf1 = cf1(2:end);
cf2 = cf2(:).*(fac(:)-2); cf2 = cf2(2:end);
j0m1d3 = surfwave.flex.even_pseval(cf1,rsus).*rm1;

fd3 = surfwave.flex.even_pseval(cf2,rsus).*rm1;
fac = fac(1:end-1);
cf1 = cf1(:).*(fac(:)-1); cf2 = cf2(:).*(fac(:)-1);
j0m1d4 = surfwave.flex.even_pseval(cf1,rsus).*rm2;
fd4 = surfwave.flex.even_pseval(cf2,rsus).*rm2;

cf1 = cf1(:).*(fac(:)-2); cf1 = cf1(2:end);
cf2 = cf2(:).*(fac(:)-2); cf2 = cf2(2:end);
j0m1d5 = surfwave.flex.even_pseval(cf1,rsus).*rm1;
fd5 = surfwave.flex.even_pseval(cf2,rsus).*rm1;

% combine to get derivative of i/4 H + log/(2*pi)
r2fac = -(1-r2logrfac)*k*k*0.25;
g0(isus) = const1*(j0m1+1+r2fac*rsus.*rsus) - o2p*(f  + logr.*j0m1);
g1(isus) = const1*(j0m1d1+2*r2fac*rsus) - o2p*(fd1 + logr.*j0m1d1 + j0m1.*rm1);
g21(isus) = const1*(j0m1d2-j0m1d1.*rm1) - o2p*(fd21 + logr.*(j0m1d2-j0m1d1.*rm1) + 2*j0m1d1.*rm1 - 2*j0m1.*rm2);
g3(isus) = const1*j0m1d3 - o2p*(fd3 + logr.*j0m1d3 + 3*j0m1d2.*rm1 - ...
    3*j0m1d1.*rm2 + 2*j0m1.*rm3);
g4(isus) = const1*j0m1d4 - o2p*(fd4 + logr.*j0m1d4 + 4*j0m1d3.*rm1 - ...
    6*j0m1d2.*rm2 + 8*j0m1d1.*rm3 - 6*j0m1.*rm4);
g5(isus) = const1*j0m1d5 - o2p*(fd5 + logr.*j0m1d5 + 5*j0m1d4.*rm1 - ...
    10*j0m1d3.*rm2 + 20*j0m1d2.*rm3 - 30*j0m1d1.*rm4 + 24*j0m1.*rm5);

end

