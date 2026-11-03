function varargout = rootsum(rts, cjs, src, targ, realtol, c3d)
%SURFWAVE.ROOTSUM   Sum over the dispersion roots of the surface wave
%   radial profiles.
%
%   [val,grad,hess,third,fourth] = SURFWAVE.ROOTSUM(rts,cjs,src,targ)
%   evaluates G(x,y) = sum_j cjs(j) K_j(|x-y|), with
%   K_j = -S_0(rho_j r) + 2i H_0^(1)(rho_j r) for a real root rho_j
%   (|angle| < realtol, default 1e-14) and S_0(-rho_j r) otherwise,
%   S_0 = H_0 - Y_0. c3d (default 0) adds c3d/(4 pi r).
%
%   src (2,ns) and targ (2,nt) are points in the plane. Outputs are
%   (nt,ns,:), ordered as in FLEX2D.HELMDIFFGREEN.
%
% See also KERNEL3D.CAPILLARY, KERNEL3D.ICEFLEX, KERNEL3D.GRAVITY

if ( nargin < 5 || isempty(realtol) ), realtol = 1e-14; end
if ( nargin < 6 ), c3d = 0; end
if ( isstruct(src) ),  src  = src.r;  end
if ( isstruct(targ) ), targ = targ.r; end
src = src(1:2,:);  targ = targ(1:2,:);

nout = max(nargout, 1);
out = cell(1, nout);
for k = 1:nout, out{k} = 0; end

for j = 1:numel(rts)
    rhoj = rts(j);
    P = cell(1, nout);
    if ( abs(angle(rhoj)) < realtol )
        S = cell(1, nout);  H = cell(1, nout);
        [S{:}] = struveKdiffgreen(rhoj, src, targ);
        [H{:}] = flex2d.helmdiffgreen(rhoj, src, targ);
        for k = 1:nout, P{k} = -S{k} + 8*H{k}; end     % 2i*(-4i H)
    else
        [P{:}] = struveKdiffgreen(-rhoj, src, targ);
    end
    for k = 1:nout, out{k} = out{k} + cjs(j)*P{k}; end
end

% the singular parts; moments that vanish up to rounding are snapped to 0
c0 = -2/pi*snap(cjs(:));
c1 = -2/pi*snap(cjs(:).*rts(:));
c3 =  2/pi*snap(cjs(:).*rts(:).^3);
dx = targ(1,:).' - src(1,:);
dy = targ(2,:).' - src(2,:);
r  = sqrt(dx.^2 + dy.^2);
[L, R] = logr_r(dx, dy, r, nout);
varargout = cell(1, nout);
varargout{1} = out{1} + c0*L{1};
for k = 2:nout
    varargout{k} = out{k} + c0*L{k} + c1*R{k};
end
% the rho_j^3 term of the third and fourth derivatives
if ( nout > 3 && c3 ~= 0 )
    ux = dx./r;  uy = dy./r;
    varargout{4} = varargout{4} + c3*cat(3, ux.^3, ux.^2.*uy, ux.*uy.^2, uy.^3);
    if ( nout > 4 )
        T1 = cat(3, ux.^4, ux.^3.*uy, ux.^2.*uy.^2, ux.*uy.^3, uy.^4);
        T2 = cat(3, 6*ux.^2, 3*ux.*uy, ux.^2 + uy.^2, 3*ux.*uy, 6*uy.^2);
        varargout{5} = varargout{5} + c3*(T2 - 7*T1)./r;
    end
end
if ( c3d ~= 0 )
    P = lap3d(dx, dy, r, nout);
    for k = 1:nout, varargout{k} = varargout{k} + c3d*P{k}; end
end
end


function [L, R] = logr_r(dx, dy, r, nout)
% log r and r with their target derivatives
x = dx;  y = dy;
L = {log(r), cat(3, x./r.^2, y./r.^2), ...
     cat(3, (r.^2 - 2*x.^2)./r.^4, -2*x.*y./r.^4, (r.^2 - 2*y.^2)./r.^4), ...
     cat(3, (8*x.^3 - 6*x.*r.^2)./r.^6, 2*y.*(4*x.^2 - r.^2)./r.^6, ...
            2*x.*(4*y.^2 - r.^2)./r.^6, (8*y.^3 - 6*y.*r.^2)./r.^6), ...
     cat(3, 6*(-8*x.^4 + 8*x.^2.*r.^2 - r.^4)./r.^8, ...
            24*x.*y.*(r.^2 - 2*x.^2)./r.^8, ...
            2*(-24*x.^2.*y.^2 + 4*(x.^2 + y.^2).*r.^2 - r.^4)./r.^8, ...
            24*x.*y.*(r.^2 - 2*y.^2)./r.^8, ...
            6*(-8*y.^4 + 8*y.^2.*r.^2 - r.^4)./r.^8)};
R = {r, cat(3, x./r, y./r), ...
     cat(3, (r.^2 - x.^2)./r.^3, -x.*y./r.^3, (r.^2 - y.^2)./r.^3), ...
     cat(3, 3*x.*(x.^2 - r.^2)./r.^5, y.*(3*x.^2 - r.^2)./r.^5, ...
            x.*(3*y.^2 - r.^2)./r.^5, 3*y.*(y.^2 - r.^2)./r.^5), ...
     cat(3, 3*(-5*x.^4 + 6*x.^2.*r.^2 - r.^4)./r.^7, ...
            3*x.*y.*(3*r.^2 - 5*x.^2)./r.^7, ...
            (-15*x.^2.*y.^2 + 3*(x.^2 + y.^2).*r.^2 - r.^4)./r.^7, ...
            3*x.*y.*(3*r.^2 - 5*y.^2)./r.^7, ...
            3*(-5*y.^4 + 6*y.^2.*r.^2 - r.^4)./r.^7)};
L = L(1:nout);  R = R(1:nout);
end

function s = snap(c)
s = sum(c);
if ( abs(s) <= 1e-12*sum(abs(c)) ), s = 0; end
end


function P = lap3d(x, y, r, nout)
% 1/(4 pi r) and its target derivatives in the plane
c = 1./(4*pi*r);
P = cell(1, nout);
P{1} = c;
if ( nout > 1 ), P{2} = -c./r.^2.*cat(3, x, y); end
if ( nout > 2 )
    P{3} = c./r.^4.*cat(3, 2*x.^2 - y.^2, 3*x.*y, 2*y.^2 - x.^2);
end
if ( nout > 3 )
    P{4} = c./r.^6.*cat(3, -3*x.*(2*x.^2 - 3*y.^2), ...
        -3*y.*(4*x.^2 - y.^2), 3*x.*(x.^2 - 4*y.^2), 3*y.*(3*x.^2 - 2*y.^2));
end
if ( nout > 4 )
    P{5} = c./r.^8.*cat(3, 3*(8*x.^4 - 24*x.^2.*y.^2 + 3*y.^4), ...
        15*x.*y.*(4*x.^2 - 3*y.^2), -3*(4*x.^4 - 27*x.^2.*y.^2 + 4*y.^4), ...
        -15*x.*y.*(3*x.^2 - 4*y.^2), 3*(3*x.^4 - 24*x.^2.*y.^2 + 8*y.^4));
end
end


function [val,grad,hess,third,fourth] = struveKdiffgreen(rhoj,src,targ,ifr2logr)
% -i R_0(rhoj r) + i H_0^(1)(rhoj r) + 2/pi log r and its target
% derivatives, with the singular parts of R_0 subtracted

if nargin < 4
    ifr2logr = false;
end

[~,ns] = size(src);
[~,nt] = size(targ);

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
    rhoj = conj(rhoj);  % upper half plane
end

zt = r*rhoj;
[cr0,cr1] = struveR(zt);

if nargout > 4
    [h0,gradh0,hessh0,thirdh0,fourthh0] = flex2d.helmdiffgreen(rhoj,src,targ,ifr2logr);
elseif nargout > 3
    [h0,gradh0,hessh0,thirdh0] = flex2d.helmdiffgreen(rhoj,src,targ,ifr2logr);
elseif nargout > 2
    [h0,gradh0,hessh0] = flex2d.helmdiffgreen(rhoj,src,targ,ifr2logr);
elseif nargout > 1
    [h0,gradh0] = flex2d.helmdiffgreen(rhoj,src,targ,ifr2logr);
else
    [h0] = flex2d.helmdiffgreen(rhoj,src,targ,ifr2logr);
end

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

end

end


function [ch0,ch1] = struveR(z)
% J_0(z) + i H_0(z) and J_1(z) + i H_1(z) (H Struve), z in the upper half
% plane

    sz = size(z);
    z = z(:);
    iloc = find(abs(z)<95);
    ifar = find(abs(z)>=95);
    ch0 = zeros(size(z));
    ch1 = zeros(size(z));
    [ch0(iloc),ch1(iloc)] = eval_struve_loc(z(iloc));
    [ch0(ifar),ch1(ifar)] = eval_struve_asymp(z(ifar));
    
    ch0 = reshape(ch0,sz);
    ch1 = reshape(ch1,sz);

end

function [ch0,ch1] = eval_struve_loc(z)

        if (numel(z) == 0)
            ch0 = zeros(size(z));
            ch1 = zeros(size(z));
            return
        end

        % 37 point quadrature
        roots = [...
     0.87186936754544766250253027960067061E-03,
     0.46139685912320967303416540407238922E-02,
     0.11403724801942921356137834302761037E-01,
     0.21273869799713696933641832975796249E-01,
     0.34193590050116829966579764433743035E-01,
     0.50066369794921376232043513262965254E-01,
     0.68744840661990788702429710889821597E-01,
     0.90048072402715190852258825130292457E-01,
     0.11377592663404126825938959160714100E0,
     0.13971978110961365171254749916299551E0,
     0.16767014236115274323005600121328955E0,
     0.19742177893584148950155333084573482E0,
     0.22877691807308378517723923313666051E0,
     0.26154695322712432159879164746999082E0,
     0.29555302753335277446131742777185285E0,
     0.33062578374754158104960283618100391E0,
     0.36660450216541491786235920888991686E0,
     0.40333578530448759280923922177159180E0,
     0.44067189263979230096015237792322106E0,
     0.47846878035451058360863654925350577E0,
     0.51658385874506676265637066649992671E0,
     0.55487344175168025198877443769071049E0,
     0.59318982695586911097452516857230420E0,
     0.63137790840299165793725647143822823E0,
     0.66927118776716681032991184652397348E0,
     0.70668701278815968450670916337028782E0,
     0.74342084125637165388666138194419335E0,
     0.77923931924677709974401783652563375E0,
     0.81387200770526215128960217429540181E0,
     0.84700176082852057052869742284312294E0,
     0.87825418155746073449521778714979575E0,
     0.90718746700691083026622346312383955E0,
     0.93328559201548523410470858115116893E0,
     0.95596033956534706418672743248507695E0,
     0.97457068671996242739973875100248174E0,
     0.98846928432620062697670165425220407E0,
     0.99708064081719044882055745173383494E0];
        weights = [...
     0.22400793162522620028337342576607655E-02,
     0.52562958545833668597087967858577443E-02,
     0.83292445719813523846209588396397472E-02,
     0.11405895332021561946770211306887681E-01,
     0.14417175501532344966893665060168281E-01,
     0.17303524427299038655530608551526178E-01,
     0.20023022646990755626055401940839219E-01,
     0.22549907244778405352983171625076791E-01,
     0.24870905511862094728629286906339806E-01,
     0.26981803593965940026041345109116380E-01,
     0.28884661448585884128859085797934419E-01,
     0.30585621981978769080658948395974632E-01,
     0.32093212833515847024838346971574105E-01,
     0.33417052689694834454821729836007566E-01,
     0.34566884458490632207669519552867110E-01,
     0.35551863517485979214669339016958224E-01,
     0.36380034954008166378022870213703710E-01,
     0.37057940602559844694600638048292195E-01,
     0.37590304080249571384582448446974135E-01,
     0.37979748759608503619743452686633684E-01,
     0.38226508779578309155122146782475474E-01,
     0.38328096296860372305947506587248322E-01,
     0.38278889190620506125746342362630497E-01,
     0.38069602854625530975733505636189940E-01,
     0.37686608908088999340201452402700335E-01,
     0.37111065688452326336389956002501603E-01,
     0.36317836731785215519238677820354814E-01,
     0.35274207268375168408050987036082763E-01,
     0.33938490583125740266040378837852324E-01,
     0.32258792010039490880959886652098818E-01,
     0.30172542060499594231165932685367417E-01,
     0.27608011550870843420481696161406415E-01,
     0.24489904826860435351837466865555922E-01,
     0.20751994230408540299014138363576147E-01,
     0.16359482915027751299740298894463263E-01,
     0.11340067141542798168025906208432430E-01,
     0.58136209987079405666670968754469649E-02];

     w0 = weights./sqrt(1-roots.*roots);
     w1 = weights.*sqrt(1-roots.*roots);

     s = sum(w0.*roots);
     sw = sum(w0);
     zz = exp(1i*z)-1;

     ch0 = 0;
     ch1 = 0;

     for nnode = 1:37
        zr_tmp = z*roots(nnode);
        ex_tmp = exp(1i*zr_tmp);

        ch0 = ch0 + ex_tmp*w0(nnode);
        ch1 = ch1 + ex_tmp*w1(nnode);

     end

     ch0 = ch0 - zz*s-sw;
     ch0 = 2/pi*(ch0 + pi/2+zz);

     ch1 = 2/pi*(ch1.*z);

end

function [ch0,ch1] = eval_struve_asymp(z)

        zv = -z;
        eps = 1.0d-15;


        ch0 = 2/pi./zv;
        ct  = ch0;

        ch1 = 2/pi;
        ct1 = ch1;

        r1 = abs(ct);

        ifdone = ones(size(zv));

        for i=1:1000
            ct =-ct./zv./zv/4/i/i.* ...
               (2*i*(2*i-1)).^2;

            ct1 =-ct1./zv./zv*4.* ...
               ((i-1)^2-1/4);

            r2 = abs(ct);

            ifdone = ifdone.*(r2<=r1);
            ifdone = ifdone.*(r2>eps);

            if(sum(ifdone) == 0)
                break
            end
            ch0 = ch0 + ct;
            ch1 = ch1 + ct1;
            r1 = r2;
        end

        h0 = besselh(0,z);
        h1 = besselh(1,z);

        indsp = find(real(zv)>0);
        indsn = find(real(zv)<=0);
        ch0(indsp) = ch0(indsp)-1i*h0(indsp);
        ch0(indsn) = ch0(indsn)+1i*h0(indsn);
        ch1(indsp) = ch1(indsp)+1i*h1(indsp);
        ch1(indsn) = ch1(indsn)-1i*h1(indsn);

        ch0 = -1i*ch0;
        ch1 =  1i*ch1;
        ch0(ifdone==1) = NaN;
        ch1(ifdone==1) = NaN;

end
