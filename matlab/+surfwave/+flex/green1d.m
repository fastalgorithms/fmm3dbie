function [val,grad,hess,third,phi] = green1d(rts,ejs,src,targ,theta)
% green1d evaluates the one-dimensional G_S centered at src pointing in the
% direction of theta
%
% All derivatives are w.r.t. the target xt

[~,ns] = size(src);
[~,nt] = size(targ);

if ns > 1
    error('source must be a single point')
end

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

dx = xt-xs;
dy = yt-ys;

dx2 = dx.*dx;
dy2 = dy.*dy;

% r2 = dx2 + dy2;
% r = sqrt(r2);
% r3 =  r2.*r;
% r4 = r2.*r2;
% r5 = r2.*r3;

val = 0;
phi = 0;
gradx = 0;
grady = 0;
hessxx = 0;
hessxy = 0;
hessyy = 0;
thirdxxx = 0;
thirdxxy = 0;
thirdxyy = 0;
thirdyyy = 0;

r = (cos(theta)*dx + sin(theta)*dy);
r2 = (cos(theta)*dx + sin(theta)*dy).^2;


for i = 1:5

    rhoj = rts(i);
    ej = ejs(i);

    val = val + 1/(2*pi)*ej*rhoj*(surfwave.flex.expeval(1i*rhoj*abs(r))+surfwave.flex.expeval(-1i*rhoj*abs(r)));

    gradx = gradx + sign(r)*cos(theta)*1i/(2*pi)*ej*rhoj^2.*(surfwave.flex.expeval(1i*rhoj*abs(r))-surfwave.flex.expeval(-1i*rhoj*abs(r)));
    grady = grady + sign(r)*sin(theta)*1i/(2*pi)*ej*rhoj^2.*(surfwave.flex.expeval(1i*rhoj*abs(r))-surfwave.flex.expeval(-1i*rhoj*abs(r)));

    hessxx = hessxx - cos(theta)^2/(2*pi)*ej*rhoj^3.*(surfwave.flex.expeval(1i*rhoj*abs(r))+surfwave.flex.expeval(-1i*rhoj*abs(r)));
    hessxy = hessxy - cos(theta)*sin(theta)/(2*pi)*ej*rhoj^3.*(surfwave.flex.expeval(1i*rhoj*abs(r))+surfwave.flex.expeval(-1i*rhoj*abs(r)));
    hessyy = hessyy - sin(theta)^2/(2*pi)*ej*rhoj^3.*(surfwave.flex.expeval(1i*rhoj*abs(r))+surfwave.flex.expeval(-1i*rhoj*abs(r)));

    thirdxxx = thirdxxx - sign(r)*cos(theta)^3*1i/(2*pi)*ej*rhoj^4.*(surfwave.flex.expeval(1i*rhoj*abs(r))-surfwave.flex.expeval(-1i*rhoj*abs(r)));
    thirdxxy = thirdxxy - sign(r)*cos(theta)^2*sin(theta)*1i/(2*pi)*ej*rhoj^4.*(surfwave.flex.expeval(1i*rhoj*abs(r))-surfwave.flex.expeval(-1i*rhoj*abs(r)));
    thirdxyy = thirdxyy - sign(r)*sin(theta)^2*cos(theta)*1i/(2*pi)*ej*rhoj^4.*(surfwave.flex.expeval(1i*rhoj*abs(r))-surfwave.flex.expeval(-1i*rhoj*abs(r)));
    thirdyyy = thirdyyy - sign(r)*sin(theta)^3*1i/(2*pi)*ej*rhoj^4.*(surfwave.flex.expeval(1i*rhoj*abs(r))-surfwave.flex.expeval(-1i*rhoj*abs(r)));

    phi = phi + 1/(4*pi)*ej*(surfwave.flex.expeval(1i*rhoj*abs(r))+surfwave.flex.expeval(-1i*rhoj*abs(r)));

    if (angle(rhoj) < 1e-8) && (angle(rhoj) >= 0) && (rhoj ~= 0)

        val = val + 1i*ej*rhoj*exp(1i*rhoj*abs(r));

        gradx = gradx - sign(r).*cos(theta)*ej*rhoj^2.*exp(1i*rhoj*abs(r));
        grady = grady - sign(r).*sin(theta)*ej*rhoj^2.*exp(1i*rhoj*abs(r));

        hessxx = hessxx - cos(theta)^2*1i*ej*rhoj^3.*exp(1i*rhoj*abs(r));
        hessxy = hessxy - cos(theta)*sin(theta)*1i*ej*rhoj^3.*exp(1i*rhoj*abs(r));
        hessyy = hessyy - sin(theta)^2*1i*ej*rhoj^3.*exp(1i*rhoj*abs(r));

        thirdxxx = thirdxxx + sign(r)*cos(theta)^3*ej*rhoj^4.*exp(1i*rhoj*abs(r));
        thirdxxy = thirdxxy + sign(r)*cos(theta)^2*sin(theta)*ej*rhoj^4.*exp(1i*rhoj*abs(r));
        thirdxyy = thirdxyy + sign(r)*cos(theta)*sin(theta)^2*ej*rhoj^4.*exp(1i*rhoj*abs(r));
        thirdyyy = thirdyyy + sign(r)*sin(theta)^3*ej*rhoj^4.*exp(1i*rhoj*abs(r));

        phi = phi + 1i/2*ej*exp(1i*rhoj*abs(r));

    end

end

grad = cat(3,gradx,grady);
hess = cat(3,hessxx,hessxy,hessyy);
third = cat(3,thirdxxx,thirdxxy,thirdxyy,thirdyyy);


end