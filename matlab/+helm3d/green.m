function [val,grad,hess] = green(k,src,targ)
%HELM3D.GREEN evaluate the helmholtz green's function
% for the given sources and targets

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);
zs = repmat(src(3,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);
zt = repmat(targ(3,:).',1,ns);

rx = xt-xs;
ry = yt-ys;
rz = zt-zs;

rx2 = rx.*rx;
ry2 = ry.*ry;
rz2 = rz.*rz;

r2 = rx2+ry2+rz2;

r = sqrt(r2);

efact = exp(1i*k*r)/(4*pi);


if nargout > 0
    val = efact./r;
end

[m,n] = size(xs);

if nargout > 1
    gd = -efact./(r.^2).*(1-1i*k*r);
    grad = zeros(m,n,3);
    grad(:,:,1) = gd.*rx./r;
    grad(:,:,2) = gd.*ry./r;
    grad(:,:,3) = gd.*rz./r;
end

if nargout > 2

    hess = zeros(m,n,3,3);
    
    gdd = efact./(r.^5);
    gdd1 = gd./r;
    gdd2 = gdd.*(3-3*1i*k*r-k^2*r2);
    
    hess(:,:,1,2) = gdd2.*rx.*ry;
    hess(:,:,2,1) = hess(:,:,1,2);
    hess(:,:,1,3) = gdd2.*rx.*rz;
    hess(:,:,3,1) = hess(:,:,1,3);
    hess(:,:,2,3) = gdd2.*ry.*rz;
    hess(:,:,3,2) = hess(:,:,2,3);
    hess(:,:,1,1) = gdd2.*rx2+gdd1;
    hess(:,:,2,2) = gdd2.*ry2+gdd1;
    hess(:,:,3,3) = gdd2.*rz2+gdd1;

end
