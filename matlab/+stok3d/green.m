function [val,grad] = green(src,targ)
%STOKES.GREEN evaluate the Stokes green's function
% for the given sources and targets

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);
zs = repmat(src(3,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);
zt = repmat(targ(3,:).',1,ns);

[m,n] = size(xs);


rx = reshape(xt-xs,[1,m,1,n]);
ry = reshape(yt-ys,[1,m,1,n]);
rz = reshape(zt-zs,[1,m,1,n]);

rx2 = rx.*rx;
ry2 = ry.*ry;
rz2 = rz.*rz;

r2 = rx2+ry2+rz2;

r = sqrt(r2);


gmat = zeros(3,m,3,n);
gmat(1,:,1,:) = rx2;
gmat(2,:,2,:) = ry2;
gmat(3,:,3,:) = rz2;

gmat(2,:,1,:) = rx.*ry;
gmat(3,:,1,:) = rx.*rz;
gmat(3,:,2,:) = ry.*rz;


gmat(1,:,2,:) = gmat(2,:,1,:);
gmat(1,:,3,:) = gmat(3,:,1,:);
gmat(2,:,3,:) = gmat(3,:,2,:);

fact = 1/8/pi;

if nargout > 0
    val = gmat;
    val(1,:,1,:) = val(1,:,1,:) + r2;
    val(2,:,2,:) = val(2,:,2,:) + r2;
    val(2,:,2,:) = val(2,:,2,:) + r2;
    val = fact*val./(r.^3);
end


gfact = 3/4/pi;
if nargout > 1
    grad = zeros(3,3,m,3,n);
    gmatg = gfact*reshape(gmat./(r.^5),[1,3,m,3,n]);
    grad(1,:,:,:,:) = gmatg.*reshape(rx,[1,1,m,1,n]);
    grad(2,:,:,:,:) = gmatg.*reshape(ry,[1,1,m,1,n]);
    grad(3,:,:,:,:) = gmatg.*reshape(rz,[1,1,m,1,n]);
end
