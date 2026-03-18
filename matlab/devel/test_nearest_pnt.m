clear all;close all
% Two torii test to get 4 digits of accuracy in the solution
nu = 1;
nv = 3;

rmajor = 1;
rminor = 0.5;
rwave = 0.0;

radii2 = [rmajor,rminor,rwave];
radii = [rminor,rmajor,rwave];
nuv = [nu,nv];
norder = 4;
nosc = 1;
scales = [1,1,1];
iptype_pat = 1;
iptyperef = 1;
tol = 1e-3;
tolref = 1e-4;

lptype = 0;

S = geometries_adap.wobblytorus_adap(radii,nosc,scales,nuv,norder,iptype_pat,tol,lptype);

%% generate targs
ntarg = 1000;
rng(1)
phitarg = rand(1,ntarg)*2*pi;
thetatarg = rand(1,ntarg)*2*pi;
radtarg = (1 + (rand(1,ntarg)*0.1-0.05)) * rminor;

xtarg = (rmajor + radtarg.*cos(phitarg)).*cos(thetatarg);
ytarg = (rmajor + radtarg.*cos(phitarg)).*sin(thetatarg);
ztarg = radtarg.*sin(phitarg);
targ = cat(1,xtarg,ytarg,ztarg);

xnear = (rmajor + rminor.*cos(phitarg)).*cos(thetatarg);
ynear = (rmajor + rminor.*cos(phitarg)).*sin(thetatarg);
znear = rminor.*sin(phitarg);

xyznear = cat(1,xnear,ynear,znear);

maxiter = 100;
tol = 1e-14;

% targ: [3, ntarg]
% Sr  : [3, npts]

[idx_init, dists_knn] = knnsearch(S.r.', targ.');

% idx  : [ntarg,1], index in Sr
% dist : [ntarg,1], Euclidean distance

[nearest_pnts, nearest_uvs,dists,flags] = ... 
    find_nearest_surf_point(S,targ,idx_init,tol,maxiter);


figure(1)
S.plot
alpha 0.4
hold on
plot3(targ(1,:),targ(2,:),targ(3,:),'r*')
plot3(nearest_pnts(1,:),nearest_pnts(2,:),nearest_pnts(3,:),'bo')
hold off

