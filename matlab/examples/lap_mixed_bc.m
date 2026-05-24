ctr1 = [0;0;0];
ctr2 = [0;0;0];


S1 = geometries.ellipsoid([1,1,1.5],[4,4,4],ctr1,6);

S2 = geometries.ellipsoid([1,1.5,1],[4,4,4],ctr2,7,11);

kerns(2,2) = kernel3d();
kerns(1,1) = 2*kernel3d('l','c',[1,1]);
kerns(2,1) = 2*kernel3d('l','cp',[1,1]);
kerns(1,2) = -2*kernel3d('l','sp');
kerns(2,2) = -2*kernel3d('l','s');
% kerns(1,1) = -2*kernel3d('l','sp');

eps = 1e-10;
srfrs = [S1,S2];
tic;
Smat = surfermat(srfrs,kerns,eps);
Smat = Smat+eye(size(Smat));
tbuild = toc
src = []; src.r = ctr1;
rhskerns(2,1) = kernel3d();
rhskerns(1) = 2*kernel3d('l','s');
rhskerns(2) = -2*kernel3d('l','sp');

rhs = zeros(size(Smat,1));
inds1 = surferids(srfrs,1);
rhs(inds1) = -2*rhskerns(1).eval(src,S1);
inds2 = surferids(srfrs,2);
rhs(inds2) = -(-2)*rhskerns(2).eval(src,S2);

tic;
sol = Smat\rhs;
tsolve = toc

kernseval = kerns(1,:);

nplot = 30;
xx = linspace(-2,6,nplot);  yy = 2*linspace(-2,2,nplot);
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).'; YY(:).'; 0*YY(:).'];

tic;
skern = kernel3d('l','s');
uin = skern.eval(src,targs);
uscat = surferkerneval(srfrs,kernseval,sol,tsolve,eps);
tplot = toc


figure(1);clf
scatter3(targs.r(1,:),targs.r(2,:),targs.r(3,:),[],real(uscat),'.');
hold on
plot(srfrs)
hold off


