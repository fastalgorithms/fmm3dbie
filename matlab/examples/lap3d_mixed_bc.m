% Solve a scattering problem with a Dirichlet obstacle and a Neumann
% obstacle

%% Setup geometries and system matrix
ctr1 = [0;0;0];
ctr2 = [3;0;0];

S1 = geometries.ellipsoid([1,1,1.5],3*[1,1,1],ctr1,6);
S2 = geometries.ellipsoid([1,1.5,1],3*[1,1,1],ctr2,6);

zk = 0.1;
kerns(2,2) = kernel3d();
kerns(1,1) = 2*kernel3d('l','c',[1,1]);
kerns(1,2) = -2*kernel3d('l','s');
kerns(2,1) = 2*kernel3d('l','cp',[1,1]);
kerns(2,2) = -2*kernel3d('l','sp');

eps = 1e-10;
srfrs = [S1,S2];
tic;
Smat = surfermat(srfrs,kerns,eps);
Smat = Smat+eye(size(Smat));
tbuild = toc

%% Get right hand side and solve

rhskerns(2,1) = kernel3d();
rhskerns(1) = kernel3d('l','s');
rhskerns(2) = kernel3d('l','sp');

src = []; src.r = ctr2;
rhs = zeros(size(Smat,1),1);
inds1 = surferids(srfrs,1);
rhs(inds1) = -rhskerns(1).eval(src,S1);
inds2 = surferids(srfrs,2);
rhs(inds2) = -rhskerns(2).eval(src,S2);

tic;
sol = Smat\rhs;
tsolve = toc

kernseval = kerns(1,:);
nplot = 100;
xx = linspace(-2,6,nplot)+1e-2;  yy = 2*linspace(-2,2,nplot)+3e-2;
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).'; YY(:).'; 0*YY(:).'];

tic;
skern = kernel3d('l','s');
uin = skern.eval(src,targs);
uscat = surferkerneval(srfrs,kernseval,sol,targs,eps);
tplot = toc

%% Plot

figure(1);clf
h = pcolor(XX,YY,reshape(log10(abs(uin+uscat)),size(XX))); h.EdgeColor = 'none';
hold on
plot(srfrs,zeros(size(rhs,1),1))
hold off
colorbar