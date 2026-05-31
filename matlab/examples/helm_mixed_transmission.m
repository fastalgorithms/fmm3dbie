% Two boundaries, one Dirichlet, one transmission. Currently fails.
ctr1 = [0;0;0];
ctr2 = [3;0;0];

S1 = geometries.ellipsoid([1,1,1.5],3*[1,1,1],ctr1,6);
S2 = geometries.ellipsoid([1,1.5,1],3*[1,1,1],ctr2,6);

zk = 0.1;
zks = [zk,0.5];
% zks = [0.5,zk];
kerns(2,2) = kernel3d();
kerns(1,1) = 2*kernel3d('h','c',zk,[1,1]);
kerns(1,2) = kernel3d('h','trans_rep',zk);
kerns(2,1) = 2*kernel3d('h','c2trans',zks);
kerns(2,2) = kernel3d('h','trans_sys_diff',zks);

eps = 1e-10;
srfrs = [S1,S2];
tic;
Smat = surfermat(srfrs,kerns,eps);
Smat = Smat+eye(size(Smat));
tbuild = toc

rhskerns(2,1) = kernel3d();
rhskerns(1) = kernel3d('h','s',zk);
rhskerns(2) = kernel3d('h','s2trans',zks);
%%
src = []; src.r = ctr1;
% src = []; src.r = ctr2;
% src = []; src.r = [1;2;0.5];
rhs = zeros(size(Smat,1),1);
inds1 = surferids(srfrs,1);
rhs(inds1) = -rhskerns(1).eval(src,S1);
inds2 = surferids(srfrs,2);
rhs(S1.npts+1:end) = -rhskerns(2).eval(src,S2);

tic;
sol = Smat\rhs;
tsolve = toc

kernseval = kerns(1,:);
% %%
nplot = 100;
xx = linspace(-2,6,nplot)+1e-2;  yy = 2*linspace(-2,2,nplot)+3e-2;
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).'; YY(:).'; 0*YY(:).'];

tic;
skern = kernel3d('h','s',zk);
uin = skern.eval(src,targs);
uscat = surferkerneval(srfrs,kernseval,sol,targs,eps);
tplot = toc
%%

figure(1);clf
h = pcolor(XX,YY,reshape(log10(abs(uin+uscat)),size(XX))); h.EdgeColor = 'none';
hold on
plot(srfrs,zeros(sum([srfrs.npts]),1))
hold off
colorbar

figure(2);clf
subplot(1,2,1)
h = pcolor(XX,YY,reshape(real(uin),size(XX))); h.EdgeColor = 'none';
colorbar
clim([0,0.1])
subplot(1,2,2)
h = pcolor(XX,YY,reshape(real(-uscat),size(XX))); h.EdgeColor = 'none';
colorbar
clim([0,0.1])
% 
% figure(3);clf
% h = pcolor(XX,YY,reshape(real(uin+uscat),size(XX))); h.EdgeColor = 'none';
% clim([0,0.1])
% colorbar
% hold on
% h = plot(srfrs,ones(size(rhs,1),1));
% h.FaceColor = 0.8*[1,1,1];
% hold off

% figure(4);clf
% % plot(srfrs,rhs);colorbar
% plot(merge(srfrs),log10(surf_fun_error(merge(srfrs),rhs)));colorbar

