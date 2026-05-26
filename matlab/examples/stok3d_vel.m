ctr1 = [0;0;0];


S = geometries.ellipsoid([1,1,1.5],[3,3,3],ctr1,6);

kerns = kernel3d('stokes','c',[1,1]);

eps = 1e-10;
tic;
Smat = surfermat(S,kerns,eps);
Smat = Smat+0.5*eye(size(Smat));
tbuild = toc

rhskerns = kernel3d('stokes','s');
%%
src = []; src.r = ctr1; src.d = [1;0;0];
src2 = []; src2.r = [-2;0;0.5]; src2.d = [1;0;0];

rhs = -[rhskerns.eval(src,S)*src.d,rhskerns.eval(src2,S)*src2.d];

tic;
sol = Smat\rhs;
tsolve = toc

kernseval = kerns;

% %%
nplot = 100;
xx = linspace(-2,6,nplot)+1e-2;  yy = 2*linspace(-2,2,nplot)+3e-2;
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).'; YY(:).'; 0*YY(:).'];

tic;
skern = kernel3d('stokes','s');
uin = skern.eval(src,targs)*src.d;
uin2 = skern.eval(src2,targs)*src2.d;
uscat = surferkerneval(S,kernseval,sol,targs,eps);
tplot = toc
%%
uscat1 = uscat(:,1); uscat2 = uscat(:,2);
uin = reshape(uin,3,[]);
uscat1 = reshape(uscat1,3,[]);
figure(1);clf
h = pcolor(XX,YY,reshape(log10(vecnorm(uin+uscat1)),size(XX))); h.EdgeColor = 'none';
hold on
plot(S,zeros(S.npts,1))
hold off
colorbar

figure(2);clf
h = pcolor(XX,YY,reshape((vecnorm(uin2+uscat2)),size(XX))); h.EdgeColor = 'none';
utots = uin2+uscat2;
utotx = reshape(utots(1,:),size(XX));
utoty = reshape(utots(2,:),size(XX));
utotz = reshape(utots(3,:),size(XX));
streamslice(XX,YY,0*YY,utotx,utoty,utotz,[],[],0)
hold on
plot(S,zeros(S.npts,1))
hold off
colorbar

%%
nplot = 40;
xx = linspace(-4,4,nplot)+1e-2;
[XX,YY,ZZ] = meshgrid(xx,xx,xx);
targs_vol = []; targs_vol.r = [XX(:).'; YY(:).'; ZZ(:).'];

tic;
skern = kernel3d('stokes','s');
uin3 = skern.eval(src2,targs_vol)*src2.d;
uscat3 = surferkerneval(S,kernseval,sol(:,2),targs_vol,eps);
tplotvol = toc

uin3 = reshape(uin3,3,[]);
uscat3 = reshape(uscat3,3,[]);

%%
figure(3);clf
utots = uin3+uscat3;
utotx = reshape(utots(1,:),size(XX));
utoty = reshape(utots(2,:),size(XX));
utotz = reshape(utots(3,:),size(XX));
zs = -2:0.6:2;
for z = zs
l = streamslice(XX,YY,ZZ,utotx,utoty,utotz,[],[],z,0.3);
set(l,'Color','k');
% streamslice(XX,YY,ZZ,utotx,utoty,utotz,[],[],3)
end
hold on
plot(S,zeros(S.npts,1))
hold off
colorbar

%%
figure(3);clf
xxstream = linspace(-4,4,15);
[startX,startY,startZ] = meshgrid([-1],xxstream,xxstream);
l = streamline(XX,YY,ZZ,utotx,utoty,utotz,startX,startY,startZ);
hold on
h = plot(S,zeros(S.npts,1));
h.FaceColor = 0.8*[1,1,1];
material('dull')
camlight('headlight')
hold off
% %%