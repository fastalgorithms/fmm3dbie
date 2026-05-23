S = geometries.ellipsoid([1,1,1.],[3,3,3],[],8);
S2 = geometries.ellipsoid([1,1.,1],[3,3,3],[4;0;0],8);


figure(1);clf
plot(S)
hold on
plot(S2)
hold off

% kerns(2,2) = kernel3d();
% % kerns(:,:) = kernel3d('h','d',0.1);
% % kerns(:,:) = kernel3d('h','c',0.1,[1,1]);
% kerns(:,:) = kernel3d('l','c',[1,1]);
kerns = kernel3d('l','c',[1,1]);
srfrs = [S,S2];
% srfrs = S;
tic;
Smat = surfermat(srfrs,kerns);
Smat = Smat + 0.5*eye(size(Smat));

syscors = surfermat(srfrs,kerns,struct('corrections',1));
syscors = syscors + 0.5*speye(size(syscors));
toc;
%%
% src = []; src.r = [2;0;1];
src = []; src.r = [0;0;0];
skern = kernel3d('l','s');
smerge = merge(srfrs);
rhs = -skern.eval(src,smerge);

sysapply = @(dens) surfermatapply(srfrs,kerns,dens,syscors);
norm(sysapply(rhs)-Smat*rhs)
tic;
% sol = Smat\rhs;
sol = gmres(Smat,rhs,[],1e-5,100);
sol2 = gmres(sysapply,rhs,[],1e-5,100);
toc;
norm(sol - sol2)/norm(sol)

%%
kernseval = kerns(1,:);

nplot = 100;
xx = linspace(-2,6,nplot);
yy = 2*linspace(-2,2,nplot);
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).';YY(:).';0*YY(:).'];

tic;
% uin = skern.eval(src,targs);
% evalmat = surferkernevalmat(srfrs,kernseval,targs);
% uscat3 = evalmat*sol;
% opts.corrections = 1;
% cors = surferkernevalmat(srfrs,kernseval,targs,opts);
% opts.corrections = cors;
% uscat2 = surferkerneval(srfrs,kernseval,sol,targs,opts);
uscat = surferkerneval(srfrs,kernseval,sol,targs);
toc

% norm(uscat - uscat2) / norm(uscat)
% norm(uscat - uscat3) / norm(uscat)

%%
figure(1);clf
plot(smerge,zeros(smerge.npts,1))
hold on
% scatter3(targs.r(1,:),targs.r(2,:),targs.r(3,:),[],uin+uscat,'.');colorbar
scatter3(targs.r(1,:),targs.r(2,:),targs.r(3,:),[],log10(abs(uin+uscat)),'.');colorbar
hold off

figure(2);clf
plot(smerge,rhs);colorbar