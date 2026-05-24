S = geometries.ellipsoid([1,1,1.1],[3,3,3],[],6);
S2 = geometries.ellipsoid([1,1.,1],[3,3,3],[4;0;0],8);


figure(1);clf
plot(S)
hold on
plot(S2)
hold off

eps = 1e-10;
% kerns(2,2) = kernel3d();
% % kerns(:,:) = kernel3d('h','d',0.1);
% % kerns(:,:) = kernel3d('h','c',0.1,[1,1]);
% kerns(:,:) = kernel3d('l','c',[1,1]);
zk = 1;
kerns = kernel3d('h','c',zk,[1,1]);
% kerns = kernel3d('h','d',zk);
% kerns = kernel3d([kerns,kerns;kerns,kerns]);
srfrs = [S,S2];
srfrs = S;
tic;
Smat = surfermat(srfrs,kerns,eps);
Smat = Smat - 0.5*eye(size(Smat));

[syscors,novers] = surfermat(srfrs,kerns,eps,struct('corrections',1));
syscors = syscors - 0.5*speye(size(syscors));

[sysnsmth,novers2] = surfermat(srfrs,kerns,eps,struct('nonsmoothonly',1));
sysnsmth = sysnsmth - 0.5*speye(size(sysnsmth));
toc;
%%
src = []; src.r = [2;0;1];
% src = []; src.r = [0;0;0.1];
skern = kernel3d('l','s');
skern = kernel3d('h','s',zk);
smerge = merge(srfrs);
rhs = -skern.eval(src,smerge);

sysapply = @(dens) surfermatapply(srfrs,kerns,dens,eps,novers,syscors);
opts = []; opts.usematlab = 0;
sysapply2 = @(dens) surfermatapply(srfrs,kerns,dens,eps,novers2,sysnsmth,opts);
matlabapply_err = norm(sysapply(rhs)-Smat*rhs);
layerapply_err = norm(sysapply2(rhs)-Smat*rhs);
[matlabapply_err,layerapply_err]
tic;
% sol = Smat\rhs;
sol = gmres(Smat,rhs,[],1e-8,100);
sol2 = gmres(sysapply,rhs,[],1e-8,100);
sol3 = gmres(sysapply2,rhs,[],1e-8,100);
toc;
[norm(sol - sol2)/norm(sol),norm(sol - sol3)/norm(sol)]

%%
kernseval = kerns(1,:);

nplot = 30;
xx = linspace(-2,6,nplot);
yy = 2*linspace(-2,2,nplot);
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).';YY(:).';0*YY(:).'];

% targs.r = [0.95;0;0];
uin = skern.eval(src,targs);


tic;
% no singular quad
uscat = surferkerneval(srfrs,kernseval,sol,targs,eps);

opts = [];opts.corrections = 1;
[evalcors,novers] = surferkernevalmat(srfrs,kernseval,targs,eps,opts);
opts.corrections = evalcors; opts.novers = novers;
uscat2 = surferkerneval(srfrs,kernseval,sol,targs,eps,opts);

evalmat = surferkernevalmat(srfrs,kernseval,targs,eps);
uscat3 = evalmat*sol;

opts2 = []; opts2.usematlab = 0;
uscat4 = surferkerneval(srfrs,kernseval,sol,targs,eps,opts2);

opts = []; opts.nonsmoothonly = 1;
[nsmth,novers] = surferkernevalmat(srfrs,kernseval,targs,eps,opts);
opts2 = []; opts2.usematlab = 0;
opts2.corrections = nsmth; opts2.novers = novers;
uscat5 = surferkerneval(srfrs,kernseval,sol,targs,eps,opts2);

toc
%%
% norm(uscat4 - uscat) / norm(uscat)
% norm(uscat4 - uscat2) / norm(uscat)
% norm(uscat4 - uscat3) / norm(uscat)
% norm(uscat4 - uscat5) / norm(uscat)
uref = uscat4;
[norm(uref - uscat), norm(uref - uscat2), ...
    norm(uref - uscat3), norm(uref - uscat4),norm(uref - uscat5)]/norm(uscat)
% [norm(uscat2 - uscat), norm(uscat3 - uscat), ...
%     norm(uscat4 - uscat), norm(uscat5 - uscat)]/norm(uscat)
% [norm(uin+uscat),norm(uin+uscat2),norm(uin+uscat3),norm(uin+uscat4),norm(uin+uscat5)]

%%
figure(1);clf
% plot(smerge,zeros(smerge.npts,1))
hold on
% scatter3(targs.r(1,:),targs.r(2,:),targs.r(3,:),[],uin+uscat,'.');colorbar
scatter3(targs.r(1,:),targs.r(2,:),targs.r(3,:),[],log10(abs(uin+uscat)),'.');colorbar
hold off

figure(2);clf
plot(S,real(rhs)); colorbar
% 
%%

targest.r = [0.9;0;0];
uintest = skern.eval(src,targest);
uscattest = surferkerneval(srfrs,kernseval,sol,targest,eps);

analytic_sol_err = abs(uintest+uscattest)