% Two boundaries, one Dirichlet, one transmission.
ctr1 = [0;0;0];
ctr2 = [3;0;0];

S1 = geometries.ellipsoid([1,1,1.5],3*[1,1,1],ctr1,6);
S2 = geometries.ellipsoid([1,1.5,1],3*[1,1,1],ctr2,6);

zk = 0.1;
zks = [zk,0.5];

% use the representation S[sigma]-D[mu] for transmission
neg2 = diag([1,-1]);

kerns(2,2) = kernel3d();
kerns(1,1) = 2*kernel3d('h','c',zk,[1,1]);
kerns(1,2) = kernel3d('h','trans_rep',zk) .* neg2;
kerns(2,1) = 2*kernel3d('h','c2trans',zk,[1,1]);
kerns(2,2) = kernel3d('h','trans_sys_diff',zks) .* neg2;

eps = 1e-10;
srfrs = [S1,S2];
tic;
Smat = surfermat(srfrs,kerns,eps);
Smat = Smat + eye(size(Smat));
tbuild = toc

rhskerns(2,1) = kernel3d();
rhskerns(1) = kernel3d('h','s',zk);
rhskerns(2) = kernel3d('h','s2trans',zk);
%%
% Build 3 right-hand sides:
%   1 - src at ctr1 to test the Neumann/Dirichlet half
%   2 - src at ctr2 to test the transmission half
%   3 - src at [1;2;0.5]
src = [];
src.r = [ctr1, ctr2, [1;2;0.5]];
nrhs = size(src.r,2);

rhs = zeros(size(Smat,1),nrhs);
inds1 = surferids(srfrs,1);
inds2 = surferids(srfrs,2);
rhs(inds1,:)          = -rhskerns(1).eval(src,S1);
rhs(S1.npts+1:end,:)  = -rhskerns(2).eval(src,S2);


tic;
sol = Smat\rhs;
tsolve = toc

% sol's S2-block second density is lambda' = -lambda.

% Evaluation kernel for targets exterior to S2
kernseval = kerns(1,:);

% Evaluation kernel for targets interior to S2
kernseval2(1,2) = kernel3d();
kernseval2(1) = kernel3d('z');
kernseval2(2) = kernel3d('h','trans_rep',zks(2)) .* neg2;

%%
nplot = 100;
xx = linspace(-2,6,nplot)+1e-2;  yy = 2*linspace(-2,2,nplot)+3e-2;
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).'; YY(:).'; 0*YY(:).'];

in1 = surferinterior(S1,targs);
in2 = surferinterior(S2,targs);
targs_out = []; targs_out.r = targs.r(:,~in2);
targs_in  = []; targs_in.r  = targs.r(:, in2);

tic;
skern = kernel3d('h','s',zk);
uin = skern.eval(src,targs);
uin(in2,:) = 0;
uin(in1,:) = NaN*(1+1i);
uscat = zeros(size(targs.r,2),nrhs);
uscat(~in2,:) = surferkerneval(srfrs,kernseval, sol,targs_out,eps);
uscat( in2,:) = surferkerneval(srfrs,kernseval2,sol,targs_in, eps);
uscat(in1,:) = NaN*(1+1i);
tplot = toc
%%

figure(1);clf
subplot(1,2,1)
h = pcolor(XX,YY,reshape(log10(abs(uin(:,1)+uscat(:,1))),size(XX))); h.EdgeColor = 'none';
colorbar
title('Neumann/Dirichlet test error (src at ctr1)')

subplot(1,2,2)
h = pcolor(XX,YY,reshape(log10(abs(uin(:,2)+uscat(:,2))),size(XX))); h.EdgeColor = 'none';
colorbar
title('Transmission test error (src at ctr2)')

figure(2);clf
h = pcolor(XX,YY,reshape(real(uin(:,3)+uscat(:,3)),size(XX))); h.EdgeColor = 'none';
colorbar
