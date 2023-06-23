%
% This file tests the Helmholtz single layer
% potenial on the sphere
%
% Additional dependencies: chebfun
%
addpath(genpath('~/git/fmm3dbie/matlab'))
S = surfer.ellipsoid([1,1,1],0.5,6,0);



tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
tic, plot(S); toc;

zk = 1.1;

zpars = complex([zk, 1.0, 0.0]);



ndeg = 2;

jn = sqrt(pi/2/zk)*besselj(ndeg+0.5,zk);
hn = sqrt(pi/2/zk)*besselh(ndeg+0.5,1,zk);
f = spherefun.sphharm(ndeg,1);
rhs = f(S.r(1,:),S.r(2,:),S.r(3,:));

rhs = rhs(:);
eps = 1e-7;
p = helm3d.dirichlet.eval(S,zpars,rhs,eps);

p_ex = rhs*jn*hn*1j*zk;

err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));

fprintf('Error in single layer potential=%d\n',err1);

%% Now test the solver + kernel evaluation routines

xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
rhs = helm3d.kern(zpars(1),src_info,S,'s');


zpars = [zk,-1j*zk,1];
sig = helm3d.dirichlet.solver(S,zpars,rhs,eps);

targ_info = [];
targ_info.r = xyz_out;

dat = helm3d.kern(zpars(1),S,targ_info,'c',zpars(2),zpars(3));



