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

zpars = complex([1e-7, 1.0, 0.0]);

ndeg = 2;
f = spherefun.sphharm(ndeg,1);
rhs = f(S.r(1,:),S.r(2,:),S.r(3,:));

rhs = rhs(:);
eps = 1e-7;
p = helm3d.dirichlet.eval(S,zpars,rhs,eps);

p_ex = rhs/(2*ndeg+1);

err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));

fprintf('Error in single layer potential=%d\n',err1);

%% Now test the solver + kernel evaluation routines


