addpath(genpath('~/git/fmm3dbie/matlab'))

% sphere(norder,nu,nref)
S = surfer.sphere(6,1,3);

ndeg = 3;
f = spherefun.sphharm(ndeg,1);
rhs = f(S.r(1,:),S.r(2,:),S.r(3,:));

surf_lap_f = get_surface_laplacian(S,rhs);
wts = cat(1,S.weights{:});


rfac = -(ndeg+0.0)*(ndeg+1.0);
errf = norm((surf_lap_f - rfac*rhs).*sqrt(wts)');
fprintf('error in surface laplacian=%d\n',errf);


