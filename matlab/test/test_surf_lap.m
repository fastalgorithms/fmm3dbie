run ../startup.m
% sphere(norder,nu,nref)
S = surfer.sphere(6,1,3);

ndeg = 1;
rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
rhs = S.r(3,:)./rr;

surf_lap_f = get_surface_laplacian(S,rhs);
wts = cat(1,S.weights{:});


rfac = -(ndeg+0.0)*(ndeg+1.0);
errf = norm((surf_lap_f - rfac*rhs).*sqrt(wts)');
fprintf('error in surface laplacian=%d\n',errf);


