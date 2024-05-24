%
% This file tests the Laplace single layer
% potenial on the sphere
%
% Additional dependencies: chebfun, surface-hps
%
run ../startup.m

% For triangles
S = surfer.ellipsoid([1,1,1],0.5,6,0);

% For quads
S = surfer.sphere(6, 1, 2, 11);

% For dan meshes

% dom = surfacemesh.sphere(7,2);
% S = surfer.surfacemesh_to_surfer(dom);

tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
rr2 = sum(wts);
rr3 = norm(S.r-S.n);

% tic, plot(S); toc;

dpars = [1.0, 0.0];
ndeg = 1;


rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
rhs = S.r(3,:)./rr;

rhs = rhs(:);
eps = 1e-7;


p = lap3d.dirichlet.eval(S,dpars,rhs,eps);

p_ex = rhs/(2*ndeg+1);

err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));

fprintf('Error in single layer potential=%d\n',err1);

%% Now test eval routine with precomputed quadrature corrections
opts_quad = [];
opts_quad.format='rsc';
Q = lap3d.dirichlet.get_quadrature_correction(S,dpars, ...
   eps,S,opts_quad);
p = lap3d.dirichlet.eval(S,dpars,rhs,eps,S,Q);
err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));

fprintf('Error in single layer potential with precomp corr=%d\n',err1);


%% Now test the solver + kernel evaluation routines

xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
rhs = lap3d.kern(src_info,S,'s');


dpars = [1,1];
sig = lap3d.dirichlet.solver(S,dpars,rhs,eps);

targ_info = [];
targ_info.r = xyz_out;

dat = lap3d.kern(S,targ_info,'c',dpars(1),dpars(2));

pot = dat*(sig.*wts);
pot_ex = lap3d.kern(src_info,targ_info,'s');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));
