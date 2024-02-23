%
% This file tests the Helmholtz single layer
% potenial on the sphere
%
% Additional dependencies: chebfun, surface-hps
%
run ../startup.m

% For triangles
S = surfer.ellipsoid([1,1,1],0.5,6,0);

% For quads
S = surfer.sphere(6,1,2, 11);

% For dan meshes

% dom = surfacemesh.sphere(7,2);
% S = surfer.surfacemesh_to_surfer(dom);

tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
rr2 = sum(wts);
rr3 = norm(S.r-S.n);

% tic, plot(S); toc;

zk = 1.1;

zpars = complex([zk, 1.0, 0.0]);
ndeg = 1;

jn = sqrt(pi/2/zk)*besselj(ndeg+0.5,zk);
hn = sqrt(pi/2/zk)*besselh(ndeg+0.5,1,zk);

rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
rhs = S.r(3,:)./rr;

rhs = rhs(:);
eps = 1e-7;


p = helm3d.dirichlet.eval(S,zpars,rhs,eps);

p_ex = rhs*jn*hn*1j*zk;

err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));

fprintf('Error in single layer potential=%d\n',err1);

%% Now test eval routine with precomputed quadrature corrections
opts_quad = [];
opts_quad.format='rsc';
Q = helm3d.dirichlet.get_quadrature_correction(S,zpars, ...
   eps,S,opts_quad);
p = helm3d.dirichlet.eval(S,zpars,rhs,eps,S,Q);
err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));

fprintf('Error in single layer potential with precomp corr=%d\n',err1);


%% Now test the solver + kernel evaluation routines

xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
rhs = helm3d.kern(zk,src_info,S,'s');


zpars = [zk,-1j*zk,1];
sig = helm3d.dirichlet.solver(S,zpars,rhs,eps);

targ_info = [];
targ_info.r = xyz_out;

dat = helm3d.kern(zk,S,targ_info,'c',zpars(2),zpars(3));

pot = dat*(sig.*wts);
pot_ex = helm3d.kern(zk,src_info,targ_info,'s');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));

% Test matrix entry generator or fast direct solver depending on size of linear system
P = zeros(S.npts,1);
opts_quad = [];
opts_quad.format='sparse';
Q = helm3d.dirichlet.get_quadrature_correction(S,zpars, ...
   eps,S,opts_quad);

A = helm3d.dirichlet.matgen(1:S.npts,1:S.npts,S,zpars,P,Q); 
sig_ds = A\rhs;
pot_ds = dat*(sig_ds.*wts);
    
fprintf('Error in direct solver=%d\n',abs(pot_ds-pot_ex)/abs(pot_ex));

