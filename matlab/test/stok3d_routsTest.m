%
% This file tests the Stokes single layer
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

% % tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
% % rr2 = sum(wts);
% % rr3 = norm(S.r-S.n);
% % 
% % % tic, plot(S); toc;
% % 
% % dpars = [1.0, 0.0];
% % ndeg = 1;
% % 
% % 
% % rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
% % rhs = S.r(3,:)./rr;
% % 
% % rhs = rhs(:);
% % rhs3 = repmat(rhs.',3)
% % eps = 1e-7;
% % 
% % 
% % p = stok3d.eval(S,dpars,rhs,eps);
% % 
% % p_ex = rhs/(2*ndeg+1);
% % 
% % err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));
% % 
% % fprintf('Error in single layer potential=%d\n',err1);
% % 
% % %% Now test eval routine with precomputed quadrature corrections
% % opts_quad = [];
% % opts_quad.format='rsc';
% % Q = stok3d.get_quadrature_correction(S,dpars, ...
% %    eps,S,opts_quad);
% % p = stok3d.eval(S,dpars,rhs,eps,S,Q);
% % err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));
% % 
% % fprintf('Error in single layer potential with precomp corr=%d\n',err1);


%% Now test the solver + kernel evaluation routines

xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
src_info.v = [1;-2;3];
rhs = sum(stok3d.kern(src_info,S,'s').*reshape(src_info.v,[1,1,3,1]),3);  


dpars = [1,1];
sig = stok3d.solver(S,dpars,rhs,eps);

targ_info = [];
targ_info.r = xyz_out;

dat = stok3d.kern(S,targ_info,'c',dpars(1),dpars(2));

pot = sum(dat.*(reshape(sig.*wts.',[1,1,3,S.npts])),[3,4]);
pot_ex = sum(stok3d.kern(src_info,targ_info,'s').*reshape(src_info.v,[1,1,3]),3);
fprintf('Error in iterative solver=%d\n',vecnorm(pot-pot_ex)/vecnorm(pot_ex));

