%
% This file tests the Stokes single layer
% potenial on the sphere
%
%
run ../startup.m

S = geometries.sphere(1, 2, [0;0;0], 4, 1);
tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;


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
targ_info.patch_id = -1;
dat = stok3d.kern(S,targ_info,'c',dpars(1),dpars(2));
pot = sum(dat.*(reshape(sig.*wts.',[1,1,3,S.npts])),[3,4]);
pot_eval = stok3d.eval(S,dpars,sig,eps,targ_info);
pot_ex = sum(stok3d.kern(src_info,targ_info,'s').*reshape(src_info.v,[1,1,3]),3);
fprintf('Error in iterative solver=%d\n',vecnorm(pot-pot_ex)/vecnorm(pot_ex));
fprintf('Error in stok3d_eval=%d\n',vecnorm(pot_eval-pot_ex)/vecnorm(pot_ex));




