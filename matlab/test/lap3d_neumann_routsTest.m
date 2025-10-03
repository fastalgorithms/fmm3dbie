%
% This file tests the Helmholtz impedance problem
%
%
%run ../startup.m
S = geometries.sphere(1, 2, [0;0;0], 9);

tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
[~, npts] = size(srcvals);

xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;

eps = 1e-10;

rhs = lap3d.kern(src_info, S, 'sprime');

[sig] = lap3d.neumann.solver(S, rhs, eps);

targ_info = [];
targ_info.r = xyz_out;

pot = lap3d.neumann.eval(S, sig, targ_info, eps);
pot_ex = lap3d.kern(src_info,targ_info,'s');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));

pot = lap3d.neumann.apply_sprime(S,sig,eps);

