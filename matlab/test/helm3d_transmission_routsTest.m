%
% This file tests the Helmholtz impedance problem
%
%
clear
run ../startup.m
S = geometries.sphere(1, 2, [0;0;0], 4, 1);

tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
[~, npts] = size(srcvals);

zks = [1.1; 2.1];
alpha0 = 1;
beta0 = 1;

alpha1 = 1;
beta1 = 1;

rep_params = [alpha0; beta0; alpha1; beta1];

eps = 1e-7;

xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info_in = [];
src_info_in.r = xyz_in;

src_info_out = [];
src_info_out.r = xyz_out;

rhs = complex(zeros(2,npts));
u0 = helm3d.kern(zks(1), src_info_in, S, 's');
dudn0 = helm3d.kern(zks(1), src_info_in, S, 'sprime');

u1 = helm3d.kern(zks(2), src_info_out, S, 's');
dudn1 = helm3d.kern(zks(2), src_info_out, S, 'sprime');

rhs(1,:) = alpha0*u0 - alpha1*u1;
rhs(2,:) = beta0*dudn0 - beta1*dudn1;


[densities] = helm3d.transmission.solver(S, rhs, eps, zks, rep_params);

targ_info = [];
targ_info.r = [xyz_out, xyz_in];

pot = helm3d.transmission.eval(S, densities, targ_info, eps, zks, rep_params);
pot_ex0 = helm3d.kern(zks(1), src_info_in,src_info_out,'s');
fprintf('Error in iterative solver=%d\n',abs(pot(1)-pot_ex0)/abs(pot_ex0));

pot_ex1 = helm3d.kern(zks(2), src_info_out, src_info_in,'s');
fprintf('Error in iterative solver=%d\n',abs(pot(2)-pot_ex1)/abs(pot_ex1));
