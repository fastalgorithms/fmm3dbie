%
% Tests surfermatgmres by solving a mixed boundary value
% problem on two spheres using the representation
%
% u = D_1[sigma] + G(t,0) * int sigma + S_2[mu]
%
run ../startup.m

S1 = geometries.sphere(1, 2, [0;0;0], 8, 1);
S2 = S1 + [4;0;0];
srfrs = [S1, S2];

eps = 1e-9;

src0 = []; src0.r = [0;0;0];
srcs = []; srcs.r = [[0.3; -0.1; 0.2],[4; 0.1; -0.3]];

skern = kernel3d('l', 's');
spkern = kernel3d('l', 'sp');

rhs1 = skern.eval(srcs, S1);
rhs2 = spkern.eval(srcs, S2);
rhs = [rhs1; rhs2];

area1 = sum(S1.wts(:));

kerns(2,2) = kernel3d();
kerns(1,1) = kernel3d('l', 'd') + (@(t) skern.eval(src0,t)) .* kernel3d.ones;
kerns(1,2) = -kernel3d('l', 's');
kerns(2,1) =  kernel3d('l', 'dp') + (@(t) spkern.eval(src0,t)) .* kernel3d.ones;
kerns(2,2) = -kernel3d('l', 'sp');

opts_corr = [];
opts_corr.corrections = 1;
[cors, objover] = surfermat(srfrs, kerns, eps, opts_corr);
cors = cors + 0.5*speye(size(cors,1));

[sigma, flag, relres, iters] = surfermatgmres(srfrs, kerns, rhs, eps, 1e-8, 100, objover, cors);

for i = 1:size(rhs,2)
fprintf('flag=%d, relres=%5.2e, iters=%d \n', flag(i), relres(i), iters(i));
end
assert(all(flag == 0), 'gmres did not converge');

targ1 = []; targ1.r = [1.7; 2.1; -0.9];
pot = surferkerneval(srfrs, kerns(1,:), sigma, targ1, eps);
pot_ex = skern.eval(srcs, targ1);

err = abs(pot - pot_ex)./abs(pot_ex);
fprintf('Error in exterior Dirichlet solve = %d\n', err(1));
assert(err(1) < 1e-6, 'src 1 error too large');
fprintf('Error in exterior Dirichlet solve = %d\n', err(2));
assert(err(2) < 1e-6, 'src 2 error too large');