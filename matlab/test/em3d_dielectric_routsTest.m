%
% This file tests the Maxwell pec 
%
%
S = geometries.sphere(1, 2, [0;0;0]);

tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
[~, npts] = size(srcvals);

om = 1.1;

ep0 = 1.4;
mu0 = 1;
ep1 = 1;
mu1 = 1;
rep_params = complex([ep0, mu0, ep1, mu1]);

zkin = om*sqrt(ep1*mu1);
zkout = om*sqrt(ep0*mu0);


eps = 1e-6;

xyz_in = [0.17; 0.23; -0.11];
xyz_out = [3.1; 1.0; 1.0];
src_info = [];
src_info.r = xyz_out;
src_info.edips = rand(3,1) + 1j*rand(3,1);
src_info.hdips = rand(3,1) + 1j*rand(3,1);


[Einc, Hinc] = em3d.incoming_sources(zkin, src_info, S, 'ehd');

%
opts = [];
opts.eps_gmres = 1e-8;
[densities, errs, rres, Q] = em3d.dielectric.solver(S, Einc, Hinc, eps, om, rep_params, opts);

%%
%
targ_info = [];
targ_info.r = xyz_in;

[E, H] = em3d.dielectric.eval(S, densities, targ_info, eps, om, rep_params);


[E_ex, H_ex] = em3d.incoming_sources(zkin, src_info, targ_info, 'ehd');

u_comp = [E; H];
u_ex = [E_ex; H_ex];
utest = u_comp + u_ex;


fprintf('Error in fields iterative solver=%d\n',norm(utest(:))/norm(u_ex(:)));
