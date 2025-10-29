%
% This file tests the Maxwell pec 
%
%
run ../startup.m
S = geometries.sphere(1, 2, [0;0;0], 4, 1);
%%
tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
[~, npts] = size(srcvals);

zk = 1.1;
alpha = 1.0;


eps = 1e-7;

xyz_in = [0.17; 0.23; -0.11];
xyz_out = [3.1; 1.0; 1.0];
src_info = [];
src_info.r = xyz_out;
src_info.edips = rand(3,1) + 1j*rand(3,1);
src_info.hdips = rand(3,1) + 1j*rand(3,1);

[Einc, Hinc] = em3d.incoming_sources(zk, src_info, S, 'ehd');

%

zpars = complex([zk, alpha]);
opts = [];
opts.eps_gmres = 1e-5;
[densities] = em3d.pec.solver(S, Einc, Hinc, eps, zk, alpha, opts);

%
targ_info = [];
targ_info.r = xyz_in;

[E, H] = em3d.pec.eval(S, densities, targ_info, eps, zk, alpha);


[E_ex, H_ex] = em3d.incoming_sources(zk, src_info, targ_info, 'ehd');

u_comp = [E; H];
u_ex = [E_ex; H_ex];
utest = u_comp + u_ex;


fprintf('Error in fields iterative solver=%d\n',norm(utest(:))/norm(u_ex(:)));

