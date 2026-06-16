%% Test script: functions used in error estimation adaptive mesh refinement

close all;clear

%%
% tolerances
tol_dens = 1e-5; % tolerance for refining patch based on error in monitor
tol_geom = 1e-10; % tolerance in resolving the geometry. Set less than tol_dens
tol_sol = 1e-10; % tolerance for reference solution
lptype = 2; % lpnorm for adaptive surface mesh generation
eta_surf = 1; % h^eta for quad tree
ifplot = true; % plot flag
compref = true; % compute reference solution
loadref = false; % load reference solution
nd_monitor = 3;
%iptype = 1;

ref_prcent = 0.25; % Percentage of number of flagged patches to refine

% % Two spheres: test to get tol_dens digits of accuracy in the solution
a = 1.1;
d = 0.001;
c0_1 = zeros(3,1);
c0_2 = [2*a + d; 0; 0];
abc = a * ones(3,1);
na = 8;
%nuv = [2,1];
norder = 8; % patch order
npols = (norder+1)*(norder+2)/2; % points per patch
iptype_pat = 1;

% generate torii
S1 = geometries.sphere(a, na, [], norder, 1);
S2 = geometries.sphere(a, na, c0_2, norder, 1);

% merge meshes
S = merge([S1, S2]);
% extract data
[srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
% assign different rhs: [0,1]
rhs1 = ones(S1.npts,1);
rhs2 = zeros(S2.npts,1);
rhs = [rhs1; rhs2];
% generate reference solution
nasol=  10;
nordersol = 12;
S1sol = geometries.sphere(a, nasol, [], nordersol, 1);
S2sol = geometries.sphere(a, nasol, c0_2, nordersol, 1);
% merge
S = merge([S1sol, S2sol]);
% extract data

% compute u and dudn due to a point source outside
sout = [2.2;4.1;3.1];
srcinfo = [];
src_info.r = sout;
uin = lap3d.kern(src_info, S, 's');
dudnin = lap3d.kern(src_info, S, 'sprime');

dpars_d = [0.0, 1.0];
Du = lap3d.dirichlet.eval(S, uin, S, tol_sol, dpars_d);

dpars_s = [1.0, 0.0];
Sdudn = lap3d.dirichlet.eval(S, dudnin, S, tol_sol, dpars_s);

err_l2 = norm((uin/2 - Sdudn + Du).*sqrt(S.wts(:)));
fprintf('err = %d\n', err_l2)
