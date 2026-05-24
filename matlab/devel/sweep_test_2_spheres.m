
% percentage of flagged patches to refine: 1 - all patches, 0.5 - 50%
ref_prcnt_in_vec = 0.05:0.05:1;
n_ref_prcnt = numel(ref_prcnt_in_vec);

sphere_rad = 1.1; % radius of spheres
sepdist = 0.001; % offset between spheres

% generate reference solution
tol_sol = 1e-10; % tolerance for reference solution
nasol=  10;
nordersol = 12;
S1sol = geometries.sphere(sphere_rad, nasol, [], nordersol, 1);
nptssol1 = S1sol.npts;
c0_2 = [2*sphere_rad + sepdist; 0; 0];
S2sol = geometries.sphere(sphere_rad, nasol, c0_2, nordersol, 1);
% merge
Ssol = merge([S1sol, S2sol]);
% extract data
[srcvalssol,srccoefssol,norderssol,ixyzssol,iptypesol,wtssol] = extract_arrays(Ssol);
% assign different rhs: [0,1]
rhs1sol = ones(S1sol.npts,1);
rhs2sol = zeros(S2sol.npts,1);
rhssol = [rhs1sol; rhs2sol];
rep_pars = [1.0; 1.0];
[layerdenssol, errs] = lap3d.solver(Ssol, 'dir', rhssol, tol_sol, rep_pars);


% allocate error vectors for outputs
err_H_Lp_vec = zeros(n_ref_prcnt,1); % error in resolving mean curvature
err_orig_vec = zeros(n_ref_prcnt,1); % error in density*x*wts before ref
err_ref_vec = zeros(n_ref_prcnt,1); % error in density*x*wts after ref
dens_diff_linf_vec = zeros(n_ref_prcnt,1); % difference in layer density in non-refined patches after new solve, in L_inf norm
dens_diff_l2_vec = zeros(n_ref_prcnt,1); % difference in layer density in non-refined patches after new solve, in L_2 norm

% loop over percentages ref_prcnt_in_vec: percentage of flagged patches to
% refine

for ii = 1:n_ref_prcnt
ref_prcnt_in = ref_prcnt_in_vec(ii);
[err_H_Lp_vec(ii),err_orig_vec(ii),err_ref_vec(ii), dens_diff_linf_vec(ii),dens_diff_l2_vec(ii)] = ...
    test_2_spheres(ref_prcnt_in,sphere_rad, sepdist,Ssol,nptssol1,layerdenssol);
end

function [err_H_Lp,err_orig,err_ref, dens_diff_linf,dens_diff_l2] = test_2_spheres(ref_prcnt_in,sphere_rad,sepdist,Ssol,nptssol1,layerdenssol)
%%
% tolerances
tol_dens = 1e-5; % tolerance for refining patch based on error in monitor
tol_solve = 1e-10; % tolerance in resolving the geometry. Set less than tol_dens
tol_sol = 1e-10; % tolerance for reference solution
lptype = 2; % lpnorm for adaptive surface mesh generation
eta_dens = 1; % scaling of h^eta_dens in error checks for refinement

ref_prcnt = ref_prcnt_in; % Percentage of number of flagged patches to refine

% % Two spheres: test to get tol_dens digits of accuracy in the solution
a = sphere_rad;
d = sepdist;
c0_2 = [2*a + d; 0; 0];
na = 8;
%nuv = [2,1];
norder = 8; % patch order
npols = (norder+1)*(norder+2)/2; % points per patch

% generate spheres
S1 = geometries.sphere(a, na, [], norder, 1);
S2 = geometries.sphere(a, na, c0_2, norder, 1);

% merge meshes
S = merge([S1, S2]);
% assign different rhs: [0,1]
rhs1 = ones(S1.npts,1);
rhs2 = zeros(S2.npts,1);
rhs = [rhs1; rhs2];

%% error estimate based on mean curvature H
Hext = 1/a; % analytic formula for H
H = S.mean_curv; % approximation
Hdiff = abs(H - Hext);
% error in Lp norm
if lptype == 0
    err_H_Lp = max(Hdiff);
elseif lptype == 1
    err_H_Lp = sum(abs(Hdiff).*S.wts)/S.area;
elseif lptype == 2
    err_H_Lp = sqrt(sum(abs(Hdiff).^2.*S.wts))/sqrt(S.area);
end

fprintf('error in resolving curvature in L%d norm: %d\n', lptype,err_H_Lp);
%% solve a BIE on torii with rhs equal to 0 or 1
rep_pars = [1.0; 1.0];
[layerdens, errs] = lap3d.solver(S, 'dir', rhs, tol_solve, rep_pars);
  
%% compute error in polarization charge
x = S1.r(1,:).';
y = S1.r(2,:).';
z = S1.r(3,:).';
c1 = sum(layerdens(1:S1.npts).*S1.wts.*x);

xsol = Ssol.r(1,1:nptssol1).';
c1_sol = sum(layerdenssol(1:nptssol1).*Ssol.wts(1:nptssol1).*xsol);

err_orig = abs(c1 - c1_sol)./(abs(c1));
fprintf('Original: error in polarization charge on obstacle 1=%d\n', err_orig);
%% Compute pointwise error in layer density at quadrature nodes S.r

fun_monitor = layerdens';
fun_monitor_sol = layerdenssol';

%nd = size(fun_monitor,1);
%[errl2abs_onSsol,fun_int2Ssol] = testfuns.interp_lp_error(S,Ssol,fun_monitor,fun_monitor_sol);
[errl2abs_onS,fun_int2S] = testfuns.interp_lp_error(Ssol,S,fun_monitor_sol,fun_monitor);

%% Find patches to refine based on interpolation error
fun_diff = fun_int2S - fun_monitor;
[errlp_patch] = testfuns.get_patch_lp_norm(S,fun_diff,[],lptype);
funlpS1 = sum(abs(fun_monitor(:,1:S1.npts)).^lptype .* S1.wts',2).^(1/lptype);
funlpS2 = sum(abs(fun_monitor(:,1+S1.npts:end)).^lptype .* S2.wts',2).^(1/lptype);

% flag patches for refinement
[refinelist1,errqlist1] = find_patch_refine_list(S1,errlp_patch(:,1:S1.npatches),funlpS1,tol_dens,[],lptype,eta_dens);
[refinelist2,errqlist2] = find_patch_refine_list(S2,errlp_patch(:,1+S1.npatches:end),funlpS2,tol_dens,[],lptype,eta_dens);

% Sort errors for patches where refinement test failed
[errq_desc1,ind_errq1] = sort(errqlist1,'descend');
[errq_desc2,ind_errq2] = sort(errqlist2,'descend');

refinelist1_orig = refinelist1;
refinelist2_orig = refinelist2;

% Pick ref_prcnt of the flagged patches
nrefine_1 = round(numel(ind_errq1) * ref_prcnt);
nrefine_2 = round(numel(ind_errq2) * ref_prcnt);

refinelist1 = refinelist1(ind_errq1(1:nrefine_1));
refinelist2 = refinelist2(ind_errq2(1:nrefine_2));


%% Refine patches
[S1refine,ipatchid_1,~,interp_1refine_mat,~,~,~,~] = split_patches(S1,refinelist1);
[S2refine,ipatchid_2,~,interp_2refine_mat,~,~,~,~] = split_patches(S2,refinelist2);
% set boundary conditions
rhs_1refine = interp_1refine_mat * rhs1; % == 1
rhs_2refine = interp_2refine_mat * rhs2; % == 0
rhs_refine = [rhs_1refine;rhs_2refine];
% merge patches
Srefine = merge([S1refine, S2refine]);

%% Solve for density on refined surface
[layerdens_refine, errs] = lap3d.solver(Srefine, 'dir', rhs_refine, tol_sol, rep_pars);

%% Compute error in layer density

layerdens_refine_1 = layerdens_refine(1:S1refine.npts);
layerdens_refine_2 = layerdens_refine(S1refine.npts+1:end);

x = S1refine.r(1,:).';

c1refine = sum(layerdens_refine(1:S1refine.npts).*S1refine.wts.*x);

err_ref = abs(c1refine - c1_sol)./(abs(c1));
fprintf('Refined: error in polarization charge on obstacle 1=%d\n', err_ref);

%% Compute update of non-refined patches after solve
% densities on non-refined patche
layerdens_refine_1_no_ref = layerdens_refine_1(~ipatchid_1(:,2));

mask_1 = ipatchid_1(~ipatchid_1(:,2),1);
no_split_patches = unique(mask_1);
n_no_split_patches = numel(no_split_patches);

no_split_pnts_idx = zeros(npols*n_no_split_patches,1);

[~,~,~,ixyzsrefine_1,~,~] = extract_arrays(S1refine);

for ipat = 1:n_no_split_patches
    ipatch = no_split_patches(ipat);
    idx_populate = ((ipat-1)*npols+1):(ipat*npols);
    no_split_pnts_idx(idx_populate) = ixyzsrefine_1(ipatch):(ixyzsrefine_1(ipatch+1)-1);
end

layerdens_no_ref_orig = layerdens(no_split_pnts_idx);
% pointwise difference in update of layer density for non-refined patches
densitydiff = abs(layerdens_refine_1_no_ref-layerdens_no_ref_orig);

dens_diff_linf = max(densitydiff(:));
dens_diff_l2 = sqrt(sum(densitydiff.^2.*S.wts(no_split_pnts_idx)))/sqrt(S.area);


return

end