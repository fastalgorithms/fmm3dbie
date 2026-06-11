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
Ssol = merge([S1sol, S2sol]);
% extract data
[srcvalssol,srccoefssol,norderssol,ixyzssol,iptypesol,wtssol] = extract_arrays(Ssol);
% assign different rhs: [0,1]
rhs1sol = ones(S1sol.npts,1);
rhs2sol = zeros(S2sol.npts,1);
rhssol = [rhs1sol; rhs2sol];
% plot spheres
if ifplot
    plot(S)
end
%% error estimate based on mean curvature H
Hext = 1/a; % analytic formula for H
H = S.mean_curv; % approximation
Hdiff = abs(H - Hext);
% error in Lp norm
if lptype == 0
    errLp = max(Hdiff);
elseif lptype == 1
    errLp = sum(abs(Hdiff).*S.wts)/S.area;
elseif lptype == 2
    errLp = sqrt(sum(abs(Hdiff).^2.*S.wts))/sqrt(S.area);
end

fprintf('error in resolving curvature in L%d norm: %d\n', lptype,errLp);
%% solve a BIE on torii with rhs equal to 0 or 1
rep_pars = [1.0; 1.0];
[layerdens, errs] = lap3d.solver(S, 'dir', rhs, 1e-8, rep_pars);
%% plot layer densities and rhs
if ifplot
    figure(2)
    clf
    plot(S, layerdens)
    title(" layer density ")
    colorbar
    figure(3)
    clf
    plot(S, rhs)
    title("right hand side")
end
%% plot tail of expansion on each patch
dcoefs = vals2coefs(S, S.r(3,:));
nppatch = S.npts/S.npatches;
dcoefs = reshape(dcoefs, [nppatch, S.npatches]);
dsub = sum(abs(dcoefs(end-norder+1:end,:)),1);
if ifplot
    figure(4);
    plot(S, log10(abs(dsub.')))
    clim([min(log10(abs(dsub)))-1 max(log10(abs(dsub)))+1])
    colorbar
    title(sprintf('Tail for tolerance = %.3e', tol_dens))
end
%% compute reference solution
% solve a BIE on torii with rhs equal to 0 or 1
if compref

[layerdenssol, errs] = lap3d.solver(Ssol, 'dir', rhssol, tol_sol, rep_pars);
nppatchsol = Ssol.npts/Ssol.npatches;

dcoefs_sol = vals2coefs(Ssol, layerdenssol);
dcoefs_sol = reshape(dcoefs_sol, [nppatchsol, Ssol.npatches]);
dsub_solsol = sum(abs(dcoefs_sol(end-nordersol+1:end,:)),1);
if ifplot
    figure(5);
    plot(Ssol, dsub_solsol.')
    colorbar
    title(sprintf('Tail for tolerance = %.3e', tol_sol))
end
elseif loadref
    load("ref_solution.m","Ssol","rhssol","tol_sol","layerdenssol", "errs")
else
    error("No reference solution")
end    
%% compute error in polarization charge
x = S1.r(1,:).';
y = S1.r(2,:).';
z = S1.r(3,:).';
c1 = sum(layerdens(1:S1.npts).*S1.wts.*x);

xsol = S1sol.r(1,:).';
%% 
ysol = S1sol.r(2,:).';
c1_sol = sum(layerdenssol(1:S1sol.npts).*S1sol.wts.*xsol);

err1 = abs(c1 - c1_sol)./(abs(c1));
fprintf('error in polarization charge on obstacle 1=%d\n', err1);
%% Compute pointwise error in layer density at quadrature nodes S.r
% interpolate on Ssol to S.r
% fun_monitor = cat(1,layerdens',srcvals(1:3,:),S.mean_curv');
% fun_monitor_sol = cat(1,layerdenssol',srcvalssol(1:3,:),Ssol.mean_curv');

%fun_monitor = cat(1,layerdens',S.mean_curv');
%fun_monitor_sol = cat(1,layerdenssol',Ssol.mean_curv');

fun_monitor = layerdens';
fun_monitor_sol = layerdenssol';

%fun_monitor = S.mean_curv';
%fun_monitor_sol = Ssol.mean_curv';

nd = size(fun_monitor,1);
[errl2abs_onSsol,fun_int2Ssol] = testfuns.interp_lp_error(S,Ssol,fun_monitor,fun_monitor_sol);
[errl2abs_onS,fun_int2S] = testfuns.interp_lp_error(Ssol,S,fun_monitor_sol,fun_monitor);
%% Get errors in polynomial approximation of densities
if ifplot
    interp_errs = S.surf_fun_error(fun_int2S - fun_monitor,lptype);
    clim_max = max(max(log10(abs(interp_errs(1,:)))),max(log10(abs(interp_errs(end,:))))) + 1;
    clim_min = min(min(log10(abs(interp_errs(1,:)))),min(log10(abs(interp_errs(end,:))))) - 1;
    figure(6)
    sgtitle("Patchwise magnitute of interpolation error for: Original mesh")
    subplot(121)
    plot(S,log10(abs(interp_errs(1,:))))
    colorbar
    clim([clim_min clim_max])
    title("Layer density")
    subplot(122)
    plot(S,log10(abs(interp_errs(end,:))))
    colorbar
    clim([clim_min clim_max])
    title("Mean curvature")
end
%% Get errors in polynomial approximation of densities
tail_errs = S.surf_fun_error(fun_monitor,lptype);

if ifplot
    clim_max = max(max(log10(abs(tail_errs(1,:)))),max(log10(abs(tail_errs(end,:))))) + 1;
    clim_min = min(min(log10(abs(tail_errs(1,:)))),min(log10(abs(tail_errs(end,:))))) - 1;
    figure(7)
    sgtitle("Patchwise magnitute of polynomial tails for: Original mesh")
    subplot(121)
    plot(S,log10(abs(tail_errs(1,:))))
    colorbar
    clim([clim_min clim_max])
    title("Layer density")
    subplot(122)
    plot(S,log10(abs(tail_errs(end,:))))
    colorbar
    clim([clim_min clim_max])
    title("Mean curvature")
end
%% Find patches to refine based on interpolation error
fun_diff = fun_int2S - fun_monitor;
[errlp_patch] = testfuns.get_patch_lp_norm(S,fun_diff,[],lptype);
funlpS1 = sum(abs(fun_monitor(:,1:S1.npts)).^lptype .* S1.wts',2).^(1/lptype);
funlpS2 = sum(abs(fun_monitor(:,1+S1.npts:end)).^lptype .* S2.wts',2).^(1/lptype);

%[list] = find_patch_refine_list(S1,errlp_patch(:,1:S1.npts),funlp1,tol_dens,[],lptype,eta_dens);
eta_dens = 1;
[refinelist1,errqlist1] = find_patch_refine_list(S1,errlp_patch(:,1:S1.npatches),funlpS1,tol_dens,[],lptype,eta_dens);
[refinelist2,errqlist2] = find_patch_refine_list(S2,errlp_patch(:,1+S1.npatches:end),funlpS2,tol_dens,[],lptype,eta_dens);

[errq_desc1,ind_errq1] = sort(errqlist1,'descend');
[errq_desc2,ind_errq2] = sort(errqlist2,'descend');

refinelist1_orig = refinelist1;
refinelist2_orig = refinelist2;


nrefine_1 = round(numel(ind_errq1) * ref_prcent);
nrefine_2 = round(numel(ind_errq2) * ref_prcent);

refinelist1 = refinelist1(ind_errq1(1:nrefine_1));
refinelist2 = refinelist2(ind_errq2(1:nrefine_2));

refinelist1 = 1:S1.npatches;
refinelist2 = 1:S2.npatches;

% plot patches to be refined
if ifplot
    plotlist1 = zeros(S1.npatches,1);
    plotlist1(refinelist1) = 1;
    plotlist2 = zeros(S2.npatches,1);
    plotlist2(refinelist2) = 1;
    
    figure(7)
    plot(S1,plotlist1)
    hold on
    plot(S2,plotlist2)
    hold off
    colorbar
    title("Patches flagged for refinement")
end

%% Refine patches
[S1refine,ipatchid_1,~,interp_1refine_mat,~,~,~,ipatsplit_1] = split_patches(S1,refinelist1);
[S2refine,ipatchid_2,~,interp_2refine_mat,~,~,~,ipatsplit_2] = split_patches(S2,refinelist2);
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

err1 = abs(c1refine - c1_sol)./(abs(c1));
fprintf('error in polarization charge on obstacle 1=%d\n', err1);

%
% densities on non-refined patches and their update
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

densitydiff = abs(layerdens_refine_1_no_ref-layerdens_no_ref_orig);

densitydiff_plot = NaN(S1.npts,1);

densitydiff_plot(no_split_pnts_idx) = densitydiff;

figure()
plot(S1,log10(densitydiff_plot))
colorbar


tail_errs = Srefine.surf_fun_error(layerdens_refine,lptype);

if ifplot
    clim_max = max(max(log10(abs(tail_errs(1,:)))),max(log10(abs(tail_errs(end,:))))) + 1;
    clim_min = min(min(log10(abs(tail_errs(1,:)))),min(log10(abs(tail_errs(end,:))))) - 1;
    figure()
    sgtitle("Patchwise magnitute of polynomial tails for: Original mesh")
    subplot(121)
    plot(Srefine,log10(abs(tail_errs(1,:))))
    colorbar
    clim([clim_min clim_max])
    title("Layer density")
    subplot(122)
    plot(Srefine,log10(abs(tail_errs(end,:))))
    colorbar
    clim([clim_min clim_max])
    title("Mean curvature")
end
return

