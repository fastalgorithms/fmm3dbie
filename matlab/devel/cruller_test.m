%% Test script: functions used in error estimation adaptive mesh refinement
% Single cruller
close all;clear

%%
% tolerances
tol_dens = 1e-6; % tolerance for refining patch based on error in monitor
tol_geom = 1e-5; % tolerance in resolving the geometry. Set less than tol_dens
tol_sol = 1e-7; % tolerance for reference solution
lptype = 2; % lpnorm for adaptive surface mesh generation
eta_surf = 1; % h^eta for quad tree
eta_dens = 0; % h^eta for refinement criterien
ifplot = true; % plot flag
nmonitor = 10; % 3, 9 or 10

dpars = [1,0]; % Single layer potential

% xyz_in(1) =  -0.12652506388863263;
% xyz_in(2) =  -0.16069536550050378;
% xyz_in(3) =  4.2853888370977941d-3;

%xyz_out = [0.1;-0.3;0.2];
xyz_out = [3.17d0,-0.03d0,6.15d0]';
radii(1) = 0.2;
radii(2) = 0.1;
radii(3) = 0.05;
radii(4) = 5.0;
radii(5) = 3.0;

scales  = [5 5 5]';

nu = 3;
nv = 2;
nuv = [nu,nv]; % split (u,v) space into nu X nv parts patches at root level. x2 for triangles
norder = 6; % patch order
iptype_pat = 1;
S = geometries_adap.cruller_adap(radii,scales,nuv,norder,iptype_pat,...
    tol_geom,lptype,eta_surf,nmonitor);
% Ssol = geometries_adap.cruller_adap(radii,scales,nuv,norder,iptype_pat,...
%     tol_sol,lptype,eta_surf,nmonitor);
[Sus,~,~,interp_mat] = split_patches(S,1:S.npatches);

% extract data
[srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);

% create rhs

src_info = [];
src_info.r = xyz_out;

rhs = lap3d.kern(src_info,S,'s');
%rhssol = lap3d.kern(src_info,Ssol,'s');
rhsus = lap3d.kern(src_info,Sus,'s');

sig = lap3d.dirichlet.solver(S,rhs,tol_dens,dpars);
%sigsol = lap3d.dirichlet.solver(Ssol,rhssol,tol_dens,dpars);
sigus = lap3d.dirichlet.solver(Sus,rhsus,tol_dens,dpars);

targ_info = S;
%%
pot = lap3d.dirichlet.eval(S,sig,S,tol_sol,dpars);
%potsol = lap3d.dirichlet.eval(Ssol,sigsol,Ssol,tol_sol,dpars);
potus = lap3d.dirichlet.eval(Sus,sigus,Sus,tol_sol,dpars);
%%
pot_ex = lap3d.kern(src_info,S,'s');
%potsol_ex = lap3d.kern(src_info,Ssol,'s');
potus_ex = lap3d.kern(src_info,Sus,'s');

%%
fprintf('Error in iterative solver=%d\n',sqrt(sum(abs(pot-pot_ex).^2.*S.wts))/sqrt(sum(abs(pot_ex).^2.*S.wts)));
%fprintf('Error in iterative solver, reference=%d\n',sqrt(sum(abs(potsol-potsol_ex).^2.*Ssol.wts))/sqrt(sum(abs(potsol_ex).^2.*Ssol.wts)));
fprintf('Error in iterative solver, upsample=%d\n',sqrt(sum(abs(potus-potus_ex).^2.*Sus.wts))/sqrt(sum(abs(potus_ex).^2.*Sus.wts)));

figure()

S.plot
alpha 0.7
hold on
plot3(src_info.r(1),src_info.r(2),src_info.r(3),'k*')
hold off


%% plot tail of expansion on each patch
patcharea = full(sparse(S.patch_id, 1, S.wts));
dcoefs = vals2coefs(S, sig);
nppatch = S.npts/S.npatches;
dcoefs = reshape(dcoefs, [nppatch, S.npatches]);
dsub = sum(abs(dcoefs(end-norder+1:end,:)),1).*patcharea.';
if ifplot
    figure(4);
    plot(S, log10(abs(dsub.')))
    clim([min(log10(abs(dsub)))-1 max(log10(abs(dsub)))+1])
    colorbar
    title(sprintf('Tail for tolerance = %.3e', tol_dens))
end


% %fun_monitor = cat(1,layerdens',S.mean_curv');
% %fun_monitor_sol = cat(1,layerdenssol',Ssol.mean_curv');
% fun_monitor = sig';
% fun_monitor_sol = sigsol';
% [errl2abs_sig,int2Ssol_sig] = testfuns.interp_lp_error(S,Ssol,fun_monitor,fun_monitor_sol);
% %[errl2abs_onS,fun_int2S] = testfuns.interp_lp_error(Ssol,S,fun_monitor_sol,fun_monitor);
% fprintf('Error in interpolating density =%d\n',errl2abs_sig/sqrt(sum(Ssol.wts.*sigsol.^2)));
% 
% fun_monitor = pot';
% fun_monitor_sol = potsol';
% [errl2abs_pot,int2Ssol_pot] = testfuns.interp_lp_error(S,Ssol,fun_monitor,fun_monitor_sol);
% %[errl2abs_onS,fun_int2S] = testfuns.interp_lp_error(Ssol,S,fun_monitor_sol,fun_monitor);
% fprintf('Error in interpolating potential =%d\n',errl2abs_pot/sqrt(sum(Ssol.wts.*potsol.^2)));
% 
% fun_monitor = S.mean_curv';
% fun_monitor_sol = Ssol.mean_curv';
% [errl2abs_hcurv,int2Ssol_hcurv] = testfuns.interp_lp_error(S,Ssol,fun_monitor,fun_monitor_sol);
% %[errl2abs_onS,fun_int2S] = testfuns.interp_lp_error(Ssol,S,fun_monitor_sol,fun_monitor);
% fprintf('Error in interpolating mean curvature =%d\n',errl2abs_hcurv/sqrt(sum(Ssol.wts.*Ssol.mean_curv.^2)));
% 

%% Using upsampling
sig0us = sig.' * interp_mat.'; 
pot0us = pot.' * interp_mat.';
fprintf('Error in upsampling density =%d\n',sqrt(sum(abs(sig0us.'- sigus).^2.*Sus.wts))/sqrt(sum(Sus.wts.*abs(sigus).^2)));
fprintf('Error in upsampling potential =%d\n',sqrt(sum(abs(pot0us.'- potus).^2.*Sus.wts))/sqrt(sum(Sus.wts.*abs(potus).^2)));
