%% Test script: functions used in error estimation adaptive mesh refinement

close all;clear

%%
% tolerances
tol_dens = 1e-5; % tolerance for refining patch based on error in monitor
tol_geom = 1e-7; % tolerance in resolving the geometry. Set less than tol_dens
tol_sol = 1e-8; % tolerance for reference solution
lptype=2; % lpnorm for adaptive surface mesh generation
eta_surf=1; % h^eta for quad tree
eta_dens = 0;
ifplot = true; % plot flag

% Two torii test to get tol_dens digits of accuracy in the solution
rmajor = 2;
rminor = 0.5;
rwave = 0.0;

radii = [rminor,rmajor,rwave];
nosc = 1; % number of oscillations (nosc = 1 default)
scales = [1,1,1]; % torus settings. Default

R = rmajor+rminor; % offset
d = -rminor * 4.5; % min distance between torii
c2 = [2*R + d; 0; 0]; % center of 2nd torus, S2
rot2 = [0,pi/2,0]; % rotate 2nd torus 90 deg around y-axis
nu = 2;
nv = 1;
nuv = [nu,nv]; % split (u,v) space into nu X nv parts patches at root level. x2 for triangles
norder = 6; % patch order
iptype_pat = 1;

% generate torii
S1 = geometries_adap.wobblytorus_adap(radii,nosc,scales,nuv,norder,iptype_pat,tol_geom,lptype,eta_surf);
S2 = S1.affine_transf(eye(3),c2);
S2 = S2.rotate(rot2);
% merge meshes
S = merge([S1, S2]);
% extract data
[srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
% assign different rhs: [0,1]
rhs1 = ones(S1.npts,1);
rhs2 = zeros(S2.npts,1);
rhs = [rhs1; rhs2];
% generate reference solution
nuvsol = [16,16];
S1sol = geometries_adap.wobblytorus_adap(radii,nosc,scales,nuvsol,norder,iptype_pat,tol_sol,lptype,eta_surf);
S2sol = S1sol.affine_transf(eye(3),c2);
S2sol = S2sol.rotate(rot2);
% merge
Ssol = merge([S1sol, S2sol]);
% extract data
[srcvalssol,srccoefssol,norderssol,ixyzssol,iptypesol,wtssol] = extract_arrays(Ssol);
% assign different rhs: [0,1]
rhs1sol = ones(S1sol.npts,1);
rhs2sol = zeros(S2sol.npts,1);
rhssol = [rhs1sol; rhs2sol];
% plot torii
if ifplot
    plot(S)
end
%% error estimate based on mean curvature H
rho = cat(2,sqrt(S1.r(1,:).^2+S1.r(2,:).^2),sqrt((S2.r(1,:)-c2(1)).^2 + (S2.r(3,:)-c2(3)).^2)).';
Hext = (2*rho-rmajor)./(2*rminor*rho); % analytic formula for H
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
dcoefs = vals2coefs(S, layerdens);
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
[layerdenssol, errs] = lap3d.solver(Ssol, 'dir', rhssol, tol_sol, rep_pars);

dcoefs_sol = vals2coefs(Ssol, layerdenssol);
dcoefs_sol = reshape(dcoefs_sol, [nppatch, Ssol.npatches]);
dsub_solsol = sum(abs(dcoefs_sol(end-norder+1:end,:)),1);
if ifplot
    figure(5);
    plot(Ssol, dsub_solsol.')
    colorbar
    title(sprintf('Tail for tolerance = %.3e', tol_sol))
end
%% compute error in polarization charge
x = S1.r(1,:).';
y = S1.r(2,:).';
z = S1.r(3,:).';
c1 = sum(layerdens(1:S1.npts).*S1.wts.*x);

xsolsol = S1sol.r(1,:).';
ysolsol = S1sol.r(2,:).';
c1_solsol = sum(layerdenssol(1:S1sol.npts).*S1sol.wts.*xsolsol);

err1 = abs(c1 - c1_solsol)./(abs(c1));
fprintf('error in polarization charge on obstacle 1=%d\n', err1);
%% Compute pointwise error in layer density at quadrature nodes S.r
% interpolate on Ssol to S.r
%fun_monitor = cat(1,layerdens',srcvals(1:3,:),S.mean_curv');
%fun_monitor_sol = cat(1,layerdenssol',srcvalssol(1:3,:),Ssol.mean_curv');

fun_monitor = cat(1,layerdens',S.mean_curv');
fun_monitor_sol = cat(1,layerdenssol',Ssol.mean_curv');

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
[refinelist1] = find_patch_refine_list(S1,errlp_patch(:,1:S1.npatches),funlpS1,tol_dens,[],lptype,eta_dens);
[refinelist2] = find_patch_refine_list(S2,errlp_patch(:,1+S1.npatches:end),funlpS2,tol_dens,[],lptype,eta_dens);

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

[S1ref,~,~,interp_1ref_mat] = split_patches(S1,refinelist1);
[S2ref,~,~,interp_2ref_mat] = split_patches(S2,refinelist2);
rhs_1ref = interp_1ref_mat * rhs1; % == 1
rhs_2ref = interp_2ref_mat * rhs2; % == 0
rhs_ref = [rhs_1ref;rhs_2ref];
Sref = merge([S1ref, S2ref]);
[densities_ref, errs] = lap3d.solver(Sref, 'dir', rhs_ref, 1e-8, rep_pars);
%%
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate error by upsampling
[Sus,~,~,interp_mat] = split_patches(S,1:S.npatches);

% interpolate
xyzus = S.r * interp_mat.';              % points xyz
Hus   = S.mean_curv.' * interp_mat.';    % mean curvature

% first half: torus centered at origin, axis along z
x = xyzus(1,1:Sus.npts/2);
y = xyzus(2,1:Sus.npts/2);
z = xyzus(3,1:Sus.npts/2);

theta = atan2(y,x);
rho   = sqrt(x.^2 + y.^2);
phi   = atan2(z, rho - rmajor);

xus_ana = (rmajor + rminor*cos(phi)).*cos(theta);
yus_ana = (rmajor + rminor*cos(phi)).*sin(theta);
zus_ana = rminor*sin(phi);

rho_safe = max(rho,1e-15);
Hext1 = (2*rho_safe - rmajor)./(2*rminor*rho_safe);

xyzus_ana = cat(1,xus_ana,yus_ana,zus_ana);

% second half: torus centered at c2, axis along z
x = xyzus(1,Sus.npts/2+1:end) - c2(1);
y = xyzus(2,Sus.npts/2+1:end) - c2(2);
z = xyzus(3,Sus.npts/2+1:end) - c2(3);

theta = atan2(y,x);
rho   = sqrt(x.^2 + y.^2);
phi   = atan2(z, rho - rmajor);

xus_ana = (rmajor + rminor*cos(phi)).*cos(theta) + c2(1);
yus_ana = (rmajor + rminor*cos(phi)).*sin(theta) + c2(2);
zus_ana = rminor*sin(phi) + c2(3);

rho_safe = max(rho,1e-15);
Hext2 = (2*rho_safe - rmajor)./(2*rminor*rho_safe);

xyzus_ana = cat(2,xyzus_ana,cat(1,xus_ana,yus_ana,zus_ana));
Hext = [Hext1, Hext2];

errL2  = norm((xyzus - xyzus_ana).*sqrt(repmat(Sus.wts',[3,1])),2);
errRMS = errL2 / sqrt(S.area);
fprintf('error in resolving the upsampled surface --- xyz = %.3e\n', errRMS);

errL2  = norm((Hus - Hext).*sqrt(Sus.wts.'));
errRMS = errL2 / sqrt(S.area);
fprintf('error in resolving the upsampled surface --- H = %.3e\n', errRMS);
%%

%%
nd = 1; % size(fun_monitor) = [S.npts,nd]
lp = 1;
fun_monitor = layerdens;

fun_monitor1 = layerdens(1:S1.npts);
fun_monitor2 = layerdens((1+S1.npts):end);

errp1 = S1.surf_fun_error(fun_monitor1,lp);
errp2 = S2.surf_fun_error(fun_monitor2,lp);

funlp1 = sum(abs(fun_monitor1).^lp .* repmat(S1.wts,[1,nd]),1);
funlpS2 = sum(abs(fun_monitor2).^lp .* repmat(S2.wts,[1,nd]),1);

[refinelist1] = find_refine_patches_list_sandbox(S1,errp1,funlp1,tol_dens);
[refinelist2] = find_refine_patches_list_sandbox(S2,errp2,funlpS2,tol_dens);

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

[S1ref,~,~,interp_1ref_mat] = split_patches(S1,refinelist1);
[S2ref,~,~,interp_2ref_mat] = split_patches(S2,refinelist2);
rhs_1ref = interp_1ref_mat * rhs1; % == 1
rhs_2ref = interp_2ref_mat * rhs2; % == 0
rhs_ref = [rhs_1ref;rhs_2ref];
Sref = merge([S1ref, S2ref]);
[densities_ref, errs] = lap3d.solver(Sref, 'dir', rhs_ref, 1e-8, rep_pars);
%%
x = S1ref.r(1,:).';
y = S1ref.r(2,:).';
z = S1ref.r(3,:).';
c1_ref = sum(densities_ref(1:S1ref.npts).*S1ref.wts.*x);
err1_ref = abs(c1_ref - c1_refsol)./(abs(c1_ref));
fprintf('Adaptively Refined: error in polarization charge on obstacle 1 =%d\n', err1_ref);

%% Refine again
funref_monitor1 = densities_ref(1:S1ref.npts);
funref_monitor2 = densities_ref((1+S1ref.npts):end);

errp1ref = S1ref.surf_fun_error(funref_monitor1,lp);
errp2ref = S2ref.surf_fun_error(funref_monitor2,lp);

funlp1ref = sum(abs(funref_monitor1).^lp .* repmat(S1ref.wts,[1,nd]),1);
funlp2ref = sum(abs(funref_monitor2).^lp .* repmat(S2ref.wts,[1,nd]),1);

[refinelist12] = find_refine_patches_list_sandbox(S1ref,errp1ref,funlp1ref,tol_dens);
[refinelist22] = find_refine_patches_list_sandbox(S2ref,errp2ref,funlp2ref,tol_dens);


if ifplot
    plotlist1 = zeros(S1ref.npatches,1);
    plotlist1(refinelist12) = 1;
    plotlist2 = zeros(S2ref.npatches,1);
    plotlist2(refinelist22) = 1;
    
    figure(8)
    plot(S1ref,plotlist1)
    hold on
    plot(S2ref,plotlist2)
    hold off
    colorbar
    title("Patches flagged for refinement")
end


[S1ref2,~,~,interp_12ref_mat] = split_patches(S1ref,refinelist12);
[S2ref2,~,~,interp_22ref_mat] = split_patches(S2ref,refinelist22);
rhs_12ref = interp_12ref_mat * rhs_1ref; % == 1
rhs_22ref = interp_22ref_mat * rhs_2ref; % == 0
rhs_ref2 = [rhs_12ref;rhs_22ref];
Sref2 = merge([S1ref2, S2ref2]);
[densities_ref2, errs2] = lap3d.solver(Sref2, 'dir', rhs_ref2, 1e-8, rep_pars);
x = S1ref2.r(1,:).';
y = S1ref2.r(2,:).';
z = S1ref2.r(3,:).';
c1_ref2 = sum(densities_ref2(1:S1ref2.npts).*S1ref2.wts.*x);
err1_ref2 = abs(c1_ref2 - c1_refsol)./(abs(c1_ref2));
fprintf('Adaptively Refined: error in polarization charge on obstacle 1 =%d\n', err1_ref2);
