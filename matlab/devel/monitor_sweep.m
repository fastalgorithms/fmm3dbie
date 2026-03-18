opts.lptype = 2; % l2 - norm
opts.iptype = 1; % triangles with VR nodes
opts.norder = 6; % not low enough to run out of memory
opts.tol_ref = 1e-7;
tols = [1e-5,1e-6,1e-7]';
etas = [1,2]'; % OBS sqrt(patch of area)
geoms = {'cruller','stellarator'};%,'trefoil'};
nd_monitors = [3,9,10]';

ntest = numel(tols)*numel(etas)*numel(geoms)*numel(nd_monitors);

errslp_strong = zeros(numel(geoms),numel(tols),numel(etas),numel(nd_monitors));
errslp_weak = zeros(numel(geoms),numel(tols),numel(etas),numel(nd_monitors));
npatches = zeros(numel(geoms),numel(tols),numel(etas),numel(nd_monitors));

for igeom = 1:numel(geoms)
    for itol = 1:numel(tols)
        for ieta = 1:numel(etas)
            for ind_monitor = 1:numel(nd_monitors)
                
                geom = geoms{igeom};
                tol = tols(itol);
                eta = etas(ieta);
                nd_monitor = nd_monitors(ind_monitor);
                fprintf("---------- \n")
                fprintf("Geometry: %s. Tol: %e eta: %d. nd_monitor: %d \n", ...
                    geom,tol,eta,nd_monitor)

                [err_strong,err_weak,npatch,~] = ...
                    run_test_adap_mesh(geom,tol,eta,nd_monitor,opts);
                errslp_strong(igeom,itol,ieta,ind_monitor) = err_strong;
                errslp_weak(igeom,itol,ieta,ind_monitor) = err_weak;
                npatches(igeom,itol,ieta,ind_monitor) = npatch;

            end
        end
    end
end



function [err_strong,err_weak,npatches,nrefpat] = run_test_adap_mesh(geom,tol,eta,nd_monitor,opts)
tol_ref = opts.tol_ref;
norder = opts.norder;
iptype_pat = opts.iptype;
lptype = opts.lptype;

% test: Single layer potential
dpars = [1,0]; % Single layer potential


switch geom
    case 'cruller'
        xyz_out = [3.17,-0.03,6.15]';
        radii(1) = 0.2;
        radii(2) = 0.1;
        radii(3) = 0.05;
        radii(4) = 5.0;
        radii(5) = 3.0;
        nu = 3;
        nv = 2;
        nuv = [nu,nv]; % split (u,v) space into nu X nv parts patches at root level. x2 for triangles

        scales  = [5 5 5]';
        S = geometries_adap.cruller_adap(radii,scales,nuv,norder,iptype_pat,...
        tol,lptype,eta,nd_monitor);
    case 'stellarator'
        nu = 1;
        nv = 2;
        nuv = [nu,nv];
        S = geometries_adap.stellarator_adap(nuv,norder,iptype_pat,tol);
        [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);       
        xyz_out = [5.5;0.5;0.1];
    case 'trefoil'
end

% build reference geometry by upsampling given geometry twice
[Sus,~,~,interp_mat] = split_patches(S,1:S.npatches);
% create rhs
src_info = [];
src_info.r = xyz_out;
rhs = lap3d.kern(src_info,S,'s');
rhsus = lap3d.kern(src_info,Sus,'s');
% solve for layer density
dens = lap3d.dirichlet.solver(S,rhs,tol,dpars);
densus = lap3d.dirichlet.solver(Sus,rhsus,tol_ref,dpars);
% compute potentials
pot = lap3d.dirichlet.eval(S,dens,S,tol_ref,dpars);
potus = lap3d.dirichlet.eval(Sus,densus,Sus,tol_ref,dpars);
% interpolate density dens and potential pot
dens0us = dens.' * interp_mat.'; 
pot0us = pot.' * interp_mat.';
% compute potentials using Green's function
pot_ex = lap3d.kern(src_info,S,'s');
potus_ex = lap3d.kern(src_info,Sus,'s');
% compute errors
err_strong_ana = sqrt(sum(abs(pot-pot_ex).^2.*S.wts))/sqrt(sum(abs(pot_ex).^2.*S.wts));
err_strong_ana_us = sqrt(sum(abs(potus-potus_ex).^2.*Sus.wts))/sqrt(sum(abs(potus_ex).^2.*Sus.wts));
err_weak = sqrt(sum(abs(dens0us.'- densus).^2.*Sus.wts))/sqrt(sum(Sus.wts.*abs(densus).^2));
err_strong = sqrt(sum(abs(pot0us.'- potus).^2.*Sus.wts))/sqrt(sum(Sus.wts.*abs(potus).^2));
% print
fprintf('Error in potential=%d\n',err_strong_ana);
fprintf('Error in potential, upsampled=%d\n',err_strong_ana_us);
fprintf('(Weak) Error in upsampling density =%d\n',err_weak);
fprintf('(Strong) Error in upsampling potential =%d\n',err_strong);
% output
npatches = S.npatches;
nrefpat = Sus.npatches;
return

end

%%
S.npatches
S.npts