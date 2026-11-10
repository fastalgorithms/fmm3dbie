% Sweep geometry monitors; strong density and weak potential errors in Lp.
% Run with the repository MATLAB/MEX dependencies available.
% Reference meshes are rebuilt from the analytic geometry, never from S.
script_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(fileparts(script_dir)));
addpath(genpath(fullfile(repo_dir,'matlab')));

lptypes = [0,1,2]; % 0 denotes Linf, matching the mesh generator
norm_names = {'Linf','L1','L2'};
opts.iptype = 1;
opts.norder = 6;
opts.tol_solve = 1e-12;       % fixed quadrature/FMM accuracy
opts.tol_gmres = 1e-12;       % estimated residual stopping tolerance
opts.tol_residual = 1e-11;    % acceptance threshold for recomputed residual
opts.maxit = 1000;
opts.rep_pars = [1,0];       % original single-layer equation

tols = [1e-5,1e-6,1e-7];
etas = [0,1,2];
geoms = {'cruller'}; % One reference shared by all cases in either mode
nd_monitors = [3,9,10];

% Successive independent reference meshes: tighten geometry tolerance,
% increase polynomial order, and double analytic chart subdivisions.
opts.ref_tol_geom = 1e-3;
opts.ref_order = 8;
opts.ref_eta = 0;
opts.ref_lptype = 0; % fixed reference generation, independent of sweep norm
opts.ref_nd_monitor = 10;
opts.ref_density_tol = min(tols)/10;
opts.ref_max_levels = 3;     % stop with an error if convergence fails

% Start with smoke mode to check the complete workflow on a small problem.
% Change to 'full' for the 81-case accuracy experiment configured above.
run_mode = 'full';
switch run_mode
    case 'smoke'
        opts.norder = 3;
        tols = 1e-2;
        etas = 0;
        nd_monitors = 3;
        % All three norms are still exercised (three sweep cases).
        opts.tol_solve = 1e-6;
        opts.tol_gmres = 1e-6;
        opts.tol_residual = 1e-5;
        opts.maxit = 200;
        opts.ref_tol_geom = 1e-3;
        opts.ref_order = 4;
        opts.ref_nd_monitor = 3;
        opts.ref_density_tol = 5e-2;
        opts.ref_max_levels = 2;
    case 'full'
opts.tol_solve = 1e-12;
opts.tol_gmres = 1e-12;
opts.tol_residual = 1e-11;
opts.norder = 6;
opts.maxit = 1000;

tols = [1e-5, 1e-6, 1e-7];
etas = [0, 1, 2];

opts.ref_nd_monitor = 10;
opts.ref_tol_geom = 1e-8;
opts.ref_order = 8;
opts.ref_density_tol = min(tols)/10;  % 1e-8
opts.ref_max_levels = 3;

        % Use the full settings above.
    otherwise
        error('monitor_sweep_density:RunMode','Unknown run mode: %s',run_mode);
end
ncase = numel(geoms)*numel(tols)*numel(etas)*numel(nd_monitors)*numel(lptypes);
fprintf('Mode: %s; %d density comparison cases.\n',run_mode,ncase);
% Smoke tolerances only check execution, not the full experiment's accuracy.

% Fixed observation surface for every case and every reference level.
% The cruller lies within radius 1.75; radius 3 is outside it.
opts.weak_radius = 3;
opts.weak_na = 4;
opts.weak_order = 6;
opts.tol_eval = min(opts.tol_solve,1e-10);
opts.ref_weak_tol = opts.ref_density_tol;
Sweak = geometries.sphere(opts.weak_radius,opts.weak_na,[0;0;0], ...
    opts.weak_order,1);
weak_targets = struct('r',Sweak.r,'patch_id',-ones(Sweak.npts,1));

dims = [numel(geoms),numel(tols),numel(etas),numel(nd_monitors),numel(lptypes)];
err_density_abs = nan(dims);
err_density_rel = nan(dims);
err_weak_abs = nan(dims);
err_weak_rel = nan(dims);
npatches = nan(dims);
npts = nan(dims);
solver_residual = nan(dims);
reference_info = cell(numel(geoms),1);
reference_densities = cell(numel(geoms),1);
reference_potentials = cell(numel(geoms),1);
% Variable-length outputs use cells with the same five axes as the errors.
density = cell(dims);
density_on_reference = cell(dims);
weak_potential = cell(dims);
gmres_iterations = nan(dims);
% One row per executed case; columns describe every printed case value.
case_columns = {'geometry_index','lptype','geometry_tol','eta','monitor', ...
    'strong_abs','strong_rel','weak_abs','weak_rel','npatches','npts', ...
    'solver_residual','gmres_iterations'};
case_data = nan(ncase,numel(case_columns));
case_index = zeros(ncase,5);
icase = 0;

for igeom = 1:numel(geoms)
    geom = geoms{igeom};
    [Sref,densref,potref,reference_info{igeom}] = ...
        build_reference(geom,opts,Sweak,weak_targets);
    reference_densities{igeom} = densref;
    reference_potentials{igeom} = potref;
    for inorm = 1:numel(lptypes)
        lptype = lptypes(inorm);
        for itol = 1:numel(tols)
            for ieta = 1:numel(etas)
                for imonitor = 1:numel(nd_monitors)
                    tol = tols(itol);
                    eta = etas(ieta);
                    nd_monitor = nd_monitors(imonitor);
                    fprintf('\n%s: %s, geometry tol %.1e, eta %g, monitor %d\n', ...
                        geom,norm_names{inorm},tol,eta,nd_monitor);
                    S = make_geometry(geom,tol,opts.norder,opts.iptype, ...
                        lptype,eta,nd_monitor,1);
                    [dens,rres,niter] = solve_density(S,geom,opts);
                    [ea,er,dens_interp] = testfuns.interp_lp_error(S,Sref,dens,densref,lptype);
                    err_density_abs(igeom,itol,ieta,imonitor,inorm) = ea;
                    err_density_rel(igeom,itol,ieta,imonitor,inorm) = er;
                    pot = lap3d.dirichlet.eval(S,dens,weak_targets, ...
                        opts.tol_eval,opts.rep_pars);
                    wa = density_norm(Sweak,pot(:)-potref(:),lptype);
                    wr = relative_error(wa,density_norm(Sweak,potref,lptype));
                    err_weak_abs(igeom,itol,ieta,imonitor,inorm) = wa;
                    err_weak_rel(igeom,itol,ieta,imonitor,inorm) = wr;
                    fprintf('Weak potential %s error: absolute %.6e, relative %.6e\n', ...
                        norm_names{inorm},wa,wr);
                    npatches(igeom,itol,ieta,imonitor,inorm) = S.npatches;
                    npts(igeom,itol,ieta,imonitor,inorm) = S.npts;
                    solver_residual(igeom,itol,ieta,imonitor,inorm) = rres;
                    density{igeom,itol,ieta,imonitor,inorm} = dens;
                    density_on_reference{igeom,itol,ieta,imonitor,inorm} = dens_interp;
                    weak_potential{igeom,itol,ieta,imonitor,inorm} = pot;
                    gmres_iterations(igeom,itol,ieta,imonitor,inorm) = niter;
                    icase = icase+1;
                    case_index(icase,:) = [igeom,itol,ieta,imonitor,inorm];
                    case_data(icase,:) = [igeom,lptype,tol,eta,nd_monitor, ...
                        ea,er,wa,wr,S.npatches,S.npts,rres,niter];
                    fprintf('Strong density %s error: absolute %.6e, relative %.6e\n', ...
                        norm_names{inorm},ea,er);
                    fprintf('Patches %d, points %d, solve residual %.3e\n', ...
                        S.npatches,S.npts,rres);
                end
            end
        end
    end % norm sweep
end
% Results, Sref, densref, Sweak, and potref remain in the workspace.
% Array axes: geometry, tol, eta, monitor, norm.
% Linf is the maximum sampled on the reference quadrature nodes.

function [Sref,densref,potref,info] = build_reference(geom,opts,Sweak,targets)
    info = struct('order',[],'geometry_tol',[],'chart_multiplier',[], ...
        'npatches',[],'npts',[],'solver_residual',[], ...
        'relative_density_change',[],'relative_weak_change',[]);
    info.level = [];
    info.gmres_iterations = [];
    info.density = cell(opts.ref_max_levels,1);
    info.potential = cell(opts.ref_max_levels,1);
    for level = 1:opts.ref_max_levels
        order = opts.ref_order + 2*(level-1);
        tol = opts.ref_tol_geom * 0.1^(level-1);
        chart_multiplier = 2^(level-1);
        Snew = make_geometry(geom,tol,order,opts.iptype,opts.ref_lptype, ...
            opts.ref_eta,opts.ref_nd_monitor,chart_multiplier);
        [dnew,rres,niter] = solve_density(Snew,geom,opts);
        pnew = lap3d.dirichlet.eval(Snew,dnew,targets,opts.tol_eval,opts.rep_pars);
        change = nan(1,3);
        weak_change = nan(1,3);
        if level > 1
            % Locate/interpolate once, then check all three density norms.
            [~,~,dprev] = testfuns.interp_lp_error(Sref,Snew,densref,dnew,2);
            delta = dprev(:)-dnew(:);
            for k = 1:3
                change(k) = relative_error(density_norm(Snew,delta,k-1), ...
                    density_norm(Snew,dnew,k-1));
                weak_change(k) = relative_error( ...
                    density_norm(Sweak,pnew(:)-potref(:),k-1), ...
                    density_norm(Sweak,pnew,k-1));
            end
        end
        info.level(level) = level;
        info.gmres_iterations(level) = niter;
        info.density{level} = dnew;
        info.potential{level} = pnew;
        info.order(level) = order;
        info.geometry_tol(level) = tol;
        info.chart_multiplier(level) = chart_multiplier;
        info.npatches(level) = Snew.npatches;
        info.npts(level) = Snew.npts;
        info.solver_residual(level) = rres;
        info.relative_density_change(level,:) = change;
        info.relative_weak_change(level,:) = weak_change;
        fprintf('%s reference %d: %d points, density change [Linf L1 L2] = [%.6e %.6e %.6e]\n', ...
            geom,level,Snew.npts,change);
        fprintf('Reference weak change [Linf L1 L2] = [%.6e %.6e %.6e]\n',weak_change);
        potref = pnew;
        Sref = Snew;
        densref = dnew;
        if level > 1 && all(isfinite(change)) && all(change <= opts.ref_density_tol) && ...
                all(isfinite(weak_change)) && all(weak_change <= opts.ref_weak_tol)
            return
        end
    end
    error('monitor_sweep_density:ReferenceNotConverged', ...
        '%s reference failed density or weak convergence in %d levels.', ...
        geom,opts.ref_max_levels);
end

function S = make_geometry(geom,tol,order,iptype,lptype,eta,monitor,mult)
    switch geom
        case 'cruller'
            S = geometries_adap.cruller_adap([0.2;0.1;0.05;5;3], ...
                [5;5;5],mult*[3,2],order,iptype,tol,lptype,eta,monitor);
        case 'stellarator'
            S = geometries_adap.stellarator_adap(mult*[1,2], ...
                order,iptype,tol,lptype,eta,monitor);
        otherwise
            error('monitor_sweep_density:UnknownGeometry','Unknown geometry: %s',geom);
    end
end

function [dens,rres,niter] = solve_density(S,geom,opts)
    switch geom
        case 'cruller'
            src.r = [3.17;-0.03;6.15];
        case 'stellarator'
            src.r = [5.5;0.5;0.1];
        otherwise
            error('monitor_sweep_density:UnknownGeometry','Unknown geometry: %s',geom);
    end
    rhs = lap3d.kern(src,S,'s');
    solve_opts.eps_gmres = opts.tol_gmres;
    solve_opts.maxit = opts.maxit;
    fprintf('Starting solve: %d patches, %d points\n', ...
    S.npatches, S.npts);
drawnow;
    [dens,history,rres] = lap3d.dirichlet.solver( ...
        S,rhs,opts.tol_solve,opts.rep_pars,solve_opts);
    niter = numel(history);
    if any(~isfinite(dens(:))) || ~isfinite(rres) || rres > opts.tol_residual
        error('monitor_sweep_density:SolveFailed', ...
            '%s density solve failed: residual %.3e (acceptance threshold %.3e).', ...
            geom,rres,opts.tol_residual);
    end
end

function value = density_norm(S,density,p)
    if p == 0
        value = max(abs(density(:)));
    else
        value = sum(S.wts(:).*abs(density(:)).^p)^(1/p);
    end
end

function value = relative_error(error_abs,reference_norm)
    if reference_norm > 0
        value = error_abs/reference_norm;
    elseif error_abs == 0
        value = 0;
    else
        value = Inf;
    end
end
