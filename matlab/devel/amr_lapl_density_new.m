%% Adaptive density refinement on two close spheres

% The sweep runs this script in a fresh function workspace for each case.
if ~exist('amr_sweep_case','var')
    clear
end

% Smoke settings are applied only when requested by the sweep.
is_smoke = exist('amr_sweep_case','var') && ...
    isfield(amr_sweep_case,'smoke') && amr_sweep_case.smoke;

%% Error metric
errtype = 'strong';       % error in density
%errtype = 'weak';         % error in potential on outer confining ellipsoid
%errtype = 'functional';  % error in user-defined functional

%% Error estimator for patch refinement
esttype = 'tails_vals';
%esttype = 'tails_coefs';
%esttype = 'upsample';
%esttype = 'nystrom';

if exist('amr_sweep_case','var')
    errtype = amr_sweep_case.errtype;
    esttype = amr_sweep_case.esttype;
end

if strcmp(errtype,'weak')
    esttype = 'weak';
end

if strcmp(errtype,'functional')
    esttype = 'functional';
end


%% Set norm and scaling
lptype = 2;   % 0 = Linf, otherwise L^p
eta_dens = 1; % scale factor for AMR, patchlength^eta

if exist('amr_sweep_case','var')
    lptype = amr_sweep_case.lptype;
    eta_dens = amr_sweep_case.eta_dens;
end

%% Set functional, if any
% Linear scalar observation of the potential on the fixed outer ellipsoid.
% L(sig) = integral over Souter of T(sig) dS, divided by area(Souter).
% T(sig) is the combined-layer potential evaluated on Souter.
functional = @(Souter,pot) sum(Souter.wts(:).*pot(:))/sum(Souter.wts(:));
fprintf("****\n")
fprintf("%s error in L%d, eta=%g. Error estimate: %s \n",...
    errtype,lptype,eta_dens,esttype)
%% Parameters
tol_amr = 1e-6; % AMR tolerance
if exist('amr_sweep_case','var')
    tol_amr = amr_sweep_case.tol_amr;
end
tol_sol = 1e-9; % tolerance for solvie BIE
refine_fraction = 1; % refine refine_fraction of flagged patches
norderhead_shift = 1; % order of tails te estimate error: 1-highest order, 2-second to highest and highest

ifplot = false;
compref = true; % compute and verify a new reference solution
ref_opts = struct('nasol',10,'nordersol',10,'prefactor',2, ...
    'max_refinements',20,'max_attempts',3,'max_npts',Inf, ...
    'maxit',200,'tol_sol',1e-12,'tol_change',min(1e-8,tol_amr/10));
reference = [];
if exist('amr_sweep_case','var')
    ref_opts.tol_change = min(1e-8,amr_sweep_case.reference_tol);
    reference = amr_sweep_case.reference;
end
if is_smoke
    tol_sol = 1e-7;
    ref_opts.nasol = 6;
    ref_opts.nordersol = 8;
    ref_opts.prefactor = 1; % fewer close-gap patches; convergence checks remain
    ref_opts.tol_sol = tol_sol;
    ref_opts.tol_change = amr_sweep_case.reference_tol;
end
reference_file = 'ref_solution01_new2.mat'; % trust ref_solution01_new2.mat to abs err 3e-9 and rel err 5e-11
%reference_file = 'ref_solution001.mat'; % trust ref_solution01_new2.mat to abs err 3e-9 and rel err 5e-11

%% Geometry and boundary data
a = 1.1; % radius
d = 0.01; % separation distance
c0_1 = [0;0;0]; % center of sphere 1
c0_2 = [2*a + d; 0; 0]; % center of sphere 2
%% Discretization parameters
na = 9; % number of initial refinements
norder = 9; % quadrature order
if is_smoke
    na = 3;
    norder = 6;
end
%% Order measure of tails in error estimate
norderhead = norder-norderhead_shift;
%% Setup geometry
S1 = geometries.sphere(a,na,c0_1,norder,1);
S2 = geometries.sphere(a,na,c0_2,norder,1);
S = merge([S1,S2]);

rhs1 = ones(S1.npts,1);
rhs2 = zeros(S2.npts,1);
rhs = [rhs1;rhs2];

rep_pars = [1.0;1.0]; % potenial = a*S + b*D, rep_pars = [a,b];

fprintf('Initial patch counts: S1 = %d, S2 = %d.\n', ...
    S1.npatches,S2.npatches)
fprintf('Initial point counts: S1 = %d, S2 = %d.\n', ...
    S1.npts,S2.npts)
%% norm name
if lptype == 0
    norm_name = 'Linf';
else
    norm_name = sprintf('L%d',lptype);
end

%% error estimate based on mean curvature H
% estimate resolution of spheres by estimating mean curvature H
Hext = 1/a; % analytic formula for H
H = S.mean_curv; % approximation
Hdiff = abs(H - Hext);

err_H = testfuns.comp_Lp(S,Hdiff(:),lptype);

fprintf('absolute error in resolving curvature in %s norm: %.16e\n', ...
    norm_name,err_H);

%% Outer target ellipsoid for weak error
if compref || strcmp(errtype,'weak') || strcmp(errtype,'functional')
    c0_out = (c0_1+c0_2)/2;

    % smallest patch legnith on S
    hsrc = sqrt(full(sparse(S.patch_id(:),1,S.wts(:), ...
        S.npatches,1)));

    h0 = 5*min(hsrc);     % minimum distance between ellipsoid and spheres
    R0 = a + h0;
    z0 = a + d/2;
    s0 = sqrt(z0^2 + 3*R0^2);
    a0 = sqrt((s0-z0)*(s0+2*z0)/3);
    b0 = a0;
    c0 = (s0+2*z0)/sqrt(3);
    %   generate geometry
    Souter = geometries.ellipsoid([c0 a0 b0], ...
        [5,5,5],c0_out,norder,1);

    targinfo_outer.r = Souter.r;
    targinfo_outer.patch_id = -ones(Souter.npts,1);

end
%% Reference solution
if compref
    reference = testfuns.two_sphere_reference(a,c0_1,c0_2,rep_pars, ...
        Souter,functional,ref_opts,reference);
    S1sol = reference.S1sol;
    S2sol = reference.S2sol;
    Ssol = reference.Ssol;
    sig_sol = reference.sig_sol;
    layerdenssol = sig_sol;
else
    % load precomputed reference solution
    refdata = load(reference_file,'S1sol','S2sol','layerdenssol','pot_sol');

    S1sol = refdata.S1sol;
    S2sol = refdata.S2sol;
    Ssol = merge([S1sol,S2sol]);
    sig_sol = refdata.layerdenssol;
    pot_sol = refdata.pot_sol;
end

%% Function to monitor for error
if strcmp(errtype,'strong')
    fun_monitor_sol = sig_sol.';
elseif strcmp(errtype,'weak') || strcmp(errtype,'functional')
    if compref
        pot_sol_outer = reference.pot_outer;
    else
        pot_sol_outer = lap3d.dirichlet.eval( ...
            Ssol,sig_sol,targinfo_outer,tol_sol,rep_pars);
    end
    fun_monitor_sol = pot_sol_outer;
    if strcmp(errtype,'functional')
        fun_monitor_sol = functional(Souter,pot_sol_outer);
    end
end

%% Initial solve of BIE for layer density sig
sig = lap3d.solver(S,'dir',rhs,tol_sol,rep_pars);
%% Function to monitor for error estimates
if strcmp(errtype,'strong')
    fun_monitor = sig;
elseif strcmp(errtype,'weak')
    fun_monitor = sig;
elseif strcmp(errtype,'functional')
    fun_monitor = sig; % density tails seed the first refinement
end
%% Initial actual error

if strcmp(errtype,'strong')
    [fun_monitor_lp_abs,fun_monitor_lp_rel,~] = ...
        testfuns.interp_lp_error(S,Ssol,fun_monitor,fun_monitor_sol,lptype);
elseif strcmp(errtype,'weak')
    pot_outer = lap3d.dirichlet.eval( ...
        S,sig,targinfo_outer,1e-12,rep_pars);
    pot_outer_diff = pot_outer-pot_sol_outer;
    % compute lp norm of difference on Souter (ellipsoid)
    fun_monitor_lp_abs = testfuns.comp_Lp(Souter,pot_outer_diff,lptype);
    % compute strong error of density
    [sig_lp_abs,sig_lp_rel,~] = ...    
        testfuns.interp_lp_error(S,Ssol,sig,sig_sol,lptype);
    % errors are relative ||sig||_p + ||f||_p
     sig_p_rhs_norm_lp = testfuns.comp_Lp(S,sig,lptype) + testfuns.comp_Lp(S,rhs,lptype);

     fprintf('Initial actual %s %s relative error in density = %.16e\n', ...
    norm_name,errtype,sig_lp_abs/sig_p_rhs_norm_lp)
elseif strcmp(errtype,'functional')
    pot_outer = lap3d.dirichlet.eval(S,sig,targinfo_outer,tol_sol,rep_pars);
    functional_coarse = functional(Souter,pot_outer);
    fun_monitor_lp_abs = abs(functional_coarse-fun_monitor_sol);
end
% compute relative scaling ||sig||_p+||rhs||_p
sig_p_rhs_norm_lp = testfuns.comp_Lp(S,sig,lptype) + testfuns.comp_Lp(S,rhs,lptype);
if strcmp(errtype,'functional')
    fun_monitor_lp_rel = fun_monitor_lp_abs/sig_p_rhs_norm_lp;
end
fprintf('Initial actual %s %s relative error = %.16e\n', ...
    norm_name,errtype,fun_monitor_lp_abs/sig_p_rhs_norm_lp)
%% Estimate errors and flag refinement candidates

% for 'upsample', 'nystrom', 'weak' we need an initial refinement. We use
% tails_vals

% flag patches to refine
if strcmp(esttype,'tails_vals')
    errlp_patch = testfuns.surf_fun_tail_val_error(S,fun_monitor,lptype, ...
        norderhead,[]);
elseif strcmp(esttype,'tails_coefs')
    errlp_patch = testfuns.surf_fun_koorn_shell_error(S,fun_monitor,...
        lptype,norderhead,[]);
elseif strcmp(esttype,'upsample') || strcmp(esttype,'nystrom') || strcmp(esttype,'weak') || strcmp(esttype,'functional')
    errlp_patch = testfuns.surf_fun_tail_val_error(S,fun_monitor,lptype, ...
        norderhead,[]);
else
    error("Estimator not implemented");
end

% get list of patches to refine
[refinelist,errqlist] = testfuns.find_patch_refine_list( ...
    S,errlp_patch,sig_p_rhs_norm_lp,tol_amr,[],lptype,eta_dens);


fprintf('First refinement: S refines %d of %d patches\n', ...
    numel(refinelist),S.npatches);

%% Adaptive refinement loop
irefine = 0; % number of refinement iterations

% variables for upsampling only
fun_monitor_coarse = fun_monitor;
fun_monitor_lp_abs_coarse = fun_monitor_lp_abs;
fun_monitor_norm_lp_coarse = sig_p_rhs_norm_lp;
S_coarse = S;
sig_coarse = sig;

% start refinement loop. Terminate when list of refinement candidates is
% empty
while ~isempty(refinelist)
    irefine = irefine + 1;
    fprintf("Refinement %d \n",irefine);

    % refine surface S->S
    [S,interp_refine_mat] = split_patches(S,refinelist);

    % get new rhs
    rhs = interp_refine_mat*rhs;

    % solve BIE on refined S
    [sig,~] = lap3d.solver(S,'dir',rhs,tol_sol,rep_pars);

    % get patchwise error estimates
    if strcmp(errtype,'strong')
        fun_monitor = sig;
    elseif strcmp(errtype,'weak')
        fun_monitor = sig;
    elseif strcmp(errtype,'functional')
        fun_monitor = sig; % density tails seed the first refinement
    end

    fun_monitor_norm_S = testfuns.comp_Lp(S,sig,lptype) ...
        + testfuns.comp_Lp(S,rhs,lptype);
    hpatch = sqrt(full(sparse(S.patch_id(:),1,S.wts(:),S.npatches,1)));

    if strcmp(esttype,'tails_vals')
        errlp_patch = testfuns.surf_fun_tail_val_error(S, fun_monitor,lptype, ...
            norderhead,[]);
        errlp_est = errlp_patch .* hpatch(:).'.^eta_dens;
    elseif strcmp(esttype,'tails_coefs')
        errlp_patch = testfuns.surf_fun_koorn_shell_error(S,fun_monitor,...
            lptype,norderhead);
        errlp_est = errlp_patch .* hpatch(:).'.^eta_dens;
    elseif strcmp(esttype,'upsample')
        monitor_fun_err = fun_monitor - interp_refine_mat*fun_monitor_coarse(:);
        errlp_patch = testfuns.get_patch_lp_norm(S,monitor_fun_err.',[],lptype);
        errlp_est = errlp_patch .* hpatch(:).'.^eta_dens;
    elseif strcmp(esttype,'nystrom')
        opts = struct( ...
            'tol',1e-12, ...
            'maxiter',100, ...
            'rfac',2);

        % Locate nodes of S_coarse on the finer S
        [~,patch_id_c2f,uvs_targ_c2f,~,flags] = ...
            get_closest_pts(S_coarse,S.r,opts);
        targ_info = [];
        targ_info.r = S.r;
        targ_info.patch_id = patch_id_c2f;
        targ_info.uvs_targ = uvs_targ_c2f;
        % evaluate coarse potential on refined surface
        pot_coarse_on_fine = lap3d.dirichlet.eval(S_coarse,sig_coarse,...
            targ_info,tol_sol,rep_pars);
        sig_coarse_nystrom = 2*(rhs - pot_coarse_on_fine);
        monitor_fun_err = sig - sig_coarse_nystrom;
        errlp_patch = testfuns.get_patch_lp_norm(S,monitor_fun_err.',[],lptype);
        errlp_est = errlp_patch .* hpatch(:).'.^eta_dens;
        % update
        S_coarse = S;
        sig_coarse = sig;

    elseif strcmp(esttype,'functional')
        pot_outer = lap3d.dirichlet.eval(S,sig,targinfo_outer,tol_sol,rep_pars);
        functional_fine = functional(Souter,pot_outer);
        % A coarse-error estimate. Each solve uses its own source quadrature.
        err_monitor_est_abs = abs(functional_fine-functional_coarse);
        % Local indicators use the density surplus on the refined surface.
        sigd = sig(:)-interp_refine_mat*sig_coarse(:);
        errlp_patch = patch_functional_error(S,sigd,Souter,rep_pars,functional);
        sig_coarse = sig;
        functional_coarse = functional_fine;

    elseif strcmp(esttype,'weak')

        [fun_monitor_lp_abs,fun_monitor_lp_rel,sig_sol_f] = ...
            testfuns.interp_lp_error(S_coarse,S,sig_coarse,sig,lptype);

        sigd = sig_sol_f - sig;
        [errlp_patch,~] = testfuns.patch_weak_norm(S,sigd,Souter,rep_pars,lptype);
        errlp_est = errlp_patch .* hpatch(:).'.^eta_dens;

        S_coarse = S;
        sig_coarse = sig;
    end

%   find patches to flag for refinement
    [refinelist,errqlist] = testfuns.find_patch_refine_list( ...
        S,errlp_patch,fun_monitor_norm_S,tol_amr,[],lptype,eta_dens);

    refinelist_orig = refinelist(:);
%   refine only refine_fraction fraction of flagged patches
    refinelist = testfuns.select_top_fraction( ...
        refinelist_orig,errqlist,refine_fraction);


%   estimate global error in lp-norm and compute actual error
    if ~strcmp(esttype,'functional')
        if lptype == 0
            err_monitor_est_abs = max(errlp_est);
        else
            err_monitor_est_abs = (sum(errlp_est.^lptype,2)).^(1/lptype);
        end
        err_monitor_est_rel = err_monitor_est_abs/fun_monitor_norm_S;
    else
        err_monitor_est_rel = err_monitor_est_abs/fun_monitor_norm_lp_coarse;
    end

    if strcmp(errtype,"strong")
        [err_monitor_lp_abs_comp,err_monitor_lp_rel_comp,~] = ...
            testfuns.interp_lp_error( ...
            S,Ssol,fun_monitor,fun_monitor_sol,lptype);
    elseif strcmp(errtype,'weak')
        pot_outer = lap3d.dirichlet.eval( ...
            S,sig,targinfo_outer,tol_sol,rep_pars);
        pot_outer_diff = pot_outer-pot_sol_outer;
       err_monitor_lp_abs_comp = testfuns.comp_Lp(Souter,pot_outer_diff,lptype);
    elseif strcmp(errtype,'functional')
        err_monitor_lp_abs_comp = abs(functional_fine-fun_monitor_sol);
    end

% print errors
    if strcmp(esttype,'functional')
        fprintf('Functional change from coarse to refined solution = %.16e\n', ...
            err_monitor_est_abs)
        fprintf('Actual absolute functional error of coarse solution = %.16e\n', ...
            fun_monitor_lp_abs_coarse)
        fprintf('Estimated relative functional error of coarse solution = %.16e\n', ...
            err_monitor_est_rel)
        fprintf('Actual relative functional error of coarse solution = %.16e\n', ...
            fun_monitor_lp_abs_coarse/fun_monitor_norm_lp_coarse)
        fun_monitor_lp_abs_coarse = err_monitor_lp_abs_comp;
        fun_monitor_norm_lp_coarse = fun_monitor_norm_S;
    else
        fprintf('Estimated absolute %s error of solution = %.16e\n',...
            norm_name,err_monitor_est_abs)
        fprintf('Estimated relative %s error of solution = %.16e\n',...
            norm_name,err_monitor_est_rel)
    end
    if strcmp(esttype,'upsample')
        fprintf('Actual absolute %s error of coarse solution    = %.16e\n',...
            norm_name,fun_monitor_lp_abs_coarse)
        fprintf('Actual relative %s error of coarse solution    = %.16e\n',...
            norm_name,fun_monitor_lp_abs_coarse/fun_monitor_norm_lp_coarse);
        fun_monitor_coarse = fun_monitor;
        fun_monitor_lp_abs_coarse = err_monitor_lp_abs_comp;
        fun_monitor_norm_lp_coarse = fun_monitor_norm_S;

    end
    if strcmp(errtype,'functional')
        fprintf('Actual absolute functional error of refined solution = %.16e\n', ...
            err_monitor_lp_abs_comp)
        fprintf('Actual relative functional error of refined solution = %.16e\n', ...
            err_monitor_lp_abs_comp/fun_monitor_norm_S)
    else
        fprintf('Actual absolute %s error of solution    = %.16e\n',...
            norm_name,err_monitor_lp_abs_comp)
        fprintf('Actual relative %s error of  solution    = %.16e\n',...
            norm_name,err_monitor_lp_abs_comp/fun_monitor_norm_S)
    end

    fprintf('Patch counts: S = %d.\n', S.npatches)
    fprintf('Point counts: S = %d.\n',S.npts)

    fprintf('S refines %d of %d flagged patches\n', ...
        numel(refinelist),numel(refinelist_orig));


    fprintf('--------\n')

end

fprintf('AMR stopped after %d refinement iterations.\n',irefine)
fprintf('Final patch counts: S = %d.\n',S.npatches)
fprintf('Final point counts: S = %d.\n',S.npts)

function errlp_patch = patch_functional_error(S,sigd,Souter,rep_pars,functional)
% Same separated-target smooth quadrature as patch_weak_norm. Apply the
% functional to each patch field before taking its magnitude for marking.
% The global estimate above keeps cancellation between patches.
    strength = S.wts(:).*sigd(:);
    targets.r = Souter.r;
    errlp_patch = zeros(1,S.npatches);
    for j = 1:S.npatches
        idx = S.ixyzs(j):S.ixyzs(j+1)-1;
        src.r = S.r(:,idx);
        src.n = S.n(:,idx);
        pot = lap3d.kern(src,targets,'c',rep_pars(1),rep_pars(2))*strength(idx);
        errlp_patch(j) = abs(functional(Souter,pot));
    end
end
