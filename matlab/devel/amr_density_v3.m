%% Adaptive density refinement on two close spheres

close all
clear

%% Error metric
errtype = 'strong';       % pointwise
%errtype = 'weak';       % integrated Lp
%errtype = 'functional'; % prescribed functional

if strcmp(errtype,'strong')
    lptype = 0;
    eta_dens = 1;
elseif strcmp(errtype,'weak')
    lptype = 0;   % user sets p
    eta_dens = 1;
elseif strcmp(errtype,'functional')
    lptype = 1;   % dummy
    eta_dens = 1;
    % polarization
    gfun = @(sigma,S,idx) sum(S.r(1,idx).' .* sigma(idx) .* S.wts(idx));
end

%% Parameters
tol_dens = 1e-7; % tolerance in density error: weak, strong, or functional
tol_sol = 1e-9; % tolerance for solvie BIE
refine_fraction = 1.0; % refine refine_fraction of flagged patches 

ifplot = false;
compref = false; % compute reference solution
reference_file = 'ref_solution001.mat';

%% Geometry and boundary data
a = 1.1; % radius
d = 0.001; % separation distance
c0_2 = [2*a + d; 0; 0];
na = 9; % number of initial refinements
norder = 9; % quadrature order

S1 = geometries.sphere(a,na,[],norder,1);
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
%% error estimate based on mean curvature H
% estimate resolution of spheres by estimating mean curvature H
Hext = 1/a; % analytic formula for H
H = S.mean_curv; % approximation
Hdiff = abs(H - Hext);

errLp = comp_ffun_lp(S,Hdiff(:).',lptype);

% relative
if lptype ~= 0
    errLp = errLp/S.area^(1/lptype);
end

fprintf('error in resolving curvature in L%d norm: %d\n', lptype,errLp);

%% Reference solution
% compute reference solution, or load one
if compref
    nasol = 10;
    nordersol = 10;

    S1sol_orig = geometries.sphere(a,nasol,[],nordersol,1);
    S2sol_orig = geometries.sphere(a,nasol,c0_2,nordersol,1);

    [S1sol,S2sol] = refine_close_patches( ...
        S1sol_orig,S2sol_orig,0.5,10);

    Ssol = merge([S1sol,S2sol]);
    rhssol = [ones(S1sol.npts,1);zeros(S2sol.npts,1)];

    [sig_sol,~] = lap3d.solver( ...
        Ssol,'dir',rhssol,tol_sol,rep_pars);
else
    refdata = load(reference_file,'S1sol','S2sol','layerdenssol');

    S1sol = refdata.S1sol;
    S2sol = refdata.S2sol;
    sig_sol = refdata.layerdenssol;
end

fun_monitor_sol = sig_sol.';

%% Initial solve
[sig,~] = lap3d.solver(S,'dir',rhs,tol_sol,rep_pars);

if ifplot
    figure(1)
    clf
    plot(S,sig)
    colorbar
    title('Layer density')

    figure(2)
    clf
    plot(S,rhs)
    colorbar
    title('Boundary data')
end

%% Reference polarization charge
if strcmp(errtype,'functional')
    pol1 = gfun(sig(1:S1.npts),S1,1:S1.npts);
    pol1_sol = gfun(sig_sol(1:S1sol.npts),S1sol,1:S1sol.npts);

    err1 = abs(pol1-pol1_sol)/abs(pol1_sol);
    fprintf('Initial polarization-charge error = %.16e\n',err1)
elseif strcmp(errtype,'weak')
    pot_sol = lap3d.dirichlet.eval(Ssol,sig_sol,S,tol_sol,rep_pars);
    pot = lap3d.dirichlet.eval(S,sig,S,tol_sol,rep_pars);
    [sig_lp_abs,sig_lp_rel,~,~,~,~,~] = ...
        density_diagnostics( ... 
        S1,S2,S1sol,S2sol,pot,pot_sol,lptype);
    fprintf('Initial actual %s error = %.16e\n', ...
        norm_name,sig_lp_abs)
elseif strcmp(errtype,'strong')
    if lptype == 0
        norm_name = 'Linf';
    else
        norm_name = sprintf('L%d',lptype);
    end

    [sig_lp_abs,sig_lp_rel,~,~,~,~,~] = ...
        density_diagnostics( ...
        S1,S2,S1sol,S2sol,sig,sig_sol,lptype);

    fprintf('Initial actual %s error = %.16e\n', ...
        norm_name,sig_lp_abs)
end
%% Initial density error and refinement candidates

% relative_scale

% compute norms
siglpS1 = comp_ffun_lp(S1,sig(1:S1.npts).',lptype);
siglpS2 = comp_ffun_lp(S2,sig(1+S1.npts:end).',lptype);

[flagged1,errqlist1] = flag_refinement_candidates( ...
    S1,sig(1:S1.npts).',siglpS1,tol_dens,lptype,eta_dens);
[flagged2,errqlist2] = flag_refinement_candidates( ...
    S2,sig(1+S1.npts:end).',siglpS2,tol_dens,lptype,eta_dens);

nflagged1 = numel(flagged1);
nflagged2 = numel(flagged2);

% Apply the percentage restriction.
refinelist1 = select_top_fraction( ...
    flagged1,errqlist1,refine_fraction);
refinelist2 = select_top_fraction( ...
    flagged2,errqlist2,refine_fraction);

fprintf('First refinement: S1 refines %d of %d flagged patches\n', ...
    numel(refinelist1),nflagged1);
fprintf('First refinement: S2 refines %d of %d flagged patches\n', ...
    numel(refinelist2),nflagged2);

plot_refinement_lists(S1,S2,refinelist1,refinelist2,ifplot,0)

%% Adaptive refinement loop
irefine = 0;

while ~isempty(refinelist1) || ~isempty(refinelist2)
    irefine = irefine + 1;
    fprintf("Refinement %d \n",irefine);
    
    sig1_old = sig(1:S1.npts);
    sig2_old = sig(S1.npts+1:end);

    % refined surface    
    [S1_new,interp_1refine_mat] = split_patches(S1,refinelist1);
    [S2_new,interp_2refine_mat] = split_patches(S2,refinelist2);

    % upsampled density    
    sig1_interp_old = interp_1refine_mat*sig1_old;
    sig2_interp_old = interp_2refine_mat*sig2_old;

    rhs1 = interp_1refine_mat*rhs1;
    rhs2 = interp_2refine_mat*rhs2;
    rhs = [rhs1;rhs2];

    S1_old = S1;
    S2_old = S2;
    S1 = S1_new;
    S2 = S2_new;
    S = merge([S1,S2]);

    % compute actual error by comparing with manufactured solution
    
if strcmp(errtype,'functional')
    pol1_interp = gfun(sig1_interp_old,S1,1:S1.npts);
    interp_actual = abs(pol1_interp-pol1_sol);
    interp_actual_rel = interp_actual/abs(pol1_sol);
else

    layerdens_interp = [sig1_interp_old;sig2_interp_old];

    [interp_actual,interp_actual_rel,~,~,~,~,~] = density_diagnostics( ...
        S1_old,S2_old,S1sol,S2sol,sig,fun_monitor_sol,lptype);

end

    % solve BIE
    [sig,~] = lap3d.solver( ...
        S,'dir',rhs,tol_sol,rep_pars);

    % estimate error. Pointwise difference
    sig_err1_old = sig(1:S1.npts) - sig1_interp_old;
    sig_err2_old = sig(S1.npts+1:end) - sig2_interp_old;

    if strcmp(errtype,'functional')

        % Actual functional error on S1.
        pol1 = gfun(sig(1:S1.npts),S1,1:S1.npts);
        pol1_sol = gfun(sig_sol(1:S1sol.npts),S1sol,1:S1sol.npts);

        % error of new solution, compared to manufactured solution
        sig_lp_abs_comp = abs(pol1-pol1_sol);
        sig_lp_rel_comp = sig_lp_abs_comp/abs(pol1_sol);
        norm_name = 'functional';

        % Patchwise estimated functional error.
        errlp1 = get_patch_functional(S1,sig_err1_old,gfun);
        errlp2 = get_patch_functional(S2,sig_err2_old,gfun);

        % Global estimated functional error.
        density_est_abs = abs(gfun(sig_err1_old,S1,1:S1.npts));
        density_lp = abs(gfun(sig(1:S1.npts),S1,1:S1.npts));

    else

        % error of new solution, compared to manufactured solution
        [sig_lp_abs_comp,sig_lp_rel_comp,~,~,siglpS1,siglpS2, ...
            norm_name] = density_diagnostics( ...
            S1,S2,S1sol,S2sol,sig,fun_monitor_sol,lptype);

        % Patchwise estimated density error.
        errlp1 = testfuns.get_patch_lp_norm( ...
            S1,sig_err1_old.',[],lptype);
        errlp2 = testfuns.get_patch_lp_norm( ...
            S2,sig_err2_old.',[],lptype);

        if lptype == 0
            density_est_abs = max([errlp1,errlp2],[],2);
            density_lp = max(siglpS1,siglpS2);
        else
            density_est_abs = ( ...
                sum(errlp1.^lptype,2) + ...
                sum(errlp2.^lptype,2)).^(1/lptype);

            density_lp = ( ...
                siglpS1.^lptype + siglpS2.^lptype).^(1/lptype);
        end

    end

    density_est_rel = density_est_abs/density_lp;

    fprintf('--------\n')

    fprintf('Actual absolute %s error of old solution    = %.16e\n',norm_name,interp_actual)
    fprintf('Estimated absolute %s error of old solution = %.16e\n',norm_name,density_est_abs)
    fprintf('Actual relative %s error of old solution    = %.16e\n',norm_name,interp_actual_rel)
    fprintf('Estimated relative %s error of old solution = %.16e\n',norm_name,density_est_rel)

    fprintf('Actual absolute %s error of new solution    = %.16e\n',norm_name,sig_lp_abs_comp)
    fprintf('Actual relative %s error of new solution    = %.16e\n',norm_name,sig_lp_rel_comp)


    fun1 = sig(1:S1.npts).';
    fun2 = sig(S1.npts+1:end).';

    if strcmp(errtype,'functional')
        siglpS1 = abs(gfun(sig(1:S1.npts),S1,1:S1.npts));
        siglpS2 = abs(gfun(sig(1+S1.npts:end),S2,1:S2.npts));
    else
        siglpS1 = comp_ffun_lp(S1,fun1,lptype);
        siglpS2 = comp_ffun_lp(S2,fun2,lptype);
    end

    [refinelist1_orig,errqlist1] = find_patch_refine_list( ...
        S1,errlp1,siglpS1,tol_dens,[],lptype,eta_dens);

    [refinelist2_orig,errqlist2] = find_patch_refine_list( ...
        S2,errlp2,siglpS2,tol_dens,[],lptype,eta_dens);

    refinelist1 = select_top_fraction( ...
        refinelist1_orig,errqlist1,refine_fraction);
    refinelist2 = select_top_fraction( ...
        refinelist2_orig,errqlist2,refine_fraction);


    fprintf('Patch counts: S1 = %d, S2 = %d.\n', ...
    S1.npatches,S2.npatches)
    fprintf('Point counts: S1 = %d, S2 = %d.\n', ...
    S1.npts,S2.npts)

    fprintf('S1 refines %d of %d flagged patches\n', ...
    numel(refinelist1),numel(refinelist1));
    fprintf('S2 refines %d of %d flagged patches\n', ...
    numel(refinelist2),numel(refinelist2));

effectivity = interp_actual ./ density_est_abs;

fprintf('Actual(old) / hierarchical increment = %.16e\n',effectivity)

end

fprintf('AMR stopped after %d refinement iterations.\n',irefine)
fprintf('Final patch counts: S1 = %d, S2 = %d.\n', ...
    S1.npatches,S2.npatches)
fprintf('Final point counts: S1 = %d, S2 = %d.\n', ...
    S1.npts,S2.npts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ffun_lp] = comp_ffun_lp(S,ffun,lptype)

if lptype == 0
    ffun_lp = max(abs(ffun),[],2);
else
    ffun_lp = sum(abs(ffun).^lptype.*S.wts(:).',2).^(1/lptype);
end

end

function [density_lp_abs,density_lp_rel,fun1,fun2,funlpS1,funlpS2, ...
    norm_name] = density_diagnostics( ...
    S1,S2,S1sol,S2sol,layerdens,fun_monitor_sol,lptype)

fun_monitor = layerdens.';

fun1 = fun_monitor(:,1:S1.npts);
fun2 = fun_monitor(:,S1.npts+1:end);

fun1sol = fun_monitor_sol(:,1:S1sol.npts);
fun2sol = fun_monitor_sol(:,S1sol.npts+1:end);

% interp_lp_error maps data from its first surface onto its second.
[~,fun_int2S1] = testfuns.interp_lp_error( ...
    S1,S1sol,fun1,fun1sol);
[~,fun_int2S2] = testfuns.interp_lp_error( ...
    S2,S2sol,fun2,fun2sol);

fun_diff1 = fun_int2S1 - fun1;
fun_diff2 = fun_int2S2 - fun2;

errlp_patch1 = testfuns.get_patch_lp_norm( ...
    S1,fun_diff1,[],lptype);
errlp_patch2 = testfuns.get_patch_lp_norm( ...
    S2,fun_diff2,[],lptype);

if lptype == 0
    density_lp_abs = max([errlp_patch1,errlp_patch2],[],2);
    density_ref_lp = max(abs([fun_int2S1,fun_int2S2]),[],2);

    funlpS1 = max(abs(fun1),[],2);
    funlpS2 = max(abs(fun2),[],2);
    norm_name = 'Linf';
else
    density_lp_abs = ( ...
        sum(errlp_patch1.^lptype,2) + ...
        sum(errlp_patch2.^lptype,2)).^(1/lptype);

    density_ref_lp = ( ...
        sum(abs(fun_int2S1).^lptype.*S1.wts(:).',2) + ...
        sum(abs(fun_int2S2).^lptype.*S2.wts(:).',2) ...
        ).^(1/lptype);

    funlpS1 = sum( ...
        abs(fun1).^lptype.*S1.wts(:).',2).^(1/lptype);
    funlpS2 = sum( ...
        abs(fun2).^lptype.*S2.wts(:).',2).^(1/lptype);
    norm_name = sprintf('L%d',lptype);
end

density_lp_rel = zeros(size(density_lp_abs));
nonzero_ref = density_ref_lp ~= 0;
density_lp_rel(nonzero_ref) = ...
    density_lp_abs(nonzero_ref)./density_ref_lp(nonzero_ref);
density_lp_rel(~nonzero_ref & density_lp_abs ~= 0) = Inf;

end

function [refinelist,errqlist] = flag_refinement_candidates( ...
    S,fun,funlp,tol_dens,lptype,eta_dens)

% if lptype == 0
%     tail_errs = S.surf_fun_error(fun);
% else
%     tail_errs = S.surf_fun_error(fun,lptype);
%     tail_errs = surf_fun_tail_val_error(S, fun, lptype, norder-3);
% end

tail_errs = surf_fun_tail_val_error(S, fun, lptype);

[refinelist,errqlist] = find_patch_refine_list( ...
    S,tail_errs,funlp,tol_dens,[],lptype,eta_dens);

refinelist = refinelist(:);

end

function selected = select_top_fraction(refinelist,errqlist,fraction)

if isempty(refinelist)
    selected = refinelist;
    return
end

% Support either one score per flagged patch or one score per mesh patch.
if numel(errqlist) == numel(refinelist)
    scores = errqlist(:);
elseif numel(errqlist) >= max(refinelist)
    scores = errqlist(refinelist);
    scores = scores(:);
else
    error(['errqlist is incompatible with refinelist: expected one ', ...
        'score per flagged patch or one score per mesh patch.'])
end

[~,order] = sort(scores,'descend');
nselect = max(1,ceil(fraction*numel(refinelist)));
selected = refinelist(order(1:nselect));

end

function print_density_error(iteration,norm_name,abs_error,rel_error)

fprintf('AMR iteration %d: absolute density %s error:', ...
    iteration,norm_name)
fprintf(' %.16e',abs_error)
fprintf('\n')

fprintf('AMR iteration %d: relative density %s error:', ...
    iteration,norm_name)
fprintf(' %.16e',rel_error)
fprintf('\n')

end

function plot_refinement_lists( ...
    S1,S2,refinelist1,refinelist2,ifplot,iteration)

if ~ifplot
    return
end

plotlist1 = zeros(S1.npatches,1);
plotlist2 = zeros(S2.npatches,1);
plotlist1(refinelist1) = 1;
plotlist2(refinelist2) = 1;

figure(3)
clf
plot(S1,plotlist1)
hold on
plot(S2,plotlist2)
hold off
colorbar
title(sprintf('Patches selected after AMR iteration %d',iteration))

end

function err = get_patch_functional(S,fun,gfun)

[~,~,~,ixyzs,~,~] = extract_arrays(S);
err = zeros(1,S.npatches);

for ip = 1:S.npatches
    idx = ixyzs(ip):(ixyzs(ip+1)-1);
    err(ip) = abs(gfun(fun,S,idx));
end

end



function [errps] = surf_fun_tail_val_error(obj, fun, p, norderhead)
%  SURF_FUN_TAIL_VAL_ERROR  Per-patch L^p norm of the high-degree tail,
%  measured in VALUE space rather than coefficient space.
%
%  Same tail selection as SURF_FUN_ERROR, but the tail coefficients are
%  mapped back to the nodes and the norm is taken over values. The result
%  is a genuine L^p(P_j) norm, directly comparable to
%  testfuns.get_patch_lp_norm.  No (|P_j|/|S|)^(1/p) scaling is applied.
%
%    p = 0 or Inf : max |tail| over the patch nodes
%    p >= 1       : ( sum |tail|^p .* wts )^(1/p)

    if nargin < 3 || isempty(p), p = 0; end
    if nargin < 4, norderhead = []; end

    npatches = obj.npatches;

    fcoefs = obj.vals2coefs(fun);
    if size(fcoefs,2) == obj.npts
        fcoefs = fcoefs.';
    end
    nfuns = size(fcoefs,2);

    errps = zeros(npatches, nfuns);

    ntmp = [obj.norders(:), obj.iptype(:)];
    [ntmp_uni,~,intmp] = unique(ntmp,'rows');

    for i = 1:size(ntmp_uni,1)
        norder = ntmp_uni(i,1);
        iptype = ntmp_uni(i,2);
        inds   = find(intmp == i);
        indv   = obj.ixyzs(inds);

        if isempty(norderhead)
            head = norder - 1;
        else
            head = norderhead;
        end

        if iptype == 1
            npols = (norder+1)*(norder+2)/2;
            ihead = 1:((head+1)*(head+2)/2);
            bmat  = koorn.coefs2vals(norder, koorn.rv_nodes(norder));
        elseif iptype == 11 || iptype == 12
            npols = (norder+1)^2;
            [N1,N2] = ndgrid(0:norder, 0:norder);
            ihead   = find(N1(:) + N2(:) <= head);
            if iptype == 11
                bmat = polytens.lege.coefs2vals(norder, polytens.lege.nodes(norder));
            else
                bmat = polytens.cheb.coefs2vals(norder, polytens.cheb.nodes(norder));
            end
        else
            error('SURF_FUN_TAIL_VAL_ERROR: unsupported iptype %d', iptype);
        end

        [I,V] = meshgrid(indv, 1:npols);
        itot  = I(:) + V(:) - 1;

        cf = reshape(fcoefs(itot,:), npols, numel(inds), nfuns);
        cf(ihead,:,:) = 0;                       % keep tail only

        v = bmat * reshape(cf, npols, []);       % tail values at nodes
        v = reshape(v, npols, numel(inds), nfuns);

        if p == 0 || isinf(p)
            e = reshape(max(abs(v),[],1), numel(inds), nfuns);
        else
            w = [obj.weights{inds}];             % (npols, ninds)
            e = reshape(sum(abs(v).^p .* w, 1), numel(inds), nfuns).^(1/p);
        end

        errps(inds,:) = e;
    end

    errps = errps.';
end
