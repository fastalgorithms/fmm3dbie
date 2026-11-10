function ref = two_sphere_reference(a,c0_1,c0_2,rep_pars,Souter,functional,opts,cached)
% Independently regenerate both analytic spheres for each reference attempt.
% Accept only when density, outer potential, and functional changes all pass.
% Relative errors use ||sig||_p + ||rhs||_p on the finer reference.
% Successive agreement is a convergence check, not a rigorous error bound.
if nargin < 8, cached = []; end
key = {a,c0_1,c0_2,rep_pars,Souter.r,Souter.wts,functions(functional),opts};
if ~isempty(cached) && isequaln(cached.key,key)
    fprintf('Reusing verified two-sphere reference (%d points).\n',cached.Ssol.npts);
    ref = cached;
    return
end
assert(opts.max_attempts >= 2,'At least two reference solves are required.');
assert(opts.tol_sol <= opts.tol_change/10, ...
    'Reference solve tolerance must be at most one tenth of the change tolerance.');
targets = struct('r',Souter.r,'patch_id',-ones(Souter.npts,1));
solve_opts = struct('eps_gmres',opts.tol_sol,'maxit',opts.maxit);
history = struct([]);
previous = [];
for attempt = 1:opts.max_attempts
    nasol = opts.nasol+attempt-1;
    nordersol = opts.nordersol+attempt-1;
    predicted_points = 24*nasol^2*(nordersol+1)*(nordersol+2)/2;
    if predicted_points > opts.max_npts
        error('two_sphere_reference:PointLimit', ...
            'Reference initial mesh needs %d points (limit %d).', ...
            predicted_points,opts.max_npts);
    end
    S1sol = geometries.sphere(a,nasol,c0_1,nordersol,1);
    S2sol = geometries.sphere(a,nasol,c0_2,nordersol,1);
    [S1sol,S2sol,refine_info] = testfuns.refine_close_patches( ...
        S1sol,S2sol,opts.prefactor,opts.max_refinements,opts.max_npts);
    if ~refine_info.converged
        error('two_sphere_reference:GeometryNotConverged', ...
            'Close-patch refinement did not converge on attempt %d.',attempt);
    end
    Ssol = merge([S1sol,S2sol]);
    rhssol = [ones(S1sol.npts,1);zeros(S2sol.npts,1)];
    fprintf('Reference attempt %d: order %d, %d patches, %d points.\n', ...
        attempt,nordersol,Ssol.npatches,Ssol.npts);
    [sig_sol,~,rres] = lap3d.solver(Ssol,'dir',rhssol,opts.tol_sol,rep_pars,solve_opts);
    if any(~isfinite(sig_sol(:))) || ~isfinite(rres) || rres > opts.tol_sol
        error('two_sphere_reference:SolveFailed', ...
            'Reference residual %.3e exceeds %.3e.',rres,opts.tol_sol);
    end
    pot_outer = lap3d.dirichlet.eval(Ssol,sig_sol,targets,opts.tol_sol,rep_pars);
    if any(~isfinite(pot_outer(:)))
        error('two_sphere_reference:EvaluationFailed','Nonfinite reference potential.');
    end
    density_change = nan(1,3);
    weak_change = nan(1,3);
    functional_change = nan(1,3);
    if ~isempty(previous)
        [~,~,sig_previous] = testfuns.interp_lp_error( ...
            previous.Ssol,Ssol,previous.sig_sol,sig_sol,2);
        delta_sig = sig_sol(:)-sig_previous(:);
        delta_pot = pot_outer(:)-previous.pot_outer(:);
        for p = 0:2
            scale = testfuns.comp_Lp(Ssol,sig_sol(:),p) + ...
                testfuns.comp_Lp(Ssol,rhssol,p);
            density_change(p+1) = testfuns.comp_Lp(Ssol,delta_sig,p)/scale;
            weak_change(p+1) = testfuns.comp_Lp(Souter,delta_pot,p)/scale;
            functional_change(p+1) = abs(functional(Souter,delta_pot))/scale;
        end
    end
    history(attempt).npts = Ssol.npts;
    history(attempt).order = nordersol;
    history(attempt).residual = rres;
    history(attempt).density_change = density_change;
    history(attempt).weak_change = weak_change;
    history(attempt).functional_change = functional_change;
    fprintf('Reference density changes [Linf L1 L2]: %.3e %.3e %.3e\n',density_change);
    fprintf('Reference weak changes [Linf L1 L2]: %.3e %.3e %.3e\n',weak_change);
    fprintf('Reference functional changes [Linf L1 L2 scaling]: %.3e %.3e %.3e\n',functional_change);
    changes = [density_change,weak_change,functional_change];
    if all(isfinite(changes)) && all(changes <= opts.tol_change)
        ref = struct('S1sol',S1sol,'S2sol',S2sol,'Ssol',Ssol, ...
            'sig_sol',sig_sol,'pot_outer',pot_outer,'history',history,'key',{key});
        return
    end
    previous = struct('Ssol',Ssol,'sig_sol',sig_sol,'pot_outer',pot_outer);
end
error('two_sphere_reference:NotConverged', ...
    'Reference did not meet %.3e in %d attempts. See the reported changes.', ...
    opts.tol_change,opts.max_attempts);
end
