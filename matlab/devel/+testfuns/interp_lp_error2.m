function [err_lp_abs,err_lp_rel,funref_on_S] = interp_lp_error2(S,Sref,fun,funref,lptype)

%% Compute pointwise error in fun at quadrature nodes on surface Sref
% map data on quadrature nodes on S to quadrature nodes on Sref.
% Compute L2 norm on Sref with funref at quadrature nodes on Sref as
% reference values


opts = struct( ...
    'tol',1e-12, ...
    'maxiter',100, ...
    'rfac',2);

% Locate nodes of Sref on the high-order reference surface
[~,patch_ref,uv_ref,~,flags] = ...
    get_closest_pts(Sref,S.r,opts);

if any(flags ~= 0)
    error('Closest-point search failed for %d targets.', ...
        nnz(flags ~= 0))
end

% Interpolate the high-order solution from Sref onto S
funref_on_S = Sref.interpolate_data( ...
    funref,patch_ref(:),uv_ref);

fun = fun(:);
funref_on_S = funref_on_S(:);
diff_fun = funref_on_S-fun;

if lptype == 0 || isinf(lptype)
    err_lp_abs = max(abs(diff_fun),[],1);
    err_ref_lp = max(abs(fun),[],1);
else
    err_lp_abs = sum( ...
        abs(diff_fun).^lptype.*Sref.wts(:),1) ...
        .^(1/lptype);
    err_ref_lp = sum( ...
        abs(fun).^lptype.*Sref.wts(:),1).^(1/lptype);
end

err_lp_rel = zeros(size(err_lp_abs));
nonzero_ref = err_ref_lp ~= 0;
err_lp_rel(nonzero_ref) = ...
    err_lp_abs(nonzero_ref)./err_ref_lp(nonzero_ref);
err_lp_rel(~nonzero_ref & err_lp_abs ~= 0) = Inf;


return

end
