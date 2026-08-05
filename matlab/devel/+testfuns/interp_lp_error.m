function [errl2abs,funref_on_S] = interp_lp_error(S,Sref,fun,funref)

%% Compute pointwise error in fun at quadrature nodes on surface Sref
% map data on quadrature nodes on S to quadrature nodes on Sref.
% Compute L2 norm on Sref with funref at quadrature nodes on Sref as
% reference values


opts = struct( ...
    'tol',1e-12, ...
    'maxiter',100, ...
    'rfac',2);

% Locate nodes of S on the high-order reference surface
[~,patch_ref,uv_ref,~,flags] = ...
    get_closest_pts(Sref,S.r,opts);

if any(flags ~= 0)
    error('Closest-point search failed for %d targets.', ...
        nnz(flags ~= 0))
end

% Interpolate the high-order solution from Sref onto S
funref_on_S = Sref.interpolate_data( ...
    funref,patch_ref(:),uv_ref);

diff_fun = funref_on_S-fun;

errl2abs = sqrt(sum( ...
    abs(diff_fun).^2.*S.wts(:).',2));

return

end
