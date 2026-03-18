function [errl2abs,fun_interp2targ] = interp_lp_error(S,Sref,fun,funref)
%% Compute pointwise error in fun at quadrature nodes on surface Sref
% map data on quadrature nodes on S to quadrature nodes on Sref.
% Compute L2 norm on Sref with funref at quadrature nodes on Sref as
% reference values

% find closest point settings
tol_newton = 1e-14;
maxiter_newton = 100;
% find closest point in S.r for each Sref.r
[idx_init, dists_knn] = knnsearch(S.r.', Sref.r.');
% for each Sref.r, find the closest point on S and (u,v) on the patch
[nearest_pnts,nearest_uvs,dists,flags] = ... 
    find_nearest_surf_point(S,Sref.r,idx_init,tol_newton,maxiter_newton);
% get nearest patch on S for each Sref.r
nearest_patch = S.patch_id(idx_init);
% interpolate function funref on S.r to points Sref.r
fun_interp2targ = S.interpolate_data(fun,nearest_patch,nearest_uvs);
% L2 error, integrate over Sref with funref as solution
diff_fun = abs(fun_interp2targ - funref);
errl2abs = sqrt(diff_fun.^2 *Sref.wts);


return

end
