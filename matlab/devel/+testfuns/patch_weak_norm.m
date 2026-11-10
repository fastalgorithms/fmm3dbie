function [errlp_patch,utot] = patch_weak_norm(S,sigd,Souter,rep_pars,lptype)
%PATCH_WEAK_NORM  Per-patch weak (potential) error indicator.
%
%   [errlp_patch,utot] = patch_weak_norm(S,sigd,Souter,rep_pars,lptype)
%
%   For every patch J of S, evaluates on Souter the combined-field potential
%   radiated by sigd restricted to patch J,
%
%       u_J(x) = rep_pars(1)*S[sigd 1_J](x) + rep_pars(2)*D[sigd 1_J](x),
%
%   and returns errlp_patch(J) = comp_Lp(Souter,u_J,lptype).
%
%   Smooth quadrature only (S.wts): Souter must lie outside the near field of
%   every patch of S.  One lap3d.kern call per patch, over that patch's nodes
%   only, so total work is S.npts*Souter.npts kernel evaluations and peak
%   memory is Souter.npts * (nodes per patch).
%
%   INPUTS
%     S         surfer, current (coarse) mesh
%     sigd      S.npts x 1, density DIFFERENCE on S
%     Souter    surfer, outer ellipsoid; Souter.r targets, Souter.wts weights
%     rep_pars  [alpha beta], potential = alpha*S + beta*D
%     lptype    0 = Linf, otherwise p
%
%   OUTPUTS
%     errlp_patch  1 x S.npatches, same orientation as the other estimators
%     utot         Souter.npts x 1, sum_J u_J.  Check once against
%                  lap3d.dirichlet.eval(S,sigd,targinfo_outer,tol_sol,rep_pars).
 
alpha = rep_pars(1);
beta  = rep_pars(2);
[srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
str = wts(:).*sigd(:);          % npts x 1, weight * density
 
targinfo   = [];
targinfo.r = Souter.r;
 
ixyzs    =ixyzs(:);
npatches = S.npatches;
 
errlp_patch = zeros(1,npatches);
utot = zeros(size(Souter.r,2),1);

for J = 1:npatches
    idx = ixyzs(J):ixyzs(J+1)-1;
 
    srcinfo   = [];
    srcinfo.r = S.r(:,idx);
    srcinfo.n = srcvals(10:12,idx);
    uJ = lap3d.kern(srcinfo,targinfo,'c',alpha,beta)*str(idx);
 
    errlp_patch(J) = testfuns.comp_Lp(Souter,uJ,lptype);
    utot = utot + uJ;
end
 
end
