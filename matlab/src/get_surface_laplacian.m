function [surf_lap_p] = get_surface_laplacian(S,p)
%
%  
%  surf_lap_p = get_surface_laplacian(S,p)
%    This subrorutine evaluates the surface laplacian of a given 
%    function p
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * p: input function on surface
%
%  Output arguments:
%    * surf_lap_p: surface laplcian of p
%

% Extract arrays
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;

    nd0 = 1;
    n2 = 2;


    dp = zeros(2,npts);
    mex_id_ = 'get_surf_grad(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c io double[xx])';
[dp] = fmm3dbie_routs(mex_id_, nd0, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, p, dp, 1, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, npts, n2, npts);

    surf_lap_p = zeros(size(p));
    mex_id_ = 'get_surf_div(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[xx], c io double[x])';
[surf_lap_p] = fmm3dbie_routs(mex_id_, nd0, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, dp, surf_lap_p, 1, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, n2, npts, npts);

    surf_lap_p = reshape(surf_lap_p,size(p));

end
%
%
%
%
%
