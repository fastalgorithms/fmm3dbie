function [cms, rads] = get_centroid_rads(npatches,norders,ixyzs,iptype,npts, ...
      srccoefs)
% get the radius of each patch
    n3 = 3;
    [n9,~] = size(srccoefs);
    cms = zeros(3,npatches);
    rads = zeros(npatches,1);
    npatp1 = npatches+1;

    mex_id_ = 'get_centroid_rads(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], io double[xx], io double[x])';
[cms, rads] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, cms, rads, 1, npatches, npatp1, npatches, 1, n9, npts, n3, npatches, npatches);
end

