function [S] = plane_a380(iref)
    npatches = double(7834*(4^iref));
    norder = 3;
    npols = int32((norder+1)*(norder+2)/2);
    npts = double(npatches*npols);
    srcvals = zeros(12,npts);
    srccoefs = zeros(9,npts);
    wts = zeros(npts,1);
    iptype = zeros(npatches,1);
    ixyzs = zeros(npatches+1,1);
    norders = norder*ones(npatches,1);
    npp1 = npatches+1;

    ndsv = 12;
    ndsc = 9;

    mex_id_ = 'get_plane_geom(i int[x], i int[x], i int[x], i int[x], io int[x], io int[x], io int[x], io double[xx], io double[xx], io double[x])';
[norders, ixyzs, iptype, srcvals, srccoefs, wts] = fmm3dbierouts(mex_id_, iref, npatches, norder, npts, norders, ixyzs, iptype, srcvals, srccoefs, wts, 1, 1, 1, 1, npatches, npp1, npatches, ndsv, npts, ndsc, npts, npts);

    S = [];
    S.norders = norders;
    S.npatches = npatches;
    S.npts = npts;
    S.srcvals = srcvals;
    S.srccoefs = srccoefs;
    S.iptype = iptype;
    S.wts = wts;
    S.ixyzs = ixyzs;
    

end

