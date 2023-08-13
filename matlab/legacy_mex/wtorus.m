function [S] = wtorus(radii,scales,nosc,nu,nv,norder)
    npatches = double(2*nu*nv);
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

    nr = 3;

    mex_id_ = 'get_wtorus_geom(i double[x], i double[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], io int[x], io int[x], io int[x], io double[xx], io double[xx], io double[x])';
[norders, ixyzs, iptype, srcvals, srccoefs, wts] = fmm3dbierouts(mex_id_, radii, scales, nosc, nu, nv, npatches, norder, npts, norders, ixyzs, iptype, srcvals, srccoefs, wts, 3, 3, 1, 1, 1, 1, 1, 1, npatches, npp1, npatches, ndsv, npts, ndsc, npts, npts);

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

