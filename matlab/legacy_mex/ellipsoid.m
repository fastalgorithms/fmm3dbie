function [S] = ellipsoid(a,b,c,rmax,ifc,norder)
%
% rmax maximum patch size
% ifc is whether mesh is conforming or not
%
%

    npatches = 0;
    mex_id_ = 'get_ellipsoid_mem(i double[x], i double[x], i double[x], i double[x], i int[x], io int[x])';
[npatches] = fmm3dbierouts(mex_id_, a, b, c, rmax, ifc, npatches, 1, 1, 1, 1, 1, 1);
    npatches = double(npatches);
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

    mex_id_ = 'get_ellipsoid_geom(i double[x], i double[x], i double[x], i double[x], i int[x], i int[x], i int[x], i int[x], io int[x], io int[x], io int[x], io double[xx], io double[xx])';
[norders, ixyzs, iptype, srcvals, srccoefs] = fmm3dbierouts(mex_id_, a, b, c, rmax, ifc, norder, npatches, npts, norders, ixyzs, iptype, srcvals, srccoefs, 1, 1, 1, 1, 1, 1, 1, 1, npatches, npp1, npatches, ndsv, npts, ndsc, npts);
    mex_id_ = 'get_qwts(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], io double[x])';
[wts] = fmm3dbierouts(mex_id_, npatches, norders, ixyzs, iptype, npts, srcvals, wts, 1, npatches, npp1, npatches, 1, ndsv, npts, npts);
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

