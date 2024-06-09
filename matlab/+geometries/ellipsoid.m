function [S] = ellipsoid(abc, nabc, c0, norder, iptype)
% ELLIPSOID Create discretized axes-aligned ellipsoid surfer object.
%
%  Syntax
%   S = geometries.ellipsoid(abc)
%   S = geometries.ellipsoid(abc, nabc)
%   S = geometries.ellipsoid(abc, nabc, c0)
%   S = geometries.ellipsoid(abc, nabc, c0, norder)
%   S = geometries.ellipsoid(abc, nabc, c0, norder, iptype)
%   If arguments nabc, c0, norder, and/or iptype are empty, defaults are used 
%
%  Input arguments:
%    * abc(3): semi-major axes in x,y, and z directions 
%    * nabc(3): (optional, [2,2,2])
%        number of patches along the coordinate directions of the
%        cube from which the ellipsoid is constructed
%    * c0(3): (optional, [0,0,0])
%        center of the ellipsoid
%    * norder: (optional, 4)
%        order of discretization
%    * iptype: (optional, 1)
%        type of patch to be used in the discretization
%        * iptype = 1, triangular patch discretized using
%                      Rokhlin-Vioreanu nodes
%        * iptype = 11, quadrangular patch discretized using tensor
%                       product Gauss-Legendre nodes
%        * iptype = 12, quadrangular patch discretized using tensor
%                       product Chebyshev 
% Example
%   S = geometries.ellipsoid([1.5,1,0.5], [5 3 2], [], 8, 11)
%
  if nargin < 2 || isempty(nabc)
    nabc = [2; 2; 2];
  end

  if nargin < 3 || isempty(c0)
    c0 = [0; 0; 0];
  end

  if nargin < 4 || isempty(norder)
    norder = 4;
  end

  if nargin < 5 || isempty(iptype)
    iptype = 1;
  end

  npatches = 0;
  npts = 0;
  
  mex_id_ = 'get_ellipsoid_npat_mem(i double[x], i int[x], i double[x], i int[x], i int[x], io int[x], io int[x])';
[npatches, npts] = fmm3dbie_routs(mex_id_, abc, nabc, c0, norder, iptype, npatches, npts, 3, 3, 3, 1, 1, 1, 1);
  

  srcvals = zeros(12,npts);
  srccoefs = zeros(9,npts);

  npatp1 = npatches + 1;
  ixyzs = zeros(npatp1,1);
  norders = zeros(npatches,1);
  iptypes = zeros(npatches,1);

  n9 = 9;
  n12 = 12;

  mex_id_ = 'get_ellipsoid_npat(i double[x], i int[x], i double[x], i int[x], i int[x], i int[x], i int[x], io int[x], io int[x], io int[x], io double[xx], io double[xx])';
[norders, ixyzs, iptypes, srccoefs, srcvals] = fmm3dbie_routs(mex_id_, abc, nabc, c0, norder, iptype, npatches, npts, norders, ixyzs, iptypes, srccoefs, srcvals, 3, 3, 3, 1, 1, 1, 1, npatches, npatp1, npatches, n9, npts, n12, npts);

  S = surfer(npatches, norder, srcvals, iptype);


end
%
%
%-------------------------------------------------------------
