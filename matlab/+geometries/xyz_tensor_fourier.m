function [S] = xyz_tensor_fourier(coefs, nfp, scales, iort, nuv, norder, iptype)
% GEOMETRIES.xyz_tensor_fourier, get toroidal double fourier surface given by
%
% hat(x) = \sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} x_{ij} b_{i} (u) b_{j} (nfp*v)
% hat(y) = \sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} y_{ij} b_{i} (u) b_{j} (nfp*v)
% hat(z) = \sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} z_{ij} b_{i} (u) b_{j} (nfp*v)
%
% x(u,v) = (hat(x) \cos(v) - hat(y) sin(v))*scales(1)
% y(u,v) = (hat(x) \sin(v) + hat(y) cos(v))*scales(2)
% z(u,v) = hat(z)*scales(3)
%
% u,v \in [0, 2\pi)^2, and
%
% b_{i}(u) = \{ 1, cos(u), \ldots cos(m u), sin(u), \ldots sin(m u) \}
%
%
%  Syntax
%   S = geometries.xyz_tensor_fourier(coefs)
%   S = geometries.xyz_tensor_fourier(coefs, nfp)
%   S = geometries.xyz_tensor_fourier(coefs, nfp, scales)
%   S = geometries.xyz_tensor_fourier(coefs, nfp, scales, iort)
%   S = geometries.xyz_tensor_fourier(coefs, nfp, scales, iort, nuv)
%   S = geometries.xyz_tensor_fourier(coefs, nfp, scales, iort, nuv, norder)
%   S = geometries.xyz_tensor_fourier(coefs, nfp, scales, iort, nuv, norder, iptype)
%   If arguments scales, iort, nuv, norder, and/or iptype are empty, defaults are used 
%
%  Input arguments:
%    * coefs(2*m+1, 2*m+1, 3): Fourier coefs
%         hat(x) = coefs(:,:,1)
%         hat(y) = coefs(:,:,2)
%         hat(z) = coefs(:,:,3)
%    * nfp: integer
%         number of field periods
%    * scales(3): (optional, [1,1,1])
%        scaling parameters for the coordinates of the surface
%    * iort: integer (optional, 1)
%        orientation flag
%        if iort = 1, then parameter space is [0, 2\pi)^2
%        if iort = -1, then parameter space is [2\pi, 0) \times [0, 2\pi)
%    * nuv(2): integer (optional, [5,5])
%        number of quad patches in u, and v direction. if iptype is 1,
%        then each quad patch is further divided into two triangular
%        patches
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
%

  if nargin < 2 || isempty(scales)
    scales = [1;1;1];
  end
 
  if nargin < 3 || isempty(iort)
    iort = 1; 
  end

  if nargin < 4 || isempty(nuv)
    nuv = [5;5];
  end

  if nargin < 5 || isempty(norder)
    norder = 4;
  end

  if nargin < 6 || isempty(iptype)
    iptype = 1;
  end

  mm = length(size(coefs));
  if mm ~= 3
    error('GEOMETRIES.XYZ_TENSOR_FOURIER: coefs must be a three dimensional tensor');
  end

  [m1, m2, m3] = size(coefs);
  if (m1 ~= m2)
    error('GEOMETRIES.XYZ_TENSOR_FOURIER: first two dimensions of coefs must be of same size');
  end

  if mod(m1, 2) ~= 1
    error('GEOMETRIES.XYZ_TENSOR_FOURIER: first two dimensions of coefs must be odd');
  end

  if m3 ~=3 
    error('GEOMETRIES.XYZ_TENSOR_FOURIER: last dimension of coefs array must be 3');
  end

  m = (m1-1)/2;

  muse = m1*m2;
  coefs_use = reshape(coefs, [muse, 3]);

  npatches = 0;
  npts = 0;
  
  mex_id_ = 'get_xyz_tensor_fourier_npat_mem(c i double[xx], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io int64_t[x], c io int64_t[x])';
[npatches, npts] = fmm3dbie_routs(mex_id_, coefs_use, m, nfp, scales, iort, nuv, norder, iptype, npatches, npts, muse, 3, 1, 1, 3, 1, 2, 1, 1, 1, 1);
  

  srcvals = zeros(12,npts);
  srccoefs = zeros(9,npts);

  npatp1 = npatches + 1;
  ixyzs = zeros(npatp1,1);
  norders = zeros(npatches,1);
  iptypes = zeros(npatches,1);

  n9 = 9;
  n12 = 12;

  mex_id_ = 'get_xyz_tensor_fourier_npat(c i double[xx], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io int64_t[x], c io int64_t[x], c io int64_t[x], c io double[xx], c io double[xx])';
[norders, ixyzs, iptypes, srccoefs, srcvals] = fmm3dbie_routs(mex_id_, coefs_use, m, nfp, scales, iort, nuv, norder, iptype, npatches, npts, norders, ixyzs, iptypes, srccoefs, srcvals, muse, 3, 1, 1, 3, 1, 2, 1, 1, 1, 1, npatches, npatp1, npatches, n9, npts, n12, npts);

  S = surfer(npatches, norder, srcvals, iptype);


end
%
%
%
%--------------------------------------------
