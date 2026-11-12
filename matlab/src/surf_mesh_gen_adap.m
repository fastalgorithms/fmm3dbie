function [S,nlevels,nboxes,ltree,rintlin,nboxesleaf] = surf_mesh_gen_adap(dpars, zpars, ipars, iort, nuv, norder, iptype, tol, lptype, eta)
% surf_tree_mem, get toroidal double fourier surface given by
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
%   S = geometries.xyz_tensor_fourier(coefs, nfp, scales, iort, nuv, norder,
%   iptype, tol, eta)
%   If arguments scales, iort, nuv, norder, and/or iptype are empty, defaults are used 
%
%  Input arguments:
%    * dpars(*): reals
%         parameters for evaluation of global chart
%    * zpars(*): complex
%         parameters for evaluation of global chart
%    * ipars(*): integers
%         parameters for evaluation of global chart
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
%    * tol: (optional, 1)
%        desired precision in resolving the surface. Patches are
%        subdivided until both nuv and tol are satisfied   
%    * lptype - (optional, 0)
%           error estimation lp-norm
%           lptype = 0 - linf
%           lptype = 1 - l1
%           lptype = 2 - l2 
%    * eta: (optional 1)
%        sets scaling parameter h^eta for error estimation, where
%        h is the length of the box. 
%
 
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

  if nargin < 7 || isempty(tol)
    tol = 1;
  end

  if nargin < 8 || isempty(lptype)
    lptype = 0;
  end

  if nargin < 9 || isempty(eta)
    eta = 0;
  end


%  standard settings
  npatches = 0;
  npts = 0;

  zk = complex(0.0);
  boxlen = 1.0;

  ifun_dummy = 1;
  
  nd = ipars(2);
  nosc = ipars(3);

% output
  nlevels = 0;
  nboxes = 0;
  ltree = 0;
  rintlin = zeros(nd,201);
  nboxesleaf = 0;

mex_id_ = 104;
[nlevels, nboxes, ltree, rintlin, nboxesleaf] = fmm3dbie_routs(mex_id_, nuv, tol, zk, boxlen, norder, lptype, eta, ifun_dummy, nd, dpars, zpars, ipars, nlevels, nboxes, ltree, rintlin, nboxesleaf, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1000, 1, 100, 1, 1, 1, nd, 201, 1);


nfvals = nd*norder^2*nboxes;
nlevelsp1 = nlevels + 1;

itree = zeros(ltree,1);
iptr = zeros(8,1);
fvals = zeros(nfvals,1);
centers = zeros(2,nboxes);
boxsizein = zeros(2,nlevelsp1);

% triangles: split leaf box into 2 triangles
if iptype == 1
   npols = (norder + 1) * (norder + 2) / 2;
   npatches = nboxesleaf * 2;
   npts = npols * npatches;
end

% quads
if (iptype == 11)||(iptype == 12)
   npols = norder^2;
   npatches = nboxesleaf;
   npts = npols * npatches;	
end



mex_id_ = 105;
[nlevels, nboxes, ltree, rintlin, itree, iptr, fvals, centers, boxsizein] = fmm3dbie_routs(mex_id_, nuv, tol, zk, boxlen, norder, lptype, eta, ifun_dummy, nd, dpars, zpars, ipars, nlevels, nboxes, ltree, rintlin, itree, iptr, fvals, centers, boxsizein, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1000, 1, 100, 1, 1, 1, nd, 201, ltree, 8, nfvals, 2, nboxes, 2, nlevelsp1);

npatp1 = npatches + 1;

srcvals = zeros(12,npts);
srccoefs = zeros(9,npts);
norders = zeros(npatches,1);
ixyzs = zeros(npatp1,1);
iptype_all = zeros(npatches,1);

mex_id_ = 106;
[srcvals, srccoefs, norders, ixyzs, iptype_all] = fmm3dbie_routs(mex_id_, norder, iptype, nd, dpars, zpars, ipars, nlevels, nboxes, ltree, itree, iptr, centers, boxsizein, npatches, npts, srcvals, srccoefs, norders, ixyzs, iptype_all, 1, 1, 1, 1000, 1, 100, 1, 1, 1, ltree, 8, 2, nboxes, 2, nlevelsp1, 1, 1, 12, npts, 9, npts, npatches, npatp1, npatches);

S = surfer(npatches, norders, srcvals, iptype_all);
     
end


%
%
%
