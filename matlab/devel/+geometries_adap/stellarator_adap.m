function [S] = stellarator_adap(nuv, norder, iptype, tol, lptype, eta, ...
    nd_monitor)
% GEOMETRIES.stellarator, get toroidal double fourier surface given by
%
% x(u,v) = hat(x)(u,v) \cos(v)
% y(u,v) = hat(x)(u,v) \sin(v)  
% z(u,v) = hat(z)(u,v)
%
% u,v \in [2\pi,0)\times [0,2\pi) , and
%
% hat(x) = \sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} x_{ij} b_{i} (u) b_{j} (v)
% hat(z) = \sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} z_{ij} b_{i} (u) b_{j} (v)
%
% with m = 2, and the only non-zero coefs are given by
%
% hat(x)(3,2) = 0.17
% hat(x)(5,4) = 0.17
% hat(x)(1,1) = 4.5
% hat(x)(2,1) = 0.75
% hat(x)(3,1) = 0.11
% hat(x)(2,2)  -0.38
% hat(x)(4,4) = -0.52
%
% hat(z)(5,2) = 0.17
% hat(z)(3,4) = -0.17
% hat(z)(5,1) = 0.11
% hat(z)(4,2) = 0.52
% hat(z)(2,4) = -0.38
%
%  Syntax
%   S = geometries.stellarator()
%   S = geometries.stellarator(nuv)
%   S = geometries.stellarator(nuv, norder)
%   S = geometries.stellarator(nuv, norder, iptype)
%   S = geometries.stellarator(nuv, norder, iptype, tol)
%   If arguments norder, iptype, and/or iort are empty, defaults are used 
%
%  Input arguments:
%    * nuv(2): integer (optional, [5,15])
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
%        subdivided until both nuv and tol are satisfied. Set tol >= 1
%        for uniform refinement
%    * lptype - (optional, 0)
%           error estimation lp-norm
%           lptype = 0 - linf
%           lptype = 1 - l1
%           lptype = 2 - l2 
%    * eta: (optional 1)
%        sets scaling parameter h^eta for error estimation, where
%        h is the length of the box. 
%    * nd_monitor: (optional 3)
%        sets the number of monitor functions for adaptive refinement for
%        the global chart r=r(u,v) = (x(u,v),y(u,v),z(u,v)).
%        * nd_monitor = 3, monitor functions are r
%        * nd_monitor = 9, monitor functions are r and derivatives r_u and
%                          r_v of r wrt u and v
%        * nd_monitor = 10, monitor functions are r,r_u,r_v and mean
%                          curvature


  if nargin < 1 || isempty(nuv)
    nuv = [5;15];
  end

  if nargin < 2 || isempty(norder)
    norder = 4;
  end

  if nargin < 3 || isempty(iptype)
    iptype = 1;
  end

  if nargin < 4 || isempty(tol)
      tol = 1; 
  end

  if nargin < 5 || isempty(lptype)
    lptype = 0;
  end

  if nargin < 6 || isempty(eta)
    eta = 0;
  end

  if nargin < 7 || isempty(nd_monitor)
    nd_monitor = 3;
  end

  iort = -1;

  m = 2;
  
  coefs = zeros(2*m+1, 2*m+1, 3);
  ncoefs = numel(coefs);
  
  coefs(3,2,1) = 0.17;
  coefs(5,4,1) = 0.17;
  coefs(1,1,1) = 4.5;
  coefs(2,1,1) = 0.75;
  coefs(3,1,1) = 0.11;
  coefs(2,2,1) = -0.38;
  coefs(4,4,1) = -0.52;

  coefs(5,2,3) = 0.17;
  coefs(3,4,3) = -0.17;
  coefs(4,1,3) = 1.25;
  coefs(5,1,3) = 0.11;
  coefs(4,2,3) = 0.52;
  coefs(2,4,3) = -0.38;

  scales = [1;1;1];
  
  nfp = 1;
  
  ndpars = 1000;
  dpars = cat(1,coefs(:),scales(:),zeros(ndpars-3-ncoefs,1));

  ipars = zeros(100,1);
  ipars(1) = 2;
  ipars(2) = nd_monitor;
  ipars(3) = nfp;
  ipars(4) = m;
  ipars(5) = 1; % ifplot = 1


  zpars = complex(0.0);
       
  S = surf_mesh_gen_adap(dpars, zpars, ipars, iort, nuv, norder, iptype,...
      tol, lptype, eta);

end
%
%
%
%
%--------------------------------------------
