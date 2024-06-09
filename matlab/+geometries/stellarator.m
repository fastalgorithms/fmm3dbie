function [S] = stellarator(nuv, norder, iptype)
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
%   If arguments norder, and/or iptype are empty, defaults are used 
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
%

  if nargin < 1 || isempty(nuv)
    nuv = [5;15];
  end

  if nargin < 2 || isempty(norder)
    norder = 4;
  end

  if nargin < 3 || isempty(iptype)
    iptype = 1;
  end

  m = 2; 
  coefs = zeros(2*m+1, 2*m+1, 3);

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


  iort = -1;

  S = geometries.xyz_tensor_fourier(coefs, scales, iort, nuv, norder, iptype);


end
%
%
%

%
%--------------------------------------------
