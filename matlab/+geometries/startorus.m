function [S] = startorus(radii, nosc, scales, nuv, norder, iptype)
% GEOMETRIES.startorus, get toroidal double fourier surface given by
%
% x(u,v) = \rho(u) \cos(v) scales(1)
% y(u,v) = \rho(u) \sin(v) scales(2)
% z(u,v) = z(u) scales(3)
%
% u,v \in [0, 2\pi)^2, and
%
% \rho(u) = Rmajor + r(u)*cos(u)
% z(u) = r(u)*sin(u)
%
% r(u) = rminor + rwave*cos(nosc*u)
%
%  Syntax
%   S = geometries.startorus(radii)
%   S = geometries.startorus(radii, nosc)
%   S = geometries.startorus(radii, nosc, scales)
%   S = geometries.startorus(radii, nosc, scales, nuv)
%   S = geometries.startorus(radii, nosc, scales, nuv, norder)
%   S = geometries.startorus(radii, nosc, scales, nuv, norder, iptype)
%
%  Input arguments:
%    * radii(3): radii of star shaped torus
%         radii(1) = rmajor
%         radii(2) = rminor
%         radii(3) = rwave
%    * nosc(1): integer (optional, 0)
%         number of oscillations in star-shaped torus
%    * scales(3): (optional, [1,1,1])
%        scaling parameters for the coordinates of the surface
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

  if nargin < 2
     nosc = 0;
  end

  if nargin < 3
    scales = [1;1;1];
  end

  if nargin < 4
    nuv = [5;5];
  end

  if nargin < 5
    norder = 4;
  end

  if nargin < 6
    iptype = 1;
  end

  m = nosc + 1;
  coefs = zeros(2*m+1, 2*m+1, 3);

  coefs(1,1,1) = radii(1);
  coefs(2,1,1) = radii(2);
  coefs(1+nosc,1,1) = radii(3)/2;
  coefs(1+abs(nosc-1),1,1) = radii(3)/2;

  coefs(m+2,1,3) = radii(2);
  coefs(2*m+1,1,3) = radii(3)/2;

  if nosc-1 > 0
    coefs(m+1+nosc-1,1,3) = -radii(3)/2;
  end

  if nosc-1 < 0
    coefs(m+1+abs(nosc-1),1,3) = radii(3)/2;
  end

  iort = 1;

  S = geometries.xyz_tensor_fourier(coefs, scales, iort, nuv, norder, iptype);


end
%
%
%
%
%
%--------------------------------------------
