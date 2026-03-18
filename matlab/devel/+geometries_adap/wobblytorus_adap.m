function S = wobblytorus_adap(radii, nosc, scales, nuv, norder, iptype,...
    tol, lptype, eta)
% STARTORUS Create discretized torus surfer deformed toroidal Fourier mode.
%
% Surface is parameterized by
%
% x(u,v) = \rho(u) \cos(v) scales(1)
% y(u,v) = \rho(u) \sin(v) scales(2)
% z(u,v) = z(u) scales(3)
%
% u,v \in [0, 2\pi)^2, and
%
% \rho(u) = Rmajor +rwave*cos(nosc*v) + r(u)*cos(u)
% z(u) = r(u)*sin(u)
%
% r(u) = rminor;
%
%  Syntax
%   S = geometries.wobbletorus(radii)
%   S = geometries.wobbletorus(radii, nosc)
%   S = geometries.wobbletorus(radii, nosc, scales)
%   S = geometries.wobbletorus(radii, nosc, scales, nuv)
%   S = geometries.wobbletorus(radii, nosc, scales, nuv, norder)
%   S = geometries.wobbletorus(radii, nosc, scales, nuv, norder, iptype)
%   If arguments nosc, scales, nuv, norder, and/or iptype are empty, defaults are used 
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
%        patches. Patches are subdivided until both nuv and tol are
%        satisfied
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
%    
%    
% Example
%   % create 8th-order GL-quad patch with 5-pointed star with surface
%   resolved to an error < 1e-5:
%   S = geometries.wobbletorus([1,0.5,0.05], 5, [], [], 8, 11,1e-5)
%
  if nargin < 2 || isempty(nosc)
     nosc = 1;
  end

  if nargin < 3 || isempty(scales)
    scales = [1;1;1];
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

  
  iort = -1;
  nd = 3;
  
  dpars = zeros(1000,1);
  
  dpars(1) = radii(1);
  dpars(2) = radii(2);
  dpars(3) = radii(3);
  dpars(4) = scales(1);
  dpars(5) = scales(2);
  dpars(6) = scales(3);
  dpars(7) = nosc;

  ipars = zeros(100,1);
  ipars(1) = 1;
  ipars(2) = nd;
  ipars(3) = nosc;

  ipars(5) = 1;


  zpars = complex(0.0);

  S = surf_mesh_gen_adap(dpars, zpars, ipars, iort, nuv, norder, iptype,...
      tol, lptype, eta);

end
%
%
%
%
%
%--------------------------------------------
