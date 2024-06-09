function [S] = sphere(a, na, c0, norder, iptype)
%GEOMETRIES.SPHERE  Create a discretized sphere surfer object.
%
%  Syntax
%   S = geometries.sphere(a)
%   S = geometries.sphere(a, na)
%   S = geometries.sphere(a, na, c0)
%   S = geometries.sphere(a, na, c0, norder)
%   S = geometries.sphere(a, na, c0, norder, iptype)
%   If any argument na, c0, norder is empty, default is used
%  
%  Input arguments:
%    * a: radius of the sphere 
%    * na: (optional, 2) 
%        number of patches along the coordinate directions of the
%        cube from which the sphere is constructed
%    * c0(3): (optional, [0,0,0])
%        center of the sphere
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
% Example
%  S = geometries.sphere(1.0, [], [], 10, 11)  % unit sphere, big 10th-ord GL-quad patches

  if nargin < 2 || isempty(na)
    na = 2;
  end

  if nargin < 3 || isempty(c0)
    c0 = [0; 0; 0];
  end

  if nargin < 4 || isempty(norder)
    norder = 4;
  end

  if nargin < 5
    iptype = 1;
  end

  abc = [a; a; a];
  nabc = [na; na; na];

  S = geometries.ellipsoid(abc, nabc, c0, norder, iptype);

end
%
%
%
%
%--------------------------------------------
