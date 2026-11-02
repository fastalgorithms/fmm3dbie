function [S] = square(npars, norder, iptype, iort)
% SQUARE Create discretized flat square [-1,1]^2 in the z = 0 plane
%
% The square is subdivided into npars(1)*npars(2) rectangles, with
% npars(1) intervals in the x direction and npars(2) in the y direction.
% If iptype is 1, each rectangle is further divided into two triangular
% patches.
%
%  Syntax
%   S = geometries.square()
%   S = geometries.square(npars)
%   S = geometries.square(npars, norder)
%   S = geometries.square(npars, norder, iptype)
%   S = geometries.square(npars, norder, iptype, iort)
%
%   If arguments npars, norder, iptype, and/or iort are empty,
%   defaults are used
%
%  Input arguments:
%    * npars(2): integer (optional, [3,3])
%        number of rectangles in the x and y directions. A scalar is
%        used for both directions.
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
%    * iort: (optional, -1)
%         orientation, the normals are (0,0,1) if iort > 0 and
%         (0,0,-1) otherwise (the default, as in geometries.disk)
%
% Example
%   % create 8th-order GL-quad patches for the square, 4 per side:
%   S = geometries.square(4, 8, 11)
%
% See also GEOMETRIES.DISK
%
  if nargin < 1 || isempty(npars)
    npars = [3;3];
  end
  if isscalar(npars)
    npars = [npars; npars];
  end

  if nargin < 2 || isempty(norder)
    norder = 4;
  end

  if nargin < 3 || isempty(iptype)
    iptype = 1;
  end

  if nargin < 4 || isempty(iort)
    iort = -1;
  end

  nx = npars(1); ny = npars(2);
  hx = 2/nx; hy = 2/ny;

  % reference nodes, and the affine maps (a, b) of each reference patch
  % to the unit cell [0,1]^2: cell coordinates = a + b*[u;v]
  if iptype == 1
    uvs  = koorn.rv_nodes(norder);
    amap = {[0;0], [1;1]};
    bmap = {eye(2), -eye(2)};
  elseif iptype == 11
    uvs  = polytens.lege.nodes(norder);
    amap = {[0.5;0.5]};
    bmap = {0.5*eye(2)};
  elseif iptype == 12
    uvs  = polytens.cheb.nodes(norder);
    amap = {[0.5;0.5]};
    bmap = {0.5*eye(2)};
  else
    error('GEOMETRIES.SQUARE: unsupported patch type iptype = %d.', iptype);
  end

  % du x dv points along +z for the maps above, swap the roles of u and v in the
  % maps to flip it
  if iort < 0
    for k = 1:numel(bmap)
      bmap{k} = bmap{k}(:,[2,1]);
    end
  end

  npols    = size(uvs, 2);
  nsub     = numel(amap);
  npatches = nsub*nx*ny;
  npts     = npatches*npols;

  srcvals = zeros(12, npts);
  ipatch = 0;
  for j = 1:ny
    for i = 1:nx
      x0 = -1 + (i-1)*hx;
      y0 = -1 + (j-1)*hy;
      for k = 1:nsub
        st = amap{k} + bmap{k}*uvs;
        ind = ipatch*npols + (1:npols);
        srcvals(1,ind) = x0 + hx*st(1,:);
        srcvals(2,ind) = y0 + hy*st(2,:);
        srcvals(4,ind) = hx*bmap{k}(1,1);
        srcvals(5,ind) = hy*bmap{k}(2,1);
        srcvals(7,ind) = hx*bmap{k}(1,2);
        srcvals(8,ind) = hy*bmap{k}(2,2);
        ipatch = ipatch + 1;
      end
    end
  end
  srcvals(12,:) = sign(iort);

  norders = norder*ones(npatches,1);
  iptype_all = iptype*ones(npatches,1);

  S = surfer(npatches, norders, srcvals, iptype_all);

end
