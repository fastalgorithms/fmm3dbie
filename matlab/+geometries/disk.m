function [S] = disk(scales, rmid, npars, norder, iptype, iort)
% DISK Create discretized scaled unit disk 
%
% Surface is the disk scaled by scales(1) in the first coordinate,
% and scales(2) in the second cooridnates. 
%
% In the unscaled version, the interior [-rmid, rmid]^2 is first
% subdivided into npars(1)*npars(1) quads. The annular region between
% the square and the boundary of the disk is subdivided into 4 sectors,
% and each sector is discretized into further npars(2)*npars(3) squares
% with npars(2) the number of intervals in the radial direction
% and npars(3) the number of intervals in the angular direction
%
%
%  Syntax
%   S = geometries.disk()
%   S = geometries.disk(scales)
%   S = geometries.disk(scales, rmid)
%   S = geometries.disk(scales, rmid, npars)
%   S = geometries.disk(scales, rmid, npars, norder)
%   S = geometries.disk(scales, rmid, npars, norder, iptype)
%   S = geometries.disk(scales, rmid, npars, norder, iptype, iort)
%
%   If arguments scales, rmid, npars, norder, iptype, and/or iort are empty,
%   defaults are used 
%
%  Input arguments:
%    * scales(3): (optional, [1,1])
%        scaling parameters for the coordinates of the disk 
%    * rmid: (optional, 1/2)
%        size of interior square in the discretization
%    * npars(3): integer (optional, [3,3,3])
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
%    * iort: (optional, 1)
%         orientation vector, normals point in the unbounded region
%         if iort = 1, and in the interior of the startorus otherwise
%    
% Example
%   % create 8th-order GL-quad patch for an ellipse with semi-major
%   axes (1,0.5):
%   S = geometries.disk([1,0.5], [], [], [], 8, 11)
%
  if nargin < 1 || isempty(scales)
    scales = [1;1];
  end

  if nargin < 2 || isempty(rmid)
    rmid = 0.5;
  end

  if nargin < 3 || isempty(npars)
    npars = [3;3;3];
  end

  if nargin < 4 || isempty(norder)
    norder = 4;
  end

  if nargin < 5 || isempty(iptype)
    iptype = 1;
  end

  if nargin < 6 || isempty(iort)
    iort = -1;
  end


  if iptype == 1
    npatches = 2*(npars(1)*npars(1) + 4*npars(2)*npars(3))
    npts = npatches*(norder+1)*(norder+2)/2
  elseif iptype == 11 || iptype == 12
    npatches = npars(1)*npars(1) + 4*npars(2)*npars(3)
    npts = npatches*(norder+1)*(norder+1)
  end

  npp1 = npatches+1;
  ixys = zeros(npp1,1);
  ptinfo = zeros(6,npts);

  mex_id_ = 'mesh_circle_pts(i double[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], io int[x], io double[xx])';
[ixys, ptinfo] = fmm3dbie_routs(mex_id_, rmid, npars, iort, norder, iptype, npatches, npts, ixys, ptinfo, 1, 3, 1, 1, 1, 1, 1, npp1, 6, npts);

  norders = norder*ones(npatches,1);
  srcvals = zeros(12,npts);
  srcvals(1,:) = scales(1)*ptinfo(1,:);
  srcvals(2,:) = scales(2)*ptinfo(2,:);
  srcvals(4,:) = scales(1)*ptinfo(3,:);
  srcvals(5,:) = scales(2)*ptinfo(4,:);
  srcvals(7,:) = scales(1)*ptinfo(5,:);
  srcvals(8,:) = scales(2)*ptinfo(6,:);
  srcvals(12,:) = iort/abs(iort);

  iptype_all = iptype*ones(npatches,1);

  S = surfer(npatches, norders, srcvals, iptype_all);



end
%
%-------------------------------------------------------------

%
%
%
%%   Common routines
%%
%
%-------------------------------------------------
