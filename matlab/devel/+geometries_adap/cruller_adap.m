function S = cruller_adap(radii, scales, nuv, norder, iptype, ...
    tol, lptype, eta, nd_monitor)
% CRULLER_ADAP Create discretized cruller-like torus.
%
% Surface is parameterized by
%
%   h(s,t) = rmean + ramp*cos(freqs*s + freqt*t)
%
%   x(s,t) = sx*(rmaj + h(s,t)*cos(t))*cos(s)
%   y(s,t) = sy*(rmaj + h(s,t)*cos(t))*sin(s)
%   z(s,t) = sz*h(s,t)*sin(t)
%
% where
%
%   radii  = [rmaj; rmean; ramp; freqs; freqt]
%   scales = [sx; sy; sz]
%
% Syntax
%   S = geometries.cruller_adap(radii)
%   S = geometries.cruller_adap(radii, scales)
%   S = geometries.cruller_adap(radii, scales, nuv)
%   S = geometries.cruller_adap(radii, scales, nuv, norder)
%   S = geometries.cruller_adap(radii, scales, nuv, norder, iptype)
%   S = geometries.cruller_adap(radii, scales, nuv, norder, iptype, tol)
%   S = geometries.cruller_adap(radii, scales, nuv, norder, iptype, tol, lptype)
%   S = geometries.cruller_adap(radii, scales, nuv, norder, iptype, tol, lptype, eta)
%
% Input arguments:
%   * radii(5):
%       radii(1) = rmaj
%       radii(2) = rmean
%       radii(3) = ramp
%       radii(4) = freqs
%       radii(5) = freqt
%   * scales(3): (optional, [1;1;1])
%       scaling parameters for x,y,z coordinates
%   * nuv(2): integer (optional, [5;5])
%       number of quad patches in s and t directions. If iptype is 1,
%       each quad patch is split into two triangular patches.
%   * norder: (optional, 4)
%       order of discretization
%   * iptype: (optional, 1)
%       patch type
%         iptype = 1  : triangular Rokhlin-Vioreanu
%         iptype = 11 : tensor-product Gauss-Legendre
%         iptype = 12 : tensor-product Chebyshev
%   * tol: (optional, 1)
%       desired precision in resolving the surface. Set tol >= 1 for
%       uniform refinement.
%   * lptype: (optional, 0)
%       error estimation lp-norm
%         lptype = 0 : linf
%         lptype = 1 : l1
%         lptype = 2 : l2
%   * eta: (optional, 0)
%       scaling parameter h^eta for error estimation
%
% Example
%   S = geometries.cruller_adap([1;0.3;0.05;4;2], [1;1;1], [], 8, 11, 1e-5);

  if nargin < 2 || isempty(scales)
    scales = [1;1;1];
  end

  if nargin < 3 || isempty(nuv)
    nuv = [5;5];
  end

  if nargin < 4 || isempty(norder)
    norder = 4;
  end

  if nargin < 5 || isempty(iptype)
    iptype = 1;
  end

  if nargin < 6 || isempty(tol)
    tol = 1;
  end

  if nargin < 7 || isempty(lptype)
    lptype = 0;
  end

  if nargin < 8 || isempty(eta)
    eta = 0;
  end

  if nargin < 9 || isempty(nd_monitor)
    nd_monitor = 3;
  end

  iort = -1;
  
  dpars = zeros(1000,1);
  dpars(1) = radii(1);   % rmaj
  dpars(2) = radii(2);   % rmean
  dpars(3) = radii(3);   % ramp
  dpars(4) = radii(4);   % freqs
  dpars(5) = radii(5);   % freqt
  dpars(6) = scales(1);  % sx
  dpars(7) = scales(2);  % sy
  dpars(8) = scales(3);  % sz

  ipars = zeros(100,1);
  ipars(1) = 5;
  ipars(2) = nd_monitor;
  ipars(5) = 1; % ifplot = 1 then plot

  zpars = complex(0.0);

  S = surf_mesh_gen_adap(dpars, zpars, ipars, iort, nuv, norder, iptype, ...
      tol, lptype, eta);

end

