function [amat] = barymat(xout, xin, w)
% BARYMAT: construct the barycentric interpolation matrix to
%  interpolate from xin nodes to xout, where w are the 
%  barycentric weights. If w is not supplied, it is assumed
%  to be weights for polynomial interpolation at chebyshev nodes
%
%  Syntax:
%    amat = barymat(xout, xin);
%    amat = barymat(xout, xin, w)


  n = length(xin);
  m = length(xout);

  if nargin < 3
    w = ones(1,n);
    w(2:2:end) = -1;
    w([1,n]) = 0.5*w([1,n]);
  else
    w = reshape(w, 1, n);
  end

  amat = xout - xin.';
  amat = w./amat;
  rnorm = 1./sum(amat, 2);
  amat = amat.*rnorm;
  amat(isnan(amat)) = 1; 

end
