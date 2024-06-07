function [x,w,v] = pts(k,opts)
%CHEB.PTS get nodes, weights, and matrices for switching between values
% and coefficients for Chebyshev nodes of the specified order
%
% IN:
%   k - the desired order (number of points)
%   opts - options structure.
%     opts.kind, default(1)
%       first (1) or second (2) kind chebyshev nodes 
% OUT: 
%   x - k Legendre nodes on [-1,1]
%   w - the corresponding integration weights
%   v - barycentric weights
% 
%

kind = 1;
if nargin > 1
   if isfield(opts, 'kind')
     kind = opts.kind;
   end
end

if kind == 1
  theta = pi*(2*(1:k)'-1)/(2*k);
  x = -cos(theta);
  if nargout > 1
    l = floor(k/2)+1;
    kk = 0:k-l;   
    wtmp = [2*exp(1i*pi*kk/k)./(1-4*kk.^2)  zeros(1,l)];
    w = real(ifft(wtmp(1:k) + conj(wtmp(k+1:-1:2))));
  end
end

if kind == 2
  K = k-1; % different from original
  theta = pi*(0:K)'/K; x = cos(theta); x = x(end:-1:1);
  theta = theta';

  if nargout > 1
     w = zeros(1,k); ii = 2:K; v = ones(1,k-2);
     if mod(K,2)==0
       w(1) = 1/(K^2-1); w(k) = w(1);
       for l=1:K/2-1
         v = v - 2*cos(2*l*theta(ii))/(4*l^2-1); 
       end
       v = v - cos(K*theta(ii))/(K^2-1);
     else
       w(1) = 1/K^2; w(k) = w(1);
       for l=1:(K-1)/2
         v = v - 2*cos(2*l*theta(ii))/(4*l^2-1); 
       end
     end
     w(ii) = 2*v/K; 
     w = w(end:-1:1); % different from original
  end
end


if nargout > 2
  xx = x - x.';
  xx(1:k+1:end) = 1;
  ww = prod(xx,1);
  v = ww(1)./ww(:);
  vmax = max(abs(v));
  v = v/vmax;
end
