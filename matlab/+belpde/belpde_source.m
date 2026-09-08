function [ut,rhs] = belpde_source(S,x_source,charge,zk,type)
%BELPDE.BELPDE_SOURCE  Exact solution and right hand side for a Beltrami PDE.
%
%   [ut,rhs] = belpde.belpde_source(S,x_source,charge,zk)
%   [ut,rhs] = belpde.belpde_source(S,x_source,charge,zk,type)
%
%   Restricts a function of the ambient space to S and returns its values ut
%   together with rhs = (\Delta_\Gamma + zk^2) ut.
%
%   zk is the wavenumber, either a scalar or a function handle evaluated as
%   zk(x) on the (3,npts) node array, giving a variable wavenumber. For a
%   wavenumber that vanishes identically, ut is normalized to have zero mean
%   on the surface.
%
%   type selects the ambient function: 1 (default) a point charge at
%   x_source, 2 a Gaussian centred at x_source, 3 a spherical harmonic
%   (assumes S is a sphere).

npts = size(S.r,2);
rhs = zeros(1,npts);

if nargin < 5, type = 1; end

if isa(zk,'function_handle')
    zkvals = reshape(zk(S.r),1,[]);
else
    zkvals = zk;
end

if type == 3
    % test on spherical harmonic right hand side. assumes a sphere
    l = 3; m = 0;
    f = spherefun.sphharm(l,m); % Chebfun spherical harmonic
    
    rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
    rhs = f(S.r(1,:)./rr,S.r(2,:)./rr,S.r(3,:)./rr);
    
    ut = rhs./(-l*(l+1)+zkvals);
    % ut = rhs;
else
    if type ==1
    fct = @(x) -1./vecnorm(x-x_source,2);
    fctgrad = @(x) (x-x_source).*(vecnorm(x-x_source,2).^(-3));
    fcthess = @(x) newton_hess(x,x_source);
    elseif type == 2
    fct = @(x) exp(-vecnorm(x-x_source,2).^2/4);
    fctgrad = @(x) -2*(x-x_source).*exp(-vecnorm(x-x_source,2).^2/4)/4;
    fcthess = @(x) gauss_hess(x,x_source);
    else % spherical harmonic
    
    end
    
    x = S.r;
    n = S.n;
    H = -S.mean_curv;
    ut = charge*fct(x);
    fgrad = fctgrad(x);
    fhess = fcthess(x);
    
    for i = 1:npts
    rhs(i) = charge*(trace(fhess(:,:,i))+2*H(i)*n(:,i)'*fgrad(:,i)-n(:,i)'*fhess(:,:,i)*n(:,i));
    end
    
    if any(abs(zkvals)>1e-14)
      rhs = rhs + zkvals.^2.*charge.*fct(x);
    else
        ut = ut - sum(ut(:).*S.wts)/sum(S.wts);
    end
    ut = ut(:);
    rhs = rhs(:);
    
end
end

function H=newton_hess(x,x_source)
  n = size(x,2);
  H = zeros(3,3,n);
  for i = 1:n
    y = x(:,i)-x_source;
    nrm = vecnorm(y);
    H(:,:,i) = eye(3)/nrm^3 - 3*(y*y')/(nrm)^5;
  end
end

function H=gauss_hess(x,x_source)
  n = size(x,2);
  H = zeros(3,3,n);
  for i = 1:n
    y = x(:,i)-x_source;
    % nrm = vecnorm(y);
    H(:,:,i) = (-2*eye(3)/4 + 4*(y*y')/16).*exp(-vecnorm(y,2).^2/4);
  end
end

