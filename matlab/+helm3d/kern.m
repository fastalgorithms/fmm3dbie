
function submat= kern(zk,srcinfo,targinfo,type,varargin)
%HELM3D.KERN standard Helmholtz layer potential kernels in 3D
% 
% Syntax: submat = helm3d.kern(zk,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points.
%  
% Kernels based on G(x,y) = exp(i*k*|x-y|)/(4*pi*|x-y|)
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
%
% Input:
%   zk - complex number, Helmholtz wave number
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (3,:) array
%                ptinfo.du - first derivative with respect to u in 
%                     underlying parameterization (3,:)
%                ptinfo.dv - first derivative with respect to u in 
%                     underlying parameterization (3,:)
%                ptinfo.n - normals (3,:)
%   targinfo - description of targets in ptinfo struct format,
%                if info not relevant (du,dv,n) it doesn't need to
%                be provided. sprime requires normal info in
%                targinfo.n.
%   type - string, determines kernel type
%                type == 'd', double layer kernel D
%                type == 's', single layer kernel S
%                type == 'sdu', derivative of single layer with respect to
%                        u parameter
%                type == 'sdv', derivative of single layer with respect to
%                        v parameter
%                type == 'sprime', normal derivative of single
%                      layer S'
%                type == 'c', combined layer kernel alph*S+beta*D
%                type == 'eval' computes matrices for both S and D
%                type == 'evalg' computes matrices for both S and D, and 
%                      their gradients
%   varargin{1} - alpha, beta in the combined layer formula, coef in eval 
%                and evalg, otherwise does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% Examples
%   s.r = [0;0;0]; s.n = [1;0;0];            % source struct
%   t.r = [1;0;0];                           % target struct
%   k = 10; S = exp(1i*k)/(4*pi); D = -(1i*k-1)*exp(1i*k)/(4*pi);
%   helm3d.kern(k,s,t,'s') - S               % gives 0
%   helm3d.kern(k,s,t,'d') - D               % gives 0
%   eta = -1i*k; cfie = eta*S + 1.0*D;
%   helm3d.kern(k,s,t,'c',eta,1.0) - cfie    % gives 0
%
% See also HELM3D.GREEN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'d')
  srcnorm = srcinfo.n;
  [~,grad] = helm3d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  nz = repmat(srcnorm(3,:),nt,1);
  submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny+grad(:,:,3).*nz);
end

if strcmpi(type,'sprime')
  targnorm = targinfo.n;
  [~,grad] = helm3d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);
  nz = repmat((targnorm(3,:)).',1,ns);
  submat = (grad(:,:,1).*nx + grad(:,:,2).*ny+grad(:,:,3).*nz);
end

if strcmpi(type,'sdu')
  du = targinfo.du;
  [~,grad] = helm3d.green(zk,src,targ);
  dx = repmat((du(1,:)).',1,ns);
  dy = repmat((du(2,:)).',1,ns);
  dz = repmat((du(3,:)).',1,ns);
  dn = sqrt(dx.*dx+dy.*dy+dz.*dz);
  submat = (grad(:,:,1).*dx./dn + grad(:,:,2).*dy+...
      grad(:,:,3).*dz)./dn;
end

if strcmpi(type,'sdv')
  dv = targinfo.dv;
  [~,grad] = helm3d.green(zk,src,targ);
  dx = repmat((dv(1,:)).',1,ns);
  dy = repmat((dv(2,:)).',1,ns);
  dz = repmat((dv(3,:)).',1,ns);
  dn = sqrt(dx.*dx+dy.*dy+dz.*dz);
  submat = (grad(:,:,1).*dx./dn + grad(:,:,2).*dy+...
      grad(:,:,3).*dz)./dn;
end

if strcmpi(type,'s')
  submat = helm3d.green(zk,src,targ);
end

if strcmpi(type,'dprime')
  disp("unsupported kernel");
  submat = 0;
end

if strcmpi(type,'c')
%%%%%%
%       .  .  .  alpha*S + beta*D
%%%%%%
  srcnorm = srcinfo.n;
  alpha = varargin{1};
  beta  = varargin{2};
  [s,grad] = helm3d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  nz = repmat(srcnorm(3,:),nt,1);
  submat = -beta*(grad(:,:,1).*nx + grad(:,:,2).*ny+grad(:,:,3).*nz)...
      +alpha*s;
end

if strcmpi(type,'all')
  disp("unsupported kernel");
  submat = 0;
end

if strcmpi(type,'eval')
  coef = varargin{1};
  
  srcnorm = srcinfo.n;
  
  submat = zeros(nt,2*ns);
  % S
  [submats,grad] = helm3d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nzsrc = repmat(srcnorm(3,:),nt,1);
  
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc+grad(:,:,3).*nzsrc);
    
  submat(:,1:2:2*ns) = coef*submatd;
  submat(:,2:2:2*ns) = submats;
end


if strcmpi(type,'evalg')
  coef = varargin{1};
  
  srcnorm = srcinfo.n;
  
  submat = zeros(nt,ns,8);
  % S
  [submats,grad,hess] = helm3d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nzsrc = repmat(srcnorm(2,:),nt,1);
  % D
  submatd  = -(grad(:,:,1).*nxsrc+grad(:,:,2).*nysrc+grad(:,:,3).*nzsrc);
    
  submat(:,:,1) = coef*submatd;
  submat(:,:,2) = submats;
  submat(:,:,3) = -coef*(hess(:,:,1,1).*nxsrc + hess(:,:,1,2).*nysrc+...
      hess(:,:,1,3).*nzsrc);
  submat(:,:,4) = grad(:,:,1);
  submat(:,:,5) = -coef*(hess(:,:,2,1).*nxsrc + hess(:,:,2,2).*nysrc+...
      hess(:,:,2,3).*nzsrc);
  submat(:,:,6) = grad(:,:,2);
    submat(:,:,7) = -coef*(hess(:,:,3,1).*nxsrc + hess(:,:,3,2).*nysrc+...
      hess(:,:,3,3).*nzsrc);
  submat(:,:,8) = grad(:,:,3);
end



if strcmpi(type,'trans1')
  disp("unsupported kernel");
  submat = 0;
end

