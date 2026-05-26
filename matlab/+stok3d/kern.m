
function submat= kern(srcinfo,targinfo,type,varargin)
%STOKES.KERN standard Stokes layer potential kernels in 3D
% 
% Syntax: submat = stokes.kern(zk,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points.
%  
% Kernels based on G_{ij}(x,y) = (r_i r_j)/(2r^3) + delta_{ij}/(2r)
%               T_{ijk}(x,y) = -3 r_i r_j r_k/ r^5 = \partiak_k G_{ij}(x,y)
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
%
% Input:
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
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
%                type == 'd', Stresslet/double layer kernel D
%                type == 's', Stokeslet/single layer kernel S
%                type == 'c', combined layer kernel alph*S+beta*D
%   varargin{1} - alpha, beta in the combined layer formula, coef in eval 
%                and evalg, otherwise does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% see also STOK3D.GREEN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'d')
  srcnorm = srcinfo.n;
  [~,grad] = stok3d.green(src,targ);
  nx = reshape(repmat(srcnorm(1,:),nt,1),[1,1,nt,1,ns]);
  ny = reshape(repmat(srcnorm(2,:),nt,1),[1,1,nt,1,ns]);
  nz = reshape(repmat(srcnorm(3,:),nt,1),[1,1,nt,1,ns]);
  submat = -(grad(1,:,:,:,:).*nx + grad(2,:,:,:,:).*ny+grad(3,:,:,:,:).*nz);
  submat = -reshape(submat,[3,nt,3,ns]);
end

if strcmpi(type,'sprime')
  % S'_{ij}(x,y) = T_{ijk}(x,y) n_k(x)
  %              = -3(r . n_x) r_i r_j / r^5 / (4pi),  r = x - y
  %
  % Identical to the DLP formula but contracting with the TARGET normal
  % n_x instead of the source normal n_y.  Uses the same double-negative
  % sign convention as the 'd' case above so that the resulting matrix
  % is sign-consistent with st3d_strac and stok_sprime_eval_addsub.
  targnorm = targinfo.n;
  [~,grad] = stok3d.green(src, targ);
  % grad(k, i, it, j, is): stok3d.green(src,targ) with r = targ - src
  nx = reshape(repmat(targnorm(1,:).', 1, ns), [1, 1, nt, 1, ns]);
  ny = reshape(repmat(targnorm(2,:).', 1, ns), [1, 1, nt, 1, ns]);
  nz = reshape(repmat(targnorm(3,:).', 1, ns), [1, 1, nt, 1, ns]);
  submat = -(grad(1,:,:,:,:).*nx + grad(2,:,:,:,:).*ny + grad(3,:,:,:,:).*nz);
  submat = -reshape(submat, [3, nt, 3, ns]);
end

if strcmpi(type,'sdu')
  disp("unsupported kernel");
  submat = 0;
end

if strcmpi(type,'sdv')
  disp("unsupported kernel");
  submat = 0;
end

if strcmpi(type,'s')
  submat = stok3d.green(src,targ);
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
  [s,grad] = stok3d.green(src,targ);
  nx = reshape(repmat(srcnorm(1,:),nt,1),[1,1,nt,1,ns]);
  ny = reshape(repmat(srcnorm(2,:),nt,1),[1,1,nt,1,ns]);
  nz = reshape(repmat(srcnorm(3,:),nt,1),[1,1,nt,1,ns]);
  submat = -(grad(1,:,:,:,:).*nx + grad(2,:,:,:,:).*ny+grad(3,:,:,:,:).*nz);
  submat = -beta*reshape(submat,[3,nt,3,ns]) + alpha*s;
end

if strcmpi(type,'all')
  disp("unsupported kernel");
  submat = 0;
end

if strcmpi(type,'eval')
  disp("unsupported kernel");
  submat = 0;
end


if strcmpi(type,'evalg')
  disp("unsupported kernel");
  submat = 0;
end



if strcmpi(type,'trans1')
  disp("unsupported kernel");
  submat = 0;
end

