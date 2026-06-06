
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
  % n_x instead of the source normal n_y.
  % grad(k, i, it, j, is) = (3/4pi) * r_k * r_i * r_j / r^5,  r = targ - src
  % so  T_{ijk} n_k^x = sum_k grad(k,...) n_k  = (3/4pi)*(r.n_x)*r_i*r_j/r^5
  % The physical kernel has a leading minus, so we negate once.
  targnorm = targinfo.n;
  [~,grad] = stok3d.green(src, targ);
  nx = reshape(repmat(targnorm(1,:).', 1, ns), [1, 1, nt, 1, ns]);
  ny = reshape(repmat(targnorm(2,:).', 1, ns), [1, 1, nt, 1, ns]);
  nz = reshape(repmat(targnorm(3,:).', 1, ns), [1, 1, nt, 1, ns]);
  submat = grad(1,:,:,:,:).*nx + grad(2,:,:,:,:).*ny + grad(3,:,:,:,:).*nz;
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

if strcmpi(type,'cprime')
  %%%%%%%
  %   cprime = traction of alpha*S + beta*D
  %
  %  Let r = x - y (targ - src), m = n^x (target normal), n = n^y (source normal).
  %
  %  The D kernel velocity is  u_i^D = T_{ijl}(x,y) n_l^y sigma_j
  %                                   = (3/4pi) r_i r_j r_l / r^5 * n_l sigma_j.
  %
  %  Traction of D (derived analytically):
  %    t_i^D = m_i (sigma.n)/(2pi r^3)
  %          + (3/4pi r^5) [
  %              r_i (m.sigma)(r.n)
  %            + r_i (r.sigma)(m.n)
  %            + (m.r) sigma_i (r.n)
  %            + (m.r)(r.sigma) n_i
  %            - 10 (m.r) r_i (r.sigma)(r.n) / r^2
  %            ]
  %
  %  Traction of S = sprime (see sprime block above).
  %  cprime = alpha*sprime + beta*dprime.
  %%%%%%
  srcnorm  = srcinfo.n;    % n^y, (3, ns)
  targnorm = targinfo.n;   % m = n^x, (3, nt)
  alpha = varargin{1};
  beta  = varargin{2};

  % ---- sprime part ----
  [~, grad] = stok3d.green(src, targ);
  % grad(k,i,it,j,is) = (3/4pi) r_k r_i r_j / r^5
  mx = reshape(repmat(targnorm(1,:).', 1, ns), [1, 1, nt, 1, ns]);
  my = reshape(repmat(targnorm(2,:).', 1, ns), [1, 1, nt, 1, ns]);
  mz = reshape(repmat(targnorm(3,:).', 1, ns), [1, 1, nt, 1, ns]);
  sprime_mat = grad(1,:,:,:,:).*mx + grad(2,:,:,:,:).*my + grad(3,:,:,:,:).*mz;
  sprime_mat = -reshape(sprime_mat, [3, nt, 3, ns]);

  % ---- dprime part ----
  % Build r, scalar dot products (all shape (nt,ns))
  xs_ = repmat(src(1,:),  nt, 1);
  ys_ = repmat(src(2,:),  nt, 1);
  zs_ = repmat(src(3,:),  nt, 1);
  xt_ = repmat(targ(1,:).', 1, ns);
  yt_ = repmat(targ(2,:).', 1, ns);
  zt_ = repmat(targ(3,:).', 1, ns);
  rxd = xt_-xs_;  ryd = yt_-ys_;  rzd = zt_-zs_;
  rd2 = rxd.^2 + ryd.^2 + rzd.^2;
  rd  = sqrt(rd2);
  ir5 = 1./rd.^5;

  % Source normals and target normals broadcast to (nt,ns)
  nx_ = repmat(srcnorm(1,:),  nt, 1);
  ny_ = repmat(srcnorm(2,:),  nt, 1);
  nz_ = repmat(srcnorm(3,:),  nt, 1);
  mx_ = repmat(targnorm(1,:).', 1, ns);
  my_ = repmat(targnorm(2,:).', 1, ns);
  mz_ = repmat(targnorm(3,:).', 1, ns);

  rdotn  = rxd.*nx_ + ryd.*ny_ + rzd.*nz_;   % r . n^y
  rdotm  = rxd.*mx_ + ryd.*my_ + rzd.*mz_;   % r . m = r . n^x

  % dprime_mat(i, it, j, is) for each (i,j) pair
  % t_i^D(sigma_j) = m_i * delta_{ij} * (sigma.n) / (2pi r^3)   <- pressure term (summed separately)
  %   + (3/4pi r^5)[r_i*(m.e_j)*(r.n) + r_i*(r.e_j)*(m.n)
  %                 + (m.r)*delta_{ij}*(r.n) + (m.r)*(r.e_j)*n_i
  %                 - 10*(m.r)*r_i*(r.e_j)*(r.n)/r^2]
  %
  % where sigma = e_j (unit vector) since we build the kernel matrix column by column.
  % Replacing sigma = e_j:  sigma.n = n_j, m.sigma = m_j, r.sigma = r_j
  %
  % t_i^D|_{sigma=e_j} = m_i n_j / (2pi r^3)
  %   + (3/4pi r^5)[r_i m_j (r.n) + r_i r_j (m.n) + (m.r) delta_{ij} (r.n)
  %                 + (m.r) r_j n_i - 10 (m.r) r_i r_j (r.n)/r^2]

  pfact = 1/2/pi;
  tfact = 3/4/pi;

  rdotm_ir5  = rdotm .* ir5;                        % (m.r)/r^5
  rdotn_ir5  = rdotn .* ir5;                        % (r.n)/r^5
  mdotn      = mx_.*nx_ + my_.*ny_ + mz_.*nz_;      % (m.n)
  ten_mr_rn  = 10*rdotm.*rdotn./rd2.*ir5;           % 10(m.r)(r.n)/r^7

  % r components as (1,nt,1,ns) for broadcasting against i,j indices
  rc = {rxd, ryd, rzd};
  mc = {mx_, my_, mz_};
  nc = {nx_, ny_, nz_};

  dprime_mat = zeros(3, nt, 3, ns);
  for ii = 1:3
    ri = reshape(rc{ii}, [1,nt,1,ns]);
    mi = reshape(mc{ii}, [1,nt,1,ns]);
    ni = reshape(nc{ii}, [1,nt,1,ns]);
    for jj = 1:3
      rj = reshape(rc{jj}, [1,nt,1,ns]);
      mj = reshape(mc{jj}, [1,nt,1,ns]);
      nj = reshape(nc{jj}, [1,nt,1,ns]);

      % Pressure term: m_i n_j / (2pi r^3)
      val = pfact * mi .* nj ./ reshape(rd.^3,[1,nt,1,ns]);

      % Symmetrized velocity gradient terms (3/4pi)/r^5 * [...]
      rdotn5 = reshape(rdotn_ir5,[1,nt,1,ns]);
      rdotm5 = reshape(rdotm_ir5,[1,nt,1,ns]);
      mn     = reshape(mdotn,    [1,nt,1,ns]);
      ten    = reshape(ten_mr_rn,[1,nt,1,ns]);

      ir5_  = reshape(ir5,  [1,nt,1,ns]);
      val = val + tfact*(ri.*mj.*rdotn5 ...           % r_i m_j (r.n)/r^5
                       + ri.*rj.*mn.*ir5_ ...         % r_i r_j (m.n)/r^5
                       + (ii==jj)*reshape(rdotm.*rdotn.*ir5,[1,nt,1,ns]) ... % (m.r)(r.n)/r^5 * delta_{ij}
                       + rdotm5.*rj.*ni ...            % (m.r)/r^5 * r_j n_i
                       - ri.*rj.*ten);                 % -10(m.r) r_i r_j (r.n)/r^7

      dprime_mat(ii,:,jj,:) = val;
    end
  end

  submat = alpha*sprime_mat + beta*dprime_mat;
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

