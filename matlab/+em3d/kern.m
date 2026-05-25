
function submat= kern(zk,srcinfo,targinfo,type,varargin)
%EM3D.KERN standard vector Helmholtz layer potential kernels 
%in 3D required for Maxwell potentials
% 
% Syntax: submat = em3d.kern(zk,srcinfo,targinfo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points.
%  
% Kernels based on S:= G(x,y) = exp(i*k*|x-y|)/(4*pi*|x-y|)
%
%
% NxDelS       = n \times \nabla S [scalar]
% NxCurlS      = n \times \nabla \times S [vector]
% NdotCurlS    = n \cdot \nabla \times S [vector]
% NdotS        = n \cdot S [vector]
% NxS          = n \times S [vector]
% NxCurlCurlS  = n \times \nabla \times \nabla \times S [vector]
% NdotCurCurlS = n \cdot \nabla \times \nabla \times S [vector]
% DeldotS      = \nabla \cdot S [vector]
%
% Input:
%   zk - complex number, Maxwell wave number
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
%                type == 'sdu', derivative of single layer with respect to
%                        u parameter
%                type == 'sdv', derivative of single layer with respect to
%                        v parameter
%                type == 'NxDelS'
%                type == 'NxCurlS'
%                type == 'NdotCurlS'
%                type == 'NdotS'
%                type == 'NxS'
%                type == 'NxCurlCurlS'
%                type == 'NdotCurlCurlS'
%                type == 'DeldotS'
%   varargin{1} - alpha, beta in the combined layer formula, coef in eval 
%                and evalg, otherwise does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

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
  [~,grad] = em3d.green(zk,src,targ);
  dx = repmat((dv(1,:)).',1,ns);
  dy = repmat((dv(2,:)).',1,ns);
  dz = repmat((dv(3,:)).',1,ns);
  dn = sqrt(dx.*dx+dy.*dy+dz.*dz);
  submat = (grad(:,:,1).*dx./dn + grad(:,:,2).*dy+...
      grad(:,:,3).*dz)./dn;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(type,'NxDelS')
  targnorm = targinfo.n;
  [~,grad] = helm3d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);
  nz = repmat((targnorm(3,:)).',1,ns);
  
  drut = targinfo.dru;
  dxut = repmat((drut(1,:)).',1,ns);
  dyut = repmat((drut(2,:)).',1,ns);
  dzut = repmat((drut(3,:)).',1,ns);
  
  drvt = targinfo.drv;
  dxvt = repmat((drvt(1,:)).',1,ns);
  dyvt = repmat((drvt(2,:)).',1,ns);
  dzvt = repmat((drvt(3,:)).',1,ns);
  
  tmatx =                 -nz.*grad(:,:,2)+ny.*grad(:,:,3);
  tmaty =  nz.*grad(:,:,1)                -nx.*grad(:,:,3);
  tmatz = -ny.*grad(:,:,1)+nx.*grad(:,:,2);
  
  submat = zeros(2,nt,ns);
  submat(1,:,:) = dxut.*tmatx+dyut.*tmaty+dzut.*tmatz;
  submat(2,:,:) = dxvt.*tmatx+dyvt.*tmaty+dzvt.*tmatz;
  submat = reshape(submat,2*nt,ns);
end

if strcmpi(type,'NxCurlS')
  targnorm = targinfo.n;
  [~,grad] = helm3d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);
  nz = repmat((targnorm(3,:)).',1,ns);
  
  drut = targinfo.dru;
  dxut = repmat((drut(1,:)).',1,ns);
  dyut = repmat((drut(2,:)).',1,ns);
  dzut = repmat((drut(3,:)).',1,ns);
  
  drvt = targinfo.drv;
  dxvt = repmat((drvt(1,:)).',1,ns);
  dyvt = repmat((drvt(2,:)).',1,ns);
  dzvt = repmat((drvt(3,:)).',1,ns);
  
  drus = srcinfo.dru;
  dxus = repmat((drus(1,:)),nt,1);
  dyus = repmat((drus(2,:)),nt,1);
  dzus = repmat((drus(3,:)),nt,1);
  
  drvs = srcinfo.drv;
  dxvs = repmat((drvs(1,:)),nt,1);
  dyvs = repmat((drvs(2,:)),nt,1);
  dzvs = repmat((drvs(3,:)),nt,1);
  
  gmat11 = -grad(:,:,3).*nz-grad(:,:,2).*ny;
  gmat12 =  grad(:,:,1).*ny;
  gmat13 =  grad(:,:,1).*nz;
  gmat21 =  grad(:,:,2).*nx;
  gmat22 = -grad(:,:,1).*nx-grad(:,:,3).*nz;
  gmat23 =  grad(:,:,2).*nz;
  gmat31 =  grad(:,:,3).*nx;
  gmat32 =  grad(:,:,3).*ny;
  gmat33 = -grad(:,:,1).*nx-grad(:,:,2).*ny;
  
  tmat11 = gmat11.*dxus + gmat12.*dyus + ...
      gmat13.*dzus;
  tmat21 = gmat21.*dxus + gmat22.*dyus + ...
      gmat23.*dzus;
  tmat31 = gmat31.*dxus + gmat32.*dyus + ...
      gmat33.*dzus;
  
  tmat12 = gmat11.*dxvs + gmat12.*dyvs + ...
      gmat13.*dzvs;
  tmat22 = gmat21.*dxvs + gmat22.*dyvs + ...
      gmat23.*dzvs;
  tmat32 = gmat31.*dxvs + gmat32.*dyvs + ...
      gmat33.*dzvs;
  
  submat = zeros(2,nt,2,ns);
  submat(1,:,1,:) = tmat11.*dxut + tmat21.*dyut + ...
      tmat31.*dzut;
  submat(1,:,2,:) = tmat12.*dxut + tmat22.*dyut + ...
      tmat32.*dzut;
  
  submat(2,:,1,:) = tmat11.*dxvt + tmat21.*dyvt + ...
      tmat31.*dzvt;
  submat(2,:,2,:) = tmat12.*dxvt + tmat22.*dyvt + ...
      tmat32.*dzvt;
  submat = reshape(submat,2*nt,2*ns);
  
end

if strcmpi(type,'NdotS')
  targnorm = targinfo.n;
  [gfun] = helm3d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);
  nz = repmat((targnorm(3,:)).',1,ns);
  
%%%% extract and copy local dr vectors

  drus = srcinfo.dru;
  dxus = repmat((drus(1,:)),nt,1);
  dyus = repmat((drus(2,:)),nt,1);
  dzus = repmat((drus(3,:)),nt,1);
  
  drvs = srcinfo.drv;
  dxvs = repmat((drvs(1,:)),nt,1);
  dyvs = repmat((drvs(2,:)),nt,1);
  dzvs = repmat((drvs(3,:)),nt,1);
  
  submat = zeros(ns,2,nt);
  submat(:,1,:) = gfun.*(dxus.*nx+dyus.*ny+dzus.*nz);
  submat(:,2,:) = gfun.*(dxvs.*nx+dyvs.*ny+dzvs.*nz);
  submat = reshape(submat,[ns,2*nt]);
  
end


if strcmpi(type,'NxS')
    
  targnorm = targinfo.n;
  [gfun] = helm3d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);
  nz = repmat((targnorm(3,:)).',1,ns);
  
  drut = targinfo.dru;
  dxut = repmat((drut(1,:)).',1,ns);
  dyut = repmat((drut(2,:)).',1,ns);
  dzut = repmat((drut(3,:)).',1,ns);
  
  drvt = targinfo.drv;
  dxvt = repmat((drvt(1,:)).',1,ns);
  dyvt = repmat((drvt(2,:)).',1,ns);
  dzvt = repmat((drvt(3,:)).',1,ns);
  
  drus = srcinfo.dru;
  dxus = repmat((drus(1,:)),nt,1);
  dyus = repmat((drus(2,:)),nt,1);
  dzus = repmat((drus(3,:)),nt,1);
  
  drvs = srcinfo.drv;
  dxvs = repmat((drvs(1,:)),nt,1);
  dyvs = repmat((drvs(2,:)),nt,1);
  dzvs = repmat((drvs(3,:)),nt,1);
  
  submat = zeros(2,nt,2,ns);
  submat(1,:,1,:) = dxut.*(ny.*dzus-nz.*dyus)+dyut.*(nz.*dxus-nx.*dzus)+...
                    dzut.*(nx.*dyus-ny.*dxus);
  submat(2,:,1,:) = dxvt.*(ny.*dzus-nz.*dyus)+dyvt.*(nz.*dxus-nx.*dzus)+...
                    dzvt.*(nx.*dyus-ny.*dxus);
  submat(1,:,2,:) = dxut.*(ny.*dzvs-nz.*dyvs)+dyut.*(nz.*dxvs-nx.*dzvs)+...
                    dzut.*(nx.*dyvs-ny.*dxvs);
  submat(2,:,2,:) = dxvt.*(ny.*dzvs-nz.*dyvs)+dyvt.*(nz.*dxvs-nx.*dzvs)+...
                    dzvt.*(nx.*dyvs-ny.*dxvs);                    
  submat = reshape(submat,[2*nt,2*ns]);
  
end

if strcmpi(type,'NxCurlCurlS')
    
  targnorm = targinfo.n;
  [gfun,~,hess] = helm3d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);
  nz = repmat((targnorm(3,:)).',1,ns);
  
  hess(:,:,1,1) = hess(:,:,1,1) + zk^2*gfun;
  hess(:,:,2,2) = hess(:,:,2,2) + zk^2*gfun;
  hess(:,:,3,3) = hess(:,:,3,3) + zk^2*gfun;
  
  drut = targinfo.dru;
  dxut = repmat((drut(1,:)).',1,ns);
  dyut = repmat((drut(2,:)).',1,ns);
  dzut = repmat((drut(3,:)).',1,ns);
  
  drvt = targinfo.drv;
  dxvt = repmat((drvt(1,:)).',1,ns);
  dyvt = repmat((drvt(2,:)).',1,ns);
  dzvt = repmat((drvt(3,:)).',1,ns);
  
  drus = srcinfo.dru;
  dxus = repmat((drus(1,:)),nt,1);
  dyus = repmat((drus(2,:)),nt,1);
  dzus = repmat((drus(3,:)),nt,1);
  
  drvs = srcinfo.drv;
  dxvs = repmat((drvs(1,:)),nt,1);
  dyvs = repmat((drvs(2,:)),nt,1);
  dzvs = repmat((drvs(3,:)),nt,1);
  
  svec = zeros(3,2,nt,ns);
  svec(1,1,:,:) = dxus;
  svec(2,1,:,:) = dyus;
  svec(3,1,:,:) = dzus;
  svec(1,2,:,:) = dxvs;
  svec(2,2,:,:) = dyvs;
  svec(3,2,:,:) = dzvs;
  
  tvec = zeros(3,2,nt,ns);
  tvec(1,1,:,:) = dyut.*nz-dzut.*ny;
  tvec(2,1,:,:) = dzut.*nx-dxut.*nz;
  tvec(3,1,:,:) = dxut.*ny-dyut.*nx;
  tvec(1,2,:,:) = dyvt.*nz-dzvt.*ny;
  tvec(2,2,:,:) = dzvt.*nx-dxvt.*nz;
  tvec(3,2,:,:) = dxvt.*ny-dyvt.*nx;
  
  submat = zeros(2,nt,2,ns);
  for ii=1:2
      for jj=1:2
          for kk=1:3
              for iii=1:3
     submat(ii,:,jj,:) = squeeze(submat(ii,:,jj,:))+...
         squeeze(svec(kk,jj,:,:)).*squeeze(tvec(iii,ii,:,:)).*squeeze(hess(:,:,iii,kk));
              end
          end
      end
  end
  
  submat = reshape(submat,[2*nt,2*ns]);
  
end

if strcmpi(type,'NdotCurlCurlS')
    
  targnorm = targinfo.n;
  [gfun,~,hess] = helm3d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);
  nz = repmat((targnorm(3,:)).',1,ns);
  
  hess(:,:,1,1) = hess(:,:,1,1) + zk^2*gfun;
  hess(:,:,2,2) = hess(:,:,2,2) + zk^2*gfun;
  hess(:,:,3,3) = hess(:,:,3,3) + zk^2*gfun;
  
  drut = targinfo.dru;
  dxut = repmat((drut(1,:)).',1,ns);
  dyut = repmat((drut(2,:)).',1,ns);
  dzut = repmat((drut(3,:)).',1,ns);
  
  drvt = targinfo.drv;
  dxvt = repmat((drvt(1,:)).',1,ns);
  dyvt = repmat((drvt(2,:)).',1,ns);
  dzvt = repmat((drvt(3,:)).',1,ns);
  
  drus = srcinfo.dru;
  dxus = repmat((drus(1,:)),nt,1);
  dyus = repmat((drus(2,:)),nt,1);
  dzus = repmat((drus(3,:)),nt,1);
  
  drvs = srcinfo.drv;
  dxvs = repmat((drvs(1,:)),nt,1);
  dyvs = repmat((drvs(2,:)),nt,1);
  dzvs = repmat((drvs(3,:)),nt,1);
  
  svec = zeros(3,2,nt,ns);
  svec(1,1,:,:) = dxus;
  svec(2,1,:,:) = dyus;
  svec(3,1,:,:) = dzus;
  svec(1,2,:,:) = dxvs;
  svec(2,2,:,:) = dyvs;
  svec(3,2,:,:) = dzvs;
  
  tvec = zeros(3,1,nt,ns);
  tvec(1,1,:,:) = nx;
  tvec(2,1,:,:) = ny;
  tvec(3,1,:,:) = nz;

  
  submat = zeros(1,nt,2,ns);
  for ii=1:1
      for jj=1:2
          for kk=1:3
              for iii=1:3
     submat(ii,:,jj,:) = squeeze(submat(ii,:,jj,:))+...
         squeeze(svec(kk,jj,:,:)).*squeeze(tvec(iii,ii,:,:)).*squeeze(hess(:,:,iii,kk));
              end
          end
      end
  end
  
  submat = reshape(submat,[nt,2*ns]);
  
end


if strcmpi(type,'DeldotS')
  [~,grad] = helm3d.green(zk,src,targ);
  
  drus = srcinfo.dru;
  dxus = repmat((drus(1,:)),nt,1);
  dyus = repmat((drus(2,:)),nt,1);
  dzus = repmat((drus(3,:)),nt,1);
  
  drvs = srcinfo.drv;
  dxvs = repmat((drvs(1,:)),nt,1);
  dyvs = repmat((drvs(2,:)),nt,1);
  dzvs = repmat((drvs(3,:)),nt,1);
  submat = zeros(nt,2,ns);
  submat(:,1,:) = dxus.*grad(:,:,1)+dyus.*grad(:,:,2)+dzus.*grad(:,:,3);
  submat(:,2,:) = dxvs.*grad(:,:,1)+dyvs.*grad(:,:,2)+dzvs.*grad(:,:,3);
  submat = reshape(submat,nt,2*ns);
  
end

if strcmpi(type,'nrccie-bc')
%
%  NRCCIE boundary-condition (system) kernel.
%
%  Density: [j_u, j_v, rho]   (3 components / source point, surface-tangent basis)
%  Output:  [pot_u, pot_v, pot_rho]  (3 equations / target point)
%  Matrix size: (3*nt) x (3*ns), component-fast row/col ordering.
%
%  Follows get_nrccie_inteq_comps_from_potgrad in em_nrccie_pec.f90.
%  For source density component d_s (= du_s or dv_s) and scalar rho:
%
%  zcurl   = grad_t G x d_s              (curl of vector single layer)
%  zvec    = n_t x zcurl                 (= -M_k[J] term, magnetic BC)
%  E       = ik*G*d_s - grad_t G * rho   (E from representation)
%  nxnxE   = (n_t . E)*n_t - E           (= n x n x E)
%  zvec3   = alpha*nxnxE - zvec
%  pot_u/v = zvec3 . du_t/dv_t
%  pot_rho = (grad_t S_k[rho] - ik*S_k[J]).n_t + alpha*(div_t S_k[J] - ik*S_k[rho])
%
%  varargin{1} = alpha  (CFIE regularization parameter)

  alpha = varargin{1};

  targnorm = targinfo.n;
  drut = targinfo.du;
  drvt = targinfo.dv;
  drus = srcinfo.du;
  drvs = srcinfo.dv;
  srcnorm = srcinfo.n;

  [gfun, grad] = helm3d.green(zk, src, targ);
  % gfun : (nt, ns)
  % grad : (nt, ns, 3)  — gradient w.r.t. target x

  gx = grad(:,:,1);  gy = grad(:,:,2);  gz = grad(:,:,3);

  % Target geometry
  nx  = repmat(targnorm(1,:).', 1, ns);
  ny  = repmat(targnorm(2,:).', 1, ns);
  nz  = repmat(targnorm(3,:).', 1, ns);
  dxut = repmat(drut(1,:).', 1, ns);  dyut = repmat(drut(2,:).', 1, ns);  dzut = repmat(drut(3,:).', 1, ns);
  dxvt = repmat(drvt(1,:).', 1, ns);  dyvt = repmat(drvt(2,:).', 1, ns);  dzvt = repmat(drvt(3,:).', 1, ns);

  % Source geometry
  dxus = repmat(drus(1,:), nt, 1);  dyus = repmat(drus(2,:), nt, 1);  dzus = repmat(drus(3,:), nt, 1);
  dxvs = repmat(drvs(1,:), nt, 1);  dyvs = repmat(drvs(2,:), nt, 1);  dzvs = repmat(drvs(3,:), nt, 1);
  nxs  = repmat(srcnorm(1,:), nt, 1);
  nys  = repmat(srcnorm(2,:), nt, 1);
  nzs  = repmat(srcnorm(3,:), nt, 1);

  % Helper: given source tangent vector ds = (dxs, dys, dzs), compute the
  % (nt, ns) matrices for each block row/col of the kernel matrix.
  %
  % zcurl = grad_t G x ds  (curl of G*ds w.r.t. target)
  %   zcurl_x = gy*dzs - gz*dys
  %   zcurl_y = gz*dxs - gx*dzs
  %   zcurl_z = gx*dys - gy*dxs
  %
  % zvec = n_t x zcurl
  %   zvec_x = ny*zcurl_z - nz*zcurl_y
  %   zvec_y = nz*zcurl_x - nx*zcurl_z
  %   zvec_z = nx*zcurl_y - ny*zcurl_x
  %
  % E = ik*G*ds  (J column; rho column adds -grad_t G term handled separately)
  %   Ex = ik*G*dxs,  Ey = ik*G*dys,  Ez = ik*G*dzs
  %
  % n_t.E = ik*G*(nx*dxs + ny*dys + nz*dzs)
  % nxnxE = (n_t.E)*n_t - E
  %   nxnxE_x = (n_t.E)*nx - ik*G*dxs
  %   etc.
  %
  % zvec3 = alpha*nxnxE - zvec
  % pot_u = zvec3 . du_t,  pot_v = zvec3 . dv_t

  % --- Column j_u ---
  curl_x_u = gy.*dzus - gz.*dyus;
  curl_y_u = gz.*dxus - gx.*dzus;
  curl_z_u = gx.*dyus - gy.*dxus;

  zvec_x_u = ny.*curl_z_u - nz.*curl_y_u;
  zvec_y_u = nz.*curl_x_u - nx.*curl_z_u;
  zvec_z_u = nx.*curl_y_u - ny.*curl_x_u;

  ndotE_u  = 1i*zk*gfun .* (nx.*dxus + ny.*dyus + nz.*dzus);
  v3x_u = alpha*(ndotE_u.*nx - 1i*zk*gfun.*dxus) - zvec_x_u;
  v3y_u = alpha*(ndotE_u.*ny - 1i*zk*gfun.*dyus) - zvec_y_u;
  v3z_u = alpha*(ndotE_u.*nz - 1i*zk*gfun.*dzus) - zvec_z_u;

  K_uu = dxut.*v3x_u + dyut.*v3y_u + dzut.*v3z_u;  % ru component, j_u col
  K_vu = dxvt.*v3x_u + dyvt.*v3y_u + dzvt.*v3z_u;  % rv component, j_u col

  % pot_rho (row 3), j_u column:
  %   (grad_t(0) - ik*G*J).n_t + alpha*(div_t(G*J) - ik*0)
  %   = -ik*G*(d_s.n_t) + alpha*(grad_t G . d_s)
  K_ru = -1i*zk*gfun.*(dxus.*nx + dyus.*ny + dzus.*nz) ...
         + alpha*(gx.*dxus + gy.*dyus + gz.*dzus);

  % --- Column j_v ---
  curl_x_v = gy.*dzvs - gz.*dyvs;
  curl_y_v = gz.*dxvs - gx.*dzvs;
  curl_z_v = gx.*dyvs - gy.*dxvs;

  zvec_x_v = ny.*curl_z_v - nz.*curl_y_v;
  zvec_y_v = nz.*curl_x_v - nx.*curl_z_v;
  zvec_z_v = nx.*curl_y_v - ny.*curl_x_v;

  ndotE_v  = 1i*zk*gfun .* (nx.*dxvs + ny.*dyvs + nz.*dzvs);
  v3x_v = alpha*(ndotE_v.*nx - 1i*zk*gfun.*dxvs) - zvec_x_v;
  v3y_v = alpha*(ndotE_v.*ny - 1i*zk*gfun.*dyvs) - zvec_y_v;
  v3z_v = alpha*(ndotE_v.*nz - 1i*zk*gfun.*dzvs) - zvec_z_v;

  K_uv = dxut.*v3x_v + dyut.*v3y_v + dzut.*v3z_v;
  K_vv = dxvt.*v3x_v + dyvt.*v3y_v + dzvt.*v3z_v;

  K_rv = -1i*zk*gfun.*(dxvs.*nx + dyvs.*ny + dzvs.*nz) ...
         + alpha*(gx.*dxvs + gy.*dyvs + gz.*dzvs);

  % --- Column rho ---
  % E_rho = -grad_t G,  J=0
  % n_t.E_rho = -(gx*nx + gy*ny + gz*nz)
  % nxnxE_rho = (n_t.E_rho)*n_t - E_rho = -(gx*nx+...)*n_t + grad_t G
  % zvec_rho = n_t x (grad_t G x 0) = 0   (no J curl term)
  % zvec3_rho = alpha*nxnxE_rho

  ndotE_r = -(gx.*nx + gy.*ny + gz.*nz);
  v3x_r = alpha*(ndotE_r.*nx + gx);
  v3y_r = alpha*(ndotE_r.*ny + gy);
  v3z_r = alpha*(ndotE_r.*nz + gz);

  K_ur = dxut.*v3x_r + dyut.*v3y_r + dzut.*v3z_r;
  K_vr = dxvt.*v3x_r + dyvt.*v3y_r + dzvt.*v3z_r;

  % pot_rho row, rho column:
  %   (grad_t G - ik*0).n_t + alpha*(0 - ik*G*rho)
  %   = (grad_t G).n_t - alpha*ik*G
  %   = S'_k[rho] - alpha*ik*G
  % where S'_k = n_s . grad_y G = -n_s . grad_x G
  Sp  = -(nxs.*gx + nys.*gy + nzs.*gz);
  K_rr = Sp - alpha*1i*zk*gfun;

  % --- Assemble (3*nt) x (3*ns) matrix, component-fast rows/cols ---
  submat = complex(zeros(3*nt, 3*ns));
  submat(1:nt,        1:ns)        = K_uu;
  submat(1:nt,        ns+1:2*ns)   = K_uv;
  submat(1:nt,        2*ns+1:3*ns) = K_ur;
  submat(nt+1:2*nt,   1:ns)        = K_vu;
  submat(nt+1:2*nt,   ns+1:2*ns)   = K_vv;
  submat(nt+1:2*nt,   2*ns+1:3*ns) = K_vr;
  submat(2*nt+1:3*nt, 1:ns)        = K_ru;
  submat(2*nt+1:3*nt, ns+1:2*ns)   = K_rv;
  submat(2*nt+1:3*nt, 2*ns+1:3*ns) = K_rr;

end

if strcmpi(type,'nrccie-eval')
%
%  NRCCIE field-evaluation kernel.
%
%  Representation:
%    E = ik * S_k[J]  -  grad(S_k)[rho]
%    H = curl(S_k[J]) = grad(S_k) x J
%
%  Density columns: [J_x, J_y, J_z, rho]  (4 Cartesian components / source)
%  Output rows:     [E_x, E_y, E_z, H_x, H_y, H_z]  (6 / target)
%  Matrix size: (6*nt) x (4*ns)
%
%  For a single (t,s) pair the 6x4 block is:
%    [ik*G    0      0      -dGdx ]
%    [0       ik*G   0      -dGdy ]
%    [0       0      ik*G   -dGdz ]
%    [0      -dGdz  +dGdy    0    ]
%    [+dGdz   0     -dGdx    0    ]
%    [-dGdy  +dGdx   0       0    ]
%  where G = S_k(t,s), dGdx/y/z = d/dx G(t,s)  (w.r.t. target).

  [G, gradG] = helm3d.green(zk, src, targ);
  % G     : (nt, ns)
  % gradG : (nt, ns, 3)

  gx = gradG(:,:,1);   % (nt, ns)
  gy = gradG(:,:,2);
  gz = gradG(:,:,3);
  ikG = 1i*zk*G;

  % Build (6, nt, 4, ns) tensor, then permute+reshape to (6*nt, 4*ns)
  T = complex(zeros(6, nt, 4, ns));

  % col 1: J_x
  T(1,:,1,:) = ikG;    % E_x
  T(5,:,1,:) = gz;     % H_y = +dGdz * J_x
  T(6,:,1,:) = -gy;    % H_z = -dGdy * J_x

  % col 2: J_y
  T(2,:,2,:) = ikG;    % E_y
  T(4,:,2,:) = -gz;    % H_x = -dGdz * J_y
  T(6,:,2,:) = gx;     % H_z = +dGdx * J_y

  % col 3: J_z
  T(3,:,3,:) = ikG;    % E_z
  T(4,:,3,:) = gy;     % H_x = +dGdy * J_z
  T(5,:,3,:) = -gx;    % H_y = -dGdx * J_z

  % col 4: rho  ->  E = -gradG * rho,  H = 0
  T(1,:,4,:) = -gx;
  T(2,:,4,:) = -gy;
  T(3,:,4,:) = -gz;

  % row index: itarg + nt*(icomp-1), col index: isrc + ns*(jcomp-1)
  % (component-fast, matching stok/lap/helm eval convention)
  % permute (6,nt,4,ns) -> (nt,6,ns,4) then reshape.
  submat = reshape(permute(T, [2,1,4,3]), [6*nt, 4*ns]);

end

if strcmpi(type,'NdotCurlS')

  targnorm = targinfo.n;
  [~,grad] = helm3d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);
  nz = repmat((targnorm(3,:)).',1,ns);
  
  drut = targinfo.dru;
  dxut = repmat((drut(1,:)).',1,ns);
  dyut = repmat((drut(2,:)).',1,ns);
  dzut = repmat((drut(3,:)).',1,ns);
  
  drvt = targinfo.drv;
  dxvt = repmat((drvt(1,:)).',1,ns);
  dyvt = repmat((drvt(2,:)).',1,ns);
  dzvt = repmat((drvt(3,:)).',1,ns);
  
  drus = srcinfo.dru;
  dxus = repmat((drus(1,:)),nt,1);
  dyus = repmat((drus(2,:)),nt,1);
  dzus = repmat((drus(3,:)),nt,1);
  
  drvs = srcinfo.drv;
  dxvs = repmat((drvs(1,:)),nt,1);
  dyvs = repmat((drvs(2,:)),nt,1);
  dzvs = repmat((drvs(3,:)),nt,1);
  
  t1 = ny.*grad(:,:,3)-nz.*grad(:,:,2);
  t2 = nz.*grad(:,:,1)-nx.*grad(:,:,3);
  t3 = nx.*grad(:,:,2)-ny.*grad(:,:,1);
  
  submat = zeros(nt,2,ns);
  submat(:,1,:) = t1.*dxus+t2.*dyus+t3.*dzus;
  submat(:,2,:) = t1.*dxvs+t2.*dyvs+t3.*dzvs;
           
  submat = reshape(submat,[nt,2*ns]);
  
end
