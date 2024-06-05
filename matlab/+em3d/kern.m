
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
