function submat= kern(rts,ejs,srcinfo,targinfo,type,varargin)
%SURFWAVE.CAPILLARY.KERN capillary surface wave kernels
%
% Syntax: submat = surfwave.capillary.kern(rts,ejs,srcinfo,targinfo,type)
%
% G_S and G_phi are the capillary Green's functions with dispersion roots
% rts and partial fraction residues ejs (SURFWAVE.CAPILLARY.FIND_ROOTS_CAPILLARY).
%
% Input:
%   rts, ejs - dispersion roots and residues
%   srcinfo  - sources in ptinfo struct format; srcinfo.n for the *_d types
%   targinfo - targets in ptinfo struct format; targinfo.n for the
%              *_sprime and *_dprime types
%   type - string, determines kernel type
%        'gs_s', 'gphi_s'           - G_S, G_phi
%        'gs_d', 'gphi_d'           - source normal derivatives
%        'gs_sprime', 'gphi_sprime' - target normal derivatives
%        'gs_dprime'                - target normal derivative of 'gs_d'
%        'lap_gs', 'lap_gphi'       - surface Laplacians
%        's3d_gphi'                 - 3D single layer applied to G_phi
%
% Output:
%   submat - (nt,ns) kernel matrix
%
% See also GREEN2D.HELM

src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'gs_d')
  srcnorm = srcinfo.n;
  [~,grad] = surfwave.capillary.gshelm(rts,ejs,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'gphi_d')
  srcnorm = srcinfo.n;
  [~,grad] = surfwave.capillary.gphihelm(rts,ejs,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

% normal derivative of single layer
if strcmpi(type,'gs_sprime')
  targnorm = targinfo.n;
  [~,grad] = surfwave.capillary.gshelm(rts,ejs,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);

  submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end

% normal derivative of single layer
if strcmpi(type,'gphi_sprime')
  targnorm = targinfo.n;
  [~,grad] = surfwave.capillary.gphihelm(rts,ejs,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);

  submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'gs_dprime')
  srcnorm = srcinfo.n;
  targnorm = targinfo.n;
  [~,~,hess] = surfwave.capillary.gshelm(rts,ejs,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  submat = -(hess(:,:,1).*nx + hess(:,:,2).*ny).*nxtarg ...
      - (hess(:,:,2).*nx + hess(:,:,3).*ny).*nytarg;
end

% single layer Gs
if strcmpi(type,'gs_s')
  submat = surfwave.capillary.gshelm(rts,ejs,src,targ);
end

% laplacian single layer
if strcmpi(type,'lap_gs')
  [~,~,hess] = surfwave.capillary.gshelm(rts,ejs,src,targ);
  submat = hess(:,:,1) + hess(:,:,3);
end

% single layer Gphi
if strcmpi(type,'gphi_s')
  submat = surfwave.capillary.gphihelm(rts,ejs,src,targ);
end

% single layer S3d Gphi
if strcmpi(type,'s3d_gphi')
  submat = surfwave.capillary.s3dgphi(rts,ejs,src,targ);
end

% single layer \Del Gphi
if strcmpi(type,'lap_gphi')
  [~,~,hess] = surfwave.capillary.gphihelm(rts,ejs,src,targ);
  submat = hess(:,:,1) + hess(:,:,3);
end

