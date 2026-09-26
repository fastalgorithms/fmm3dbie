function submat= kern(srcinfo,targinfo,type,varargin)
%SURFWAVE.FLEX.KERN flexural surface wave kernels
%
% Syntax: submat = surfwave.flex.kern(srcinfo,targinfo,type,nu,rts,ejs)
%
% G_S and G_phi are the flexural Green's functions with dispersion roots
% rts and partial fraction residues ejs (SURFWAVE.FLEX.FIND_ROOTS_FLEX);
% nu is the Poisson ratio.
%
% Input:
%   srcinfo, targinfo - ptinfo structs; boundary points carry n, d, d2
%   type - string, determines kernel type
%     volume kernels:
%        'gs_s', 'gphi_s'         - G_S, G_phi
%        's3d_gphi'               - 3D single layer applied to G_phi
%        'gphi_bilap'             - bilaplacian of G_phi
%     volume to boundary:
%        'gs_v2b', 'gphi_v2b'     - free plate conditions applied to
%                                   G_S, G_phi, (2*nt,ns)
%     boundary to volume:
%        'free_plate_gs_eval', 'free_plate_gphi_eval'
%                                 - free plate representation, (nt,3*ns)
%     boundary to boundary:
%        'free plate first part', 'free plate first part bh',
%        'free plate hilbert', 'free plate hilbert bh'
%                                 - free plate system, split into the
%                                   parts acting on the density and on its
%                                   Hilbert transform
%        'hilb'                   - Hilbert transform kernel
%
% Output:
%   submat - kernel matrix

src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

nu = varargin{1};
rts = varargin{2};
ejs = varargin{3};

% free plate kernel for the modified biharmonic problem (K11 with no
% hilbert transform subtraction, K12 kernel, K22 with no curvature part) K21 is
% handled in a separate type.
if strcmpi(type, 'free plate first part')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;

   [~, ~, hess, third] = surfwave.flex.gsflex(rts,ejs,src,targ);     
   [~, ~, ~, ~, fourth] = surfwave.flex.gsflex(rts,ejs,src,targ,true);   
   %[~, ~, ~, ~, fourthbh] = surfwave.flex.bhgreen(src, targ);    
   % fourth = fourth + 2*zk^2*fourthbh;

   zk = 1/sqrt(2);
   hess = 2*zk^2*hess;
   third = 2*zk^2*third;
   fourth = 2*zk^2*fourth;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy); 

   taux = dx ./ ds;
   tauy = dy ./ ds;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; % target curvature

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   
   K11 = -(1/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 4).*(nytarg.*nytarg.*ny))) - ...
       nu./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*ny)) ;  % first kernel with no hilbert transforms (G_{nx nx ny + nu G_{taux taux ny}).

   K12 =  1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)+...
           nu/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg) ;    % G_{nx nx} + nu G_{taux taux}
   
   K21 = kappax./(2*zk^2).*(1-nu).*((third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
            third(:, :, 4).*(nytarg.*nytarg.*ny)) - ...
           (third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*ny)) ) ...
        - 1/(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          fourth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          fourth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          fourth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny)) - ...
          ((2-nu)/(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          fourth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          fourth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          fourth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny)) ) ; % - ...          
          % (1+nu)/(4*pi).*((taux.*tauxtarg + tauy.*tauytarg)./(r2) - 2*(rx.*tauxtarg + ry.*tauytarg).*(rx.*taux + ry.*tauy)./(r2.^2));

   K22 = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))  +...
        (2-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg)) + ... % G_{nx nx nx} + (2-nu) G_{taux taux nx}
        + kappax.*(1-nu).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg));
          
    
  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;

end

% free plate kernel for the modified biharmonic problem (K11 with no
% hilbert transform subtraction, K12 kernel, K22 with no curvature part) K21 is
% handled in a separate type.
if strcmpi(type, 'free plate first part bh')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;

   alpha = 1./sum(ejs.*rts.^4);

   % [~, ~, hess, third, ~] = surfwave.flex.hkdiffgreen(zk, src, targ);     
   % [~, ~, ~, ~, fourth] = surfwave.flex.hkdiffgreen(zk, src, targ, true);   
   [~, ~, ~, ~, fourthbh] = surfwave.flex.bhgreen(src, targ);    

   zk = 1/sqrt(2); % get rid of zk (they currently cancel anyway)
   fourth = 2*zk^2*fourthbh;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy); 

   taux = dx ./ ds;
   tauy = dy ./ ds;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; % target curvature

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   
   K11 = 0;

   K12 =  0;

   K21 = - 1/(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          fourth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          fourth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          fourth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny)) - ...
          ((2-nu)/(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          fourth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          fourth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          fourth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny)) ) - ...          
          (1+nu)/(4*pi).*((taux.*tauxtarg + tauy.*tauytarg)./(r2) - 2*(rx.*tauxtarg + ry.*tauytarg).*(rx.*taux + ry.*tauy)./(r2.^2));

   K21 = K21*2/alpha;

   K22 = 0;
    
  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;

end

if strcmpi(type,'hilb')
    srcnorm = [srcinfo.d(2,:); -srcinfo.d(1,:)]./vecnorm(srcinfo.d(1:2,:));
    [~,grad] = green2d.lap(src,targ,true);
    nx = repmat((srcnorm(1,:)),nt,1);
    ny = repmat((srcnorm(2,:)),nt,1);

    submat = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);
end

% kernels in K11 with hilbert transform subtractions. 
% (i.e. beta*(G_{nx nx tauy} + 1/4 H + nu*(G_{taux taux tauy} + 1/4 H))

% kernels in K21 that are coupled with Hilbert transforms. 
if strcmpi(type, 'free plate hilbert')                                  
   srctang = srcinfo.d;
   srcnorm = srcinfo.n;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;
    
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
  
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; % target curvature

   [~, ~,~, third, ~] = surfwave.flex.gsflex(rts,ejs, src, targ, true);   
   zk = 1/sqrt(2);
   third = 2*zk^2*third;

   % [~, ~,~, thirdbh, ~] = surfwave.flex.bhgreen(src, targ);            % Hankel part
   % third = third + 2*zk^2*thirdbh;

   [~, ~,~, ~, fourth] = surfwave.flex.gsflex(rts,ejs, src, targ, false);
   fourth = 2*zk^2*fourth;

   K11 =  -(1+ nu)/2*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) )  - ...
       (1+ nu)/2*nu.*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy))) ;

   K21 = kappax.*(1-nu).*(-((1+ nu)/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + ...
        ((1+ nu)/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)))) ...
        -(1+ nu)/2.*(1/(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*taux) + ...
          fourth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*tauy + 3*nxtarg.*nxtarg.*nytarg.*taux) + ...
          fourth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*tauy + 3*nxtarg.*nytarg.*nytarg.*taux) +...
          fourth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*tauy + nytarg.*nytarg.*nytarg.*taux) +...
          fourth(:, :, 5).*(nytarg.*nytarg.*nytarg.*tauy)) ) - ...
          ((2-nu)/2)*(1+nu).*(1/(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*taux) + ...
          fourth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*tauy + tauxtarg.*tauxtarg.*nytarg.*taux + 2*tauxtarg.*tauytarg.*nxtarg.*taux) + ...
          fourth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*tauy + 2*tauxtarg.*tauytarg.*nxtarg.*tauy + tauytarg.*tauytarg.*nxtarg.*taux + 2*tauxtarg.*tauytarg.*nytarg.*taux) +...
          fourth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*tauy + 2*tauxtarg.*tauytarg.*nytarg.*tauy + tauytarg.*tauytarg.*nytarg.*taux) +...
         fourth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*tauy)) ) ;

  K12 = 0;

  K22 = 0;

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;
    
end

% kernels in K21 that are coupled with Hilbert transforms. 
if strcmpi(type, 'free plate hilbert bh')                                  
   srctang = srcinfo.d;
   srcnorm = srcinfo.n;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;

   alpha = 1./sum(ejs.*rts.^4);
    
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
  
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; % target curvature

   [~,grad] = green2d.lap(src,targ,true); 
   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);

   zk = 1/sqrt(2); % get rid of zk (they currently cancel anyway)
   [~, ~,~, thirdbh, ~] = surfwave.flex.bhgreen(src, targ);            % Hankel part
   third = 2*zk^2*thirdbh;

   K11 =  -(1+ nu)/2*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) )  - ...
       (1+ nu)/2*nu.*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + (1+ nu)/2*(1+nu).*0.25*hilb;

   K21 = kappax.*(1-nu).*(-((1+ nu)/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + ...
        ((1+ nu)/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy))))  ;

  K11 = K11*2/alpha;
  K21 = K21*2/alpha;

  K12 = 0;

  K22 = 0;

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;
    
end

% Updated part in K21 that is not coupled with Hilbert transform. 
%(i.e. (1-nu)*kappa*(-G_{nx nx ny} + G_{taux taux ny}). 

% Updated part in K21 that is coupled with Hilbert transform. 
%(i.e. (1-nu)*(beta*(G_{taux taux tauy} + 1/4 H) - beta*(G_{nx nx tauy} + 1/4 H)). 

% Updated part in K22. (i.e. (1-nu)*(G_{taux taux} - G_{nx nx})

% free plate kernel K21 for the interior modified biharmonic problem. This part
% handles the singularity subtraction and swap the evaluation to its
% asymptotic expansions if the targets and sources are close. 

% Updated part in K21 (interior) that is not coupled with Hilbert transform. 
%(i.e. (1-nu)*(G_{nx nx ny} - G_{taux taux ny}). 

if strcmpi(type, 'gs_s')                                          % G = 1/(2k^2) (i/4 H_0^{1} - 1/2pi K_0)

   submat = surfwave.flex.gsflex(rts,ejs,src,targ); 

end

if strcmpi(type, 'gphi_s')                                          % G = 1/(2k^2) (i/4 H_0^{1} - 1/2pi K_0)

   submat = surfwave.flex.gphiflex(rts,ejs,src,targ); 

end

if strcmpi(type, 's3d_gphi')                                          % G = 1/(2k^2) (i/4 H_0^{1} - 1/2pi K_0)

   submat = surfwave.flex.s3dgphiflex(rts,ejs,src,targ); 

end

if strcmpi(type, 'gphi_bilap')                                          % G = 1/(2k^2) (i/4 H_0^{1} - 1/2pi K_0)

   [~,~,~,~,submat] = surfwave.flex.gphiflex(rts,ejs,src,targ); 

end

if strcmpi(type, 'free_plate_gs_eval')                                               % G_{ny}

   submat = zeros(nt,3*ns);

   srcnorm = srcinfo.n;
   srctang = srcinfo.d;

   [val,grad] = surfwave.flex.gsflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   gsn = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 
   gstau = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*taux + grad(:, :, 2).*tauy));                    % G_{tauy}
   gs = 1/(2*zk^2).*val ;

   submat(:,1:3:end) = gsn;
   submat(:,2:3:end) = gstau;
   submat(:,3:3:end) = gs;

end

if strcmpi(type, 'free_plate_gphi_eval')                                               % G_{ny}

   submat = zeros(nt,3*ns);

   srcnorm = srcinfo.n;
   srctang = srcinfo.d;

   [val,grad] = surfwave.flex.gphiflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   gphin = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 
   gphitau = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*taux + grad(:, :, 2).*tauy));                    % G_{tauy}
   gphi = 1/(2*zk^2).*val ;

   submat(:,1:3:end) = gphin;
   submat(:,2:3:end) = gphitau;
   submat(:,3:3:end) = gphi;

end

if strcmpi(type, 'gs_v2b')

   submat = zeros(2*nt,ns);

   % srcnorm = srcinfo.n;
   % srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;

   [~, ~, hess, third] = surfwave.flex.gsflex(rts,ejs,src,targ);   

   zk = 1/sqrt(2);
   hess = 2*zk^2*hess;
   third = 2*zk^2*third;
   % fourth = 2*zk^2*fourth;

   % nx = repmat(srcnorm(1,:),nt,1);
   % ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   % dx = repmat(srctang(1,:),nt,1);
   % dy = repmat(srctang(2,:),nt,1);

   % ds = sqrt(dx.*dx+dy.*dy); 

   % taux = dx ./ ds;
   % tauy = dy ./ ds;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; % target curvature

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   K12 =  1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)+...
           nu/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg) ;    % G_{nx nx} + nu G_{taux taux}

   K22 = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))  +...
        (2-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg)) + ... % G_{nx nx nx} + (2-nu) G_{taux taux nx}
        + kappax.*(1-nu).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg));
    
   submat(1:2:end,:) = K12;
   submat(2:2:end,:) = K22;
    
end

if strcmpi(type, 'gphi_v2b')

   submat = zeros(2*nt,ns);

   % srcnorm = srcinfo.n;
   % srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;

   [~, ~, hess, third] = surfwave.flex.gphiflex(rts,ejs,src,targ);   

   zk = 1/sqrt(2);
   hess = 2*zk^2*hess;
   third = 2*zk^2*third;
   % fourth = 2*zk^2*fourth;

   % nx = repmat(srcnorm(1,:),nt,1);
   % ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   % dx = repmat(srctang(1,:),nt,1);
   % dy = repmat(srctang(2,:),nt,1);

   % ds = sqrt(dx.*dx+dy.*dy); 

   % taux = dx ./ ds;
   % tauy = dy ./ ds;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; % target curvature

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   K12 =  1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)+...
           nu/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg) ;    % G_{nx nx} + nu G_{taux taux}

   K22 = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))  +...
        (2-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg)) + ... % G_{nx nx nx} + (2-nu) G_{taux taux nx}
        + kappax.*(1-nu).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg));
    
   submat(1:2:end,:) = K12;
   submat(2:2:end,:) = K22;
    
end

end

