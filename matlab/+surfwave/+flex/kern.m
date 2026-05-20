function submat= kern(srcinfo,targinfo,type,varargin)
%FLEX2D.KERN standard Modified biharmonic layer potential kernels in 2D
% 
% Syntax: submat = surfwave.flex.kern(zk,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%  
% Kernels based on G(x,y) = i/4 H_0^{(1)}(zk |x-y|)
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
% D'(x,y) = \nabla_{n_x} \nabla_{n_y} G(x,y)
%
% Input:
%   zk - complex number, Helmholtz wave number
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   targinfo - description of targets in ptinfo struct format,
%                if info not relevant (d/d2) it doesn't need to
%                be provided. sprime requires tangent info in
%                targinfo.d
%   type - string, determines kernel type
%                type == 'd', double layer kernel D
%                type == 's', single layer kernel S
%                type == 'sprime', normal derivative of single
%                      layer S'
%                type == 'dprime', normal derivative of double layer D'
%                type == 'c', combined layer kernel coef(1) D + coef(2) S
%                type == 'stau', tangential derivative of single layer
%                type == 'all', returns all four layer potentials, 
%                       [coef(1,1)*D coef(1,2)*S; coef(2,1)*D' coef(2,2)*S']
%                type == 'c2trans' returns the combined field, and the 
%                          normal derivative of the combined field
%                        [coef(1)*D + coef(2)*S; coef(1)*D' + coef(2)*S']
%                type == 'trans_rep' returns the potential corresponding
%                           to the transmission representation
%                        [coef(1)*D coef(2)*S]
%                type == 'trans_rep_prime' returns the normal derivative
%                          corresponding to the transmission representation
%                        [coef(1)*D' coef(2)*S']
%                type == 'trans_rep_grad' returns the gradient corresponding
%                         to the transmission representation
%                        [coef(1)*d_x D coef(2)*d_x S;
%                         coef(1)*d_y D coef(2)*d_y S]
%   varargin{1} - nu 
%   varargin{2} - rts
%   varargin{3} - ejs
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% see also CHNK.HELM2D.GREEN

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
    srcnorm = chnk.normal2d(srcinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((srcnorm(1,:)),nt,1);
    ny = repmat((srcnorm(2,:)),nt,1);

    submat = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);
end

% free plate kernel K21 first part (biharmonic part only) 
if strcmpi(type, 'free plate K21 first part bh')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;

   [~, ~,~, ~, fourth] = surfwave.flex.bhgreen(src, targ);  % biharmonic part

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

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   tauxtargnsrc = tauxtarg.*nx + tauytarg.*ny;
   tauxtargntarg = tauxtarg.*nxtarg + tauytarg.*nytarg;


   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   normaldot = nx.*nxtarg + ny.*nytarg;

   rn = rx.*nx + ry.*ny;

   rntarg = rx.*nxtarg + ry.*nytarg;

   
   rtautarg = rx.*tauxtarg + ry.*tauytarg;


   submat  = -(fourth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          fourth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          fourth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          fourth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny))  - ...
          (2-nu).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          fourth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          fourth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          fourth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny))  - ...           
          (1+nu)/(4*pi).*((taux.*tauxtarg + tauy.*tauytarg)./(r2) - 2*(rx.*tauxtarg + ry.*tauytarg).*(rx.*taux + ry.*tauy)./(r2.^2));
   % third kernel with singularity subtraction and H[delta']
   % 
 
end

% kernels in K11 with hilbert transform subtractions. 
% (i.e. beta*(G_{nx nx tauy} + 1/4 H + nu*(G_{taux taux tauy} + 1/4 H))

if strcmpi(type, 'free plate hilbert subtract')                                 
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;

   [~,~,~, third, ~] = surfwave.flex.hkdiffgreen(zk, src, targ, true);            % Hankel part

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   [~,grad] = chnk.lap2d.green(src,targ,true); 

   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);
  
  

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    


   submat = -((1+ nu)/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*taux) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*tauy))) + ((1+ nu)/2).*0.25*hilb -...
       ((1+ nu)/2)*nu.*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*tauy)));   % hilbert subtraction 

end




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

   [~,grad] = chnk.lap2d.green(src,targ,true); 
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

if strcmpi(type, 'free plate K21 second part')                                                         
   srcnorm = srcinfo.n;
   targnorm = targinfo.n;
   targtang = targinfo.d;
    
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);


   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


  
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 


   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
   
   zkimag = (1i)*zk;

   [~,~,~, third, ~] = surfwave.flex.helmdiffgreen(zk, src, targ);            % Hankel part
   [~,~,~, thirdK, ~] = surfwave.flex.helmdiffgreen(zkimag, src, targ);     % modified bessel K part

   submat = ((1-nu)/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 4).*(nytarg.*nytarg.*ny)) - ...
        (1-nu)/(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*nx) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*ny))) - ...
       ((1-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*ny)) - ...
        (1-nu)/(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*ny)));                                 % (1-nu)*(- G_{nx nx ny} + G_{ny taux taux})  
end

% Updated part in K21 that is coupled with Hilbert transform. 
%(i.e. (1-nu)*(beta*(G_{taux taux tauy} + 1/4 H) - beta*(G_{nx nx tauy} + 1/4 H)). 

if strcmpi(type, 'free plate K21 hilbert part bh')                                                         
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;

   [~, ~, ~,third, ~] = surfwave.flex.bhgreen(zk, src, targ);            % Hankel part
   third = 2*zk^2*third; 

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   [~,grad] = chnk.lap2d.green(src,targ,true); 
   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);
  

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    


   submat =  (1-nu).*(-((1+ nu)/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + ...
        ((1+ nu)/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*taux) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*tauy))))  ;  % + ((1+ nu)/2).*0.25*hilb - ((1+ nu)/2).*0.25*hilb)



end

if strcmpi(type, 'free plate K21 hilbert part')                                                         
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;

   [~, ~, ~,third, ~] = surfwave.flex.hkdiffgreen(zk, src, targ, true);            % Hankel part

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

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    

   submat =  (1-nu).*(-((1+ nu)/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + ...
        ((1+ nu)/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*taux) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*tauy))))  ;  % + ((1+ nu)/2).*0.25*hilb - ((1+ nu)/2).*0.25*hilb)



end



% Updated part in K22. (i.e. (1-nu)*(G_{taux taux} - G_{nx nx})

if strcmpi(type, 'free plate K22 second part')
   targnorm = targinfo.n;
   targtang = targinfo.d;



   [~, ~, hess, ~, ~] = surfwave.flex.helmdiffgreen(zk, src, targ);            % Hankel part
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, ~, hessK, ~, ~] = surfwave.flex.helmdiffgreen(zkimag, src, targ);     % modified bessel K part
   
   

  

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


  
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 


   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    
   
   submat =  (1-nu).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hessK(:, :, 1).*tauxtarg.*tauxtarg + hessK(:, :, 2).*(2*tauxtarg.*tauytarg) + ...
           hessK(:, :, 3).*tauytarg.*tauytarg)-...
           (1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)-...
           1/(2*zk^2).*(hessK(:, :, 1).*nxtarg.*nxtarg + hessK(:, :, 2).*(2*nxtarg.*nytarg) + hessK(:, :, 3).*nytarg.*nytarg)));


end

if strcmpi(type, 'free plate hilbert unsubtract')                                 
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;

   [~,~,~, third, ~] = surfwave.flex.helmdiffgreen(zk, src, targ);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~,~, ~, thirdK, ~] = surfwave.flex.helmdiffgreen(zkimag, src, targ);     % modified bessel K part

  
  

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    


   submat = -((1+ nu)/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*taux) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*tauy))) -...
       ((1+ nu)/2)*nu.*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*tauy))); % + ((1+ nu)/2)*nu.*0.25*hilb;   % hilbert unsubtraction 

end


% free plate kernel for the interior modified biharmonic problem (K11 with no
% hilbert transform subtraction, K12 kernel, K22 with no curvature part) K21 is
% handled in a separate type.
if strcmpi(type, 'free plate first part interior')
   srcnorm = srcinfo.n;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   



   [~, ~, hess, third, ~] = surfwave.flex.helmdiffgreen(zk, src, targ);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, ~, hessK, thirdK, ~] = surfwave.flex.helmdiffgreen(zkimag, src, targ);     % modified bessel K part
   


   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   
   
   
   K11 = (1/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 4).*(nytarg.*nytarg.*ny)) - ...
        1/(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*nx) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*ny))) + ...
       nu.*(1/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*ny)) - ...
        1/(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*ny)));  % first kernel with no hilbert transforms (G_{nx nx ny + nu G_{taux taux ny}).
       
    
       

   K12 =  1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)-...
           1/(2*zk^2).*(hessK(:, :, 1).*nxtarg.*nxtarg + hessK(:, :, 2).*(2*nxtarg.*nytarg) + hessK(:, :, 3).*nytarg.*nytarg)+...
           nu/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           nu/(2*zk^2).*(hessK(:, :, 1).*tauxtarg.*tauxtarg + hessK(:, :, 2).*(2*tauxtarg.*tauytarg) + ...
           hessK(:, :, 3).*tauytarg.*tauytarg) ;    % G_{nx nx} + nu G_{taux taux}
   
   K21 = 0;




   K22 = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + thirdK(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg)+...
        thirdK(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + thirdK(:, :, 4).*(nytarg.*nytarg.*nytarg)) +...
        (2-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg)) - ...
        (2-nu)/(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + thirdK(:, :, 4).*(tauytarg.*tauytarg.*nytarg)); % G_{nx nx nx} + (2-nu) G_{taux taux nx}
    
  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;

end

% free plate kernel K21 for the interior modified biharmonic problem. This part
% handles the singularity subtraction and swap the evaluation to its
% asymptotic expansions if the targets and sources are close. 

if strcmpi(type, 'free plate K21 first part interior')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;

   [~, ~,~, ~, fourth] = surfwave.flex.helmdiffgreen(zk, src, targ, true);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, ~, ~,~, fourthK] = surfwave.flex.helmdiffgreen(zkimag, src, targ, true);     % modified bessel K part 

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   tauxtargnsrc = tauxtarg.*nx + tauytarg.*ny;
   tauxtargntarg = tauxtarg.*nxtarg + tauytarg.*nytarg;


   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   normaldot = nx.*nxtarg + ny.*nytarg;

   rn = rx.*nx + ry.*ny;

   rntarg = rx.*nxtarg + ry.*nytarg;

   
   rtautarg = rx.*tauxtarg + ry.*tauytarg;


   submat  = (1/(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          fourth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          fourth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          fourth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny)) -...
          1/(2*zk^2).*(fourthK(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + fourthK(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          fourthK(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          fourthK(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          fourthK(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny))) + ...
          ((2-nu)/(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          fourth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          fourth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          fourth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny)) -...
          (2-nu)/(2*zk^2).*(fourthK(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + fourthK(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          fourthK(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          fourthK(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          fourthK(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny))) - ...
          (3/(2*pi).*(rn.*rntarg./(r2.^2)) + 3/(2*pi).*(rntarg.^2).*(normaldot)./(r2.^2) -...
            (2/pi).*(rn.*(rntarg.^3))./(r2.^3)) - ...
            (2-nu).*(-1/(2*pi).*(tauxtargnsrc).*(tauxtargntarg)./(r2) +...
            1/(pi).*(rn.*rtautarg.*tauxtargntarg)./(r2.^2) + ...
            1/(pi).*(rtautarg.*tauxtargnsrc.*rntarg)./(r2.^2) -....
            2/(pi).*(rn.*rntarg.*(rtautarg.*rtautarg))./(r2.^3) + ...
            1/(2*pi).*(rn.*rntarg)./(r2.^2)) +...           
          ((1+nu)/(4*pi)).*((taux.*tauxtarg + tauy.*tauytarg)./(r2) - 2*(rx.*tauxtarg + ry.*tauytarg).*(rx.*taux + ry.*tauy)./(r2.^2))+...
          (3/(4*pi)).*(nxtarg.*nx + nytarg.*ny)./r2  + (2-nu)/(4*pi).*(nxtarg.*nx + nytarg.*ny)./r2 - ...
          (2-nu)/(2*pi).*(nxtarg.*nx + nytarg.*ny).*(rx.*tauxtarg + ry.*tauytarg).*(rx.*tauxtarg + ry.*tauytarg)./(r2.^2);
   % third kernel with singularity subtraction and H[delta']
   % 
 
end

% Updated part in K21 (interior) that is not coupled with Hilbert transform. 
%(i.e. (1-nu)*(G_{nx nx ny} - G_{taux taux ny}). 

if strcmpi(type, 'free plate K21 second part interior')                                                         
   srcnorm = srcinfo.n;
   targnorm = targinfo.n;
   targtang = targinfo.d;
    
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);
  
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
   
   zkimag = (1i)*zk;

   [~,~,~, third, ~] = surfwave.flex.helmdiffgreen(zk, src, targ);            % Hankel part
   [~,~,~, thirdK, ~] = surfwave.flex.helmdiffgreen(zkimag, src, targ);     % modified bessel K part

   submat = -((1-nu)/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 4).*(nytarg.*nytarg.*ny)) - ...
        (1-nu)/(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*nx) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*ny))) + ...
       ((1-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*ny)) - ...
        (1-nu)/(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*ny)));                                 % (1-nu)*(- G_{nx nx ny} + G_{ny taux taux})  
end


if strcmpi(type, 'free plate eval first')                                               % G_{ny}
   srcnorm = srcinfo.n;
   [~,grad] = surfwave.flex.gsflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   submat = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 

end

if strcmpi(type, 'free plate eval first hilbert')                                               % G_{tauy} (supposed to coupled with the hilbert transform)
   srctang = srcinfo.d;

   [~,grad] = surfwave.flex.gsflex(rts,ejs,src,targ);       
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);
 
   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   submat = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*taux + grad(:, :, 2).*tauy));                    % G_{tauy}
end


if strcmpi(type, 'free plate eval second')                                          % G = 1/(2k^2) (i/4 H_0^{1} - 1/2pi K_0)

   [val,~] = surfwave.flex.gsflex(rts,ejs,src,targ); 
   zk = 1/sqrt(2);
   val = 2*zk^2*val;

   submat = 1/(2*zk^2).*val ;

end


if strcmpi(type, 'free plate eval first int eq')                                               % G_{ny}
   srcnorm = srcinfo.n;
   [~,gradgs] = surfwave.flex.gsflex(rts,ejs,src,targ);
   [~,gradgphi] = surfwave.flex.gphiflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   gradgs = 2*zk^2*gradgs;
   gradgphi = 2*zk^2*gradgphi;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   submat = 1/2*(-1/(2*zk^2).*(gradgs(:, :, 1).*(nx) + gradgs(:, :, 2).*ny)) ...
                - (-1/(2*zk^2).*(gradgphi(:, :, 1).*(nx) + gradgphi(:, :, 2).*ny)); 

end

if strcmpi(type, 'free plate eval first hilbert int eq')                                               % G_{tauy} (supposed to coupled with the hilbert transform)
   srctang = srcinfo.d;

   [~,gradgs] = surfwave.flex.gsflex(rts,ejs,src,targ); 
   [~,gradgphi] = surfwave.flex.gphiflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   gradgs = 2*zk^2*gradgs;
   gradgphi = 2*zk^2*gradgphi;

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);
 
   taux = dx./ds; 
   tauy = dy./ds;

   submat = 1/2*((1 + nu)/2).*(-1/(2*zk^2).*(gradgs(:, :, 1).*taux + gradgs(:, :, 2).*tauy)) ...
            - ((1 + nu)/2).*(-1/(2*zk^2).*(gradgphi(:, :, 1).*taux + gradgphi(:, :, 2).*tauy));
end


if strcmpi(type, 'free plate eval second int eq')                                          % G = 1/(2k^2) (i/4 H_0^{1} - 1/2pi K_0)

   [valgs,~] = surfwave.flex.gsflex(rts,ejs,src,targ); 
   [valgphi,~] = surfwave.flex.gphiflex(rts,ejs,src,targ); 
   zk = 1/sqrt(2);
   valgs = 2*zk^2*valgs;
   valgphi = 2*zk^2*valgphi;

   submat = 1/2*1/(2*zk^2).*valgs - 1/(2*zk^2).*valgphi ;

end


if strcmpi(type, 'free plate eval first gs')                                               % G_{ny}
   srcnorm = srcinfo.n;
   [~,grad] = surfwave.flex.gsflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   submat = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 

end

if strcmpi(type, 'free plate eval first hilbert gs')                                               % G_{tauy} (supposed to coupled with the hilbert transform)
   srctang = srcinfo.d;

   [~,grad] = surfwave.flex.gsflex(rts,ejs,src,targ);       
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);
 
   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   submat = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*taux + grad(:, :, 2).*tauy));                    % G_{tauy}
end


if strcmpi(type, 'free plate eval second gs')                                          % G = 1/(2k^2) (i/4 H_0^{1} - 1/2pi K_0)

   [val] = surfwave.flex.gsflex(rts,ejs,src,targ); 
   zk = 1/sqrt(2);
   val = 2*zk^2*val;

   submat = 1/(2*zk^2).*val ;

end

if strcmpi(type, 'free plate eval first gphi')                                               % G_{ny}
   srcnorm = srcinfo.n;
   [~,grad] = surfwave.flex.gphiflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   submat = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 

end

if strcmpi(type, 'free plate eval first hilbert gphi')                                               % G_{tauy} (supposed to coupled with the hilbert transform)
   srctang = srcinfo.d;

   [~,grad] = surfwave.flex.gphiflex(rts,ejs,src,targ);       
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);
 
   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   submat = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*taux + grad(:, :, 2).*tauy));                    % G_{tauy}
end

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


if strcmpi(type, 'free plate eval second gphi')                                          % G = 1/(2k^2) (i/4 H_0^{1} - 1/2pi K_0)

   [val,~] = surfwave.flex.gphiflex(rts,ejs,src,targ); 
   zk = 1/sqrt(2);
   val = 2*zk^2*val;

   submat = 1/(2*zk^2).*val ;

end

if strcmpi(type, 'free plate eval kerns gs')                                               % G_{ny}
   
   submat = zeros(3*nt,ns);

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

   submat(1:3:end,:) = gsn;
   submat(2:3:end,:) = gstau;
   submat(3:3:end,:) = gs;

end

if strcmpi(type, 'free_plate_eval_gs_1')                                               % G_{ny}
   
   srcnorm = srcinfo.n;

   [~,grad] = surfwave.flex.gsflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   submat = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 

end

if strcmpi(type, 'free_plate_eval_gs_2')                                               % G_{ny}
   
   srctang = srcinfo.d;

   [~,grad] = surfwave.flex.gsflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);
 
   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;
   
   submat = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*taux + grad(:, :, 2).*tauy));                    % G_{tauy}

end

if strcmpi(type, 'free_plate_eval_gs_3')                                               % G_{ny}

   val = surfwave.flex.gsflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   
   submat = 1/(2*zk^2).*val ;

end

if strcmpi(type, 'free_plate_eval_gp_1')                                               % G_{ny}
   
   srcnorm = srcinfo.n;

   [~,grad] = surfwave.flex.gphiflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   submat = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 

end

if strcmpi(type, 'free_plate_eval_gp_2')                                               % G_{ny}
   
   srctang = srcinfo.d;

   [~,grad] = surfwave.flex.gphiflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   grad = 2*zk^2*grad;

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);
 
   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;
   
   submat = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*taux + grad(:, :, 2).*tauy));                    % G_{tauy}

end

if strcmpi(type, 'free_plate_eval_gp_3')                                               % G_{ny}

   val = surfwave.flex.gphiflex(rts,ejs,src,targ);
   zk = 1/sqrt(2);
   
   submat = 1/(2*zk^2).*val ;

end

if strcmpi(type, 'free plate eval kerns gphi')                                               % G_{ny}
   
   submat = zeros(3*nt,ns);

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

   submat(1:3:end,:) = gphin;
   submat(2:3:end,:) = gphitau;
   submat(3:3:end,:) = gphi;

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


% free plate kernel for the modified biharmonic problem (K11 with no
% hilbert transform subtraction, K12 kernel, K22 with no curvature part) K21 is
% handled in a separate type.
if strcmpi(type, 'v2b_bc_1')
   % srcnorm = srcinfo.n;
   % srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;

   [~, ~, hess] = surfwave.flex.gsflex(rts,ejs,src,targ);   

   zk = 1/sqrt(2);
   hess = 2*zk^2*hess;

   % nx = repmat(srcnorm(1,:),nt,1);
   % ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);


   % dx = repmat(srctang(1,:),nt,1);
   % dy = repmat(srctang(2,:),nt,1);

   % ds = sqrt(dx.*dx+dy.*dy); 
   % 
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

    
  submat = K12;

end

% free plate kernel for the modified biharmonic problem (K11 with no
% hilbert transform subtraction, K12 kernel, K22 with no curvature part) K21 is
% handled in a separate type.
if strcmpi(type, 'v2b_bc_2')
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


   K22 = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))  +...
        (2-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg)) + ... % G_{nx nx nx} + (2-nu) G_{taux taux nx}
        + kappax.*(1-nu).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg));
          
    
  submat = K22;

end

if strcmpi(type, 'gs_v2b_bc_1')
   targnorm = targinfo.n;
   targtang = targinfo.d;

   [~, ~, hess] = surfwave.flex.gsflex(rts,ejs,src,targ);   

   zk = 1/sqrt(2);
   hess = 2*zk^2*hess;

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   submat =  1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)+...
           nu/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg) ;    % G_{nx nx} + nu G_{taux taux}

end

if strcmpi(type, 'gs_v2b_bc_2')
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;

   [~, ~, hess, third] = surfwave.flex.gsflex(rts,ejs,src,targ);   

   zk = 1/sqrt(2);
   hess = 2*zk^2*hess;
   third = 2*zk^2*third;

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; 

   submat = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))  +...
        (2-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg)) + ... % G_{nx nx nx} + (2-nu) G_{taux taux nx}
        + kappax.*(1-nu).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg));

end

if strcmpi(type, 'gp_v2b_bc_1')
   targnorm = targinfo.n;
   targtang = targinfo.d;

   [~, ~, hess] = surfwave.flex.gphiflex(rts,ejs,src,targ);   

   zk = 1/sqrt(2);
   hess = 2*zk^2*hess;

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   submat =  1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)+...
           nu/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg) ;    % G_{nx nx} + nu G_{taux taux}

end

if strcmpi(type, 'gp_v2b_bc_2')
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;

   [~, ~, hess, third] = surfwave.flex.gphiflex(rts,ejs,src,targ);   

   zk = 1/sqrt(2);
   hess = 2*zk^2*hess;
   third = 2*zk^2*third;

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; 

   submat = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))  +...
        (2-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg)) + ... % G_{nx nx nx} + (2-nu) G_{taux taux nx}
        + kappax.*(1-nu).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg));

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


