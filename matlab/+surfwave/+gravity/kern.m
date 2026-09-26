function submat= kern(g,srcinfo,targinfo,type,coefs)
%SURFWAVE.GRAVITY.KERN gravity surface wave kernels
%
% Syntax: submat = surfwave.gravity.kern(g,srcinfo,targinfo,type,coefs)
%
% Input:
%   g        - gravity parameter, dispersion root rho = g/2
%   srcinfo  - sources in ptinfo struct format
%   targinfo - targets in ptinfo struct format; targinfo.n, .d and .d2
%              for 'free_plate_bcs'
%   type - string, determines kernel type
%        'gs_s', 'gphi_s'           - G_S, G_phi = G_S/g
%        'gsgrad_s', 'gphigrad_s'   - target gradients, (2*nt,ns) with
%                                     d/dx and d/dy interleaved by row
%        'free_plate_bcs'           - free plate conditions applied to
%                                     G_S, (2*nt,ns); coefs = nu
%        'vol'                      - (a*bilaplacian + b) G_S - g G_phi;
%                                     coefs = [a b]
%   coefs    - see type
%
% Output:
%   submat - kernel matrix
%
% See also SURFWAVE.GRAVITY.GSGRAV, SURFWAVE.GRAVITY.GSGRAV_DER

src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'gs_s')
  submat = surfwave.gravity.gsgrav(g/2, 1, src, targ);
end

if strcmpi(type,'gphi_s')
  submat = surfwave.gravity.gsgrav(g/2, 1/g, src, targ);
end

if strcmpi(type,'gsgrad_s')
[~,grad] = surfwave.gravity.gsgrav_der(g/2, 1, src, targ);

  gx = grad(:,:,1);
  gy = grad(:,:,2);

  submat = zeros(2*nt,ns);
  submat(1:2:end,:) = gx;
  submat(2:2:end,:) = gy;

elseif strcmpi(type,'gphigrad_s')
[~,grad] = surfwave.gravity.gsgrav_der(g/2, 1/g, src, targ);

  gx = grad(:,:,1);
  gy = grad(:,:,2);

  submat = zeros(2*nt,ns);
  submat(1:2:end,:) = gx;
  submat(2:2:end,:) = gy;

end

% bending moment and Kirchhoff shear
if strcmpi(type,'free_plate_bcs')
  nu = coefs(1);
  [~,~,hess,third] = surfwave.gravity.gsgrav_der(g/2, 1, src, targ);

  nxt = repmat(targinfo.n(1,:).',1,ns);  nyt = repmat(targinfo.n(2,:).',1,ns);
  dx1 = repmat(targinfo.d(1,:).',1,ns);  dy1 = repmat(targinfo.d(2,:).',1,ns);
  ds1 = sqrt(dx1.^2+dy1.^2);
  taux = dx1./ds1;  tauy = dy1./ds1;
  d2x = repmat(targinfo.d2(1,:).',1,ns); d2y = repmat(targinfo.d2(2,:).',1,ns);
  kappat = (dx1.*d2y - d2x.*dy1) ./ (ds1.^3);

  Hxx = hess(:,:,1); Hxy = hess(:,:,2); Hyy = hess(:,:,3);
  Txxx= third(:,:,1); Txxy=third(:,:,2); Txyy=third(:,:,3); Tyyy=third(:,:,4);

  firstbc = (Hxx.*(nxt.*nxt) + Hxy.*(2*nxt.*nyt) + Hyy.*(nyt.*nyt)) + ...
        nu.*(Hxx.*(taux.*taux) + Hxy.*(2*taux.*tauy) + Hyy.*(tauy.*tauy));

  secondbc = (Txxx.*(nxt.^3) + Txxy.*(3*nxt.^2.*nyt) + Txyy.*(3*nxt.*nyt.^2) + Tyyy.*(nyt.^3)) + ...
       (2-nu).*(Txxx.*(taux.^2.*nxt) + Txxy.*(taux.^2.*nyt + 2*taux.*tauy.*nxt) + ...
                Txyy.*(2*taux.*tauy.*nyt + tauy.^2.*nxt) + Tyyy.*(tauy.^2.*nyt)) + ...
       (1-nu).*kappat.*((Hxx.*taux.^2 + Hxy.*(2*taux.*tauy) + Hyy.*tauy.^2) - ...
                        (Hxx.*nxt.^2  + Hxy.*(2*nxt.*nyt)   + Hyy.*nyt.^2));

  submat = zeros(2*nt,ns);
  submat(1:2:end,:) = firstbc;
  submat(2:2:end,:) = secondbc;
end

if strcmpi(type,'vol')
  a = coefs(1);  b = coefs(2);
  [sig,~,~,~,fourth] = surfwave.gravity.gsgrav_der(g/2, 1, src, targ);
  Ssig = surfwave.gravity.gsgrav(g/2, 1/g, src, targ);
  bilap= fourth(:,:,1) + 2*fourth(:,:,3) + fourth(:,:,5);    % G_xxxx+2G_xxyy+G_yyyy
  submat = a*bilap + b*sig - g*Ssig;
end
