function submat= kern(zk,srcinfo,targinfo,type,varargin)
%FLEX2D.KERN free space flexural wave kernels in 2D
%
% Syntax: submat = flex2d.kern(zk,srcinfo,targinfo,type,varargin)
%
% Kernels based on the flexural wave Green's function
%         G(x,y) = 1/zk^2*(i/4 H_0^(1)(k|x-y|) - 1/(2 pi)*K_0(k|x-y|))
% with x the targets and y the sources.
%
% Input:
%   zk - wave number: scalar zk (flexural, or biharmonic if |zk| < 1e-6),
%        or [zk1, zk2]
%   srcinfo - sources in ptinfo struct format, srcinfo.r (2,:)
%   targinfo - targets in ptinfo struct format; the *_bcs kernels need
%        targinfo.n, and targinfo.d and targinfo.d2 where noted
%   type - string, determines kernel type
%        type == 's', the Green's function G
%        type == 'clamped_plate_bcs', clamped plate boundary conditions
%               applied to G. Requires n.
%        type == 'free_plate_bcs', free plate boundary conditions applied
%               to G. Requires n, d and d2.
%        type == 'supported_plate_bcs', supported plate boundary
%               conditions applied to G. Requires n and d.
%   varargin{1} - nu: Poisson's ratio, needed for the free and supported
%        plate kernels
%
% Output:
%   submat - kernel matrix; the *_bcs kernels return (2*nt,ns) with the
%            two conditions interleaved by row
%
% Reference: Nekrasov, P., Su, Z., Askham, T., & Hoskins, J. G. (2024).
% Boundary Integral Formulations for Flexural Wave Scattering in Thin
% Plates. arXiv:2409.19160.

src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

if isscalar(zk)
    if abs(zk) < 1e-6 % biharmonic (\Delta^2 u) 
        gfunc1 = @(s,t) surfwave.flex.bhgreen(s,t);
        gfunc2 = @(s,t) dumfunc(s,t,ns,nt);

        zk1 = 1;
        zk2 = 0;
    else % flexural (\Delta^2 u - k^4 u) 
        zk1 = zk;
        zk2 = 1i*zk;

        gfunc1 = @(s,t) surfwave.flex.helmdiffgreen(zk1,s,t);
        gfunc2 = @(s,t) surfwave.flex.helmdiffgreen(zk2,s,t);
    end
elseif length(zk) == 2
    if any(abs(zk) < 1e-6) % Stokes (\Delta^2 u + k^2 \Delta u) 
        zk1 = zk(abs(zk) >= 1e-6);
        zk2 = 0;

        gfunc1 = @(s,t) surfwave.flex.helmdiffgreen(zk1,s,t);
        gfunc2 = @(s,t) dumfunc(s,t,ns,nt);
    else % modified flexural (a \Delta^2 u - b \Delta u - c u) 
        zk1 = zk(1);
        zk2 = zk(2);
    
        gfunc1 = @(s,t) surfwave.flex.helmdiffgreen(zk1,s,t);
        gfunc2 = @(s,t) surfwave.flex.helmdiffgreen(zk2,s,t);
    end
end

switch lower(type)
case {'s', 'single'}

   val1 = gfunc1(src,targ);
   val2 = gfunc2(src,targ);

   val = 1/(zk1^2-zk2^2)*(val1-val2);

   submat = val;

case {'clamped_plate_bcs'}
    nxtarg = targinfo.n(1,:).'; 
    nytarg = targinfo.n(2,:).';  
    submat = zeros(2*nt,ns);
    
    [val1, grad1] = gfunc1(src,targ);
    [val2, grad2] = gfunc2(src,targ);

    val = 1/(zk1^2-zk2^2)*(val1-val2);
    grad = 1/(zk1^2-zk2^2)*(grad1-grad2);

    firstbc = val ;
    secondbc = grad(:, :, 1).*nxtarg + grad(:, :, 2).*nytarg;
   
    submat(1:2:end,:) = firstbc;
    submat(2:2:end,:) = secondbc;

case {'free_plate_bcs'}
    targnorm = targinfo.n;
    targtang = targinfo.d;
    targd2 = targinfo.d2;
    nu = varargin{1};

    [~, ~, hess1, third1] = gfunc1(src,targ);
    [~, ~, hess2, third2] = gfunc2(src,targ);

    hess = 1/(zk1^2-zk2^2)*(hess1 - hess2);
    third = 1/(zk1^2-zk2^2)*(third1 - third2);

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
    
    kappatarg = numer ./ denom; % target curvature
    
    firstbc = (hess(:, :, 1).*(nxtarg.*nxtarg) + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*(nytarg.*nytarg))+...
    nu.*(hess(:, :, 1).*(tauxtarg.*tauxtarg) + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*(tauytarg.*tauytarg));
    
    secondbc = (third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
    third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))+...
    (2-nu).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
    third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg+ tauytarg.*tauytarg.*nxtarg) +...
    + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg))+...
    (1-nu).*kappatarg.*((hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
    ((hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)));

    submat = zeros(2*nt,ns);
    submat(1:2:end,:) = firstbc;
    submat(2:2:end,:) = secondbc;

case {'supported_plate_bcs'}
    nxtarg = targinfo.n(1,:).'; 
    nytarg = targinfo.n(2,:).';  
    dx = targinfo.d(1,:).';
    dy = targinfo.d(2,:).';
    ds = sqrt(dx.*dx+dy.*dy);
    tauxtarg = (dx./ds);                                                                       % normalization
    tauytarg = (dy./ds);

    nu = varargin{1};

    [val1, ~, hess1] = gfunc1(src,targ);    
    [val2, ~, hess2] = gfunc2(src,targ);  

    val = 1/(zk1^2-zk2^2)*(val1 - val2);
    hess = 1/(zk1^2-zk2^2)*(hess1 - hess2);
    
    firstbc = val ;
    
    secondbc = (hess(:, :, 1).*(nxtarg.*nxtarg) + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*(nytarg.*nytarg))+...
               nu.*(hess(:, :, 1).*(tauxtarg.*tauxtarg) + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*(tauytarg.*tauytarg));
    
    submat = zeros(2*nt,ns);

    submat(1:2:end,:) = firstbc;
    submat(2:2:end,:) = secondbc;

otherwise
    error('FLEX2D.KERN: unknown kernel type ''%s''.', type);
end

function varargout = dumfunc(src,targ,ns,nt)

    out1 = zeros(nt,ns,1);
    out2 = zeros(nt,ns,2);
    out3 = zeros(nt,ns,3);
    out4 = zeros(nt,ns,4);
    out5 = zeros(nt,ns,5);
    out6 = zeros(nt,ns,6);

    switch nargout
        case 0
            return
        case 1 
            varargout = {out1};
        case 2
            varargout = {out1, out2};
        case 3
            varargout = {out1, out2, out3};
        case 4
            varargout = {out1, out2, out3, out4};
        case 5
            varargout = {out1, out2, out3, out4, out5};
        case 6
            varargout = {out1, out2, out3, out4, out5, out6};
    end
end

end

