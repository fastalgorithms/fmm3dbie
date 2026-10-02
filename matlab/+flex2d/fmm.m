function varargout = fmm(eps, zk, srcinfo, targinfo, type, sigma, nu, pfun)
%FLEX2D.FMM   Fast multipole evaluation of the 2D flexural kernels.
%
% Syntax: pot = flex2d.fmm(eps, zk, srcinfo, targinfo, type, sigma, nu)
%
% Uses G = (G_{zk1} - G_{zk2})/(zk1^2 - zk2^2), with G_k the Helmholtz
% Green's function (the Laplace one if zk2 = 0), so each call is two
% FMM2D calls (hfmm2d, rfmm2d). The biharmonic case (zk = 0) uses
%
%   |x-y|^2 log|x-y| = (|x|^2 - 2 x.y + |y|^2) log|x-y|,
%
% i.e. four Laplace FMMs with charges sigma, sigma*y1, sigma*y2 and
% sigma*|y|^2, done in one rfmm2d call, with the origin shifted to the
% centroid of the sources to limit cancellation. The free plate
% boundary condition (third derivatives) is not supported.
%
% Input:
%   eps - precision requested
%   zk - wavenumbers [zk1, zk2], as in FLEX2D.KERN
%   srcinfo - sources, srcinfo.r (2 or 3,ns)
%   targinfo - targets, targinfo.r (2 or 3,nt); the plate bcs kernels
%                need targinfo.n, and 'supported_plate_bcs' targinfo.d
%   type - string, kernel type as in FLEX2D.KERN
%                type == 's'
%                type == 'clamped_plate_bcs'
%                type == 'supported_plate_bcs'
%                type == 'varcoef', needs nu and pfun, see FLEX2D.KERN
%   sigma - (ns,1) density, already multiplied by quadrature weights
%   nu - Poisson's ratio, needed for 'supported_plate_bcs'
%
% Output:
%   pot - (nt,1) potential, or (2*nt,1) with the two conditions
%         interleaved for the plate bcs kernels
%
% see also FLEX2D.KERN, FLEX2D.GREEN

if ( nargout == 0 )
    warning('FLEX2D:fmm:empty', ...
        'Nothing to compute in FLEX2D.FMM. Returning empty array.');
    return
end

zk = zk(:).';
if ( isscalar(zk) ), zk = [zk, 1i*zk]; end
if ( any(abs(zk) < 1e-6) )
    zk1 = zk(abs(zk) >= 1e-6);
    zk2 = 0;
    isbh = isempty(zk1);
else
    zk1 = zk(1);
    zk2 = zk(2);
    isbh = false;
end

switch lower(type)
    case {'s', 'single'}
        pgt = 1;
    case {'clamped_plate_bcs'}
        pgt = 2;
    case {'supported_plate_bcs'}
        pgt = 3;
    case {'varcoef'}
        pgt = 3;
    otherwise
        error('FLEX2D:fmm:type', 'Unsupported kernel type ''%s''.', type);
end

try
    targuse = targinfo.r(1:2,:);
catch
    targuse = targinfo(1:2,:);
end

grad = []; hess = [];
if ( isbh )
    [val, grad, hess] = bh_fmm(eps, srcinfo.r(1:2,:), targuse, sigma, pgt);
else
    srcuse = [];
    srcuse.sources = srcinfo.r(1:2,:);
    srcuse.charges = sigma(:).'/(zk1^2 - zk2^2);

    % G_{zk1} part
    U = hfmm2d(eps, zk1, srcuse, 0, targuse, pgt);

    % minus G_{zk2} part
    if ( zk2 ~= 0 )
        U2 = hfmm2d(eps, zk2, srcuse, 0, targuse, pgt);
    else
        % Laplace, G_0 = -log|x-y|/(2 pi), real charges only
        c = srcuse.charges;
        srcuse.charges = -real(c)/(2*pi);
        U2 = rfmm2d(eps, srcuse, 0, targuse, pgt);
        srcuse.charges = -imag(c)/(2*pi);
        U2i = rfmm2d(eps, srcuse, 0, targuse, pgt);
        U2.pottarg = U2.pottarg + 1i*U2i.pottarg;
        if ( pgt > 1 ), U2.gradtarg = U2.gradtarg + 1i*U2i.gradtarg; end
        if ( pgt > 2 ), U2.hesstarg = U2.hesstarg + 1i*U2i.hesstarg; end
    end

    val = U.pottarg - U2.pottarg;
    if ( pgt > 1 ), grad = U.gradtarg - U2.gradtarg; end
    if ( pgt > 2 ), hess = U.hesstarg - U2.hesstarg; end

    % grad Lap G, using Lap G_k = -k^2 G_k away from the source
    if ( strcmpi(type, 'varcoef') )
        gradlap = -zk1^2*U.gradtarg + zk2^2*U2.gradtarg;
    end
end

if ( strcmpi(type, 'varcoef') )
    if ( isbh )
        % Lap G = (log r + 1)/(2 pi), so grad Lap G = grad log r/(2 pi)
        srcuse = [];
        srcuse.sources = srcinfo.r(1:2,:);
        srcuse.nd = 2;
        srcuse.charges = [real(sigma(:).'); imag(sigma(:).')]/(2*pi);
        U = rfmm2d(eps, srcuse, 0, targuse, 2);
        nt = size(targuse, 2);
        gl = reshape(U.gradtarg, 2, 2, nt);
        gradlap = reshape(gl(1,:,:) + 1i*gl(2,:,:), 2, nt);
    end
    t = targinfo;
    if ( isnumeric(t) ), t = []; t.r = targinfo; end
    cf = flex2d.plate_coefs(pfun(t), nu, zk);
    pot = cf(1,:).*gradlap(1,:) + cf(2,:).*gradlap(2,:) + ...
        cf(3,:).*(hess(1,:) + hess(3,:)) + cf(4,:).*hess(3,:) + ...
        cf(5,:).*hess(1,:) + cf(6,:).*hess(2,:) + cf(7,:).*val(:).';
    varargout{1} = pot(:);
    return
end

switch lower(type)
    case {'s', 'single'}
        varargout{1} = val(:);
        return

    case {'clamped_plate_bcs'}
        nx = targinfo.n(1,:); ny = targinfo.n(2,:);
        row2 = grad(1,:).*nx + grad(2,:).*ny;

    case {'supported_plate_bcs'}
        nx = targinfo.n(1,:); ny = targinfo.n(2,:);
        ds = sqrt(targinfo.d(1,:).^2 + targinfo.d(2,:).^2);
        tx = targinfo.d(1,:)./ds; ty = targinfo.d(2,:)./ds;
        row2 = hess(1,:).*nx.*nx + 2*hess(2,:).*nx.*ny + hess(3,:).*ny.*ny + ...
            nu*(hess(1,:).*tx.*tx + 2*hess(2,:).*tx.*ty + hess(3,:).*ty.*ty);
end

pot = zeros(2*numel(val), 1, 'like', val);
pot(1:2:end) = val;
pot(2:2:end) = row2;
varargout{1} = pot;

end


function [val, grad, hess] = bh_fmm(eps, src, targ, sigma, pgt)
%BH_FMM  sum_j sigma_j |x-y_j|^2 log|x-y_j| / (8 pi) and its gradient and
% Hessian (xx, xy, yy) at the targets, from Laplace FMMs.
%
% With A, B = (B1,B2), C the log potentials of the charges sigma,
% sigma*y, sigma*|y|^2,
%
%   8 pi u = |x|^2 A - 2 x.B + C

nt = size(targ, 2);

% shift the origin to limit cancellation
cen = mean(src, 2);
y = src - cen;
x = targ - cen;
x1 = x(1,:); x2 = x(2,:); x2n = x1.^2 + x2.^2;

sig = sigma(:).';
q = [sig; sig.*y(1,:); sig.*y(2,:); sig.*(y(1,:).^2 + y(2,:).^2)];

% rfmm2d takes real charges, kernel log|x-y|
srcuse = [];
srcuse.sources = src;
srcuse.nd = 8;
srcuse.charges = [real(q); imag(q)];
U = rfmm2d(eps, srcuse, 0, targ, pgt);

P = reshape(U.pottarg, 8, nt);
P = P(1:4,:) + 1i*P(5:8,:);
A = P(1,:); B1 = P(2,:); B2 = P(3,:); C = P(4,:);
val = (x2n.*A - 2*(x1.*B1 + x2.*B2) + C)/(8*pi);

grad = []; hess = [];
if ( pgt > 1 )
    G = reshape(U.gradtarg, 8, 2, nt);
    G = G(1:4,:,:) + 1i*G(5:8,:,:);
    Ax  = reshape(G(1,1,:),1,nt); Ay  = reshape(G(1,2,:),1,nt);
    B1x = reshape(G(2,1,:),1,nt); B1y = reshape(G(2,2,:),1,nt);
    B2x = reshape(G(3,1,:),1,nt); B2y = reshape(G(3,2,:),1,nt);
    Cx  = reshape(G(4,1,:),1,nt); Cy  = reshape(G(4,2,:),1,nt);
    grad = zeros(2, nt, 'like', val);
    grad(1,:) = 2*x1.*A + x2n.*Ax - 2*B1 - 2*(x1.*B1x + x2.*B2x) + Cx;
    grad(2,:) = 2*x2.*A + x2n.*Ay - 2*B2 - 2*(x1.*B1y + x2.*B2y) + Cy;
    grad = grad/(8*pi);
end
if ( pgt > 2 )
    H = reshape(U.hesstarg, 8, 3, nt);
    H = H(1:4,:,:) + 1i*H(5:8,:,:);
    h = @(k,l) reshape(H(k,l,:),1,nt);
    hess = zeros(3, nt, 'like', val);
    hess(1,:) = 2*A + 4*x1.*Ax + x2n.*h(1,1) - 4*B1x ...
        - 2*(x1.*h(2,1) + x2.*h(3,1)) + h(4,1);
    hess(2,:) = 2*x1.*Ay + 2*x2.*Ax + x2n.*h(1,2) - 2*B1y - 2*B2x ...
        - 2*(x1.*h(2,2) + x2.*h(3,2)) + h(4,2);
    hess(3,:) = 2*A + 4*x2.*Ay + x2n.*h(1,3) - 4*B2y ...
        - 2*(x1.*h(2,3) + x2.*h(3,3)) + h(4,3);
    hess = hess/(8*pi);
end

end
