function varargout = fmm(eps, zk, srcinfo, targinfo, type, sigma, nu)
%FLEX2D.FMM   Fast multipole evaluation of the 2D flexural kernels.
%
% Syntax: pot = flex2d.fmm(eps, zk, srcinfo, targinfo, type, sigma, nu)
%
% Uses G = (G_{zk1} - G_{zk2})/(zk1^2 - zk2^2), with G_k the Helmholtz
% Green's function (the Laplace one if zk2 = 0), so each call is two
% FMM2D calls (hfmm2d, rfmm2d). The biharmonic case and the free plate
% kernel (third derivatives) are not supported.
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
    if ( isempty(zk1) )
        error('FLEX2D:fmm:biharmonic', ...
            'The biharmonic kernel is not supported by FLEX2D.FMM.');
    end
else
    zk1 = zk(1);
    zk2 = zk(2);
end

switch lower(type)
    case {'s', 'single'}
        pgt = 1;
    case {'clamped_plate_bcs'}
        pgt = 2;
    case {'supported_plate_bcs'}
        pgt = 3;
    otherwise
        error('FLEX2D:fmm:type', 'Unsupported kernel type ''%s''.', type);
end

try
    targuse = targinfo.r(1:2,:);
catch
    targuse = targinfo(1:2,:);
end

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

switch lower(type)
    case {'s', 'single'}
        varargout{1} = val(:);
        return

    case {'clamped_plate_bcs'}
        grad = U.gradtarg - U2.gradtarg;
        nx = targinfo.n(1,:); ny = targinfo.n(2,:);
        row2 = grad(1,:).*nx + grad(2,:).*ny;

    case {'supported_plate_bcs'}
        hess = U.hesstarg - U2.hesstarg;
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
