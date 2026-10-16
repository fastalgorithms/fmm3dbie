function varargout = fmm(eps, srcinfo, targinfo, type, sigma)
%LAP2D.FMM   Fast multipole evaluation of the 2D Laplace kernels.
%
% Syntax: [pot, grad, hess] = lap2d.fmm(eps, srcinfo, targinfo, type, sigma)
%
% Only the first two coordinates of the sources and targets are used.
% Requires the FMM2D MATLAB interface (rfmm2d).
%
% Input:
%   eps - precision requested
%   srcinfo - sources, srcinfo.r (2 or 3,ns)
%   targinfo - targets, targinfo.r (2 or 3,nt), or a (2 or 3,nt)
%                array; 'sprime' requires normals targinfo.n
%   type - string, kernel type as in LAP2D.KERN
%                type == 's', single layer S
%                type == 'sprime', normal derivative of S at the target
%   sigma - (ns,1) density, already multiplied by quadrature weights
%
% Output:
%   pot  - (nt,1) potential
%   grad - (2,nt) gradient at the targets, type 's' only
%   hess - (3,nt) Hessian (xx, xy, yy) at the targets, type 's' only
%
% see also LAP2D.KERN

if ( nargout == 0 )
    warning('LAP2D:fmm:empty', ...
        'Nothing to compute in LAP2D.FMM. Returning empty array.');
    return
end

srcuse = [];
srcuse.sources = srcinfo.r(1:2,:);
srcuse.charges = -1/(2*pi)*real(sigma(:).');

try
    targuse = targinfo.r(1:2,:);
catch
    targuse = targinfo(1:2,:);
end

switch lower(type)
    case {'s', 'single'}
        pgt = min(nargout, 3);
    case {'sp', 'sprime'}
        pgt = 2;
    otherwise
        error('LAP2D:fmm:type', 'Unknown kernel type ''%s''.', type);
end

% rfmm2d takes real charges; split a complex density
if ( isreal(sigma) )
    U = rfmm2d(eps, srcuse, 0, targuse, pgt);
else
    U  = rfmm2d(eps, srcuse, 0, targuse, pgt);
    srcuse.charges = -1/(2*pi)*imag(sigma(:).');
    Ui = rfmm2d(eps, srcuse, 0, targuse, pgt);
    U.pottarg = U.pottarg + 1i*Ui.pottarg;
    if ( pgt > 1 ), U.gradtarg = U.gradtarg + 1i*Ui.gradtarg; end
    if ( pgt > 2 ), U.hesstarg = U.hesstarg + 1i*Ui.hesstarg; end
end

switch lower(type)
    case {'s', 'single'}
        varargout{1} = U.pottarg(:);
    case {'sp', 'sprime'}
        if ( ~isfield(targinfo, 'n') && ~isprop(targinfo, 'n') )
            error('LAP2D:fmm:normals', ...
                'targinfo.n required for kernel type ''%s''.', type);
        end
        n = targinfo.n;
        varargout{1} = ( U.gradtarg(1,:).*n(1,:) + ...
                         U.gradtarg(2,:).*n(2,:) ).';
end

if ( nargout > 1 )
    switch lower(type)
        case {'s', 'single'}
            varargout{2} = U.gradtarg;
        otherwise
            error('LAP2D:fmm:grad', ...
                'Gradients not supported for kernel type ''%s''.', type);
    end
end

if ( nargout > 2 )
    switch lower(type)
        case {'s', 'single'}
            varargout{3} = U.hesstarg;
        otherwise
            error('LAP2D:fmm:hess', ...
                'Hessians not supported for kernel type ''%s''.', type);
    end
end

end
