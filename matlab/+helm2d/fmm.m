function varargout = fmm(eps, zk, srcinfo, targinfo, type, sigma)
%HELM2D.FMM   Fast multipole evaluation of the 2D Helmholtz kernels.
%
% Syntax: [pot, grad, hess] = helm2d.fmm(eps, zk, srcinfo, targinfo, type, sigma)
%
% Only the first two coordinates of the sources and targets are used.
% Requires the FMM2D MATLAB interface (hfmm2d).
%
% Input:
%   eps - precision requested
%   zk - complex number, Helmholtz wave number
%   srcinfo - sources, srcinfo.r (2 or 3,ns)
%   targinfo - targets, targinfo.r (2 or 3,nt), or a (2 or 3,nt)
%                array; 'sprime' requires normals targinfo.n
%   type - string, kernel type as in HELM2D.KERN
%                type == 's', single layer S
%                type == 'sprime', normal derivative of S at the target
%   sigma - (ns,1) density, already multiplied by quadrature weights
%
% Output:
%   pot  - (nt,1) potential
%   grad - (2,nt) gradient at the targets, type 's' only
%   hess - (3,nt) Hessian (xx, xy, yy) at the targets, type 's' only
%
% see also HELM2D.KERN

if ( nargout == 0 )
    warning('HELM2D:fmm:empty', ...
        'Nothing to compute in HELM2D.FMM. Returning empty array.');
    return
end

srcuse = [];
srcuse.sources = srcinfo.r(1:2,:);
srcuse.charges = sigma(:).';

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
        error('HELM2D:fmm:type', 'Unknown kernel type ''%s''.', type);
end

U = hfmm2d(eps, zk, srcuse, 0, targuse, pgt);

switch lower(type)
    case {'s', 'single'}
        varargout{1} = U.pottarg(:);
    case {'sp', 'sprime'}
        if ( ~isfield(targinfo, 'n') && ~isprop(targinfo, 'n') )
            error('HELM2D:fmm:normals', ...
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
            error('HELM2D:fmm:grad', ...
                'Gradients not supported for kernel type ''%s''.', type);
    end
end

if ( nargout > 2 )
    switch lower(type)
        case {'s', 'single'}
            varargout{3} = U.hesstarg;
        otherwise
            error('HELM2D:fmm:hess', ...
                'Hessians not supported for kernel type ''%s''.', type);
    end
end

end
