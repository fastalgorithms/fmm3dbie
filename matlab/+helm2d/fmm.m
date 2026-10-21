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
%                type == 'sg', gradient of S at the target, (2*nt,1)
%                    with d/dx and d/dy interleaved by target
%   sigma - (ns,ndens) densities, already multiplied by quadrature
%           weights
%
% Output:
%   pot  - (nt,ndens) potential, (2*nt,ndens) for 'sg'
%   grad - (2,nt) gradient at the targets, type 's' only, or
%          (ndens,2,nt)
%   hess - (3,nt) Hessian (xx, xy, yy) at the targets, type 's' only, or
%          (ndens,3,nt)
%
% see also HELM2D.KERN

if ( nargout == 0 )
    warning('HELM2D:fmm:empty', ...
        'Nothing to compute in HELM2D.FMM. Returning empty array.');
    return
end

ns = size(srcinfo.r(:,:), 2);
ndens = numel(sigma)/ns;
if ( ndens > 1 )
    sigma = reshape(sigma, ns, ndens);
    out = cell(1, max(nargout, 1));
    varargout = cell(1, numel(out));
    for j = 1:ndens
        [out{:}] = helm2d.fmm(eps, zk, srcinfo, targinfo, type, sigma(:,j));
        if ( j == 1 )
            varargout{1} = zeros(numel(out{1}), ndens, 'like', out{1});
            for k = 2:numel(out)
                varargout{k} = zeros([ndens, size(out{k})], 'like', out{k});
            end
        end
        varargout{1}(:,j) = out{1}(:);
        for k = 2:numel(out)
            varargout{k}(j,:,:) = reshape(out{k}, [1, size(out{k})]);
        end
    end
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
    case {'sp', 'sprime', 'sg', 'sgrad'}
        pgt = 2;
    otherwise
        error('HELM2D:fmm:type', 'Unknown kernel type ''%s''.', type);
end

U = hfmm2d(eps, zk, srcuse, 0, targuse, pgt);

switch lower(type)
    case {'s', 'single'}
        varargout{1} = U.pottarg(:);
    case {'sg', 'sgrad'}
        varargout{1} = U.gradtarg(:);
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
