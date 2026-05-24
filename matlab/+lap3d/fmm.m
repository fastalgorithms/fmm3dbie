function varargout = fmm(eps, srcinfo, targinfo, type, sigma, varargin)
%LAP3D.FMM   Fast multipole methods for evaluating Laplace layer
%   potentials, their gradients, and Hessians in 3D.
%
% Syntax:
%   pot          = lap3d.fmm(eps, srcinfo, targinfo, type, sigma)
%   [pot, grad]  = lap3d.fmm(eps, srcinfo, targinfo, type, sigma)
%   [pot, grad, hess] = lap3d.fmm(eps, srcinfo, targinfo, type, sigma)
%   ...          = lap3d.fmm(eps, srcinfo, targinfo, 'c', sigma, coefs)
%
% Let x be targets and y be sources, with n_x and n_y the corresponding
% unit normals.  Kernels based on G(x,y) = 1/(4*pi*|x-y|).
%
% Input:
%   eps      - precision requested
%   srcinfo  - ptinfo struct: .r (3,:), .n (3,:)
%   targinfo - ptinfo struct: .r (3,:), .n (3,:) (normals needed for sprime)
%   type     - 's'      single layer S
%              'd'      double layer D
%              'sp'/'sprime'  normal deriv of S at target, S'
%              'dp'/'dprime'  normal deriv of D at target, D'
%              'c'      combined layer coefs(1)*S + coefs(2)*D
%              'cp'/'cprime'  combined prime coefs(1)*S' + coefs(2)*D'
%   sigma    - density (ns x 1), already scaled by quadrature weights
%   varargin{1} - coefs [alpha; beta] for combined layer
%
% Output:
%   pot  - potential (nt x 1)
%   grad - gradient (3 x nt)
%   hess - Hessian (6 x nt), stored as [xx xy xz yy yz zz]

if ( nargout == 0 )
    warning('LAP3D:fmm:empty', ...
        'Nothing to compute in LAP3D.FMM. Returning empty array.');
    return
end

srcuse = [];
srcuse.sources = srcinfo.r;

switch lower(type)
    case {'s', 'sprime', 'sp'}
        srcuse.charges = sigma(:).';
    case {'d', 'dprime', 'dp'}
        srcuse.dipoles = sigma(:).'.*srcinfo.n;
    case {'c', 'combined'}
        coefs = varargin{1};
        srcuse.charges = coefs(1) * sigma(:).';
        srcuse.dipoles = coefs(2) * sigma(:).'.*srcinfo.n;
    case {'cp', 'cprime'}
        coefs = varargin{1};
        srcuse.charges = coefs(1) * sigma(:).';
        srcuse.dipoles = coefs(2) * sigma(:).'.*srcinfo.n;
    otherwise
        error('LAP3D:fmm:type', 'Unknown kernel type ''%s''.', type);
end

try
    targuse = targinfo.r(:,:);
catch
    targuse = targinfo;
end

pg  = 0;
pgt = min(nargout, 3);
switch lower(type)
    case {'sprime', 'sp', 'dprime', 'dp', 'cprime', 'cp'}
        pgt = max(min(nargout + 1, 3), 2);
end

U = lfmm3d(eps, srcuse, pg, targuse, pgt);

% Potentials
if ( nargout > 0 )
    switch lower(type)
        case {'s', 'd', 'c', 'combined'}
            varargout{1} = U.pottarg(:);
        case {'sprime', 'sp', 'dprime', 'dp', 'cprime', 'cp'}
            if ( ~isfield(targinfo, 'n') )
                error('LAP3D:fmm:normals', ...
                    'targinfo.n required for kernel type ''%s''.', type);
            end
            n = targinfo.n;
            U.gradtarg = reshape(U.gradtarg,3,[]);
            varargout{1} = ( U.gradtarg(1,:).*n(1,:) + ...
                             U.gradtarg(2,:).*n(2,:) + ...
                             U.gradtarg(3,:).*n(3,:) ).';
    end
end

% Gradients
if ( nargout > 1 )
    switch lower(type)
        case {'s', 'd', 'c', 'combined'}
            varargout{2} = U.gradtarg;
        otherwise
            error('LAP3D:fmm:grad', ...
                'Gradients not supported for kernel type ''%s''.', type);
    end
end

% Hessians
if ( nargout > 2 )
    switch lower(type)
        case {'s', 'd', 'c', 'combined'}
            varargout{3} = U.hesstarg;
        otherwise
            error('LAP3D:fmm:hess', ...
                'Hessians not supported for kernel type ''%s''.', type);
    end
end

end
