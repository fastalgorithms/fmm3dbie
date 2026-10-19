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
%   targinfo - ptinfo struct: .r (3,:), .n (3,:) 
%   type     - 's'      single layer S
%              'd'      double layer D
%              'sp'/'sprime'  normal deriv of S at target, S'
%              'dp'/'dprime'  normal deriv of D at target, D'
%              'c'      combined layer coefs(1)*S + coefs(2)*D
%              'cp'/'cprime'  combined prime coefs(1)*S' + coefs(2)*D'
%              'sg', 'dg', 'cg'  target gradients (3x1) of 's', 'd', 'c'
%   sigma    - density (ns x ndens), already scaled by quadrature weights
%   coefs    - (optional) [alpha; beta] for combined layer
%
% Output:
%   pot  - potential (nt x ndens)
%   grad - gradient (3 x nt), or (ndens x 3 x nt) as in LFMM3D
%   hess - Hessian (6 x nt), ordered [xx yy zz xy xz yz] as in LFMM3D,
%          or (ndens x 6 x nt)

if ( nargout == 0 )
    warning('LAP3D:fmm:empty', ...
        'Nothing to compute in LAP3D.FMM. Returning empty array.');
    return
end

srcuse = [];
srcuse.sources = srcinfo.r;
ns = size(srcinfo.r, 2);

% the densities, (ndens, ns), and the dipoles x_j n_j, (ndens, 3, ns)
ndens = numel(sigma)/ns;
q = reshape(sigma, ns, ndens).';
nvec = @(x) reshape(x, ndens, 1, ns).*reshape(srcinfo.n(:,:), 1, 3, ns);

switch lower(type)
    case {'s', 'sprime', 'sp', 'sg'}
        srcuse.charges = q;
    case {'d', 'dprime', 'dp', 'dg'}
        srcuse.dipoles = nvec(q);
    case {'c', 'combined', 'cp', 'cprime', 'cg'}
        coefs = varargin{1};
        srcuse.charges = coefs(1) * q;
        srcuse.dipoles = coefs(2) * nvec(q);
    otherwise
        error('LAP3D:fmm:type', 'Unknown kernel type ''%s''.', type);
end
if ndens == 1 && isfield(srcuse, 'dipoles')
    srcuse.dipoles = reshape(srcuse.dipoles, 3, ns);
end
srcuse.nd = ndens;

try
    targuse = targinfo.r(:,:);
catch
    targuse = targinfo;
end
nt = size(targuse, 2);

pg  = 0;
pgt = min(nargout, 3);
switch lower(type)
    case {'sprime', 'sp', 'dprime', 'dp', 'cprime', 'cp', 'sg', 'dg', 'cg'}
        pgt = max(min(nargout + 1, 3), 2);
end

U = lfmm3d(eps, srcuse, pg, targuse, pgt);

% the potentials, (nt, ndens)
if ( nargout > 0 )
    switch lower(type)
        case {'s', 'd', 'c', 'combined'}
            varargout{1} = reshape(U.pottarg, ndens, nt).';
        case {'sprime', 'sp', 'dprime', 'dp', 'cprime', 'cp'}
            if ( ~isfield(targinfo, 'n') )
                error('LAP3D:fmm:normals', ...
                    'targinfo.n required for kernel type ''%s''.', type);
            end
            g = reshape(U.gradtarg, ndens, 3, nt);
            n = reshape(targinfo.n(:,:), 1, 3, nt);
            varargout{1} = reshape(sum(g.*n, 2), ndens, nt).';
        case {'sg', 'dg', 'cg'}
            g = reshape(U.gradtarg, ndens, 3*nt);
            varargout{1} = g.';
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
