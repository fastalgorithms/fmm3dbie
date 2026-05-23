function varargout = fmm(eps, zk, srcinfo, targinfo, type, sigma, varargin)
%HELM3D.FMM   Fast multipole method for evaluating Helmholtz layer
%   potentials and their gradients in 3D.
%
% Syntax:
%   pot         = helm3d.fmm(eps, zk, srcinfo, targinfo, type, sigma)
%   [pot, grad] = helm3d.fmm(eps, zk, srcinfo, targinfo, type, sigma)
%   ...         = helm3d.fmm(eps, zk, srcinfo, targinfo, 'c', sigma, coefs)
%
% Let x be targets and y be sources, with n_x and n_y the corresponding
% unit normals.  Kernels based on G(x,y) = exp(ik|x-y|)/(4*pi*|x-y|).
%
% Input:
%   eps      - precision requested
%   zk       - Helmholtz wavenumber (complex)
%   srcinfo  - ptinfo struct: .r (3,:), .n (3,:)
%   targinfo - ptinfo struct: .r (3,:), .n (3,:) (normals needed for sprime)
%   type     - 's'           single layer S
%              'd'           double layer D
%              'sp'/'sprime' normal deriv of S at target
%              'c'           combined layer coefs(1)*S + coefs(2)*D
%   sigma    - density (ns x 1), already scaled by quadrature weights
%   varargin{1} - coefs [alpha; beta] for combined layer (type 'c' only)
%
% Output:
%   pot  - potential (nt x 1)
%   grad - gradient (3 x nt)
%
% Note: Helmholtz FMM does not support Hessians.

if ( nargout == 0 )
    warning('HELM3D:fmm:empty', ...
        'Nothing to compute in HELM3D.FMM. Returning empty array.');
    return
end

srcuse = [];
srcuse.sources = srcinfo.r;

switch lower(type)
    case {'s', 'sprime', 'sp'}
        srcuse.charges = sigma(:).';
    case {'d'}
        srcuse.dipoles = sigma(:).'.*srcinfo.n;
    case {'c', 'combined'}
        coefs = varargin{1};
        srcuse.charges = coefs(1) * sigma(:).';
        srcuse.dipoles = coefs(2) * sigma(:).'.*srcinfo.n;
    otherwise
        error('HELM3D:fmm:type', 'Unknown kernel type ''%s''.', type);
end

if ( isstruct(targinfo) )
    targuse = targinfo.r;
else
    targuse = targinfo;
end

pg  = 0;
pgt = min(nargout, 2);   % hfmm3d supports pgt = 1 (pot) or 2 (pot+grad)
switch lower(type)
    case {'sprime', 'sp'}
        pgt = max(min(nargout + 1, 2), 2);
end

U = hfmm3d(eps, zk, srcuse, pg, targuse, pgt);

% Potentials
if ( nargout > 0 )
    switch lower(type)
        case {'s', 'd', 'c', 'combined'}
            varargout{1} = U.pottarg(:);
        case {'sprime', 'sp'}
            if ( ~isfield(targinfo, 'n') )
                error('HELM3D:fmm:normals', ...
                    'targinfo.n required for kernel type ''%s''.', type);
            end
            n = targinfo.n;
            U.gradtarg = reshape(U.gradtarg, 3, []);
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
            error('HELM3D:fmm:grad', ...
                'Gradients not supported for kernel type ''%s''.', type);
    end
end

end
