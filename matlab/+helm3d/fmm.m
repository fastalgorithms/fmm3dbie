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
%   targinfo - ptinfo struct: .r (3,:), .n (3,:)
%   type     - 's'           single layer S
%              'd'           double layer D
%              'sp'/'sprime' normal deriv of S at target
%              'c'           combined layer coefs(1)*S + coefs(2)*D
%              'dp'/'dprime' normal deriv of D at target
%              'cp'/'cprime' coefs(1)*S' + coefs(2)*D'
%              's2trans'     [S; S'],  'd2trans' [D; D'],
%              'c2trans'     [C; C'] with coefs [alpha; beta]
%              'trans_rep'   [coefs(1) D, coefs(2) S] on [rho; sigma]
%              'trans_sys'   [c11 D, c12 S; c21 D', c22 S'], coefs 2x2
%   sigma    - density (ns x 1), already scaled by quadrature weights
%   coefs    - [alpha; beta] for 'c', 'cp', 'c2trans'; see above for
%              'trans_rep' and 'trans_sys'
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
ns = size(srcinfo.r, 2);
type = lower(type);
coefs = [];
if ~isempty(varargin), coefs = varargin{1}; end

% sources: charges and/or dipoles, and which targets outputs are needed
needgrad = any(strcmp(type, {'sp', 'sprime', 'dp', 'dprime', 'cp', ...
    'cprime', 's2trans', 'd2trans', 'c2trans', 'trans_sys'})) || nargout > 1;
nvec = @(x) x(:).'.*srcinfo.n;              % dipoles x_j n_j
switch type
    case {'s', 'sp', 'sprime', 's2trans'}
        srcuse.charges = sigma(:).';
    case {'d', 'dp', 'dprime', 'd2trans'}
        srcuse.dipoles = nvec(sigma);
    case {'c', 'combined', 'cp', 'cprime', 'c2trans'}
        srcuse.charges = coefs(1)*sigma(:).';
        srcuse.dipoles = coefs(2)*nvec(sigma);
    case 'trans_rep'
        % [coefs(1) D, coefs(2) S] on [rho; sigma] interleaved by source
        sg = reshape(sigma, 2, ns);
        srcuse.charges = coefs(2)*sg(2,:);
        srcuse.dipoles = coefs(1)*nvec(sg(1,:));
    case 'trans_sys'
        % rows [c11 D + c12 S; c21 D' + c22 S'] on [rho; sigma]: two
        % density sets in one call, the first for its potential and the
        % second for its normal derivative
        sg = reshape(sigma, 2, ns);
        srcuse.nd = 2;
        srcuse.charges = [coefs(1,2)*sg(2,:); coefs(2,2)*sg(2,:)];
        dip = zeros(2, 3, ns);
        dip(1,:,:) = reshape(coefs(1,1)*nvec(sg(1,:)), 1, 3, ns);
        dip(2,:,:) = reshape(coefs(2,1)*nvec(sg(1,:)), 1, 3, ns);
        srcuse.dipoles = dip;
    otherwise
        error('HELM3D:fmm:type', 'Unknown kernel type ''%s''.', type);
end

try
    targuse = targinfo.r(:,:);
catch
    targuse = targinfo;
end
nt = size(targuse, 2);

pgt = 1 + needgrad;
U = hfmm3d(eps, zk, srcuse, 0, targuse, pgt);

nd = 1;  if isfield(srcuse, 'nd'), nd = srcuse.nd; end
pot = reshape(U.pottarg, nd, nt);
if needgrad
    grad = reshape(U.gradtarg, nd, 3, nt);
end
dn = @(k) ndot(squeeze(grad(k,:,:)), targinfo, type);   % n . grad

switch type
    case {'s', 'd', 'c', 'combined', 'trans_rep'}
        out = pot;
    case {'sp', 'sprime', 'dp', 'dprime', 'cp', 'cprime'}
        out = dn(1);
    case {'s2trans', 'd2trans', 'c2trans'}
        out = [pot; dn(1)];
    case 'trans_sys'
        out = [pot(1,:); dn(2)];
end
varargout{1} = out(:);

% Gradients
if ( nargout > 1 )
    switch type
        case {'s', 'd', 'c', 'combined'}
            varargout{2} = reshape(grad, 3, nt);
        otherwise
            error('HELM3D:fmm:grad', ...
                'Gradients not supported for kernel type ''%s''.', type);
    end
end

end


function v = ndot(g, targinfo, type)
% normal derivative from the gradient (3, nt) and targinfo.n
if ~isfield(targinfo, 'n') && ~isprop(targinfo, 'n')
    error('HELM3D:fmm:normals', ...
        'targinfo.n required for kernel type ''%s''.', type);
end
n = targinfo.n(:,:);
g = reshape(g, 3, []);
v = g(1,:).*n(1,:) + g(2,:).*n(2,:) + g(3,:).*n(3,:);
end
