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
%              'trep_prime'  [coefs(1) D', coefs(2) S'] on [rho; sigma]
%              'sg', 'dg', 'cg', 'trep_grad'  the target gradients of 's',
%                            'd', 'c' and 'trans_rep', (3x1) and (3x2)
%              'trans_sys'   [c11 D, c12 S; c21 D', c22 S'], coefs 2x2
%   sigma    - density (n*ns x ndens) for n components per source,
%              already scaled by quadrature weights
%   coefs    - [alpha; beta] for 'c', 'cp', 'c2trans'; see above for
%              'trans_rep', 'trep_prime' and 'trans_sys'
%
% Output:
%   pot  - potential (m*nt x ndens) for m rows of the kernel
%   grad - gradient (3 x nt), or (ndens x 3 x nt)
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

% the densities, (n, ns, ndens) for n components per source
n = 1;
if any(strcmp(type, {'trans_rep', 'trep_prime', 'trep_grad', 'trans_sys'})), n = 2; end
ndens = numel(sigma)/(n*ns);
sg = reshape(sigma, n, ns, ndens);
q = @(k) reshape(sg(k,:,:), ns, ndens).';                  % (ndens, ns)
nvec = @(x) reshape(x, [], 1, ns).*reshape(srcinfo.n(:,:), 1, 3, ns);

% sources: charges and/or dipoles, and which targets outputs are needed
needgrad = any(strcmp(type, {'sp', 'sprime', 'dp', 'dprime', 'cp', ...
    'cprime', 's2trans', 'd2trans', 'c2trans', 'trans_sys', 'trep_prime', ...
    'sg', 'dg', 'cg', 'trep_grad'})) ...
    || nargout > 1;
switch type
    case {'s', 'sp', 'sprime', 's2trans', 'sg'}
        srcuse.charges = q(1);
    case {'d', 'dp', 'dprime', 'd2trans', 'dg'}
        srcuse.dipoles = nvec(q(1));
    case {'c', 'combined', 'cp', 'cprime', 'c2trans', 'cg'}
        srcuse.charges = coefs(1)*q(1);
        srcuse.dipoles = coefs(2)*nvec(q(1));
    case {'trans_rep', 'trep_prime', 'trep_grad'}
        % [coefs(1) D, coefs(2) S] on [rho; sigma] interleaved by source
        srcuse.charges = coefs(2)*q(2);
        srcuse.dipoles = coefs(1)*nvec(q(1));
    case 'trans_sys'
        % rows [c11 D + c12 S; c21 D' + c22 S'] on [rho; sigma]: two
        % density sets per density in one call, the first for its potential
        % and the second for its normal derivative
        srcuse.charges = [coefs(1,2)*q(2); coefs(2,2)*q(2)];
        srcuse.dipoles = [coefs(1,1)*nvec(q(1)); coefs(2,1)*nvec(q(1))];
    otherwise
        error('HELM3D:fmm:type', 'Unknown kernel type ''%s''.', type);
end
nd = ndens;
if strcmp(type, 'trans_sys'), nd = 2*ndens; end
if nd == 1 && isfield(srcuse, 'dipoles')
    srcuse.dipoles = reshape(srcuse.dipoles, 3, ns);
end
srcuse.nd = nd;

try
    targuse = targinfo.r(:,:);
catch
    targuse = targinfo;
end
nt = size(targuse, 2);

pgt = 1 + needgrad;
U = hfmm3d(eps, zk, srcuse, 0, targuse, pgt);

pot = reshape(U.pottarg, nd, nt);
if needgrad
    grad = reshape(U.gradtarg, nd, 3, nt);
end
dn = @(k) ndot(grad(k,:,:), targinfo, type);         % n . grad, (numel(k), nt)

% the potentials, (m*nt, ndens) for m rows of the kernel
d = 1:ndens;
switch type
    case {'s', 'd', 'c', 'combined', 'trans_rep'}
        out = pot;
    case {'sp', 'sprime', 'dp', 'dprime', 'cp', 'cprime', 'trep_prime'}
        out = dn(d);
    case {'s2trans', 'd2trans', 'c2trans'}
        out = [pot; dn(d)];
    case 'trans_sys'
        out = [pot(d,:); dn(ndens+d)];
    case {'sg', 'dg', 'cg', 'trep_grad'}
        out = reshape(grad, 3*ndens, nt);
end
m = size(out, 1)/ndens;
varargout{1} = reshape(permute(reshape(out, ndens, m, nt), [2 3 1]), m*nt, ndens);

% Gradients
if ( nargout > 1 )
    switch type
        case {'s', 'd', 'c', 'combined'}
            varargout{2} = reshape(grad, 3*ndens, nt);
            if ndens > 1, varargout{2} = grad; end
        otherwise
            error('HELM3D:fmm:grad', ...
                'Gradients not supported for kernel type ''%s''.', type);
    end
end

end


function v = ndot(g, targinfo, type)
% normal derivative from the gradients (k, 3, nt) and targinfo.n
if ~isfield(targinfo, 'n') && ~isprop(targinfo, 'n')
    error('HELM3D:fmm:normals', ...
        'targinfo.n required for kernel type ''%s''.', type);
end
n = reshape(targinfo.n(:,:), 1, 3, []);
v = reshape(sum(g.*n, 2), size(g, 1), []);
end
