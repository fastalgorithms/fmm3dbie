function pot = fmm(eps, zk, srcinfo, targinfo, type, sigma, varargin)
%EM3D.FMM   Fast multipole method for Maxwell (em3d) layer potentials in 3D.
%
% Syntax:
%   pot = em3d.fmm(eps, zk, srcinfo, targinfo, type, sigma)
%   pot = em3d.fmm(eps, zk, srcinfo, targinfo, 'nrccie-bc', sigma, alpha)
%
% Evaluates the far-field contribution of the Maxwell vector-Helmholtz
% representation using emfmm3d.  The representation is
%
%   E = ik S_k[J]  +  curl S_k[J]  -  grad S_k[rho]
%   H = curl S_k[J]
%
% where the density components depend on the kernel type (see below).
%
% Input:
%   eps      - precision requested
%   zk       - Maxwell wavenumber (complex)
%   srcinfo  - ptinfo struct with fields:
%                .r  (3, ns)  source positions
%                .n  (3, ns)  unit outward normals
%                .du (3, ns)  first tangent vector (dxyz/du)
%                .dv (3, ns)  second tangent vector (dxyz/dv)
%   targinfo - ptinfo struct with fields:
%                .r  (3, nt)  target positions
%                .n  (3, nt)  unit outward normals   [needed for nrccie-bc]
%                .du (3, nt)  first tangent vector   [needed for nrccie-bc]
%                .dv (3, nt)  second tangent vector  [needed for nrccie-bc]
%   type     - kernel type string:
%                'nrccie-bc'    Non-resonant CFIE boundary-condition operator.
%                               sigma (3*ns x 1): [j_u(1:ns); j_v(1:ns); rho(1:ns)]
%                               where j_u, j_v are the tangential components of J
%                               in the local orthonormal frame (ru, rv), and rho
%                               is the surface charge.
%                               pot (3*nt x 1): [M_u(1:nt); M_v(1:nt); E_n(1:nt)]
%                               i.e. the tangential components of H and the normal
%                               component of E at the target points.
%
%                'nrccie-eval'  NRCCIE field evaluation.
%                               sigma (4*ns x 1): [J_x; J_y; J_z; rho] Cartesian
%                               pot   (6*nt x 1): [E_x; E_y; E_z; H_x; H_y; H_z]
%
%   sigma    - density, flat column vector (already scaled by quadrature weights)
%   varargin{1} - alpha: CFIE regularisation parameter (only for 'nrccie-bc')
%
% Output:
%   pot - potential column vector:
%         nrccie-bc:   (3*nt x 1)
%         nrccie-eval: (6*nt x 1)
%
% Note: This routine returns the FMM (far-field) contribution only.
%       Near-field corrections must be added separately.
%
% See also: EM3D.KERN, KERNEL3D.EM3D, EMFMM3D

% -------------------------------------------------------------------------
% Extract geometry
% -------------------------------------------------------------------------
src  = srcinfo.r;
[~, ns] = size(src);

try
    targ = targinfo.r(:,:);
catch
    targ = targinfo;
end
nt = size(targ, 2);

% -------------------------------------------------------------------------
% Compute local orthonormal frames at source points (ru, rv)
% ru  = du / |du|
% rv  = (n x du) / |n x du|
% -------------------------------------------------------------------------
du = srcinfo.du;
n_src = srcinfo.n;

ru = zeros(3, ns);
rv = zeros(3, ns);
for i = 1:ns
    u = du(:,i);
    u = u / norm(u);
    ru(:,i) = u;
    tmp = cross(n_src(:,i), u);
    rv(:,i) = tmp / norm(tmp);
end

% -------------------------------------------------------------------------
% Dispatch on kernel type
% -------------------------------------------------------------------------
switch lower(type)

    % ======================================================================
    case 'nrccie-bc'
    % ======================================================================
    %
    %  Representation:
    %    E = ik S_k[J]  -  grad S_k[rho]
    %    H = curl S_k[J]
    %
    %  sigma = [j_u(1:ns); j_v(1:ns); rho(1:ns)]  (already wt-scaled)
    %
    %  We reconstruct the Cartesian current J = j_u * ru + j_v * rv
    %  and call emfmm3d with:
    %    h_current = J          -> contributes curl S_k[J]   to E (not used here directly)
    %    e_current = ik * J     -> contributes S_k[ik*J]     to E
    %    e_charge  = -rho       -> contributes grad S_k[-rho] = -grad S_k[rho] to E
    %
    %  We request both E (= ik S_k[J] - grad S_k[rho]) and curlE (= H).
    %  Then project to local components [M_u; M_v; E_n].
    %
    %  varargin{1} = alpha (scalar, required but not used in FMM itself —
    %  the projection below is independent of alpha; alpha only appears in
    %  the near-quadrature and the integral-equation assembler).

    sigma = sigma(:);
    j_u = sigma(1:ns);
    j_v = sigma(ns+1:2*ns);
    rho = sigma(2*ns+1:3*ns);

    % Cartesian current (3 x ns)
    J = bsxfun(@times, ru, j_u.') + bsxfun(@times, rv, j_v.');

    % Pack srcinfo for emfmm3d
    src_em = [];
    src_em.sources   = src;
    src_em.h_current = J;                      % (3, ns), nd=1
    src_em.e_current = (1i*zk) * J;            % (3, ns)
    src_em.e_charge  = -rho(:).';              % (1, ns)

    % Evaluate E and curlE (= H) at targets
    ifE     = 1;
    ifcurlE = 1;
    ifdivE  = 0;
    U = emfmm3d(eps, zk, src_em, targ, ifE, ifcurlE, ifdivE);

    % E  (3 x nt), H = curlE  (3 x nt)
    E = U.E;      % 3 x nt
    H = U.curlE;  % 3 x nt

    % Project onto local target frame
    du_t  = targinfo.du;
    dv_t  = targinfo.dv;
    n_t   = targinfo.n;

    % Compute target orthonormal frame (ru_t, rv_t)
    ru_t = zeros(3, nt);
    rv_t = zeros(3, nt);
    for i = 1:nt
        u = du_t(:,i);
        u = u / norm(u);
        ru_t(:,i) = u;
        tmp = cross(n_t(:,i), u);
        rv_t(:,i) = tmp / norm(tmp);
    end

    % M_u = ru_t . H,  M_v = rv_t . H,  E_n = n_t . E
    M_u = sum(ru_t .* H, 1).';   % (nt x 1)
    M_v = sum(rv_t .* H, 1).';   % (nt x 1)
    E_n = sum(n_t  .* E, 1).';   % (nt x 1)

    pot = [M_u; M_v; E_n];       % (3*nt x 1)

    % ======================================================================
    case 'nrccie-eval'
    % ======================================================================
    %
    %  Representation:
    %    E = ik S_k[J]  -  grad S_k[rho]
    %    H = curl S_k[J]
    %
    %  sigma = [J_x(1:ns); J_y(1:ns); J_z(1:ns); rho(1:ns)]  (wt-scaled)
    %  pot   = [E_x; E_y; E_z; H_x; H_y; H_z]                (6*nt x 1)

    sigma = sigma(:);
    J   = reshape(sigma(1:3*ns), 3, ns);    % Cartesian current (3 x ns)
    rho = sigma(3*ns+1:4*ns).';             % (1 x ns)

    src_em = [];
    src_em.sources   = src;
    src_em.h_current = J;            % curl S_k[J]   -> H
    src_em.e_current = (1i*zk) * J;  % ik S_k[J]    -> E
    src_em.e_charge  = -rho;         % -grad S_k[rho] -> E

    ifE     = 1;
    ifcurlE = 1;
    ifdivE  = 0;
    U = emfmm3d(eps, zk, src_em, targ, ifE, ifcurlE, ifdivE);

    E = U.E;      % (3 x nt)
    H = U.curlE;  % (3 x nt)

    pot = [E(:); H(:)];   % (6*nt x 1)  [E_x E_y E_z then H_x H_y H_z]

    % Reshape to column-major (component varies slowest, matches kern.m convention):
    %   pot(1:nt)        = E_x
    %   pot(nt+1:2*nt)   = E_y
    %   pot(2*nt+1:3*nt) = E_z
    %   pot(3*nt+1:4*nt) = H_x
    %   pot(4*nt+1:5*nt) = H_y
    %   pot(5*nt+1:6*nt) = H_z
    pot = [E(1,:).'; E(2,:).'; E(3,:).'; ...
           H(1,:).'; H(2,:).'; H(3,:).'];

    otherwise
        error('EM3D:fmm:type', 'Unknown em3d kernel type ''%s''.', type);
end

end
