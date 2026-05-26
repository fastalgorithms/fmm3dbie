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
%   srcinfo  - source ptinfo struct 
%   targinfo - targ ptinfo struct
%   type     - kernel type string:
%                'nrccie-bc'    Non-resonant CFIE boundary-condition operator.
%                               sigma (3*ns x 1): interleaved [j_ru;j_rv;rho] per point
%                                 j_ru = coeff of ru_s = du_s/|du_s| (orthonormal frame)
%                                 j_rv = coeff of rv_s = n_s x ru_s
%                               pot (3*nt x 1): interleaved [pot_ru;pot_rv;pot_rho] per point
%                                 pot_ru/rv projected onto orthonormal target frame
%
%                'nrccie-eval'  NRCCIE field evaluation.
%                               sigma (4*ns x 1): [J_x; J_y; J_z; rho] Cartesian
%                               pot   (6*nt x 1): [E_x; E_y; E_z; H_x; H_y; H_z]
%
%   sigma    - density
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
% Dispatch on kernel type
% -------------------------------------------------------------------------
switch lower(type)

    % ======================================================================
    case 'nrccie-bc'
    % ======================================================================
    %
    %  Orthonormal-frame convention matching em3d.kern 'nrccie-bc' and Fortran.
    %
    %  sigma (3*ns x 1): interleaved [j_ru(1);j_rv(1);rho(1); j_ru(2);...] (wt-scaled)
    %    j_ru = coefficient of ru_s = du_s/|du_s|
    %    j_rv = coefficient of rv_s = n_s x ru_s
    %  pot   (3*nt x 1): interleaved [pot_ru(1);pot_rv(1);pot_rho(1); pot_ru(2);...]
    %    pot_ru = zvec3 . ru_t,  pot_rv = zvec3 . rv_t,  ru_t = du_t/|du_t|
    %
    %  zvec3 = alpha*nxnxE - n_t x (curl S_k[J])
    %  pot_rho = (grad_x S_k[rho]).n_t - ik*(S_k[J].n_t)
    %            + alpha*(div_x S_k[J] - ik*S_k[rho])

    alpha = varargin{1};

    sigma3 = reshape(sigma, 3, ns);   % rows: j_ru, j_rv, rho
    j_ru = sigma3(1,:);   % (1 x ns)
    j_rv = sigma3(2,:);
    rho  = sigma3(3,:);

    % Build orthonormal source frame and Cartesian current J = j_ru*ru_s + j_rv*rv_s
    du_s = srcinfo.du;   % (3, ns)
    n_s  = srcinfo.n;
    ru_s = du_s ./ vecnorm(du_s);              % (3, ns)
    rv_s = cross(n_s ./ vecnorm(n_s), ru_s, 1); % (3, ns)
    J = bsxfun(@times, ru_s, j_ru) + bsxfun(@times, rv_s, j_rv);   % (3, ns)

    % FMM A: H = curl S_k[J]  (used for n_t x H term)
    src_A = struct('sources', src, 'h_current', J);
    U_A   = emfmm3d(eps, zk, src_A, targ, 1, 0, 0);
    H     = U_A.E;   % (3, nt)

    % FMM B: S_k[J] (vector) and div S_k[J]
    src_B  = struct('sources', src, 'e_current', J);
    U_B    = emfmm3d(eps, zk, src_B, targ, 1, 0, 1);
    SkJ    = U_B.E;      % (3, nt)
    divSkJ = U_B.divE;   % (1, nt)

    % FMM C: S_k[rho] and grad_x S_k[rho]
    src_C     = struct('sources', src, 'charges', rho);
    U_C       = hfmm3d(eps, zk, src_C, 0, targ, 2);
    SkRho     = U_C.pottarg(:).';    % (1, nt)
    gradSkRho = U_C.gradtarg;        % (3, nt)

    % Build orthonormal target frame: ru_t = du_t/|du_t|,  rv_t = n_t x ru_t
    du_t = targinfo.du;   % (3, nt)
    n_t  = targinfo.n;
    ru_t = du_t ./ vecnorm(du_t);               % (3, nt)
    rv_t = cross(n_t ./ vecnorm(n_t), ru_t, 1); % (3, nt)

    % zvec3 = alpha*nxnxE - n_t x H,  where E = ik*S_k[J] - grad S_k[rho]
    nxH_x = n_t(2,:).*H(3,:) - n_t(3,:).*H(2,:);
    nxH_y = n_t(3,:).*H(1,:) - n_t(1,:).*H(3,:);
    nxH_z = n_t(1,:).*H(2,:) - n_t(2,:).*H(1,:);

    E_x = (1i*zk)*SkJ(1,:) - gradSkRho(1,:);
    E_y = (1i*zk)*SkJ(2,:) - gradSkRho(2,:);
    E_z = (1i*zk)*SkJ(3,:) - gradSkRho(3,:);
    ndotE    = n_t(1,:).*E_x + n_t(2,:).*E_y + n_t(3,:).*E_z;
    anx2E_x  = alpha*(ndotE.*n_t(1,:) - E_x);
    anx2E_y  = alpha*(ndotE.*n_t(2,:) - E_y);
    anx2E_z  = alpha*(ndotE.*n_t(3,:) - E_z);

    zvec3_x = anx2E_x - nxH_x;
    zvec3_y = anx2E_y - nxH_y;
    zvec3_z = anx2E_z - nxH_z;

    % Project onto orthonormal target frame
    pot_ru = ru_t(1,:).*zvec3_x + ru_t(2,:).*zvec3_y + ru_t(3,:).*zvec3_z;   % (1, nt)
    pot_rv = rv_t(1,:).*zvec3_x + rv_t(2,:).*zvec3_y + rv_t(3,:).*zvec3_z;

    % pot_rho = (grad_x S_k[rho]).n_t - ik*(S_k[J].n_t) + alpha*(div S_k[J] - ik*S_k[rho])
    gradSkRho_nt = n_t(1,:).*gradSkRho(1,:) + n_t(2,:).*gradSkRho(2,:) + n_t(3,:).*gradSkRho(3,:);
    SkJ_nt       = n_t(1,:).*SkJ(1,:)       + n_t(2,:).*SkJ(2,:)       + n_t(3,:).*SkJ(3,:);
    pot_rho = gradSkRho_nt - 1i*zk*SkJ_nt + alpha*(divSkJ(:).' - 1i*zk*SkRho);

    % Interleaved output: [pot_ru(1);pot_rv(1);pot_rho(1); pot_ru(2);...]
    pot = reshape([pot_ru; pot_rv; pot_rho], 3*nt, 1);

    case 'nrccie-eval'
    % ======================================================================
    %
    %  sigma (4*ns x 1): interleaved [Jx(1);Jy(1);Jz(1);rho(1); Jx(2);...] (wt-scaled)
    %  pot   (6*nt x 1): interleaved [Ex(1);Ey(1);Ez(1);Hx(1);Hy(1);Hz(1); Ex(2);...]
    %
    %  H = curl S_k[J],  E = ik S_k[J] - grad S_k[rho]

    sigma4 = reshape(sigma, 4, ns);   % rows: Jx, Jy, Jz, rho  (interleaved input)
    J   = sigma4(1:3,:);              % (3, ns)
    rho = sigma4(4,:);                % (1, ns)

    % Two separate emfmm3d calls
    src_H = struct('sources', src, 'h_current', J);
    U_H   = emfmm3d(eps, zk, src_H, targ, 1, 0, 0);
    H     = U_H.E;   % (3, nt)

    src_E = struct('sources', src, 'e_current', (1i*zk)*J, 'e_charge', -rho);
    U_E   = emfmm3d(eps, zk, src_E, targ, 1, 0, 0);
    E     = U_E.E;   % (3, nt)

    % Interleaved output: [Ex(1);Ey(1);Ez(1);Hx(1);Hy(1);Hz(1); Ex(2);...]
    EH  = [E; H];    % (6, nt)
    pot = EH(:);

    otherwise
        error('EM3D:fmm:type', 'Unknown em3d kernel type ''%s''.', type);
end

end
