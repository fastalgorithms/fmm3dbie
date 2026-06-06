function U = fmm(eps, src, targ, type, sigma, coefs)
%STOK3D.FMM   Fast multipole method for evaluating Stokes layer
%   potentials in 3D.
%
% Syntax:
%   U  = stok3d.fmm(eps, zk, srcinfo, targinfo, type, sigma, coefs)
%
% Input:
%   eps      - precision requested
%   src      - ptinfo struct: .r (3,:), .n (3,:)
%   targ     - ptinfo struct: .r (3,:), .n (3,:)
%   type     - 's'           single layer S
%              'd'           double layer D
%              'sp','sprime' traction of S at target (S')
%              'c'           combined layer coefs(1)*S + coefs(2)*D
%   sigma    - density (3 x ns)
%   coefs    - [alpha; beta] for combined layer (type 'c' only)
%
% Output:
%   pot  - potential (3 x nt)

sigma = reshape(sigma,3,[]);
switch lower(type)
    case {'s', 'single'}
        srcinfo = struct();
        srcinfo.sources = src.r;
        srcinfo.stoklet = sigma;
        U_out = stfmm3d(eps, srcinfo, 0, targ.r, 1);
        U = reshape(U_out.pottarg, 3, []);
    case {'d', 'double', 'traction'}
        srcinfo = struct();
        srcinfo.sources = src.r;
        srcinfo.strslet = sigma;
        srcinfo.strsvec = src.n;
        U_out = stfmm3d(eps, srcinfo, 0, targ.r, 1);
        U = reshape(U_out.pottarg, 3, []);
    case {'sp', 'sprime','strac'}
        % Traction of S at target: t(S)_{ij} n_j(x) = T_{ijk}(x,y) n_k(x)
        % Use stfmm3d with ifppregtarg=3 to get velocity + pressure + gradient,
        % then assemble traction using target normals.
        srcinfo = struct();
        srcinfo.sources = src.r;
        srcinfo.stoklet = sigma;
        % ifppregtarg=3: request velocity, pressure, and velocity gradient at targs
        U_out = stfmm3d(eps, srcinfo, 0, targ.r, 3);
        % gradtarg(i,k,it) = du_i/dx_k at target it
        % traction: t_i = n_k * (-p delta_{ik} + du_i/dx_k + du_k/dx_i)
        nt = size(targ.r, 2);
        n = targ.n;                          % (3, nt)
        pre = reshape(U_out.pretarg, 1, nt); % (1, nt)
        % gradtarg is (3, 3, nt): grad(i,k,it) = du_i/dx_k
        % Squeeze each (1,1,nt) slice to (1,nt) for broadcasting
        g11 = reshape(U_out.gradtarg(1,1,:), 1, nt);
        g12 = reshape(U_out.gradtarg(1,2,:), 1, nt);
        g13 = reshape(U_out.gradtarg(1,3,:), 1, nt);
        g21 = reshape(U_out.gradtarg(2,1,:), 1, nt);
        g22 = reshape(U_out.gradtarg(2,2,:), 1, nt);
        g23 = reshape(U_out.gradtarg(2,3,:), 1, nt);
        g31 = reshape(U_out.gradtarg(3,1,:), 1, nt);
        g32 = reshape(U_out.gradtarg(3,2,:), 1, nt);
        g33 = reshape(U_out.gradtarg(3,3,:), 1, nt);
        % Physical traction: t_i = n_k(-p delta_ik + du_i/dx_k + du_k/dx_i)
        U = zeros(3, nt);
        U(1,:) = n(1,:).*(-pre + 2*g11) + n(2,:).*(g12 + g21) + n(3,:).*(g13 + g31);
        U(2,:) = n(1,:).*(g21 + g12)    + n(2,:).*(-pre + 2*g22) + n(3,:).*(g23 + g32);
        U(3,:) = n(1,:).*(g31 + g13)    + n(2,:).*(g32 + g23) + n(3,:).*(-pre + 2*g33);
    case {'c', 'combined'}
        srcinfo = struct();
        srcinfo.sources = src.r;
        if coefs(1) ~= 0
            srcinfo.stoklet = coefs(1) * sigma;
        end
        if coefs(2) ~= 0
            srcinfo.strslet = coefs(2) * sigma;
            srcinfo.strsvec = src.n;
        end
        U_out = stfmm3d(eps, srcinfo, 0, targ.r, 1);
        U = reshape(U_out.pottarg, 3, []);
    case {'cp', 'cprime','ctrac'}
        % Traction of C at target
        % Use stfmm3d with ifppregtarg=3 to get velocity + pressure + gradient,
        % then assemble traction using target normals.
        srcinfo = struct();
        srcinfo.sources = src.r;
        if coefs(1) ~= 0
            srcinfo.stoklet = coefs(1) * sigma;
        end
        if coefs(2) ~= 0
            srcinfo.strslet = coefs(2) * sigma;
            srcinfo.strsvec = src.n;
        end
        % ifppregtarg=3: request velocity, pressure, and velocity gradient at targs
        U_out = stfmm3d(eps, srcinfo, 0, targ.r, 3);
        % gradtarg(i,k,it) = du_i/dx_k at target it
        % traction: t_i = n_k * (-p delta_{ik} + du_i/dx_k + du_k/dx_i)
        nt = size(targ.r, 2);
        n = targ.n;                          % (3, nt)
        pre = reshape(U_out.pretarg, 1, nt); % (1, nt)
        % gradtarg is (3, 3, nt): grad(i,k,it) = du_i/dx_k
        % Squeeze each (1,1,nt) slice to (1,nt) for broadcasting
        g11 = reshape(U_out.gradtarg(1,1,:), 1, nt);
        g12 = reshape(U_out.gradtarg(1,2,:), 1, nt);
        g13 = reshape(U_out.gradtarg(1,3,:), 1, nt);
        g21 = reshape(U_out.gradtarg(2,1,:), 1, nt);
        g22 = reshape(U_out.gradtarg(2,2,:), 1, nt);
        g23 = reshape(U_out.gradtarg(2,3,:), 1, nt);
        g31 = reshape(U_out.gradtarg(3,1,:), 1, nt);
        g32 = reshape(U_out.gradtarg(3,2,:), 1, nt);
        g33 = reshape(U_out.gradtarg(3,3,:), 1, nt);
        % Physical traction: t_i = n_k(-p delta_ik + du_i/dx_k + du_k/dx_i)
        U = zeros(3, nt);
        U(1,:) = n(1,:).*(-pre + 2*g11) + n(2,:).*(g12 + g21) + n(3,:).*(g13 + g31);
        U(2,:) = n(1,:).*(g21 + g12)    + n(2,:).*(-pre + 2*g22) + n(3,:).*(g23 + g32);
        U(3,:) = n(1,:).*(g31 + g13)    + n(2,:).*(g32 + g23) + n(3,:).*(-pre + 2*g33);
end
end