function comps = radpp_tensors(dx,dy,dz,breaks,coefs,nderiv,dim)
%RADPP_TENSORS derivative tensors of a radial function from its on-axis
%  profiles, the matlab twin of radpp_tensors in src/kernels/radpp.f90.
%
%  comps = radpp_tensors(dx,dy,dz,breaks,coefs,nderiv)
%  comps = radpp_tensors(dx,dy,dz,breaks,coefs,nderiv,dim)
%
%  dx, dy, dz are displacements target - source (any shape, dz may be 0).
%  breaks, coefs are a radpp_fit of the on-axis profiles
%    [G | G_x | G_xx G_yy | G_xxx G_xyy]  (up to nderiv)
%  comps is (numel(dx), ncomp), derivatives with respect to the target:
%    dim = 3 (default), ncomp = 1, 4, 10, 20
%      [G | G_x G_y G_z | H_xx H_xy H_xz H_yy H_yz H_zz |
%       T_xxx T_xxy T_xxz T_xyy T_xyz T_xzz T_yyy T_yyz T_yzz T_zzz]
%    dim = 2, ncomp = 1, 3, 6, 10, in chunkie order (in-plane only)
%      [G | G_x G_y | H_xx H_xy H_yy | T_xxx T_xxy T_xyy T_yyy]
%
%  See also RADPP_FIT, KERNEL3D.RADPP

if nargin < 7 || isempty(dim), dim = 3; end

dx = dx(:);  dy = dy(:);  dz = dz(:) + 0*dx;
r = sqrt(dx.^2 + dy.^2 + dz.^2);
prof = radpp_eval(r,breaks,coefs);
comps = prof(:,1);

if dim == 3
    u = [dx dy dz]./r;
    u(r == 0,:) = repmat([1 0 0], nnz(r == 0), 1);
    if nderiv >= 1
        comps = [comps, prof(:,2).*u];
    end
    if nderiv >= 2
        % H_ij = (p-q) u_i u_j + q d_ij
        p = prof(:,3);  q = prof(:,4);
        for i = 1:3, for j = i:3
            comps(:,end+1) = (p-q).*u(:,i).*u(:,j) + q*(i==j);
        end, end
    end
    if nderiv >= 3
        % T_ijk = (a-3b) u_i u_j u_k + b (d_ij u_k + d_ik u_j + d_jk u_i)
        a = prof(:,5);  b = prof(:,6);
        for i = 1:3, for j = i:3, for k = j:3
            comps(:,end+1) = (a-3*b).*u(:,i).*u(:,j).*u(:,k) + ...
                b.*((i==j)*u(:,k) + (i==k)*u(:,j) + (j==k)*u(:,i));
        end, end, end
    end
    return
end

c = dx./r;   s = dy./r;
c(r == 0) = 1;  s(r == 0) = 0;

if nderiv >= 1
    gr = prof(:,2);
    comps = [comps, gr.*c, gr.*s];
end
if nderiv >= 2
    p = prof(:,3);  q = prof(:,4);
    comps = [comps, p.*c.^2 + q.*s.^2, (p-q).*c.*s, p.*s.^2 + q.*c.^2];
end
if nderiv >= 3
    a = prof(:,5);  b = prof(:,6);
    comps = [comps, a.*c.^3 + 3*b.*c.*s.^2, ...
                    a.*c.^2.*s + b.*(s.^3 - 2*c.^2.*s), ...
                    a.*c.*s.^2 + b.*(c.^3 - 2*c.*s.^2), ...
                    a.*s.^3 + 3*b.*s.*c.^2];
end

end
