function [val, grad, hess] = helmdiff_green(zk, coefs, src, targ)
%HELM3D.HELMDIFF_GREEN difference of Helmholtz Green's functions,
%   G = a0 g_{k0} - a1 g_{k1}, g_k = exp(ik|x-y|)/(4 pi |x-y|), with its
%   gradient and Hessian, in the helm3d.green format.
%
%   [val, grad, hess] = helm3d.helmdiff_green([k0; k1], [a0; a1], src, targ)
%
%   The singular terms are combined analytically, so for a0 = a1 the
%   bounded G and the 1/r Hessian are free of cancellation at small r.
%   coefs defaults to [1; 1].
%
% See also HELM3D.GREEN

if isempty(coefs), coefs = [1; 1]; end
c = [coefs(1); -coefs(2)];
k = zk(:);

x = targ(1,:).' - src(1,:);
y = targ(2,:).' - src(2,:);
z = targ(3,:).' - src(3,:);
r = sqrt(x.^2 + y.^2 + z.^2);

f = 1/(4*pi);
val = f*expsum(c, k, r, 1, 1);                   % G
if nargout < 2, return; end

g1 = f*expsum(c, k, r, [-1 1], 2);               % G'
U = cat(3, x, y, z)./r;
grad = g1.*U;
if nargout < 3, return; end

% hess = G'' u u' + (G'/r)(I - u u'), with G'' and G'/r each summed
% directly: for a0 = a1, G'' is bounded while G'/r ~ 1/r, so G'' must not
% be formed by cancelling two 1/r terms
g2 = f*expsum(c, k, r, [2 -2 1], 3);             % G''
gr = f*expsum(c, k, r, [-1 1], 3);               % G'/r
hess = zeros([size(r) 3 3]);
for i = 1:3
    for j = 1:3
        uu = U(:,:,i).*U(:,:,j);
        hess(:,:,i,j) = g2.*uu + gr.*((i==j) - uu);
    end
end

end


function s = expsum(c, k, r, p, m)
% sum_j c_j exp(z_j) p(z_j) / r^m, z_j = i k_j r, p ascending coefficients.
% exp(z) p(z) = T(z) + R(z), T its Taylor polynomial of degree m-1: the
% singular part sum_n t_n (sum_j c_j (i k_j)^n) r^(n-m) uses exact
% coefficient sums, and R(z) = O(z^m) comes from a series for |z| < 1
N = 40;
t = zeros(1, N+1);                               % Taylor coefficients
for n = 0:N
    for q = 0:min(n, numel(p) - 1)
        t(n+1) = t(n+1) + p(q+1)/factorial(n - q);
    end
end

s = zeros(size(r));
for n = 0:m-1
    Sn = sum(c.*(1i*k).^n);
    if Sn ~= 0
        s = s + t(n+1)*Sn*r.^(n - m);
    end
end

for j = 1:numel(c)
    zz = 1i*k(j)*r;
    R = zeros(size(r));
    sm = abs(zz) < 1;
    zs = zz(sm);
    Rs = zeros(size(zs));
    zn = zs.^m;
    for n = m:N
        Rs = Rs + t(n+1)*zn;
        zn = zn.*zs;
    end
    R(sm) = Rs;
    zb = zz(~sm);
    T = zeros(size(zb));
    for n = m-1:-1:0
        T = T.*zb + t(n+1);
    end
    R(~sm) = exp(zb).*polyval(fliplr(p), zb) - T;
    s = s + c(j)*R./r.^m;
end

end
