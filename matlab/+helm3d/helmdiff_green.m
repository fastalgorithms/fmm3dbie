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

% Taylor coefficients t_n = sum_q p_q/(n-q)!, n = 0..N. The |z| < 1 series
% tail past N is bounded by sum|p|/(N-numel(p)+1)!, which at N = 30 is
% O(1e-31) -- far below double precision.
N  = 30;
np = numel(p);
invf = ones(1, N+1);                             % 1/n! by recurrence
for n = 1:N, invf(n+1) = invf(n)/n; end
t = zeros(1, N+1);
for n = 0:N
    for q = 0:min(n, np-1)
        t(n+1) = t(n+1) + p(q+1)*invf(n-q+1);
    end
end

% r^(-1) .. r^(-m) by repeated multiplication rather than m calls to power
ri = 1./r;
rin = cell(1, m);  rin{1} = ri;
for q = 2:m, rin{q} = rin{q-1}.*ri; end

s = zeros(size(r));
for n = 0:m-1
    Sn = sum(c.*(1i*k).^n);
    if Sn ~= 0
        s = s + (t(n+1)*Sn)*rin{m-n};            % r^(n-m)
    end
end

ph = fliplr(p(:).');                             % descending, for Horner
for j = 1:numel(c)
    ik = 1i*k(j);
    % |z| < 1 <=> r < 1/|k|, since r >= 0: no complex abs over the array
    if k(j) == 0
        sm = true(size(r));
    else
        sm = r < 1/abs(k(j));
    end
    R = zeros(size(r));

    % near: R = z^m * sum_{i=0}^{N-m} t_{m+i} z^i, by Horner
    zs = ik*r(sm);
    Rs = t(N+1);
    for n = N-1:-1:m, Rs = Rs.*zs + t(n+1); end
    zm = zs;
    for q = 2:m, zm = zm.*zs; end
    R(sm) = Rs.*zm;

    % far: R = exp(z) p(z) - T(z), both by Horner
    zb = ik*r(~sm);
    P = ph(1);
    for q = 2:numel(ph), P = P.*zb + ph(q); end
    T = t(m);
    for n = m-2:-1:0, T = T.*zb + t(n+1); end
    R(~sm) = exp(zb).*P - T;

    s = s + c(j)*(R.*rin{m});
end

end
