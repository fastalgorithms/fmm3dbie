% kernel3d(@green, type) with a tensor green: the Maxwell dyadic
% G_ij = (d_ij + d_i d_j / k^2) g, g = exp(ikr)/(4 pi r). A template for
% vector problems whose green is not a scalar.

% run ../startup.m

rng(42);

test_dyad_checks();
test_dyad_eval();
test_dyad_mfie_vs_scalar();
test_dyad_ipv();
test_dyad_getquad();


function test_dyad_checks()
% tensor greens are detected and checked at build time

zk = 1.3;
K = kernel3d(@(s,t) dyadg(zk, s, t), 's');
assert(K.params.green.tensor && K.params.dim == 13, 'dyad: not a tensor green');
assert(isequal(K.opdims, [3 3]), 'dyad: opdims');

bad = {@() kernel3d(@(s,t) skewg(s, t), 's'), ...            % not isotropic
       @() kernel3d(@(s,t) dyadg(zk, s, t), 's', 'dim', 2), ...
       @() kernel3d(@(s,t) wronggrad(zk, s, t), 'curl')};
for k = 1:numel(bad)
    caught = false;
    try
        bad{k}();
    catch
        caught = true;
    end
    assert(caught, 'dyad: bad input %d not rejected', k);
end

end


function test_dyad_eval()
% eval matches the formulas written out from dyadg

zk = 1.3;
src.r  = randn(3, 20);
targ.r = randn(3, 15) + [3;0;0];
targ.n = randn(3, 15);  targ.n = targ.n./vecnorm(targ.n);
ns = 20;  nt = 15;
[G, dG] = dyadg(zk, src.r, targ.r);              % (nt,ns,3,3), (..,3)

E = @(i,j,k) (j-i)*(k-i)*(k-j)/2;
ref = struct('s', zeros(3,nt,3,ns), 'efie', 0, 'curl', 0, 'mfie', 0);
[ref.efie, ref.curl, ref.mfie] = deal(zeros(3,nt,3,ns));
n = targ.n.';
for i = 1:3, for j = 1:3
    ref.s(i,:,j,:) = reshape(G(:,:,i,j), 1, nt, 1, ns);
    ef = 0;  cu = 0;  mf = 0;
    for a = 1:3, for b = 1:3
        ef = ef + E(i,a,b)*n(:,a).*G(:,:,b,j);
        cu = cu + E(i,a,b)*dG(:,:,b,j,a);
        for c = 1:3, for d = 1:3
            mf = mf + E(i,c,d)*E(d,a,b)*n(:,c).*dG(:,:,b,j,a);
        end, end
    end, end
    ref.efie(i,:,j,:) = reshape(ef + 0*G(:,:,1,1), 1, nt, 1, ns);
    ref.curl(i,:,j,:) = reshape(cu + 0*G(:,:,1,1), 1, nt, 1, ns);
    ref.mfie(i,:,j,:) = reshape(mf + 0*G(:,:,1,1), 1, nt, 1, ns);
end, end

for type = fieldnames(ref).'
    K = kernel3d(@(s,t) dyadg(zk, s, t), type{1});
    B = reshape(ref.(type{1}), 3*nt, 3*ns);
    err = norm(K.eval(src, targ) - B, 'fro')/norm(B, 'fro');
    assert(err < 1e-12, 'dyad %s eval: %.2e', type{1}, err);
end

end


function test_dyad_mfie_vs_scalar()
% n x curl (G J) from the dyadic equals the scalar MFIE: curl kills the
% grad grad part

zk = 1.3;
src.r  = randn(3, 20);
targ.r = randn(3, 15) + [3;0;0];
targ.n = randn(3, 15);  targ.n = targ.n./vecnorm(targ.n);
Kd = kernel3d(@(s,t) dyadg(zk, s, t), 'mfie');
Ks = kernel3d(@(s,t) helm3d.green(zk, s, t), 'mfie');
A = Kd.eval(src, targ);  B = Ks.eval(src, targ);
err = norm(A - B, 'fro')/norm(B, 'fro');
assert(err < 1e-10, 'dyad mfie vs scalar: %.2e', err);

end


function test_dyad_ipv()
% the surface probe: G and n x G are hypersingular (1/r^3), curl and
% n x curl are principal value (1/r^2), zero entries get 0

zk = 1.3;
g = @(s,t) dyadg(zk, s, t);
assert(all(kernel3d(g, 's').params.ipv == 2), 'dyad ipv s');
assert(all(kernel3d(g, 'efie').params.ipv == 2), 'dyad ipv efie');
ipv = kernel3d(g, 'curl').params.ipv;
assert(all(ipv([1 5 9]) == 0) && all(ipv([2 3 4 6 7 8]) == 1), ...
       'dyad ipv curl: %s', mat2str(ipv));
assert(all(kernel3d(g, 'mfie').params.ipv == 1), 'dyad ipv mfie');

end


function test_dyad_getquad()
% near quadrature on a curved patch against the smooth integral

zk = 1.3;
g = @(s,t) dyadg(zk, s, t);
for type = {'efie', 'curl', 'mfie'}
    err = kernel3d.check_quad(kernel3d(g, type{1}), 1e-7);
    assert(err < 1e-6, 'dyad %s check: %.2e', type{1}, err);
end

end


function [G, dG] = dyadg(zk, s, t)
% Maxwell dyadic G_ij = (d_ij + d_i d_j / zk^2) g and d_k G_ij, from the
% radial derivatives of g = exp(i zk r)/(4 pi r)
x = t(1,:).' - s(1,:);  y = t(2,:).' - s(2,:);  z = t(3,:).' - s(3,:);
r = sqrt(x.^2 + y.^2 + z.^2);
U = cat(3, x, y, z)./r;
f  = exp(1i*zk*r)./(4*pi*r);  h = 1i*zk - 1./r;
f1 = f.*h;  f2 = f.*(h.^2 + 1./r.^2);  f3 = f.*(h.^3 + 3*h./r.^2 - 2./r.^3);
A  = f2 - f1./r;                    % hess g = A u u + (f1/r) I
Bt = f3 - f2./r + f1./r.^2 - 2*A./r;
G = zeros([size(x) 3 3]);  dG = zeros([size(x) 3 3 3]);
for i = 1:3, for j = 1:3
    G(:,:,i,j) = f*(i==j) + (A.*U(:,:,i).*U(:,:,j) + (i==j)*f1./r)/zk^2;
    for k = 1:3
        T = Bt.*U(:,:,i).*U(:,:,j).*U(:,:,k) + A./r.*((i==j)*U(:,:,k) ...
            + (i==k)*U(:,:,j) + (j==k)*U(:,:,i));
        dG(:,:,i,j,k) = f1.*U(:,:,k)*(i==j) + T/zk^2;
    end
end, end
end


function [G, dG] = wronggrad(zk, s, t)
[G, dG] = dyadg(zk, s, t);
dG = 1.01*dG;
end


function G = skewg(s, t)
% a tensor that is not isotropic: fixed weight on the x axis
x = t(1,:).' - s(1,:);  y = t(2,:).' - s(2,:);  z = t(3,:).' - s(3,:);
r = sqrt(x.^2 + y.^2 + z.^2);
G = zeros([size(x) 3 3]);
G(:,:,1,1) = 1./r;  G(:,:,2,2) = 2./r;  G(:,:,3,3) = 2./r;
end
