% Verify kernel3d.radpp: the fit, the tensor evaluation, the built-in
% kernels and the near quadrature.

% run ../startup.m

rng(42);

test_fit_limits();
% 2D (in-plane derivatives, dim = 2)
test_tensors_2d();
test_builtins_2d(@gauss, 1e-11);
test_builtins_2d(@(s,t) helmgreen(1.1, s, t), 1e-10);
test_helm_sp_eval_2d();
test_getquad_rows(2);

% 3D (full derivatives, dim = 3, the default)
test_tensors_3d();
test_builtins_3d();
test_getquad_sp_3d();
test_getquad_rows(3);

% value kernels and interfaces
test_getquad_lap();
test_kernel3d_eval();
test_kernel3d_getquad_helm();
test_radcheb_interface();
test_bad_inputs();


function test_fit_limits()
% nlevmax gives a warning and keeps the panels, maxsub is an error.

f = @(r) 1./r;
lastwarn('');
warning('off', 'RADPP_FIT:nlevmax');
[~, ~, info] = radpp_fit(f, 1e-6, 1, struct('nlevmax', 3));
warning('on', 'RADPP_FIT:nlevmax');
[~, id] = lastwarn();
assert(strcmp(id, 'RADPP_FIT:nlevmax'), 'radpp_fit: no nlevmax warning');
assert(info.nbin == 4, 'radpp_fit: unresolved panels not kept');

caught = false;
try
    radpp_fit(f, 1e-6, 1, struct('maxsub', 10));
catch
    caught = true;
end
assert(caught, 'radpp_fit: maxsub not an error');

end


function test_tensors_2d()
% In-plane tensor components from the on-axis profiles match a radial
% function with known derivatives, G = exp(-r^2/2).

K = kernel3d.radpp(@gauss, [1e-6 5], @(t) ones(size(t.r,2), 10), ...
                   'nderiv', 3, 'dim', 2);
p = K.params;

ns = 20;  nt = 30;
s = randn(2, ns);  t = randn(2, nt);
dx = t(1,:).' - s(1,:);  dy = t(2,:).' - s(2,:);
C = radpp_tensors(dx, dy, 0*dx, p.breaks, p.coefs, 3, 2);

[v, g, h, th] = gauss(s, t);
D = [v(:), reshape(g, [], 2), reshape(h, [], 3), reshape(th, [], 4)];
err = max(abs(C(:) - D(:)));
assert(err < 1e-11, 'radpp tensors 2d: %.2e', err);

end


function test_tensors_3d()
% Full 3D tensor components match G = exp(-r^2/2) off the plane:
%   G_i = -x_i G,  H_ij = (x_i x_j - d_ij) G,
%   T_ijk = (d_ij x_k + d_ik x_j + d_jk x_i - x_i x_j x_k) G

K = kernel3d.radpp(@gauss, [1e-6 6], @(t) ones(size(t.r,2), 20), ...
                   'nderiv', 3);
p = K.params;
assert(p.ncomp == 20 && p.dim == 3, 'radpp: 3d default not used');

n = 500;
X = 1.5*randn(n, 3);
C = radpp_tensors(X(:,1), X(:,2), X(:,3), p.breaks, p.coefs, 3);

G = exp(-sum(X.^2, 2)/2);
D = [G, -X.*G];
for i = 1:3, for j = i:3
    D(:,end+1) = (X(:,i).*X(:,j) - (i==j)).*G;
end, end
for i = 1:3, for j = i:3, for k = j:3
    D(:,end+1) = ((i==j)*X(:,k) + (i==k)*X(:,j) + (j==k)*X(:,i) ...
                  - X(:,i).*X(:,j).*X(:,k)).*G;
end, end, end
err = max(abs(C(:) - D(:)));
assert(err < 1e-11, 'radpp tensors 3d: %.2e', err);

% the plane z = 0 agrees with the 2D tensors
K2 = kernel3d.radpp(@gauss, [1e-6 6], @(t) ones(size(t.r,2), 10), ...
                    'nderiv', 3, 'dim', 2);
C3 = radpp_tensors(X(:,1), X(:,2), 0*X(:,1), p.breaks, p.coefs, 3, 3);
C2 = radpp_tensors(X(:,1), X(:,2), 0*X(:,1), K2.params.breaks, ...
                   K2.params.coefs, 3, 2);
i2 = [1 2 3 5 6 8 11 12 14 17];      % x,y-only components of the 3D list
err = max(abs(C3(:,i2) - C2), [], 'all');
assert(err < 1e-11, 'radpp tensors 3d vs 2d: %.2e', err);

end


function test_builtins_2d(green, tol)
% Built-in kernels, including the flexural boundary conditions, match
% directional derivatives of green contracted with full 2x2 and 2x2x2
% tensors.

nu = 0.3;
ns = 7;  nt = 9;
src.r = 2*rand(2, ns) - 1;
th = 2*pi*rand(1, nt);
targ.r  = 1.7*[cos(th); sin(th)];
targ.d  = 1.7*[-sin(th); cos(th)];
targ.d2 = -targ.r;
targ.n  = [cos(th); sin(th)];
kappa = 1/1.7;

[v, g, h, t3] = green(src.r, targ.r);
G1 = @(u) g(:,:,1).*u(1) + g(:,:,2).*u(2);
H  = {h(:,:,1), h(:,:,2); h(:,:,2), h(:,:,3)};
T  = cell(2,2,2);
T{1,1,1} = t3(:,:,1);
[T{1,1,2}, T{1,2,1}, T{2,1,1}] = deal(t3(:,:,2));
[T{1,2,2}, T{2,1,2}, T{2,2,1}] = deal(t3(:,:,3));
T{2,2,2} = t3(:,:,4);

ref = struct('s', zeros(nt, ns), 'sp', 0, 'clamped_bc', 0, ...
             'supported_bc', 0, 'free_bc', 0);
ref.sp = zeros(nt, ns);
ref.clamped_bc = zeros(2*nt, ns);
ref.supported_bc = zeros(2*nt, ns);
ref.free_bc = zeros(2*nt, ns);
for i = 1:nt
    n = targ.n(:,i);  tau = targ.d(:,i)/norm(targ.d(:,i));
    hnn = H2(H, n, n, i);  htt = H2(H, tau, tau, i);
    ref.s(i,:) = v(i,:);
    gn = G1(n);  ref.sp(i,:) = gn(i,:);
    ref.clamped_bc(2*i-1,:) = v(i,:);
    ref.clamped_bc(2*i,:)   = gn(i,:);
    ref.supported_bc(2*i-1,:) = v(i,:);
    ref.supported_bc(2*i,:)   = hnn + nu*htt;
    ref.free_bc(2*i-1,:) = hnn + nu*htt;
    ref.free_bc(2*i,:)   = T3(T, n, n, n, i) + (2-nu)*T3(T, n, tau, tau, i) ...
                         + (1-nu)*kappa*(htt - hnn);
end

types = fieldnames(ref);
for k = 1:numel(types)
    K = kernel3d.radpp(green, [1e-6 6], types{k}, 'nu', nu, 'dim', 2);
    A = K.eval(src, targ);
    B = ref.(types{k});
    err = norm(A - B, 'fro')/norm(B, 'fro');
    assert(err < tol, 'radpp %s: %.2e', types{k}, err);
end

end


function val = H2(H, u, w, i)
val = 0;
for a = 1:2, for b = 1:2
    val = val + H{a,b}(i,:)*u(a)*w(b);
end, end
end


function val = T3(T, u, w, z, i)
val = 0;
for a = 1:2, for b = 1:2, for c = 1:2
    val = val + T{a,b,c}(i,:)*u(a)*w(b)*z(c);
end, end, end
end


function test_getquad_lap()
% 's' with the 3D Laplace Green's function matches lap3d on a sphere.

eps = 1e-7;
S = geometries.sphere(1, 2, [0;0;0], 4, 1);

lap = @(s,t) 1./(4*pi*sqrt((t(1,:).' - s(1,:)).^2 + (t(2,:).' - s(2,:)).^2));
K = kernel3d.radpp(lap, [1e-9 4], 's');

Q  = K.getquad(S, eps);
Qe = kernel3d.lap3d('s').getquad(S, eps);

assert(~K.ifcomplex, 'radpp: laplace kernel flagged complex');
err = norm(Q - Qe, 'fro')/norm(Qe, 'fro');
assert(err < 1e-6, 'radpp getquad lap: %.2e', err);

end


function test_getquad_rows(dim)
% The rows of a multi-row kernel are quadratures of their own weights:
% row 1 of 'clamped_bc' is 's', and a custom kernel with the weights of
% 'clamped_bc' reproduces it exactly. (Scaling the weights is not exact:
% the adaptive integration refines on an absolute tolerance.)

eps = 1e-7;
S = geometries.sphere(1, 2, [0;0;0], 4, 1);
targ.r = 1.5*randn(3, 40);
targ.r = targ.r./vecnorm(targ.r)*1.05;
targ.n = targ.r./vecnorm(targ.r);

zk = 1.1;
helm = @(s,t) helmgreen(zk, s, t);

Ks = kernel3d.radpp(helm, [1e-9 4], 's', 'dim', dim);
Kc = kernel3d.radpp(helm, [1e-9 4], 'clamped_bc', 'ipv', [0 0], ...
                    'dim', dim);
Kw = kernel3d.radpp(helm, [1e-9 4], ...
        @(t) cat(3, [1+0*t.r(1,:).', zeros(size(t.r,2),dim)], ...
                    [zeros(size(t.r,2),1), t.n(1:dim,:).']), ...
        'nderiv', 1, 'opdims', [2 1], 'targ_fields', {'n'}, 'dim', dim);

Qs = Ks.getquad(S, eps, targ);
Qc = Kc.getquad(S, eps, targ);
Qw = Kw.getquad(S, eps, targ);

err1 = norm(Qc(1:2:end,:) - Qs, 'fro')/norm(Qs, 'fro');
err2 = norm(Qw - Qc, 'fro')/norm(Qc, 'fro');
assert(err1 < 1e-12, 'radpp getquad rows (dim %d): %.2e', dim, err1);
assert(err2 == 0, 'radpp getquad custom (dim %d): %.2e', dim, err2);

end


function test_kernel3d_eval()
% A kernel3d green (Helmholtz single layer): type defaults to 's', zk is
% inherited, and K.eval matches the kernel for points in full 3D.

zk = 0.8;
fk = kernel3d.helm3d('s', zk);
K = kernel3d.radpp(fk, [1e-2 6]);

assert(K.zk == zk, 'radpp: zk not inherited from kernel3d');
assert(K.ifcomplex, 'radpp: helmholtz kernel flagged real');

ns = 200;  nt = 150;
src.r  = randn(3, ns);  src.r  = src.r./vecnorm(src.r);
targ.r = randn(3, nt);  targ.r = targ.r./vecnorm(targ.r)*2.5;

A  = K.eval(src, targ);
Ae = fk.eval(src, targ);
err = norm(A - Ae, 'fro')/norm(Ae, 'fro');
assert(isequal(size(A), [nt ns]), 'radpp eval: wrong size');
assert(err < 1e-10, 'radpp kernel3d eval: %.2e', err);

% laplace is real
Kl = kernel3d.radpp(kernel3d.lap3d('s'), [1e-2 6]);
assert(~Kl.ifcomplex, 'radpp: laplace kernel3d flagged complex');

end


function test_kernel3d_getquad_helm()
% Near quadrature of the radpp form of helm3d('s') matches helm3d, on
% surface and for nearby off-surface targets.

zk = 1.0;
eps = 1e-7;
S = geometries.sphere(1, 2, [0;0;0], 4, 1);
fk = kernel3d.helm3d('s', zk);
K = kernel3d.radpp(fk, [1e-9 4], 's', 'opts', struct('eps', 1e-13));

Q  = K.getquad(S, eps);
Qe = fk.getquad(S, eps);
err = norm(Q - Qe, 'fro')/norm(Qe, 'fro');
assert(err < 1e-6, 'radpp getquad helm: %.2e', err);

targ.r = randn(3, 50);
targ.r = targ.r./vecnorm(targ.r)*1.05;
Q  = K.getquad(S, eps, targ);
Qe = fk.getquad(S, eps, targ);
err = norm(Q - Qe, 'fro')/norm(Qe, 'fro');
assert(err < 1e-6, 'radpp getquad helm off-surface: %.2e', err);

end


function test_radcheb_interface()
% The old kernel3d.radcheb(fkern, rmax, opts) call still works.

zk = 1.3;
fk = kernel3d.helm3d('s', zk);
K = kernel3d.radpp(fk, 4, struct('eps', 1e-13, 'rmin', 1e-6, ...
                                 'norder', 16, 'zk', 2.0));
assert(K.params.rmin == 1e-6 && K.params.rmax == 4, ...
       'radpp: rmin from opts not used');
assert(K.zk == 2.0, 'radpp: zk from opts not used');
assert(K.params.fitopts.norder == 16, 'radpp: fit opts not passed');

r = linspace(1e-6, 4, 500).';
val = radpp_eval(r, K.params.breaks, K.params.coefs);
ex  = exp(1i*zk*r)./(4*pi*r);
err = max(abs(val(:,1) - ex)./abs(ex));
assert(err < 1e-10, 'radpp radcheb interface fit: %.2e', err);

assert(all(isnan(radpp_eval([1e-7; 5], K.params.breaks, K.params.coefs))), ...
       'radpp eval: out of range radii should be NaN');

end


function test_helm_sp_eval_2d()
% A derivative kernel with a standard Green's function: 'sp' from the 3D
% Helmholtz Green's function in the plane z = 0 is its in-plane normal
% derivative.

zk = 1.1;
K = kernel3d.radpp(@(s,t) helmgreen(zk, s, t), [1e-6 6], 'sp', 'dim', 2);

ns = 30;  nt = 40;
src.r = 2*rand(2, ns) - 1;
th = 2*pi*rand(1, nt);
targ.r = 1.7*[cos(th); sin(th)];
targ.n = [cos(th); sin(th)];

[~, g] = helmgreen(zk, src.r, targ.r);
B = g(:,:,1).*targ.n(1,:).' + g(:,:,2).*targ.n(2,:).';
A = K.eval(src, targ);
err = norm(A - B, 'fro')/norm(B, 'fro');
assert(err < 1e-10, 'radpp helm sp 2d: %.2e', err);

end


function test_builtins_3d()
% 's', 'sp' and 'clamped_bc' from the Helmholtz Green's function match
% helm3d for points in full 3D.

zk = 1.1;
helm = @(s,t) helmgreen(zk, s, t);
ns = 50;  nt = 60;
src.r  = randn(3, ns);  src.r  = src.r./vecnorm(src.r);
targ.r = randn(3, nt);  targ.r = targ.r./vecnorm(targ.r)*1.7;
targ.n = randn(3, nt);  targ.n = targ.n./vecnorm(targ.n);

S  = kernel3d.helm3d('s', zk).eval(src, targ);
SP = kernel3d.helm3d('sp', zk).eval(src, targ);
CL = zeros(2*nt, ns);
CL(1:2:end,:) = S;  CL(2:2:end,:) = SP;

ref = struct('s', S, 'sp', SP, 'clamped_bc', CL);
types = fieldnames(ref);
for k = 1:numel(types)
    K = kernel3d.radpp(helm, [1e-6 6], types{k});
    A = K.eval(src, targ);
    B = ref.(types{k});
    err = norm(A - B, 'fro')/norm(B, 'fro');
    assert(err < 1e-10, 'radpp %s 3d: %.2e', types{k}, err);
end

end


function test_getquad_sp_3d()
% Near quadrature of 'sp' in 3D matches helm3d('sp') on a sphere.

zk = 1.1;
eps = 1e-7;
S = geometries.sphere(1, 2, [0;0;0], 4, 1);
K = kernel3d.radpp(@(s,t) helmgreen(zk, s, t), [1e-9 4], 'sp', ...
                   'opts', struct('eps', 1e-13));

Q  = K.getquad(S, eps);
Qe = kernel3d.helm3d('sp', zk).getquad(S, eps);
err = norm(Q - Qe, 'fro')/norm(Qe, 'fro');
assert(err < 1e-6, 'radpp getquad sp 3d: %.2e', err);

end


function test_bad_inputs()
% Vector valued kernel3d greens, ones with field dependence, and
% derivative kernels from a value-only kernel3d are rejected.

bad = {@() kernel3d.radpp(kernel3d.stok3d('s'), 1), ...
       @() kernel3d.radpp(kernel3d.lap3d('d'), 1), ...
       @() kernel3d.radpp(kernel3d.lap3d('s'), 1, 'sp'), ...
       @() kernel3d.radpp(@gauss, 1, 'free_bc', 'nu', 0.3)};
for k = 1:numel(bad)
    caught = false;
    try
        bad{k}();
    catch
        caught = true;
    end
    assert(caught, 'radpp: bad input %d not rejected', k);
end

end


function varargout = gauss(s, t)
% G = exp(-r^2/2) and its target derivatives, chunkie format
x = t(1,:).' - s(1,:);  y = t(2,:).' - s(2,:);
G = exp(-(x.^2 + y.^2)/2);
varargout{1} = G;
varargout{2} = cat(3, -x.*G, -y.*G);
varargout{3} = cat(3, (x.^2-1).*G, x.*y.*G, (y.^2-1).*G);
varargout{4} = cat(3, (3*x - x.^3).*G, (y - x.^2.*y).*G, ...
                      (x - x.*y.^2).*G, (3*y - y.^3).*G);
end


function varargout = helmgreen(zk, s, t)
% 3D Helmholtz Green's function in the plane z = 0 and its target
% derivatives up to third order, chunkie format
x = t(1,:).' - s(1,:);  y = t(2,:).' - s(2,:);
r = sqrt(x.^2 + y.^2);
ux = x./r;  uy = y./r;

% radial derivatives of f = exp(i zk r)/(4 pi r), with h = i zk - 1/r
f  = exp(1i*zk*r)./(4*pi*r);
h  = 1i*zk - 1./r;
f1 = f.*h;
f2 = f.*(h.^2 + 1./r.^2);
f3 = f.*(h.^3 + 3*h./r.^2 - 2./r.^3);

% H_ij = A u_i u_j + (f1/r) d_ij
% T_ijk = (A' - 2A/r) u_i u_j u_k + (A/r)(d_ij u_k + d_ik u_j + d_jk u_i)
A  = f2 - f1./r;
B  = (f3 - f2./r + f1./r.^2) - 2*A./r;
C  = A./r;

varargout{1} = f;
varargout{2} = cat(3, f1.*ux, f1.*uy);
varargout{3} = cat(3, A.*ux.^2 + f1./r, A.*ux.*uy, A.*uy.^2 + f1./r);
varargout{4} = cat(3, B.*ux.^3 + 3*C.*ux, ...
                      B.*ux.^2.*uy + C.*uy, ...
                      B.*ux.*uy.^2 + C.*ux, ...
                      B.*uy.^3 + 3*C.*uy);
end
