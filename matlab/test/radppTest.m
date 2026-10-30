% Verify kernel3d.radpp: the fit, the tensor evaluation, the built-in
% kernels and the near quadrature.

% run ../startup.m

rng(42);

test_fit_limits();
test_fit_battery();
test_radpp_eval_derivs();
% 2D (in-plane derivatives, dim = 2)
test_tensors_2d();
test_builtins_2d(@gauss, 1e-11);
test_builtins_2d(@(s,t) helmgreen(1.1, s, t), 1e-10);
test_full_2d_format();
test_helm_sp_eval_2d();
test_src_eval_2d();
test_getquad_rows(2);

% 3D (full derivatives, dim = 3, the default)
test_tensors_3d();
test_builtins_3d();
test_getquad_sp_3d();
test_src_eval_3d();
test_getquad_src_3d();
test_grad_3d();
test_hess_3d();
test_trans_3d();
test_getquad_rows(3);

% value kernels and interfaces
test_getquad_lap();
test_kernel3d_eval();
test_kernel3d_getquad_helm();
test_radcheb_interface();
test_aliases();
test_bad_inputs();

% kernel3d(@green, type): exact eval, radpp quadrature
test_radgreen_checks();
test_radgreen_ipv();
test_radgreen_order();
test_radgreen_eval();
test_radgreen_arith();
test_radgreen_getquad();
test_radgreen_rlb_getquad();
test_radgreen_mfie();
test_radgreen_fmm();
test_radgreen_2d();
test_fmm_option();
test_grad_ntrace();
test_radpp_sum();

% beyond helmholtz on one surface: laplace, yukawa, two surfers, custom
% source weights
test_real_and_decaying();
test_two_surfers();
test_srcwfun();



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

% noise in f below the tolerance: a warning, not maxsub
f = @(r) exp(1i*100*r)./r;
lastwarn('');
warning('off', 'RADPP_FIT:noise');
[~, ~, info] = radpp_fit(f, 1e-3, 10, struct('eps', 1e-15));
warning('on', 'RADPP_FIT:noise');
[~, id] = lastwarn();
assert(strcmp(id, 'RADPP_FIT:noise'), 'radpp_fit: no noise warning');
assert(info.err < 1e-11, 'radpp_fit: noise level %.1e', info.err);

% a green plus 1e-6 noise stalls at the noise level
s0 = rng(7);
f = @(r) [exp(1i*r)./(1 + r), r.^2.*log(r), r.*log(r)] + 1e-6*rand(numel(r), 3);
lastwarn('');
warning('off', 'RADPP_FIT:noise');
[~, ~, info] = radpp_fit(f, 1e-9, 3, struct('eps', 1e-13));
warning('on', 'RADPP_FIT:noise');
rng(s0);
[~, id] = lastwarn();
assert(strcmp(id, 'RADPP_FIT:noise'), 'radpp_fit: noisy green, no noise warning');
assert(info.err > 1e-8 && info.err < 1e-4, 'radpp_fit: noisy green, error %.1e', info.err);
assert(info.nbin < 300, 'radpp_fit: noisy green, %d panels', info.nbin);

f = @(r) flexprof([1.7, 0.5i], r);
lastwarn('');
[~, ~, info] = radpp_fit(f, 1e-9, 3, struct('eps', 1e-13));
assert(isempty(lastwarn()), 'radpp_fit: r log r warned');
assert(info.err < 1e-12, 'radpp_fit: r log r, error %.1e', info.err);

end


function test_radpp_eval_derivs()
% radpp_eval_derivs vs radpp_rotate(radpp_eval), dim = 2, 3, 13.

f = @(r) exp(1i*1.3*r(:))./r(:).*(1:6) + cos((1:6).*r(:));
[ipars, dpars, info] = radpp_fit(f, 1e-3, 5);
n = 500;
d = randn(3, n);
d = d./vecnorm(d).*(1e-3 + 4*rand(1, n));
d(:,1) = [0.5; 0; 0];                        % on the x axis
for dim = [2 3 13]
    nfits = 0:3;
    if dim == 13, nfits = 0:1; end
    dz = d(3,:);
    if dim == 2, dz = 0*dz; end
    r = sqrt(d(1,:).^2 + d(2,:).^2 + dz.^2);
    prof = radpp_eval(r, info.breaks, info.coefs);
    for nfit = nfits
        ip = [ipars(:); nfit; 1; dim; 0];
        A = radpp_rotate(d(1,:), d(2,:), dz, prof, nfit, dim);
        B = radpp_eval_derivs(d(1,:), d(2,:), dz, ip, dpars);
        assert(isequal(size(A), size(B)), ...
            'radpp_eval_derivs: size, dim %d, nfit %d', dim, nfit);
        err = max(abs(A(:) - B(:)))/max(abs(A(:)));
        assert(err < 1e-13, 'radpp_eval_derivs: dim %d, nfit %d, %.2e', ...
            dim, nfit, err);
    end
end

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
C = tens(dx, dy, 0*dx, p, 3, 2);

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
C = tens(X(:,1), X(:,2), X(:,3), p, 3);

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
C3 = tens(X(:,1), X(:,2), 0*X(:,1), p, 3, 3);
C2 = tens(X(:,1), X(:,2), 0*X(:,1), K2.params, 3, 2);
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
S = neighbours(geometries.sphere(1, 1, [0;0;0], 4, 1), 3);

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
S = neighbours(geometries.sphere(1, 1, [0;0;0], 4, 1), 3);
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

% a green in the helm3d.green format, full hess
K3 = kernel3d.radpp(@(s,t) helm3d.green(zk, s, t), [1e-2 6], 'cp');
src.n = randn(3, ns);  src.n = src.n./vecnorm(src.n);
targ.n = randn(3, nt);  targ.n = targ.n./vecnorm(targ.n);
B = kernel3d.helm3d('cp', zk, [1; 1]).eval(src, targ);
err = norm(K3.eval(src, targ) - B, 'fro')/norm(B, 'fro');
assert(err < 1e-10, 'radpp helm3d.green format: %.2e', err);

% laplace is real
Kl = kernel3d.radpp(kernel3d.lap3d('s'), [1e-2 6]);
assert(~Kl.ifcomplex, 'radpp: laplace kernel3d flagged complex');

end


function test_kernel3d_getquad_helm()
% Near quadrature of the radpp form of helm3d('s') matches helm3d, on
% surface and for nearby off-surface targets.

zk = 1.0;
eps = 1e-7;
S = neighbours(geometries.sphere(1, 1, [0;0;0], 4, 1), 3);
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


function test_src_eval_2d()
% Source derivative kernels in the plane: 'd' = G_{n_y} = -grad G . n_y
% and 'dp' = G_{n_x n_y} = -n_x' H n_y, from helmgreen.

zk = 1.1;
helm = @(s,t) helmgreen(zk, s, t);
ns = 30;  nt = 40;
ths = 2*pi*rand(1, ns);  tht = 2*pi*rand(1, nt);
src.r  = 0.5*[cos(ths); sin(ths)];  src.n  = [cos(ths); sin(ths)];
targ.r = 1.7*[cos(tht); sin(tht)];  targ.n = [cos(tht); sin(tht)];

[~, g, h] = helm(src.r, targ.r);
nx = targ.n;  ny = src.n;
D  = -(g(:,:,1).*ny(1,:) + g(:,:,2).*ny(2,:));
DP = -(h(:,:,1).*nx(1,:).'.*ny(1,:) + h(:,:,2).*(nx(1,:).'.*ny(2,:) ...
       + nx(2,:).'.*ny(1,:)) + h(:,:,3).*nx(2,:).'.*ny(2,:));

ref = struct('d', D, 'dp', DP);
types = fieldnames(ref);
for k = 1:numel(types)
    K = kernel3d.radpp(helm, [1e-6 6], types{k}, 'dim', 2);
    A = K.eval(src, targ);
    B = ref.(types{k});
    err = norm(A - B, 'fro')/norm(B, 'fro');
    assert(err < 1e-10, 'radpp %s 2d: %.2e', types{k}, err);
end

end


function test_src_eval_3d()
% Source derivative kernels match helm3d for points in full 3D.

zk = 1.1;  coefs = [0.7 -1.3];
helm = @(s,t) helmgreen(zk, s, t);
ns = 50;  nt = 60;
src.r  = randn(3, ns);  src.r  = src.r./vecnorm(src.r);
src.n  = randn(3, ns);  src.n  = src.n./vecnorm(src.n);
targ.r = randn(3, nt);  targ.r = targ.r./vecnorm(targ.r)*1.7;
targ.n = randn(3, nt);  targ.n = targ.n./vecnorm(targ.n);

types = {'d', 'c', 'dp', 'cp', 's2trans', 'd2trans', 'c2trans'};
for k = 1:numel(types)
    ca = cfs(types{k}, coefs);
    K  = kernel3d.radpp(helm, [1e-6 6], types{k}, ca{:});
    Ke = kernel3d.helm3d(types{k}, zk, coefs(:));
    A = K.eval(src, targ);
    B = Ke.eval(src, targ);
    err = norm(A - B, 'fro')/norm(B, 'fro');
    assert(err < 1e-10, 'radpp %s 3d: %.2e', types{k}, err);
end

end


function test_getquad_src_3d()
% Near quadrature of the source derivative kernels matches helm3d on a
% sphere, including the hypersingular 'dp' (ipv = 2).

zk = 1.1;  coefs = [0.7 -1.3];
eps = 1e-7;
S = neighbours(geometries.sphere(1, 1, [0;0;0], 4, 1), 3);
helm = @(s,t) helmgreen(zk, s, t);

types = {'d', 'c', 'dp', 'cp', 'c2trans'};
for k = 1:numel(types)
    ca = cfs(types{k}, coefs);
    K  = kernel3d.radpp(helm, [1e-9 4], types{k}, ca{:}, ...
                        'opts', struct('eps', 1e-13));
    Ke = kernel3d.helm3d(types{k}, zk, coefs(:));
    Q  = K.getquad(S, eps);
    Qe = Ke.getquad(S, eps);
    err = norm(Q - Qe, 'fro')/norm(Qe, 'fro');
    assert(err < 1e-6, 'radpp getquad %s 3d: %.2e', types{k}, err);
end

end


function test_grad_3d()
% 'grad' is grad_x G, and 'dgrad' is grad_x of 'd' (checked by finite
% differences in the target), and 'cgrad' and 'trep_grad' are built from
% them.

zk = 1.1;
helm = @(s,t) helmgreen(zk, s, t);
ns = 20;  nt = 25;
src.r  = randn(3, ns);  src.r  = src.r./vecnorm(src.r);
src.n  = randn(3, ns);  src.n  = src.n./vecnorm(src.n);
targ.r = randn(3, nt);  targ.r = targ.r./vecnorm(targ.r)*1.7;

Kg = kernel3d.radpp(helm, [1e-6 6], 'grad');
Kd = kernel3d.radpp(helm, [1e-6 6], 'd');
Kdg = kernel3d.radpp(helm, [1e-6 6], 'dgrad');
assert(isequal(Kg.opdims, [3 1]), 'radpp grad: opdims');

Ks = kernel3d.helm3d('s', zk);
G = Kg.eval(src, targ);
DG = Kdg.eval(src, targ);
h = 1e-5;
for j = 1:3
    tp = targ;  tm = targ;
    tp.r(j,:) = tp.r(j,:) + h;  tm.r(j,:) = tm.r(j,:) - h;
    fs = (Ks.eval(src, tp) - Ks.eval(src, tm))/(2*h);
    fd = (Kd.eval(src, tp) - Kd.eval(src, tm))/(2*h);
    e1 = norm(G(j:3:end,:) - fs, 'fro')/norm(fs, 'fro');
    e2 = norm(DG(j:3:end,:) - fd, 'fro')/norm(fd, 'fro');
    assert(e1 < 1e-8, 'radpp grad 3d: %.2e', e1);
    assert(e2 < 1e-8, 'radpp dgrad 3d: %.2e', e2);

% 'cgrad' and 'trep_grad' from 'grad' and 'dgrad'
c = [0.7 2i];
A = kernel3d.radpp(helm, [1e-6 6], 'cgrad', c).eval(src, targ);
err = norm(A - c(1)*G - c(2)*DG, 'fro')/norm(A, 'fro');
assert(err < 1e-12, 'radpp cgrad 3d: %.2e', err);
A = kernel3d.radpp(helm, [1e-6 6], 'trep_grad', c).eval(src, targ);
err = max(norm(A(:,1:2:end) - c(1)*DG, 'fro'), norm(A(:,2:2:end) - c(2)*G, 'fro'))/norm(A, 'fro');
assert(err < 1e-12, 'radpp trep_grad 3d: %.2e', err);
end

end


function test_hess_3d()
% 'hess' is the Hessian of G: finite differences of 'grad'.

helm = @(s,t) helmgreen(1.1, s, t);
src.r  = randn(3, 20);  src.r  = src.r./vecnorm(src.r);
targ.r = randn(3, 25);  targ.r = targ.r./vecnorm(targ.r)*1.7;

Kg = kernel3d.radpp(helm, [1e-6 6], 'grad');
Kh = kernel3d.radpp(helm, [1e-6 6], 'hess');
assert(isequal(Kh.opdims, [6 1]), 'radpp hess: opdims');
H = Kh.eval(src, targ);

h = 1e-5;
ih = [1 2 3; 2 4 5; 3 5 6];        % (i,j) -> xx xy xz yy yz zz
for j = 1:3
    tp = targ;  tm = targ;
    tp.r(j,:) = tp.r(j,:) + h;  tm.r(j,:) = tm.r(j,:) - h;
    fg = (Kg.eval(src, tp) - Kg.eval(src, tm))/(2*h);
    for i = 1:3
        f = fg(i:3:end,:);
        err = norm(H(ih(i,j):6:end,:) - f, 'fro')/norm(f, 'fro');
        assert(err < 1e-8, 'radpp hess 3d: %.2e', err);
    end
end

end


function test_trans_3d()
% 'trans_sys' (2x2), 'trans_rep' and 'trep_prime' (1x2) match helm3d, eval
% and near quadrature on a sphere.

zk = 1.1;
C = [0.7 -1.3; 0.4 2.1];  c = [0.9 -0.6];
helm = @(s,t) helmgreen(zk, s, t);
fopts = struct('eps', 1e-13);

src.r  = randn(3, 30);  src.r  = src.r./vecnorm(src.r);
src.n  = randn(3, 30);  src.n  = src.n./vecnorm(src.n);
targ.r = randn(3, 35);  targ.r = targ.r./vecnorm(targ.r)*1.7;
targ.n = randn(3, 35);  targ.n = targ.n./vecnorm(targ.n);

S = neighbours(geometries.sphere(1, 1, [0;0;0], 4, 1), 3);
cases = {{'trans_sys', C}, {'trans_rep', c}, {'trep_prime', c}};
for k = 1:numel(cases)
    [type, cf] = cases{k}{:};
    K  = kernel3d.radpp(helm, [1e-9 4], type, 'coefs', cf, 'opts', fopts);
    Ke = kernel3d.helm3d(type, zk, cf);
    assert(isequal(K.opdims, Ke.opdims), 'radpp %s: opdims', type);

    A = K.eval(src, targ);  B = Ke.eval(src, targ);
    err = norm(A - B, 'fro')/norm(B, 'fro');
    assert(err < 1e-10, 'radpp %s eval: %.2e', type, err);

    Q = K.getquad(S, 1e-7);  Qe = Ke.getquad(S, 1e-7);
    err = norm(Q - Qe, 'fro')/norm(Qe, 'fro');
    assert(err < 1e-6, 'radpp getquad %s: %.2e', type, err);
end

end


function test_radgreen_checks()
% Declare time errors: too few derivatives, not radial, wrong derivative.

zk = 1.1;
g = @(s,t) helm3d.green(zk, s, t);
bad = {@() kernel3d(@(s,t) helm3d.green(zk,s,t) + 0, 'sp'), ...   % 1 output
       @() kernel3d(@(s,t) notradial(s, t), 's'), ...
       @() kernel3d(@(s,t) wronghess(zk, s, t), 'dp')};
for k = 1:numel(bad)
    caught = false;
    try
        bad{k}();
    catch
        caught = true;
    end
    assert(caught, 'radgreen: bad green %d not rejected', k);
end

% kernel_order = max(ipv)-1
assert(kernel3d(g, 's').kernel_order == -1, 'radgreen: kernel_order s');
assert(kernel3d(g, 'cp').kernel_order == 1, 'radgreen: kernel_order cp');
assert(kernel3d(g, 'cp', 'ipv', 0).kernel_order == -1, ...
       'radgreen: kernel_order from ipv');

end


function test_radgreen_ipv()
% ipv measured from each kernel: helmholtz gives the helm3d choices, log r
% is weakly singular, an ipv given wins.

zk = 1.1;
g = @(s,t) helm3d.green(zk, s, t);
cases = {'s', 0; 'sp', 1; 'd', 1; 'c', 1; 'dp', 2; 'cp', 2; ...
         'clamped_bc', [0 1]; 'trans_sys', [1 2 0 1]; 'hess', 2};
for k = 1:size(cases, 1)
    K = kernel3d(g, cases{k,1});
    assert(isequal(K.params.ipv, repmat(cases{k,2}, 1, ...
           numel(K.params.ipv)/numel(cases{k,2}))), ...
           'radgreen ipv %s: %s', cases{k,1}, mat2str(K.params.ipv));
end

% log r: weakly singular
assert(kernel3d(@logg, 's').params.ipv == 0, 'radgreen: log s');
assert(kernel3d(@logg, 'sp').params.ipv == 0, 'radgreen: log sp');
assert(kernel3d(@logg, 'dp').params.ipv == 1, 'radgreen: log dp');
assert(kernel3d(@logg, 'dp', 'dim', 2).params.ipv == 1, 'radgreen: log dp 2d');

% 1/r^2, smooth, and a two-wavenumber difference (bounded at 0)
inv2 = @(s,t) 1./((t(1,:).'-s(1,:)).^2 + (t(2,:).'-s(2,:)).^2 + ...
                  (t(3,:).'-s(3,:)).^2);
assert(kernel3d(inv2, 's').params.ipv == 1, 'radgreen: 1/r^2 s');
assert(kernel3d(@gaussg, 'dp').params.ipv == 0, 'radgreen: smooth dp');
assert(kernel3d(@gaussg, 's').params.ipv == 0, 'radgreen: smooth s');
gd = @(s,t) helm3d.green(1.1, s, t) - helm3d.green(2.3, s, t);
assert(kernel3d(gd, 's').params.ipv == 0, 'radgreen: diff s');

% ipv given
assert(kernel3d(g, 'sp', 'ipv', 0).params.ipv == 0, 'radgreen: ipv given');
assert(kernel3d(@logg, 'dp', 'ipv', 2, 'dim', 2).params.ipv == 2, ...
       'radgreen: ipv given 2d');

end


function test_radgreen_order()
% ipv at build time for differences of
% helmholtz kernels (bounded at 0) and the Laplace- and Helmholtz-
% Beltrami kernels (log), with eval checked against helm3d and belpde.

k0 = 1.1;  k1 = 2.3;
src.r  = randn(3, 30);  src.n  = randn(3, 30);  src.n  = src.n./vecnorm(src.n);
targ.r = randn(3, 20) + [3;0;0];
targ.n = randn(3, 20);  targ.n = targ.n./vecnorm(targ.n);
targ.mean_curv = randn(20, 1);

% differences: G_k0 - G_k1 is bounded, so one order weaker than helm3d
gd = @(s,t) diffg(k0, k1, s, t);
cases = {'s', 0; 'sp', 0; 'd', 0; 'dp', 0};    % dp_diff ~ 1/r, isotropic
for k = 1:size(cases, 1)
    type = cases{k,1};
    K = kernel3d(gd, type);
    assert(K.params.ipv == cases{k,2}, 'radgreen ipv %s_diff: %d', type, ...
           K.params.ipv);
    B = kernel3d.helm3d([type '_diff'], [k0; k1]).eval(src, targ);
    err = norm(K.eval(src, targ) - B, 'fro')/norm(B, 'fro');
    assert(err < 1e-12, 'radgreen %s_diff eval: %.2e', type, err);
end

% klb and khb: log singular, weakly singular for 's'
kb = 1.3;
klb = @(s,t) belpde.kern(0, struct('r', s), struct('r', t), 'klb');
khb = @(s,t) belpde.kern(kb^2, struct('r', s), struct('r', t), 'khb');
for c = {{klb, kernel3d.belpde('klb')}, {khb, kernel3d.belpde('khb', kb)}}
    [g, Kb] = c{1}{:};
    K = kernel3d(g, 's');
    assert(K.params.ipv == 0, 'radgreen: %s not weakly singular', Kb.type);
    B = Kb.eval(src, targ);
    err = norm(K.eval(src, targ) - B, 'fro')/norm(B, 'fro');
    assert(err < 1e-13, 'radgreen %s eval: %.2e', Kb.type, err);
end

% rlb = (Laplacian - d_nn - 2 H d_n) klb, as custom weights on klb
K = kernel3d(@klbg, @w_rlb, 'nderiv', 2, 'targ_fields', {'n', 'mean_curv'});
assert(K.params.ipv == 0, 'radgreen: rlb ipv %d', K.params.ipv);     % bounded on the surface
assert(K.kernel_order == -1, 'radgreen: rlb kernel_order');
B = kernel3d.belpde('rlb').eval(src, targ);
err = norm(K.eval(src, targ) - B, 'fro')/norm(B, 'fro');
assert(err < 1e-12, 'radgreen rlb eval: %.2e', err);

end


function varargout = diffg(k0, k1, s, t)
% G_k0 - G_k1 and its derivatives, helm3d.green format
a = cell(1, nargout);  b = a;
[a{:}] = helm3d.green(k0, s, t);
[b{:}] = helm3d.green(k1, s, t);
varargout = cellfun(@minus, a, b, 'UniformOutput', false);
end


function varargout = klbg(s, t)
% klb = log(r^2)/(4 pi) with derivatives, helm3d.green format
x = t(1,:).' - s(1,:);  y = t(2,:).' - s(2,:);  z = t(3,:).' - s(3,:);
r2 = x.^2 + y.^2 + z.^2;
X = cat(3, x, y, z);
varargout{1} = log(r2)/(4*pi);
varargout{2} = X./r2/(2*pi);
H = zeros([size(x) 3 3]);
for i = 1:3, for j = 1:3
    H(:,:,i,j) = ((i==j)./r2 - 2*X(:,:,i).*X(:,:,j)./r2.^2)/(2*pi);
end, end
varargout{3} = H;
end


function W = w_rlb(t)
% weights of (Laplacian - d_nn - 2 H d_n) on [G | grad | xx xy xz yy yz zz]
n = t.n.';  H = t.mean_curv(:);
W = [0*H, -2*H.*n, 1 - n(:,1).^2, -2*n(:,1).*n(:,2), -2*n(:,1).*n(:,3), ...
     1 - n(:,2).^2, -2*n(:,2).*n(:,3), 1 - n(:,3).^2];
end


function varargout = logg(s, t)
% G = log r with target derivatives, points (2,n) or (3,n)
d = size(t, 1);
X = t(1,:).' - s(1,:);
for i = 2:d, X = cat(3, X, t(i,:).' - s(i,:)); end
x = X(:,:,1);
r2 = sum(X.^2, 3);
varargout{1} = log(r2)/2;
varargout{2} = X./r2;
H = zeros([size(x) d d]);
for i = 1:d, for j = 1:d
    H(:,:,i,j) = (i==j)./r2 - 2*X(:,:,i).*X(:,:,j)./r2.^2;
end, end
varargout{3} = H;
end


function varargout = gaussg(s, t)
% G = exp(-r^2/2) with target derivatives, helm3d.green format
x = t(1,:).' - s(1,:);  y = t(2,:).' - s(2,:);  z = t(3,:).' - s(3,:);
X = cat(3, x, y, z);
G = exp(-(x.^2 + y.^2 + z.^2)/2);
varargout{1} = G;
varargout{2} = -X.*G;
H = zeros([size(x) 3 3]);
for i = 1:3, for j = 1:3
    H(:,:,i,j) = (X(:,:,i).*X(:,:,j) - (i==j)).*G;
end, end
varargout{3} = H;
end


function test_radgreen_eval()
% Eval is exact: matches helm3d for every type to rounding.

zk = 1.1;  coefs = [0.7 -1.3];
g = @(s,t) helm3d.green(zk, s, t);
src.r  = randn(3, 40);  src.n  = randn(3, 40);  src.n  = src.n./vecnorm(src.n);
targ.r = randn(3, 30) + [3;0;0];
targ.n = randn(3, 30);  targ.n = targ.n./vecnorm(targ.n);

types = {'s', 'd', 'c', 'sp', 'dp', 'cp'};
for k = 1:numel(types)
    ca = cfs(types{k}, coefs);
    K  = kernel3d(g, types{k}, ca{:});
    Ke = kernel3d.helm3d(types{k}, zk, coefs(:));
    A = K.eval(src, targ);  B = Ke.eval(src, targ);
    err = norm(A - B, 'fro')/norm(B, 'fro');
    assert(strcmp(K.name, 'radgreen'), 'radgreen: name');
    assert(err < 1e-13, 'radgreen eval %s: %.2e', types{k}, err);
end

end


function test_radgreen_arith()
% Kernels with the same green combine into one green kernel: scaling, sums
% (S + D = c) and interleave (trans_sys). The same type on different greens
% adds the expansions; otherwise different greens fall back to the generic
% arithmetic.

zk = 1.1;
g = @(s,t) helm3d.green(zk, s, t);
src.r  = randn(3, 40);  src.n  = randn(3, 40);  src.n  = src.n./vecnorm(src.n);
targ.r = randn(3, 30) + [3;0;0];
targ.n = randn(3, 30);  targ.n = targ.n./vecnorm(targ.n);

Ks = kernel3d(g, 's');  Kd = kernel3d(g, 'd');
Ksp = kernel3d(g, 'sp');  Kdp = kernel3d(g, 'dp');

K = 2*Ks - 3*Kd;
assert(strcmp(K.name, 'radgreen'), 'radgreen: sum not merged');
B = kernel3d.helm3d('c', zk, [2; -3]).eval(src, targ);
err = norm(K.eval(src, targ) - B, 'fro')/norm(B, 'fro');
assert(err < 1e-13, 'radgreen 2s - 3d: %.2e', err);
assert(K.params.isrc == 1 && K.params.nderiv == 0, 'radgreen: c spec');

K = Ksp + 0.5*Kdp;                 % mixed derivative orders
B = kernel3d.helm3d('cp', zk, [1; 0.5]).eval(src, targ);
err = norm(K.eval(src, targ) - B, 'fro')/norm(B, 'fro');
assert(strcmp(K.name, 'radgreen') && err < 1e-13, 'radgreen cp sum: %.2e', err);

K = -Ks;
err = norm(K.eval(src, targ) + Ks.eval(src, targ), 'fro');
assert(strcmp(K.name, 'radgreen') && err == 0, 'radgreen: uminus');

C = [0.7 -1.3; 0.4 2.1];
K = kernel3d.interleave([C(1,1)*Kd, C(1,2)*Ks; C(2,1)*Kdp, C(2,2)*Ksp]);
assert(strcmp(K.name, 'radgreen') && isequal(K.opdims, [2 2]), ...
       'radgreen: interleave not merged');
Kt = kernel3d(g, 'trans_sys', 'coefs', C);
B = Kt.eval(src, targ);
err = norm(K.eval(src, targ) - B, 'fro')/norm(B, 'fro');
assert(err < 1e-13, 'radgreen interleave: %.2e', err);
assert(isequal(K.params.ipv, Kt.params.ipv), 'radgreen: interleave ipv');

% fmms are carried through the merges (stand-in fmm: dense eval)
F = @(K) @(eps, s, t, sig) K.eval(s, t)*sig;
Ks.fmm = F(Ks);  Kd.fmm = F(Kd);  Ksp.fmm = F(Ksp);  Kdp.fmm = F(Kdp);
for K = {2*Ks - 3*Kd, ...
         kernel3d.interleave([C(1,1)*Kd, C(1,2)*Ks; C(2,1)*Kdp, C(2,2)*Ksp])}
    K = K{1};
    sig = randn(K.opdims(2)*size(src.r, 2), 1);
    err = norm(K.fmm(1e-7, src, targ, sig) - K.eval(src, targ)*sig);
    assert(strcmp(K.name, 'radgreen') && err < 1e-12, 'radgreen merged fmm: %.2e', err);
end

k2 = 2.3;
g2 = @(s,t) helm3d.green(k2, s, t);            % a different green
h = @(type, k) kernel3d.helm3d(type, k).eval(src, targ);
K = 2*Ks + kernel3d(g2, 's');
assert(strcmp(K.name, 'radgreen') && strcmp(K.type, 's') && ...
       numel(K.params.green.parts) == 2, 'radgreen: s + s not merged');
B = 2*h('s', zk) + h('s', k2);
err = norm(K.eval(src, targ) - B, 'fro')/norm(B, 'fro');
assert(err < 1e-13, 'radgreen 2 s(k1) + s(k2): %.2e', err);

% ipv and kernel_order of a sum are the worse of the two
K = Ksp - kernel3d(g2, 'sp');
assert(strcmp(K.name, 'radgreen') && K.params.ipv == 1 && ...
       K.kernel_order == 0, 'radgreen: sp - sp ipv');
B = h('sp', zk) - h('sp', k2);
err = norm(K.eval(src, targ) - B, 'fro')/norm(B, 'fro');
assert(err < 1e-13, 'radgreen sp(k1) - sp(k2): %.2e', err);

K = Ks + kernel3d(g2, 'd');
assert(~strcmp(K.name, 'radgreen'), 'radgreen: s(g1) + d(g2) merged');

% separately created handles: merged when the text and captured values
% agree, not otherwise
z = zk;   K1 = kernel3d(@(s,t) helm3d.green(z, s, t), 's');
K2 = kernel3d(@(s,t) helm3d.green(z, s, t), 'd');
K = K1 + K2;
assert(strcmp(K.name, 'radgreen') && numel(K.params.green.parts) == 1, ...
       'radgreen: equal anonymous greens not merged');
z = k2;   K3 = kernel3d(@(s,t) helm3d.green(z, s, t), 'd');
K = K1 + K3;
assert(~strcmp(K.name, 'radgreen'), 'radgreen: different captured zk merged');
B = h('s', zk) + h('d', k2);
err = norm(K.eval(src, targ) - B, 'fro')/norm(B, 'fro');
assert(err < 1e-13, 'radgreen s(zk) + d(k2): %.2e', err);

end


function test_radgreen_getquad()
% The automatic range quadrature matches helm3d on and off the sphere,
% for a single kernel and merged sums, on one green or two.

zk = 1.1;  eps = 1e-7;
g = @(s,t) helm3d.green(zk, s, t);
S = neighbours(geometries.sphere(1, 1, [0;0;0], 4, 1), 3);
targ.r = randn(3, 40);  targ.r = targ.r./vecnorm(targ.r)*1.05;
targ.n = targ.r./vecnorm(targ.r);

cases = {{kernel3d(g, 's'),  kernel3d.helm3d('s', zk)}, ...
         {kernel3d(g, 'cp'), kernel3d.helm3d('cp', zk, [1; 1])}, ...
         {kernel3d(g, 's') + 2*kernel3d(g, 'd'), ...
          kernel3d.helm3d('c', zk, [1; 2])}, ...
         {kernel3d(g, 'dp') - kernel3d(@(s,t) helm3d.green(2.3, s, t), 'dp'), ...
          kernel3d.helm3d('dp', zk) - kernel3d.helm3d('dp', 2.3)}, ...
         {2*kernel3d(g, 'sp'), 2*kernel3d.helm3d('sp', zk)}};
for k = 1:numel(cases)
    [K, Ke] = cases{k}{:};
    for t = {S, targ}
        Q = K.getquad(S, eps, t{1});  Qe = Ke.getquad(S, eps, t{1});
        err = norm(Q - Qe, 'fro')/norm(Qe, 'fro');
        assert(err < 1e-6, 'radgreen getquad case %d: %.2e', k, err);
    end
end

end


function test_radgreen_rlb_getquad()
% rlb as (Laplacian - d_nn - 2 H d_n) klb with weights on the target
% fields n and mean_curv: near quadrature matches belpde on a sphere.

eps = 1e-7;
S = neighbours(geometries.sphere(1, 1, [0;0;0], 4, 1), 3);
K  = kernel3d(@klbg, @w_rlb, 'nderiv', 2, 'targ_fields', {'n', 'mean_curv'});
Kb = kernel3d.belpde('rlb');

Q  = K.getquad(S, eps);
Qb = Kb.getquad(S, eps);
err = norm(Q - Qb, 'fro')/norm(Qb, 'fro');
assert(err < 1e-6, 'radgreen rlb getquad: %.2e', err);

% one curved patch: getquad against the smooth integral
for d = [0.5 1]
    err = kernel3d.check_quad(K, eps, struct('dist', d));
    assert(err < 1e-6, 'radgreen rlb check, dist %g: %.2e', d, err);
end
g = @(s,t) helm3d.green(1.1, s, t);
for type = {'s', 'sp', 'dp', 'cp', 'trans_sys'}
    err = kernel3d.check_quad(kernel3d(g, type{1}), eps);
    assert(err < 1e-6, 'radgreen %s check: %.2e', type{1}, err);
end

end


function test_radgreen_mfie()
% 'mfie' = n x curl S_k[J] and 'curls' = curl S_k[J] on Cartesian
% densities match the formulas from helm3d.green, and the MFIE near
% quadrature matches a smooth integral on a curved patch.

zk = 1.5;
g = @(s,t) helm3d.green(zk, s, t);
src.r  = randn(3, 20);
targ.r = randn(3, 15) + [3;0;0];
targ.n = randn(3, 15);  targ.n = targ.n./vecnorm(targ.n);
[~, gr] = helm3d.green(zk, src.r, targ.r);         % (nt, ns, 3)

ns = 20;  nt = 15;
B = zeros(3, nt, 3, ns);  C = B;
% (n x (g x J))_i = g_i (n.J) - J_i (n.g)
ng = zeros(nt, ns);
for b = 1:3, ng = ng + targ.n(b,:).'.*gr(:,:,b); end
for i = 1:3
    for c = 1:3
        B(i,:,c,:) = reshape(targ.n(c,:).'.*gr(:,:,i) - (i==c)*ng, 1, nt, 1, ns);
        k = 6 - i - c;
        if i ~= c
            C(i,:,c,:) = reshape((k-i)*(c-i)*(c-k)/2*gr(:,:,k), 1, nt, 1, ns);
        end
    end
end
B = reshape(B, 3*nt, 3*ns);  C = reshape(C, 3*nt, 3*ns);

Km = kernel3d(g, 'nxcurls');
Kc = kernel3d(g, 'curls');
assert(isequal(Km.opdims, [3 3]) && strcmp(Km.type, 'mfie'), 'radgreen: mfie');
err = norm(Km.eval(src, targ) - B, 'fro')/norm(B, 'fro');
assert(err < 1e-13, 'radgreen mfie eval: %.2e', err);
err = norm(Kc.eval(src, targ) - C, 'fro')/norm(C, 'fro');
assert(err < 1e-13, 'radgreen curls eval: %.2e', err);

err = kernel3d.check_quad(Km, 1e-7);
assert(err < 1e-6, 'radgreen mfie check: %.2e', err);

end


function test_radgreen_fmm()
% A radgreen kernel with an fmm attached: em3d.fmm 'mfie' and 'curls'
% match eval, and surfermatapply and surferkerneval with the fmm match
% the dense matrix and the fmm-free evaluation.

zk = 1.5;  eps = 1e-8;
g = @(s,t) helm3d.green(zk, s, t);
Km = kernel3d(g, 'mfie');  Kc = kernel3d(g, 'curls');

% fmm types against eval, well-separated points
ns = 10000;
src.r = randn(3, ns);
targ.r = randn(3, 200) + [5;0;0];
targ.n = randn(3, 200);  targ.n = targ.n./vecnorm(targ.n);
sig = complex(randn(3*ns, 1), randn(3*ns, 1));
for c = {{Km, 'mfie'}, {Kc, 'curls'}}
    [K, type] = c{1}{:};
    A = K.eval(src, targ)*sig;
    B = em3d.fmm(eps, zk, src, targ, type, sig);
    err = norm(A - B)/norm(A);
    assert(err < 1e-6, 'em3d.fmm %s: %.2e', type, err);
end

Km.fmm = @(eps, s, t, sigma) em3d.fmm(eps, zk, s, t, 'mfie', sigma);
Kc.fmm = @(eps, s, t, sigma) em3d.fmm(eps, zk, s, t, 'curls', sigma);
S = neighbours(geometries.sphere(1, 1, [0;0;0], 3, 1), 3);       % small: 3x3 kernel
J = complex(randn(3*S.npts, 1), randn(3*S.npts, 1));

% surfermatapply with the fmm against the dense matrix
K = kernel3d.eye(0.5*eye(3)) - Km;
A = surfermat(S, K, eps)*J;
B = surfermatapply(S, K, J, eps);
err = norm(A - B)/norm(A);
assert(err < 1e-6, 'radgreen surfermatapply fmm: %.2e', err);

% surferkerneval with and without the fmm, including near targets
t2.r = randn(3, 10000);  t2.r = t2.r./vecnorm(t2.r).*(1.02 + rand(1, 10000));
A = surferkerneval(S, Kc, J, t2, eps, struct('usefmm', 0));
B = surferkerneval(S, Kc, J, t2, eps);
err = norm(A - B)/norm(A);
assert(err < 1e-6, 'radgreen surferkerneval fmm: %.2e', err);

end


function test_fmm_option()
% The 'fmm' option of radpp and radgreen (RADPP_FMM) against eval, with
% two densities, the generic mfie and curls against em3d.fmm, and the fmm
% of sums and interleaved kernels.

zk = 1.3;  ndens = 2;
g = @(s,t) helm3d.green(zk, s, t);
fmm = @(eps,s,pg,t,pgt) hfmm3d(eps, zk, s, pg, t, pgt);
src.r = randn(3, 200);  src.n = randn(3, 200);  src.n = src.n./vecnorm(src.n);
targ.r = randn(3, 150) + [5;0;0];
targ.n = randn(3, 150);  targ.n = targ.n./vecnorm(targ.n);

types = {'s','d','c','sp','dp','cp','grad','dgrad','trans_sys', ...
    'trans_rep','clamped_bc','d2trans','c2trans','mfie','curls'};
for i = 1:numel(types)
    a = targs3(types{i});
    Ks = {kernel3d.radpp(g, [1e-6 20], types{i}, a{:}, 'fmm', fmm), ...
          kernel3d(g, types{i}, a{:}, 'fmm', fmm)};
    for k = 1:2
        A = Ks{k}.eval(src, targ);
        sig = complex(randn(size(A,2), ndens), randn(size(A,2), ndens));
        u = Ks{k}.fmm(1e-12, src, targ, sig);
        err = norm(u - A*sig, 'fro')/norm(A*sig, 'fro');
        assert(err < 1e-10, 'fmm option %s (%d): %.2e', types{i}, k, err);
    end
end

% in-plane derivatives, dim = 2
s2.r = src.r(1:2,:);  s2.n = src.n(1:2,:)./vecnorm(src.n(1:2,:));
t2.r = targ.r(1:2,:);  t2.n = targ.n(1:2,:)./vecnorm(targ.n(1:2,:));
g2 = @(s,t) helmgreen(zk, s, t);
for type = {'s','dp','cp','grad','clamped_bc','trans_sys'}
    a = targs3(type{1});
    K = kernel3d.radpp(g2, [1e-6 20], type{1}, a{:}, 'dim', 2, 'fmm', fmm);
    A = K.eval(s2, t2);
    sig = randn(size(A,2), ndens);
    err = norm(K.fmm(1e-12, s2, t2, sig) - A*sig, 'fro')/norm(A*sig, 'fro');
    assert(err < 1e-10, 'fmm option dim 2 %s: %.2e', type{1}, err);
end

% the generic maxwell fmms against em3d.fmm
sig = complex(randn(600, 1), randn(600, 1));
for type = {'mfie', 'curls'}
    K = kernel3d(g, type{1}, 'fmm', fmm);
    A = em3d.fmm(1e-12, zk, src, targ, type{1}, sig);
    err = norm(K.fmm(1e-12, src, targ, sig) - A)/norm(A);
    assert(err < 1e-10, 'fmm option vs em3d.fmm %s: %.2e', type{1}, err);
end

% hessians, with the laplace fmm
gl = @(s,t) lapg(s, t);
K = kernel3d.radpp(gl, [1e-6 20], 'hess', 'fmm', ...
        @(eps,s,pg,t,pgt) lfmm3d(eps, s, pg, t, pgt));
A = K.eval(src, targ);
sig = randn(size(A,2), ndens);
err = norm(K.fmm(1e-12, src, targ, sig) - A*sig, 'fro')/norm(A*sig, 'fro');
assert(err < 1e-10, 'fmm option hess: %.2e', err);

% sums and interleaved kernels of radgreen kernels keep the fmm
Ks = kernel3d(g, 's', 'fmm', fmm);
Kd = kernel3d(g, 'd', 'fmm', fmm);
Ksp = kernel3d(g, 'sp', 'fmm', fmm);
for K = {Ks + 2*Kd, Ks - Ksp, kernel3d.interleave([Ks Kd; Ksp Ks])}
    A = K{1}.eval(src, targ);
    sig = complex(randn(size(A,2), ndens), randn(size(A,2), ndens));
    u = K{1}.fmm(1e-12, src, targ, sig);
    err = norm(u - A*sig, 'fro')/norm(A*sig, 'fro');
    assert(err < 1e-10, 'fmm option arithmetic: %.2e', err);
end

end


function test_radpp_sum()
% Sums and multiples of radpp kernels of one type are radpp kernels whose
% fit is the sum of the fits, on the common distances, with the fmms
% summed.

g1 = @(s,t) helm3d.green(1.1, s, t);
g2 = @(s,t) helm3d.green(2.3, s, t);
f1 = @(eps,s,pg,t,pgt) hfmm3d(eps, 1.1, s, pg, t, pgt);
f2 = @(eps,s,pg,t,pgt) hfmm3d(eps, 2.3, s, pg, t, pgt);
src.r = 0.5*randn(3, 30);  src.n = randn(3, 30);  src.n = src.n./vecnorm(src.n);
targ.r = 0.5*randn(3, 20) + [1;0;0];
targ.n = randn(3, 20);  targ.n = targ.n./vecnorm(targ.n);

K1 = kernel3d.radpp(g1, [1e-6 6], 'cp', [1 2i], 'fmm', f1);
K2 = kernel3d.radpp(g2, [1e-5 8], 'cp', [1 2i], 'fmm', f2);
A1 = K1.eval(src, targ);  A2 = K2.eval(src, targ);
tab = {2*K1 - K2, 2*A1 - A2; K1 + K2/3, A1 + A2/3; -K1, -A1; 1i*K1, 1i*A1};
for i = 1:size(tab, 1)
    K = tab{i,1};
    assert(strcmp(K.name, 'radpp'), 'radpp sum %d: not radpp', i);
    err = norm(K.eval(src, targ) - tab{i,2}, 'fro')/norm(tab{i,2}, 'fro');
    assert(err < 1e-13, 'radpp sum %d: %.2e', i, err);
    sig = randn(size(A1, 2), 2);
    err = norm(K.fmm(1e-12, src, targ, sig) - tab{i,2}*sig, 'fro')/ ...
          norm(tab{i,2}*sig, 'fro');
    assert(err < 1e-10, 'radpp sum %d fmm: %.2e', i, err);
end
K = K1 + K2;
assert(K.params.rmin == 1e-5 && K.params.rmax == 6, 'radpp sum: range');
K = K1 + kernel3d.radpp(g1, 6, 'sp');
assert(~strcmp(K.name, 'radpp'), 'radpp sum: different types merged');

end


function test_grad_ntrace()
% grad and normaltrace (ntrace) against the gradient and prime types, with
% their fmms, for radpp, radgreen and the built-in kernels

zk = 1.3;
g = @(s,t) helm3d.green(zk, s, t);
fmm = @(eps,s,pg,t,pgt) hfmm3d(eps, zk, s, pg, t, pgt);
src.r = randn(3, 30);  src.n = randn(3, 30);  src.n = src.n./vecnorm(src.n);
targ.r = randn(3, 20) + [5;0;0];
targ.n = randn(3, 20);  targ.n = targ.n./vecnorm(targ.n);

tab = {'s', 'sp', 'grad'; 'd', 'dp', 'dgrad'; 'c', 'cp', 'cgrad'; ...
    'trans_rep', 'trep_prime', 'trep_grad'};
mks = {@(type, a) kernel3d.radpp(g, [1e-6 20], type, a{:}, 'fmm', fmm), ...
       @(type, a) kernel3d(g, type, a{:}, 'fmm', fmm)};
for i = 1:size(tab, 1)
    a = targs3(tab{i,1});
    for k = 1:2
        K = mks{k}(tab{i,1}, a);
        E = mks{k}(tab{i,2}, a);
        Kn = normaltrace(K);
        A = E.eval(src, targ);
        sig = complex(randn(size(A,2), 2), randn(size(A,2), 2));
        assert(norm(Kn.eval(src, targ) - A, 'fro') < 1e-14*norm(A, 'fro'), ...
            'normaltrace %s (%d)', tab{i,1}, k);
        assert(norm(Kn.fmm(1e-12, src, targ, sig) - A*sig, 'fro') < ...
            1e-10*norm(A*sig, 'fro'), 'normaltrace fmm %s (%d)', tab{i,1}, k);
        assert(isequal(Kn.params.ipv, E.params.ipv), 'normaltrace ipv %s (%d)', ...
            tab{i,1}, k);
        assert(isequal(ntrace(K).eval(src, targ), Kn.eval(src, targ)), ...
            'ntrace %s (%d)', tab{i,1}, k);
        A = mks{k}(tab{i,3}, a).eval(src, targ);
        assert(norm(grad(K).eval(src, targ) - A, 'fro') < ...
            1e-14*norm(A, 'fro'), 'grad %s (%d)', tab{i,1}, k);
    end
end

% the built-in kernels
Ks = {kernel3d.lap3d('d'), kernel3d.lap3d('dp'); ...
      kernel3d.helm3d('c', zk, [1 2i]), kernel3d.helm3d('cp', zk, [1 2i]); ...
      kernel3d.helm3d('s_diff', [zk; 2.1], [1; 0.5]), ...
      kernel3d.helm3d('sp_diff', [zk; 2.1], [1; 0.5]); ...
      kernel3d.helm3d('c_diff', [zk; 2.1], [1; 2i; 1; 0.5]), ...
      kernel3d.helm3d('cp_diff', [zk; 2.1], [1; 2i; 1; 0.5]); ...
      kernel3d.helm3d('trans_rep_diff', [zk; 2.1], [0.7; 2i; 1; 0.5]), ...
      kernel3d.helm3d('trep_prime_diff', [zk; 2.1], [0.7; 2i; 1; 0.5])};
for i = 1:size(Ks, 1)
    A = Ks{i,2}.eval(src, targ);
    assert(norm(normaltrace(Ks{i,1}).eval(src, targ) - A, 'fro') < ...
        1e-14*norm(A, 'fro'), 'built-in normaltrace %d', i);
end

% trans_rep of helm3d, an interleaved kernel
cf = [0.7; 2i];
P = kernel3d(@(s,t) helm3d.green(zk, s, t), 'trep_prime', cf.').eval(src, targ);
Kn = normaltrace(kernel3d.helm3d('trans_rep', zk, cf));
assert(strcmp(Kn.type, 'trep_prime') && ...
    norm(Kn.eval(src, targ) - P, 'fro') < 1e-12*norm(P, 'fro'), ...
    'helm3d trans_rep normaltrace');

% the gradients of lap3d and helm3d, as grad of the layer potentials,
% against the radgreen kernels, with their fmms
gl = @(s,t) lap3d.green(s, t);
tab = {kernel3d.lap3d('s'), gl, 'grad', {}; ...
       kernel3d.lap3d('c', [0.7; 2]), gl, 'cgrad', {[0.7 2]}; ...
       kernel3d.helm3d('d', zk), g, 'dgrad', {}; ...
       kernel3d.helm3d('c', zk, cf), g, 'cgrad', {cf.'}; ...
       kernel3d.helm3d('trans_rep', zk, cf), g, 'trep_grad', {cf.'}};
for i = 1:size(tab, 1)
    Kg = grad(tab{i,1});
    A = kernel3d(tab{i,2}, tab{i,3}, tab{i,4}{:}).eval(src, targ);
    sig = randn(size(A,2), 2);
    assert(norm(Kg.eval(src, targ) - A, 'fro') < 1e-12*norm(A, 'fro'), ...
        'built-in grad %d', i);
    assert(norm(Kg.fmm(1e-12, src, targ, sig) - A*sig, 'fro') < ...
        1e-10*norm(A*sig, 'fro'), 'built-in grad fmm %d', i);
end

% sums of radpp and radgreen kernels keep the fmm of their greens
g2 = @(s,t) helm3d.green(2.1, s, t);
fmm2 = @(eps,s,pg,t,pgt) hfmm3d(eps, 2.1, s, pg, t, pgt);
for k = 1:2
    if k == 1
        K1 = kernel3d.radpp(g, [1e-6 20], 'c', cf.', 'fmm', fmm);
        K2 = kernel3d.radpp(g2, [1e-6 20], 'c', cf.', 'fmm', fmm2);
    else
        K1 = kernel3d(g, 'c', cf.', 'fmm', fmm);
        K2 = kernel3d(g2, 'c', cf.', 'fmm', fmm2);
    end
    Kn1 = normaltrace(K1);
    Kn2 = normaltrace(K2);
    Kn = normaltrace(2*K1 - K2);
    A = 2*Kn1.eval(src, targ) - Kn2.eval(src, targ);
    sig = randn(size(A,2), 2);
    assert(norm(Kn.eval(src, targ) - A, 'fro') < 1e-12*norm(A, 'fro'), ...
        'sum normaltrace (%d)', k);
    assert(norm(Kn.fmm(1e-12, src, targ, sig) - A*sig, 'fro') < ...
        1e-10*norm(A*sig, 'fro'), 'sum normaltrace fmm (%d)', k);
end

% 2D (in-plane) radpp kernels keep their dimension
g2 = @(s,t) helmgreen(1.1, s, t);
src2.r = randn(2, 30);  src2.n = randn(2, 30);  src2.n = src2.n./vecnorm(src2.n);
targ2.r = randn(2, 20) + [5; 0];
targ2.n = randn(2, 20);  targ2.n = targ2.n./vecnorm(targ2.n);
tab = {'s', {}, 'grad', 'sp'; 'c', {[1 2i]}, 'cgrad', 'cp'};
for i = 1:size(tab, 1)
    a = tab{i,2};
    K = kernel3d.radpp(g2, [1e-6 20], tab{i,1}, a{:}, 'dim', 2);
    Kg = grad(K);
    Kn = normaltrace(K);
    A = kernel3d.radpp(g2, [1e-6 20], tab{i,3}, a{:}, 'dim', 2).eval(src2, targ2);
    B = kernel3d.radpp(g2, [1e-6 20], tab{i,4}, a{:}, 'dim', 2).eval(src2, targ2);
    assert(Kg.params.dim == 2 && norm(Kg.eval(src2, targ2) - A, 'fro') < ...
        1e-12*norm(A, 'fro'), '2D grad %s', tab{i,1});
    assert(Kn.params.dim == 2 && norm(Kn.eval(src2, targ2) - B, 'fro') < ...
        1e-12*norm(B, 'fro'), '2D normaltrace %s', tab{i,1});
end

% kernels that grad and normaltrace reject: a prime type, a custom type,
% a sum of different types and an interleave of different kernels
Krej = {kernel3d.radpp(g, 20, 'sp'), ...
        kernel3d.radpp(g, [1e-6 20], @(t) ones(size(t.r,2), 1), 'nderiv', 0), ...
        kernel3d(g, 's') + kernel3d(g, 'dp'), ...
        kernel3d.interleave([kernel3d(g, 's'), kernel3d.lap3d('s')])};
for i = 1:numel(Krej)
    for op = {@grad, @normaltrace}
        caught = false;
        try
            op{1}(Krej{i});
        catch
            caught = true;
        end
        assert(caught, '%s: kernel %d not rejected', func2str(op{1}), i);
    end
end
end


function a = targs3(type)
% the parameter of the types that take one
a = {};
if any(strcmp(type, {'c','cp','cgrad','c2trans','trans_rep', ...
        'trep_prime','trep_grad','trans_rep_prime'})), a = {[1 2i]}; end
if strcmp(type, 'trans_sys'), a = {[1 2; 3i 4]}; end
end


function [val, grad, hess] = lapg(s, t)
% 1/(4 pi r) and its target derivatives
d = permute(t, [2 3 1]) - permute(s, [3 2 1]);
r = sqrt(sum(d.^2, 3));
val = 1./(4*pi*r);
grad = -d./(4*pi*r.^3);
hess = zeros([size(r) 3 3]);
for i = 1:3
    for j = 1:3
        hess(:,:,i,j) = 3*d(:,:,i).*d(:,:,j)./(4*pi*r.^5) - (i == j)./(4*pi*r.^3);
    end
end
end


function v = notradial(s, t)
v = exp(-(t(1,:).' - s(1,:)).^2 - 2*(t(2,:).' - s(2,:)).^2);
end


function [v, g, h] = wronghess(zk, s, t)
[v, g, h] = helm3d.green(zk, s, t);
h = 1.01*h;
end


function test_aliases()
% Every name of a type gives the same kernel.

helm = @(s,t) helmgreen(1.1, s, t);
names = {{'s', 'single', 'slp'}, {'d', 'double', 'dlp'}, ...
         {'c', 'comb', 'combined'}, {'sp', 'sprime'}, ...
         {'dp', 'dprime'}, {'cp', 'combprime', 'cprime'}, ...
         {'grad', 'sgrad', 'sg'}, {'dgrad', 'dg'}, {'cgrad', 'cg'}, ...
         {'hess', 'shess'}, ...
         {'trans_sys', 'transmission_sys', 'trans', 'transmission', 'tsys'}, ...
         {'trans_rep', 'transmission_rep', 'trep'}, ...
         {'trep_prime', 'trepp', 'trep_p', 'trans_rep_prime', 'trans_rep_p'}, ...
         {'trep_grad', 'trepg', 'trep_g', 'trans_rep_grad', 'trans_rep_g'}, ...
         {'clamped_bc', 'clamped', 'clamped_plate_bcs'}, {'s2trans', 's2t'}, ...
         {'d2trans', 'd2t'}, {'c2trans', 'c2t', 'c2tr'}};
src.r = randn(3, 5);  src.n = src.r./vecnorm(src.r);
targ.r = 2*randn(3, 6);  targ.n = targ.r./vecnorm(targ.r);
for k = 1:numel(names)
    A = kernel3d.radpp(helm, [1e-3 20], names{k}{1}).eval(src, targ);
    for j = 2:numel(names{k})
        K = kernel3d.radpp(helm, [1e-3 20], upper(names{k}{j}));
        assert(isequal(K.eval(src, targ), A), ...
               'radpp: alias ''%s'' differs', names{k}{j});
    end
end

% positional coefs and nu, as in helm3d
c = [1 1i];
A = kernel3d.radpp(helm, [1e-3 20], 'cp', 'coefs', c).eval(src, targ);
B = kernel3d.radpp(helm, [1e-3 20], 'cp', c).eval(src, targ);
assert(isequal(A, B), 'radpp: positional coefs');
g = @(s,t) helm3d.green(1.1, s, t);
A = kernel3d(g, 'c', 'coefs', c).eval(src, targ);
B = kernel3d(g, 'c', c).eval(src, targ);
assert(isequal(A, B), 'radgreen: positional coefs');
t2 = struct('r', targ.r(1:2,:), 'n', targ.n(1:2,:), ...
            'd', targ.n(1:2,:), 'd2', targ.n(1:2,:));
s2 = struct('r', src.r(1:2,:), 'n', src.n(1:2,:));
A = kernel3d.radpp(helm, [1e-3 20], 'supported', 'nu', 0.3, 'dim', 2).eval(s2, t2);
B = kernel3d.radpp(helm, [1e-3 20], 'supported', 0.3, 'dim', 2).eval(s2, t2);
assert(isequal(A, B), 'radpp: positional nu');
for bad = {{'s', [1 2]}, {'s', 'coefs', [1 2]}, {'c', [1 2], 'coefs', [3 4]}}
    try
        kernel3d(g, bad{1}{:});
        error('radpp:test', 'radgreen: bad positional coefs accepted');
    catch ME
        assert(~strcmp(ME.identifier, 'radpp:test'), ME.message);
    end
end

end


function test_getquad_sp_3d()
% Near quadrature of 'sp' in 3D matches helm3d('sp') on a sphere.

zk = 1.1;
eps = 1e-7;
S = neighbours(geometries.sphere(1, 1, [0;0;0], 4, 1), 3);
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


function test_full_2d_format()
% A green on points in the plane with full derivative arrays, hess
% (nt,ns,2,2) and third (nt,ns,2,2,2) as from chunkie, gives the same
% kernels as with packed arrays.

src.r = 2*rand(2, 7) - 1;
th = 2*pi*rand(1, 9);
targ.r  = 1.7*[cos(th); sin(th)];
targ.d  = 1.7*[-sin(th); cos(th)];
targ.d2 = -targ.r;
targ.n  = [cos(th); sin(th)];
for type = {'s', 'sp', 'clamped_bc', 'supported_bc', 'free_bc'}
    A = kernel3d.radpp(@gauss, [1e-6 6], type{1}, 'nu', 0.3, 'dim', 2);
    B = kernel3d.radpp(@gaussfull2, [1e-6 6], type{1}, 'nu', 0.3, 'dim', 2);
    err = norm(A.eval(src, targ) - B.eval(src, targ), 'fro')/ ...
          norm(A.eval(src, targ), 'fro');
    assert(err < 1e-14, 'radpp full 2d format %s: %.2e', type{1}, err);
end

end


function varargout = gaussfull2(s, t)
% gauss with full derivative arrays on points in the plane: hess
% (nt,ns,2,2), third (nt,ns,2,2,2)
out = cell(1, max(nargout, 1));
[out{:}] = gauss(s(1:2,:), t(1:2,:));
varargout = out;
if nargout >= 3
    h = out{3};
    H = cat(4, cat(3, h(:,:,1), h(:,:,2)), cat(3, h(:,:,2), h(:,:,3)));
    varargout{3} = H;
end
if nargout >= 4
    t3 = out{4};                       % xxx xxy xyy yyy
    T = zeros([size(t3, 1), size(t3, 2), 2, 2, 2]);
    for i = 1:2, for j = 1:2, for k = 1:2
        T(:,:,i,j,k) = t3(:,:,1 + (i==2) + (j==2) + (k==2));
    end, end, end
    varargout{4} = T;
end
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


function test_real_and_decaying()
% laplace (real, so real near quadrature) and yukawa (imaginary
% wavenumber) kernels against lap3d and helm3d on a sphere

eps = 1e-7;
S = neighbours(geometries.sphere(1, 2, [0;0;0], 4, 1), 3);

for type = {'s', 'd', 'sp', 'dp'}
    K = kernel3d(@lap3d.green, type{1});
    Kl = kernel3d.lap3d(type{1});
    assert(~K.ifcomplex, 'radgreen lap %s: flagged complex', type{1});
    assert(isequal(K.params.ipv, kernel3d(@(s,t) helm3d.green(1, s, t), ...
           type{1}).params.ipv), 'radgreen lap %s: ipv', type{1});
    Q = K.getquad(S, eps);  Ql = Kl.getquad(S, eps);
    assert(isreal(Q), 'radgreen lap %s: complex quadrature', type{1});
    err = norm(Q - Ql, 'fro')/norm(Ql, 'fro');
    assert(err < 1e-6, 'radgreen lap %s getquad: %.2e', type{1}, err);
end

zk = 3i;
for type = {'s', 'dp'}
    K = kernel3d(@(s,t) helm3d.green(zk, s, t), type{1});
    Q = K.getquad(S, eps);  Qh = kernel3d.helm3d(type{1}, zk).getquad(S, eps);
    err = norm(Q - Qh, 'fro')/norm(Qh, 'fro');
    assert(err < 1e-6, 'radgreen yukawa %s getquad: %.2e', type{1}, err);
end

end


function test_two_surfers()
% two surfers, each the patch facing the other across a 0.3 gap and the
% patch on the far side: near quadrature between surfaces, through
% surfermat

eps = 1e-7;  zk = 1.1;
c1 = [0;0;0];  c2 = [2.0;0.3;0];
S1 = facing(geometries.sphere(1, 2, c1, 4, 1), c2);
S2 = facing(geometries.sphere(0.7, 2, c2, 4, 1), c1);
g = @(s,t) helm3d.green(zk, s, t);
for type = {'s', 'sp', 'dp'}
    A = surfermat({S1, S2}, kernel3d(g, type{1}), eps);
    B = surfermat({S1, S2}, kernel3d.helm3d(type{1}, zk), eps);
    err = norm(A - B, 'fro')/norm(B, 'fro');
    assert(err < 1e-6, 'radgreen two surfers %s: %.2e', type{1}, err);
end

end


function test_srcwfun()
% custom source weights reproduce the built-in double layer, eval and
% near quadrature

zk = 1.1;  eps = 1e-7;
g = @(s,t) helm3d.green(zk, s, t);
one = @(t) ones(size(t.r, 2), 1);
Kc = kernel3d(g, @(t) 0*one(t), 'nderiv', 0, 'srcwfun', one);
Kd = kernel3d(g, 'd');
assert(isequal(Kc.src_fields, {'n'}) && Kc.params.ipv == Kd.params.ipv, ...
       'radgreen srcwfun: fields or ipv');

src.r = randn(3, 20);  src.n = randn(3, 20);  src.n = src.n./vecnorm(src.n);
targ.r = randn(3, 15) + [3;0;0];
err = norm(Kc.eval(src, targ) - Kd.eval(src, targ), 'fro')/norm(Kd.eval(src, targ), 'fro');
assert(err < 1e-14, 'radgreen srcwfun eval: %.2e', err);

S = neighbours(geometries.sphere(1, 2, [0;0;0], 4, 1), 3);
Q = Kc.getquad(S, eps);  Qd = Kd.getquad(S, eps);
err = norm(Q - Qd, 'fro')/norm(Qd, 'fro');
assert(err < 1e-12, 'radgreen srcwfun getquad: %.2e', err);

end


function S = facing(S, x)
% the patch of S nearest to the point x and the one farthest from it
d = vecnorm(S.cms - x);
[~, inear] = min(d);  [~, ifar] = max(d);
S = slicesurfer(S, [inear, ifar]);
end


function S = neighbours(S, n)
% patch 1 of S and its n-1 nearest patches, by centroid: self, adjacent
% and near pairs at a fraction of the cost of the whole surface
[~, idx] = sort(vecnorm(S.cms - S.cms(:,1)));
S = slicesurfer(S, idx(1:n));
end


function C = tens(dx, dy, dz, p, varargin)
% tensors from the fit in the params p of a radpp kernel
r = sqrt(dx(:).^2 + dy(:).^2 + dz(:).^2);
C = radpp_rotate(dx, dy, dz, radpp_eval(r, p.breaks, p.coefs), varargin{:});
end


function a = cfs(type, coefs)
% the 'coefs' option, for the types that take it
a = {};
if any(strcmp(type, {'c', 'cp', 'c2trans', 'trans_sys', 'trans_rep'}))
    a = {'coefs', coefs};
end
end


function f = flexprof(zk, r)
% on-axis profiles [G, G_x] of flex2d.green
[v, g] = flex2d.green(zk, [0; 0], [r(:).'; 0*r(:).']);
f = [v(:), g(:,1,1)];
end


function test_radgreen_2d()
% radgreen with dim 2 greens: flex2d.green against its fit, packed vs full derivatives

nu = 0.3;
src = []; src.r = [randn(2, 7); zeros(1, 7)];
th = 2*pi*rand(1, 5);
t = []; t.r = [2*cos(th); 2*sin(th); 0*th];
t.n = [cos(th); sin(th); 0*th];
t.d = [-2*sin(th); 2*cos(th); 0*th];
t.d2 = -t.r;
types = {'s', 'sp', 'grad', 'clamped_bc', 'supported_bc', 'free_bc'};
for zk = {1.7, [1.7, 0.5i], [1.7, 0], 0}
    z = zk{1};
    if isscalar(z) && z ~= 0, z = [z, 1i*z]; end
    g = @(s,t) flex2d.green(z, s, t);
    for i = 1:numel(types)
        args = {'dim', 2, 'nu', nu, 'zk', z};
        K = kernel3d(g, types{i}, args{:});
        A = K.eval(src, t);
        B = kernel3d.radpp(g, [1e-4 10], types{i}, args{:}).eval(src, t);
        err = norm(A - B)/norm(B);
        assert(err < 1e-10, 'radgreen 2d %s, zk = %s: %.1e', types{i}, ...
               mat2str(z), err);
        assert(all(K.params.ipv == 0), 'radgreen 2d %s: ipv %s', types{i}, ...
               mat2str(K.params.ipv));
    end
end

% packed and full derivatives give the same kernel
A = kernel3d(@gauss, 'free', 'nu', nu, 'dim', 2).eval(src, t);
B = kernel3d(@gaussfull2, 'free', 'nu', nu, 'dim', 2).eval(src, t);
assert(norm(A - B) <= 1e-14*norm(A), 'radgreen 2d: packed vs full');

end


function test_fit_battery()
% radpp_fit resolves singular, smooth, oscillatory and flexural profiles
% on short and long ranges, at any scale, without warnings

prof = { ...
  'r^m log r', @(r) [log(r), r.*log(r), r.^2.*log(r), r.^3.*log(r)]; ...
  '1/r^p',     @(r) [1./r, 1./r.^2, 1./r.^3]; ...
  'smooth',    @(r) [exp(-r.^2), exp(-r)]; ...
  'helm 3d',   @(r) [exp(1i*r)./r, exp(30i*r)./r]; ...
  'helm 2d',   @(r) [besselh(0, 1, 2*r), besselh(0, 1, 30*r)]; ...
  'flex',      @(r) [flexprof(1.7, r), flexprof([1.7, 0.5i], r), ...
                     flexprof([1.7, 0], r), flexprof(0, r)]};
ranges = [1e-9 1; 1e-4 10; 1e-6 100];

for i = 1:size(prof, 1)
    for j = 1:size(ranges, 1)
        for sc = [1 1e6 1e-6]
            g = @(r) sc*prof{i,2}(r);
            a = ranges(j,1);  b = ranges(j,2);
            lbl = sprintf('%s on [%g, %g], scale %g', prof{i,1}, a, b, sc);
            lastwarn('');
            [~, ~, info] = radpp_fit(g, a, b, struct('eps', 1e-12));
            assert(isempty(lastwarn()), 'radpp_fit %s: %s', lbl, lastwarn());
            assert(info.err < 1e-11, 'radpp_fit %s: error %.1e', lbl, info.err);
            r = min(max(logspace(log10(a), log10(b), 2000).', a), b);
            fv = g(r);
            pv = radpp_eval(r, info.breaks, info.coefs);
            fsc = max(abs(fv(r >= (a + b)/2, :)), [], 1);
            fsc(fsc == 0) = 1;
            den = max(movmax(abs(fv), 101), fsc);
            err = max(abs(pv - fv)./den, [], 'all');
            assert(err < 1e-10, 'radpp_fit %s: eval error %.1e', lbl, err);
        end
    end
end

end

