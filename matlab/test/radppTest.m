% Verify kernel3d.radpp: the fit, the tensor evaluation, the built-in
% kernels and the near quadrature.

% run ../startup.m

rng(42);

test_fit_limits();
test_tensors();
test_builtins();
test_getquad_lap();
test_getquad_rows();


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


function test_tensors()
% Tensor components from the on-axis profiles match a radial function
% with known derivatives, G = exp(-r^2/2).

K = kernel3d.radpp(@gauss, [1e-6 5], @(t) ones(size(t.r,2), 10), ...
                   'nderiv', 3);
p = K.params;

ns = 20;  nt = 30;
s = randn(2, ns);  t = randn(2, nt);
dx = t(1,:).' - s(1,:);  dy = t(2,:).' - s(2,:);
C = radpp_tensors(dx, dy, 0*dx, p.breaks, p.coefs, 3);

[v, g, h, th] = gauss(s, t);
D = [v(:), reshape(g, [], 2), reshape(h, [], 3), reshape(th, [], 4)];
err = max(abs(C(:) - D(:)));
assert(err < 1e-11, 'radpp tensors: %.2e', err);

end


function test_builtins()
% Built-in kernels match directional derivatives of G = exp(-r^2/2),
% contracted here with full 2x2 and 2x2x2 tensors.

nu = 0.3;
ns = 7;  nt = 9;
src.r = 2*rand(2, ns) - 1;
th = 2*pi*rand(1, nt);
targ.r  = 1.7*[cos(th); sin(th)];
targ.d  = 1.7*[-sin(th); cos(th)];
targ.d2 = -targ.r;
targ.n  = [cos(th); sin(th)];
kappa = 1/1.7;

[v, g, h, t3] = gauss(src.r, targ.r);
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
    K = kernel3d.radpp(@gauss, [1e-6 6], types{k}, 'nu', nu);
    A = K.eval(src, targ);
    B = ref.(types{k});
    err = norm(A - B, 'fro')/norm(B, 'fro');
    assert(err < 1e-11, 'radpp %s: %.2e', types{k}, err);
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


function test_getquad_rows()
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

Ks = kernel3d.radpp(helm, [1e-9 4], 's');
Kc = kernel3d.radpp(helm, [1e-9 4], 'clamped_bc', 'ipv', [0 0]);
Kw = kernel3d.radpp(helm, [1e-9 4], ...
        @(t) cat(3, [1+0*t.r(1,:).', zeros(size(t.r,2),2)], ...
                    [zeros(size(t.r,2),1), t.n(1:2,:).']), ...
        'nderiv', 1, 'opdims', [2 1], 'targ_fields', {'n'});

Qs = Ks.getquad(S, eps, targ);
Qc = Kc.getquad(S, eps, targ);
Qw = Kw.getquad(S, eps, targ);

err1 = norm(Qc(1:2:end,:) - Qs, 'fro')/norm(Qs, 'fro');
err2 = norm(Qw - Qc, 'fro')/norm(Qc, 'fro');
assert(err1 < 1e-12, 'radpp getquad rows: %.2e', err1);
assert(err2 == 0, 'radpp getquad custom: %.2e', err2);

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
% 3D Helmholtz Green's function in the plane z = 0 and its gradient
x = t(1,:).' - s(1,:);  y = t(2,:).' - s(2,:);
r = sqrt(x.^2 + y.^2);
G = exp(1i*zk*r)./(4*pi*r);
dG = G.*(1i*zk - 1./r)./r;
varargout{1} = G;
varargout{2} = cat(3, dG.*x, dG.*y);
end
