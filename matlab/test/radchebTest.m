% Verify the radcheb adaptive Chebyshev interpolant of a radial kernel:
% build accuracy, vectorized eval, and the near quadrature correction.

% run ../startup.m

rng(42);

%% Now run the tests

test_build_real();
test_build_complex();
test_eval_matrix();
test_eval_reference();
test_getquad_lap();
test_getquad_helm();
test_bad_inputs();


function test_build_real()
% Interpolant of a real radial kernel matches the exact kernel.

rmin = 1e-3;  rmax = 4;  tol = 1e-12;
fr = @(r) exp(-r)./r;

opts = struct('eps', tol, 'norder', 16);
kern = kernel3d.radcheb(radkern(fr), [rmin, rmax], opts);

r = rmin + (rmax-rmin)*rand(2000,1);
val = radcheb_eval(r, kern.params.breaks, kern.params.coefs);
err = max(abs(val - fr(r)))/max(abs(fr(r)));
assert(~kern.ifcomplex, 'radcheb: real kernel flagged complex');
assert(err < 1e-10, 'radcheb real build: %.2e', err);

end


function test_build_complex()
% Complex kernel, checked against the exact Helmholtz single layer.

zk = 1.1 + 0.2i;
rmin = 1e-3;  rmax = 3;
fr = @(r) exp(1i*zk*r)./(4*pi*r);

opts = struct('eps', 1e-12, 'norder', 16);
kern = kernel3d.radcheb(radkern(fr), [rmin, rmax], opts);

r = rmin + (rmax-rmin)*rand(2000,1);
val = radcheb_eval(r, kern.params.breaks, kern.params.coefs);
err = max(abs(val - fr(r)))/max(abs(fr(r)));
assert(kern.ifcomplex, 'radcheb: complex kernel flagged real');
assert(kern.params.ipars(5) == 2, 'radcheb: nv should be 2');
assert(err < 1e-10, 'radcheb complex build: %.2e', err);

end


function test_eval_matrix()
% kern.eval reproduces the underlying kernel matrix.

zk = 0.8;
rmin = 1e-2;  rmax = 6;
fk = kernel3d.helm3d('s', zk);
kern = kernel3d.radcheb(fk, [rmin, rmax], struct('eps', 1e-12));

ns = 200;  nt = 150;
src.r  = randn(3,ns);  src.r  = src.r ./vecnorm(src.r);
targ.r = randn(3,nt);  targ.r = targ.r./vecnorm(targ.r)*2.5;

A  = kern.eval(src, targ);
Ae = fk.eval(src, targ);
err = norm(A - Ae, 'fro')/norm(Ae, 'fro');
assert(isequal(size(A), [nt ns]), 'radcheb eval: wrong size');
assert(err < 1e-10, 'radcheb eval matrix: %.2e', err);

end


function test_eval_reference()
% The vectorized eval agrees with a scalar Horner reference.

kern = kernel3d.radcheb(radkern(@(r) besselj(0, 3*r)), [1e-2, 5], ...
           struct('eps', 1e-12));

breaks = kern.params.breaks;
coefs  = kern.params.coefs;
n      = size(coefs, 2);
nbin   = size(coefs, 1);

r = [1e-2; linspace(1e-2, 5, 200).'; 5];
val = radcheb_eval(r, breaks, coefs);

ref = zeros(size(r));
for k = 1:numel(r)
    ib = 1;
    for j = 1:nbin
        if ( r(k) >= breaks(j) ), ib = j; end
    end
    x = r(k) - breaks(ib);
    v = coefs(ib,1);
    for i = 2:n
        v = coefs(ib,i) + x*v;
    end
    ref(k) = v;
end

err = max(abs(val - ref));

% radii outside the build interval return NaN
rout = [1e-3; 0; -1; 5.0001; 7];
assert(all(isnan(radcheb_eval(rout, breaks, coefs))), ...
       'radcheb eval: out of range radii should be NaN');
assert(~any(isnan(val)), 'radcheb eval: in range radii returned NaN');

assert(err == 0, 'radcheb eval reference: %.2e', err);

end


function test_getquad_lap()
% Near quadrature for the radcheb form of the Laplace single layer
% matches the quadrature for the kernel itself.

eps = 1e-7;

S = geometries.sphere(1, 2, [0;0;0], 4, 1);

fk = kernel3d.lap3d('s');
kern = kernel3d.radcheb(fk, [1e-9, 4], struct('eps', 1e-13));

Q  = kern.getquad(S, eps);
Qe = fk.getquad(S, eps);

assert(~kern.ifcomplex, 'radcheb: laplace kernel flagged complex');
err = norm(Q - Qe, 'fro')/norm(Qe, 'fro');
assert(err < 1e-6, 'radcheb getquad lap: %.2e', err);

end


function test_getquad_helm()
% Near quadrature for the radcheb form of the Helmholtz single layer
% matches the quadrature for the kernel itself.

zk = 1.0;
eps = 1e-7;

S = geometries.sphere(1, 2, [0;0;0], 4, 1);

fk = kernel3d.helm3d('s', zk);
kern = kernel3d.radcheb(fk, [1e-9, 4], struct('eps', 1e-13, 'zk', zk));

Q  = kern.getquad(S, eps);
Qe = fk.getquad(S, eps);

err = norm(Q - Qe, 'fro')/norm(Qe, 'fro');
assert(err < 1e-6, 'radcheb getquad: %.2e', err);

end


function test_bad_inputs()
% Vector valued kernels and kernels with field dependence are rejected.

caught = 0;
try
    kernel3d.radcheb(kernel3d.stok3d('s'), 1);
catch
    caught = caught + 1;
end
try
    kernel3d.radcheb(kernel3d.lap3d('d'), 1);
catch
    caught = caught + 1;
end
assert(caught == 2, 'radcheb: bad inputs not rejected');

end


function f = radkern(fr)
% Wrap a function of r in the standard kernel calling sequence.

f = @(s,t) fr(sqrt((t.r(1,:).' - s.r(1,:)).^2 + (t.r(2,:).' - s.r(2,:)).^2 ...
                 + (t.r(3,:).' - s.r(3,:)).^2));

end
