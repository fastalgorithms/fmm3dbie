%TEST_CPRIME_QUAD  Targeted test of the cprime near-field quadrature.
%
%  Tests three things, in increasing proximity to the surface:
%
%  (1) OFF-SURFACE (r=5): layer_eval vs smooth quadrature — same as
%      test_layer_pot3d but focused on prime kernels only.
%      Near-field correction is zero; exercises FMM path.
%
%  (2) NEAR-SURFACE (r=1+h for small h): layer_eval vs smooth quadrature.
%      Near-field quadrature corrections are active; this stresses the
%      Fortran getquad / eval_addsub path.
%
%  (3) JUMP RELATION for S' (Laplace):
%      The normal derivative of the single layer potential satisfies
%          lim_{h->0+} S'[sigma](x + h n(x)) = -sigma(x)/2 + S'[sigma]_PV(x)
%          lim_{h->0-} S'[sigma](x + h n(x)) = +sigma(x)/2 + S'[sigma]_PV(x)
%      so  u^+ - u^- = -sigma  (exterior minus interior).
%      We verify this numerically.

run ../../startup.m
rng(42);

eps  = 1e-7;
tol  = 1e-4;

fprintf('=== test_cprime_quad ===\n\n');

S    = geometries.sphere(1, 4, [0;0;0], 5, 1);
wts  = S.wts(:);
src.r = S.r;  src.n = S.n;  src.du = S.du;  src.dv = S.dv;

sigma_r = cos(S.r(1,:).' + 0.3*S.r(2,:).') .* exp(-S.r(3,:).'.^2);
sigma_r = sigma_r / norm(sigma_r);

zk     = 1.1 + 0.4i;
sigma_c = (S.r(3,:).^2 + 1i*S.r(1,:)).';
sigma_c = sigma_c / norm(sigma_c);

coefs_lap = [0.7; 1.3];
coefs_h   = [1i*zk; 1.0];

nfail = 0;

% =========================================================================
fprintf('(1) Off-surface targets (r=5) — FMM path\n');

targ_far.r = randn(3,40);
targ_far.r = targ_far.r ./ vecnorm(targ_far.r) * 5.0;
targ_far.n = targ_far.r ./ vecnorm(targ_far.r);

prime_tests = {
  'lap sp',    kernel3d('l','sp'),               sigma_r;
  'lap dp',    kernel3d('l','dp'),               sigma_r;
  'lap cp',    kernel3d('l','cp',coefs_lap),     sigma_r;
  'helm sp',   kernel3d('h','sp',zk),            sigma_c;
  'helm dp',   kernel3d('h','dp',zk),            sigma_c;
  'helm cp',   kernel3d('h','cp',zk,coefs_h),   sigma_c;
};

for k = 1:size(prime_tests,1)
    label = prime_tests{k,1};
    K     = prime_tests{k,2};
    sigma = prime_tests{k,3};

    ref = K.eval(src, targ_far) * (wts .* sigma);
    lp  = K.layer_eval(S, sigma, targ_far, eps);

    err = norm(lp - ref) / max(norm(ref), 1e-30);
    if err <= tol, status = 'PASS'; else, status = 'FAIL'; nfail = nfail+1; end
    fprintf('  %-14s  err = %.2e  [%s]\n', label, err, status);
end
%%
% =========================================================================
fprintf('\n(2) Near-surface targets (r=1.5) — near-field quadrature active\n');

h_near = 0.5;
targ_near = [];
targ_near.r = S.r * (1 + h_near);   % just outside unit sphere
targ_near.n = S.r ./ vecnorm(S.r);   % outward normals

% Reference: high-order oversampled smooth quadrature (valid since targets
% are off-surface at distance h_near, and the smooth rule converges for
% Calderon-type kernels with the singularity removed by the gap).
nover_ref = max(S.norders)*2;
Sover = oversample(S, nover_ref * ones(S.npatches, 1));
[~, xinterp] = S.oversample(nover_ref * ones(S.npatches, 1));
src_over.r = Sover.r;  src_over.n = Sover.n;
src_over.du = Sover.du;  src_over.dv = Sover.dv;
wts_over = Sover.wts(:);

for k = 1:size(prime_tests,1)
    label = prime_tests{k,1};
    K     = prime_tests{k,2};
    sigma = prime_tests{k,3};
    sigma_ov = xinterp * sigma(:);

    ref = K.eval(src_over, targ_near) * (wts_over .* sigma_ov);
    lp  = K.layer_eval(S, sigma, targ_near, eps);

    err = norm(lp - ref) / max(norm(ref), 1e-30);
    if err <= tol*10, status = 'PASS'; else, status = 'FAIL'; nfail = nfail+1; end
    fprintf('  %-14s  err = %.2e  [%s]\n', label, err, status);
end
%%
% =========================================================================
fprintf('\n(3) Jump relation for lap S'' (exterior - interior = -sigma)\n');

h_jump = 1e-4;
targ_ext.r = S.r * (1 + h_jump);
targ_ext.n = S.r ./ vecnorm(S.r);
targ_int.r = S.r * (1 - h_jump);
targ_int.n = (S.r ./ vecnorm(S.r));   % inward for interior limit

K_sp = kernel3d('l','sp');


u_ext = K_sp.layer_eval(S, sigma_r, targ_ext, eps);
u_int = K_sp.layer_eval(S, sigma_r, targ_int, eps);

jump = u_ext - u_int;   % should be -sigma_r (up to O(h) error)
err_jump = norm(jump + sigma_r) / norm(sigma_r);
if err_jump <= 0.05   % allow O(h) error
    status = 'PASS';
else
    status = 'FAIL';
    nfail = nfail + 1;
end
fprintf('  lap sp  u^+ - u^- vs -sigma  err = %.2e  [%s]\n', err_jump, status);
% =========================================================================
fprintf('\n(3) Jump relation for Helm S'' (exterior - interior = -sigma)\n');
K_sp = kernel3d('h','sp',zk);

u_ext = K_sp.layer_eval(S, sigma_r, targ_ext, eps);
u_int = K_sp.layer_eval(S, sigma_r, targ_int, eps);

jump = u_ext - u_int;   % should be -sigma_r (up to O(h) error)
err_jump = norm(jump + sigma_r) / norm(sigma_r);
if err_jump <= 0.05   % allow O(h) error
    status = 'PASS';
else
    status = 'FAIL';
    nfail = nfail + 1;
end
fprintf('  lap sp  u^+ - u^- vs -sigma  err = %.2e  [%s]\n', err_jump, status);

%%
h_jump = 1e-4;
targ_ext.r = S.r * (1 + h_jump);
targ_int.r = S.r * (1 - h_jump);

K_sp = kernel3d('l','d');

u_ext = K_sp.layer_eval(S, sigma_r, targ_ext, eps);
u_int = K_sp.layer_eval(S, sigma_r, targ_int, eps);

jump = u_ext - u_int;   % should be -sigma_r (up to O(h) error)
err_jump = norm(jump - sigma_r) / norm(sigma_r);
if err_jump <= 0.05   % allow O(h) error
    status = 'PASS';
else
    status = 'FAIL';
    nfail = nfail + 1;
end
fprintf('  lap d  u^+ - u^- vs sigma  err = %.2e  [%s]\n', err_jump, status);


% =========================================================================
fprintf('\nDone. %d FAILED.\n', nfail);
if nfail > 0
    error('test_cprime_quad: %d failure(s)', nfail);
end
