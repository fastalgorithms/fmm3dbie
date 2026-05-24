%TEST_LAYER_POT3D   Verify that kernel3d.layer_eval agrees with smooth
%   (Gauss-Legendre) quadrature for each kernel, for a well-separated target.
%
%   For each kernel the layer potential
%
%       u(x) = integral_Gamma K(x,y) sigma(y) dS_y
%
%   is computed two ways:
%
%     (1) SMOOTH QUADRATURE – direct sum over all surface quadrature nodes
%         with weights, using the pointwise kernel evaluator K.eval.
%         This is exact (to machine precision) when x is well separated
%         from Gamma, so no near-field correction is needed.
%
%     (2) LAYER_EVAL – K.layer_eval, which uses the FMM plus near-field
%         quadrature corrections (the corrections are zero for well-separated
%         targets, but are still exercised).
%
%   The relative l2 error between (1) and (2) is printed for every kernel.
%   A FAIL is reported if the error exceeds the tolerance.

run ../../startup.m

% ---- global parameters ----
zk   = 1.1 + 0.4i;   % Helmholtz wavenumber
eps  = 1e-9;          % requested accuracy for layer_eval
tol  = 1e-5;          % pass/fail threshold

rng(17);

% ---- surface: unit sphere ----
S = geometries.sphere(1, 3, [0;0;0], 5, 1);

% smooth quadrature weights
wts = S.wts(:);   % (npts,1)

% random smooth densities
sigma_r = (S.r(3,:)).';
sigma_c = (S.r(3,:).^2 + 1i*S.r(1,:)).';

% ---- target: single well-separated point (|x| = 5) ----
targ.r = [3; 3; 3] / norm([3;3;3]) * 5;
targ.n = targ.r / norm(targ.r);   % outward normal (needed by cprime)
targ.du = randn(3,1); targ.dv = randn(3,1);
% srcinfo struct built from the surfer object
src.r  = S.r;
src.n  = S.n;
src.du = S.du;
src.dv = S.dv;

nfail = 0;

coefs_lap = [0.7; 1.3];
coefs_h   = [1i*zk; 1.0];

% ---- list of tests: {label, kernel, density} ----
tests = { ...
    'lap  s',                  kernel3d('l','s'),                sigma_r; ...
    'lap  d',                  kernel3d('l','d'),                sigma_r; ...
    'lap  sp',                 kernel3d('l','sp'),               sigma_r; ...
    'lap  dp',                 kernel3d('l','dp'),               sigma_r; ...
    'lap  c',                  kernel3d('l','c',coefs_lap),      sigma_r; ...
    'lap  cp',                 kernel3d('l','cp',coefs_lap),     sigma_r; ...
    'helm s',                  kernel3d('h','s', zk),            sigma_c; ...
    'helm d',                  kernel3d('h','d', zk),            sigma_c; ...
    'helm sp',                 kernel3d('h','sp',zk),            sigma_c; ...
    'helm dp',                 kernel3d('h','dp',zk),            sigma_c; ...
    'helm c',                  kernel3d('h','c', zk,coefs_h),   sigma_c; ...
    'helm cp',                 kernel3d('h','cp',zk,coefs_h),   sigma_c; ...
};

fprintf('\n=== All prime and non-prime kernels (well-separated target) ===\n\n');
for k = 1:size(tests,1)
    label = tests{k,1};
    K     = tests{k,2};
    sigma = tests{k,3};

    % (1) smooth quadrature reference: exact for well-separated target
    ref = K.eval(src, targ) * (sigma .* wts);

    % (2) layer_eval (FMM + near-field quadrature)
    lp  = K.layer_eval(S, sigma, targ, eps);

    err = norm(lp - ref) / max(norm(ref), 1e-30);
    if err <= tol
        status = 'PASS';
    else
        status = 'FAIL';
        nfail  = nfail + 1;
    end
    fprintf('  %-30s rel err = %.2e   [%s]\n', label, err, status);
end

fprintf('\nDone. %d test(s) run, %d FAILED.\n', size(tests,1), nfail);
if nfail > 0
    error('test_layer_pot3d: %d test(s) FAILED.', nfail);
end
