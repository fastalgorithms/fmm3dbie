%TEST_TANGENTKERN  Test kernel3d.tangent_kern wrapping stok3d('s').
%
% Checks:
%   1. eval    - against manually projected Cartesian kernel matrix
%   2. fmm     - against wrapped eval
%   3. getquad - quadrature correction has correct size and is consistent
%                with layer_eval (smooth-quad reference)
%   4. layer_eval with precomp_quadrature - against layer_eval without

run ../../startup.m

rng(42);
eps = 1e-6;

% --- geometry ---
S    = geometries.sphere(1, 4, [0;0;0], 8, 1);
S    = slicesurfer(S, 1);
ns   = S.npts;

targ.r  = randn(3, 200);
targ.r  = targ.r ./ vecnorm(targ.r) * 3.0;
targ.n  = randn(3, 200); targ.n  = targ.n  ./ vecnorm(targ.n);
targ.du = randn(3, 200); targ.dv = randn(3, 200);
nt = size(targ.r, 2);

% --- kernels ---
K_cart = kernel3d('stok', 's');   % opdims [3 3]
ids    = [1;2;3];
jds    = [1;2;3];
K_tan  = kernel3d.tangent_kern(K_cart, ids, jds);  % opdims [2 2]

% tangent density: 2 components per point
sigma_tan = randn(2*ns, 1);

failures = {};

% =========================================================================
% 1. eval: compare against manually projected Cartesian matrix
% =========================================================================

src_s.r  = S.r(:,:);
src_s.n  = S.n(:,:);
src_s.du = S.du(:,:);
src_s.dv = S.dv(:,:);

Ps_exp = tangent_block(src_s, ids, zeros(0,1), 2, false);  % (3*ns) x (2*ns)
Pt     = tangent_block(targ,  jds, zeros(0,1), 2, true);   % (2*nt) x (3*nt)

K_cart_mat = K_cart.eval(src_s, targ);
pot_ref    = Pt * K_cart_mat * Ps_exp * sigma_tan;
pot_tan    = K_tan.eval(src_s, targ) * sigma_tan;

err = norm(pot_tan - pot_ref) / norm(pot_ref);
failures = report(failures, 'eval vs manual projection', err, 1e-13);

% =========================================================================
% 2. fmm: compare against wrapped eval
% =========================================================================

pot_eval = K_tan.eval(src_s, targ) * sigma_tan;
pot_fmm  = K_tan.fmm(eps, src_s, targ, sigma_tan);

err = norm(pot_fmm(:) - pot_eval(:)) / norm(pot_eval(:));
failures = report(failures, 'fmm vs eval', err, 1e-5);

% =========================================================================
% 3. layer_eval: compare against smooth-quadrature reference
%    (off-surface targets, smooth sigma -> no near singularity)
% =========================================================================

wts         = S.wts(:).';
sigma_wt    = repmat(wts, 2, 1) .* reshape(sigma_tan, 2, ns);
pot_le      = K_tan.layer_eval(S, reshape(sigma_tan, 2, ns), targ, eps);
pot_le_ref  = K_tan.eval(src_s, targ) * sigma_wt(:);

err = norm(pot_le(:) - pot_le_ref(:)) / norm(pot_le_ref(:));
failures = report(failures, 'layer_eval vs smooth quad', err, 1e-4);

% =========================================================================
% 4. layer_eval with precomp_quadrature matches layer_eval without
% =========================================================================

Q    = K_tan.getquad(S, eps, targ);
opts = struct('precomp_quadrature', Q);

pot_precomp = K_tan.layer_eval(S, reshape(sigma_tan, 2, ns), targ, eps, opts);

err = norm(pot_precomp(:) - pot_le(:)) / norm(pot_le(:));
failures = report(failures, 'layer_eval precomp vs plain', err, 1e-10);

% =========================================================================
fprintf('\n=== Summary ===\n');
nfail = numel(failures);
if nfail == 0
    fprintf('All tests passed.\n');
else
    fprintf('%d test(s) FAILED:\n', nfail);
    for k = 1:nfail
        fprintf('  [%d] %s\n', k, failures{k});
    end
    assert(false, 'test_tangentkern: %d failure(s).', nfail);
end

% =========================================================================
% Helpers
% =========================================================================

function failures = report(failures, label, err, tol)
    if err > tol, s = 'FAIL'; else, s = 'PASS'; end
    fprintf('  %-34s  err=%.2e  [%s]\n', label, err, s);
    if strcmp(s, 'FAIL')
        failures{end+1} = sprintf('%s: err=%.2e > tol=%.2e', label, err, tol);
    end
end
