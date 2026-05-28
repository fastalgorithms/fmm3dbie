%TEST_TANGENTKERN  Test kernel3d.tangent_kern wrapping stok3d('s').
%
% Checks (stok s, full 3->2 projection on both sides):
%   1. eval    - against manually projected Cartesian kernel matrix
%   2. fmm     - against wrapped eval
%   3. layer_eval (off-surface) - against smooth-quad reference
%   4. layer_eval with precomp_quadrature - against layer_eval without
%   5. on-surface layer_eval  - one-sided limit via jump relation
%
% Checks (nrccie-eval, partial ids/jds: J->tangent + rho passthrough,
%         both E and H projected into tangent = multi-column jds):
%   6. eval    - against manually projected Cartesian matrix
%   7. layer_eval (off-surface) - against smooth-quad reference

run ../../startup.m

rng(42);
eps  = 1e-6;
zk   = 1.3;

% -------------------------------------------------------------------------
% Geometry: one patch of a sphere, order 8
% -------------------------------------------------------------------------
S    = geometries.sphere(1, 4, [0;0;0], 8, 1);
S    = slicesurfer(S, 1);
ns   = S.npts;

% Off-surface targets
targ.r  = randn(3, 200);
targ.r  = targ.r ./ vecnorm(targ.r) * 3.0;
targ.n  = randn(3, 200); targ.n  = targ.n  ./ vecnorm(targ.n);
targ.du = randn(3, 200); targ.dv = randn(3, 200);

src_s.r  = S.r(:,:);
src_s.n  = S.n(:,:);
src_s.du = S.du(:,:);
src_s.dv = S.dv(:,:);

% -------------------------------------------------------------------------
% Kernel A: stok s, full 3->2 projection on both src and targ
% -------------------------------------------------------------------------
K_cart = kernel3d('stok', 's');   % opdims [3 3]
ids    = [1;2;3];
jds    = [1;2;3];
K_tan  = kernel3d.tangent_kern(K_cart, ids, jds);  % opdims [2 2]

sigma_tan = randn(2*ns, 1);
wts       = S.wts(:).';

failures = {};

fprintf('=== test_tangentkern ===\n\n');
fprintf('--- stok s, full tangent projection ---\n');

% =========================================================================
% 1. eval: compare against manually projected Cartesian matrix
% =========================================================================

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
% 3. layer_eval (off-surface): compare against smooth-quad reference
% =========================================================================

sigma_wt   = repmat(wts, 2, 1) .* reshape(sigma_tan, 2, ns);
pot_le     = K_tan.layer_eval(S, reshape(sigma_tan, 2, ns), targ, eps);
pot_le_ref = K_tan.eval(src_s, targ) * sigma_wt(:);

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
% 5. On-surface layer_eval via jump relation
%
%   stok s has no jump, so the on-surface PV value should equal the
%   average of the one-sided limits:
%     u_pv(x) = (u_ext(x) + u_int(x)) / 2
%   We verify:  surfermatapply(at ic) ≈ (u_ext + u_int) / 2
% =========================================================================

h        = 1e-4;
tol_jump = 0.05;

rcen = mean(S.r, 2);
[~, ic] = min(vecnorm(S.r - rcen));
n0  = S.n(:, ic);  n0 = n0 / norm(n0);
r0  = S.r(:, ic);
du0 = S.du(:, ic);
dv0 = S.dv(:, ic);

targ_ext.r = r0 + h*n0; targ_ext.n = n0; targ_ext.du = du0; targ_ext.dv = dv0;
targ_int.r = r0 - h*n0; targ_int.n = n0; targ_int.du = du0; targ_int.dv = dv0;

sigma_tan2 = reshape(sigma_tan, 2, ns);
u_ext_on = K_tan.layer_eval(S, sigma_tan2, targ_ext, eps);
u_int_on = K_tan.layer_eval(S, sigma_tan2, targ_int, eps);
u_avg    = (u_ext_on + u_int_on) / 2;

pot_on_full = surfermatapply(S, K_tan, sigma_tan, eps);
od = K_tan.opdims(1);
u_on_pv = pot_on_full((ic-1)*od+1 : ic*od);

err = norm(u_on_pv(:) - u_avg(:)) / (norm(u_avg(:)) + 1e-30);
failures = report(failures, 'on-surface PV vs avg(ext,int)', err, tol_jump);

% =========================================================================
% --- nrccie-eval kernel: partial ids + multi-column jds ---
%
%   nrccie-eval: opdims [6, 4]
%     src density: [J_x, J_y, J_z, rho]  (components 1:3 = J, 4 = rho)
%     tgt output:  [E_x, E_y, E_z, H_x, H_y, H_z]  (rows 1:3 = E, 4:6 = H)
%
%   ids = [1;2;3]: project J into tangent basis, rho (row 4) passes through
%     => src opdim: 2 (tangent) + 1 (scalar rho) = 3
%
%   jds = [[1;2;3], [4;5;6]]: project both E and H into tangent (multi-col)
%     => tgt opdim: 2 + 2 = 4  (no scalar passthrough rows)
%
%   K_em_tan: opdims [4, 3]
% =========================================================================

fprintf('\n--- nrccie-eval, partial src ids + multi-column tgt jds ---\n');

K_em   = kernel3d('em', 'nrccie-eval', zk);   % opdims [6, 4]

ids_em = [1;2;3];             % project J block; rho (row 4) passes through
jds_em = [[1;2;3], [4;5;6]]; % two 3-blocks on target, both projected

K_em_tan = kernel3d.tangent_kern(K_em, ids_em, jds_em);  % opdims [4, 3]

sigma_em = randn(3*ns, 1);   % 2 tangent + 1 scalar (rho) per source point

% =========================================================================
% 6. eval: against manually projected Cartesian matrix
% =========================================================================

scalar_src = 4;                 % rho passes through on source side
scalar_tgt = zeros(0,1);        % no passthrough on target side (all projected)

Ps_exp_em = tangent_block(src_s, ids_em, scalar_src, 3, false);  % (4*ns) x (3*ns)
Pt_em     = tangent_block(targ,  jds_em, scalar_tgt, 4, true);   % (4*nt) x (6*nt)

K_em_mat   = K_em.eval(src_s, targ);
pot_em_ref = Pt_em * K_em_mat * Ps_exp_em * sigma_em;
pot_em_tan = K_em_tan.eval(src_s, targ) * sigma_em;

err = norm(pot_em_tan(:) - pot_em_ref(:)) / norm(pot_em_ref(:));
failures = report(failures, 'nrccie eval vs manual', err, 1e-13);

% =========================================================================
% 7. layer_eval (off-surface): against smooth-quad reference
% =========================================================================

sigma_em_wt = repmat(wts, 3, 1) .* reshape(sigma_em, 3, ns);
pot_em_le     = K_em_tan.layer_eval(S, reshape(sigma_em, 3, ns), targ, eps);
pot_em_le_ref = K_em_tan.eval(src_s, targ) * sigma_em_wt(:);

err = norm(pot_em_le(:) - pot_em_le_ref(:)) / norm(pot_em_le_ref(:));
failures = report(failures, 'nrccie layer_eval vs smooth', err, 1e-4);

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
    fprintf('  %-36s  err=%.2e  [%s]\n', label, err, s);
    if strcmp(s, 'FAIL')
        failures{end+1} = sprintf('%s: err=%.2e > tol=%.2e', label, err, tol);
    end
end
