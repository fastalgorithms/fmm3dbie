%TEST_TANGENTKERN  Test kernel3d.tangent_kern and kernel3d/times matrix-fn multiply.
%
% Checks (stok s, full 3->2 tangent projection on both sides):
%   1. eval    - against manually projected Cartesian kernel matrix
%   2. direct times eval - same, using f_left .* K .* f_right directly
%   3. fmm     - against wrapped eval
%   4. layer_eval (off-surface) - against smooth-quad reference
%   5. layer_eval with precomp_quadrature - against layer_eval without
%   6. on-surface layer_eval  - one-sided limit via jump relation
%
% Checks (nrccie-eval, partial ids/jds):
%   7. eval    - against manually projected Cartesian matrix
%   8. layer_eval (off-surface) - against smooth-quad reference
%   9. layer_eval with precomp_quadrature - against layer_eval without
%
% Checks (n-dot left multiply: f(t) = n(t)', stok s -> scalar output):
%  10. eval    - against manual n.'*K_cart
%  11. fmm     - against eval
%  12. layer_eval (off-surface) - against smooth-quad reference
%  13. layer_eval with precomp_quadrature - against layer_eval without

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
targ = [];
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

% Also construct directly via times to verify the new path
src_scalar = zeros(0,1);
tgt_scalar = zeros(0,1);
f_right = @(s) local_tangent_pages(s, ids, src_scalar, 2, false);  % (3 x 2 x ns)
f_left  = @(t) local_tangent_pages(t, jds, tgt_scalar, 2, true);   % (2 x 3 x nt)
K_tan_direct = f_left .* K_cart .* f_right;

sigma_tan = randn(2*ns, 1);
wts       = S.wts(:).';

failures = {};

fprintf('=== test_tangentkern ===\n\n');
fprintf('--- stok s, full tangent projection ---\n');

% =========================================================================
% 1. eval: compare against manually projected Cartesian matrix
% =========================================================================

Ps_exp = tangent_block(src_s, ids, zeros(0,1), 2, false);  % (3 x 2 x ns)
Pt     = tangent_block(targ,  jds, zeros(0,1), 2, true);   % (2 x 3 x nt)

nt_s = size(targ.r, 2);
K_cart_mat = K_cart.eval(src_s, targ);  % (3*nt x 3*ns)
% Apply Ps_exp on right: K_cart_mat*(Ps_exp*sigma_tan)
%   sigma_tan (2*ns x 1) -> (2 x 1 x ns); Ps_exp*(2x1) = (3x1) per page -> (3*ns x 1)
sig_pages  = reshape(sigma_tan, 2, 1, ns);
Ps_sig     = reshape(pagemtimes(Ps_exp, sig_pages), 3*ns, 1);   % (3*ns x 1)
Kps_sig    = K_cart_mat * Ps_sig;                                % (3*nt x 1)
% Apply Pt on left: (2x3 x nt) * (3x1 x nt) -> (2x1 x nt) -> (2*nt x 1)
pot_ref    = reshape(pagemtimes(Pt, reshape(Kps_sig, 3, 1, nt_s)), 2*nt_s, 1);
pot_tan    = K_tan.eval(src_s, targ) * sigma_tan;

err = norm(pot_tan - pot_ref) / norm(pot_ref);
failures = report(failures, 'eval vs manual projection', err, 1e-13);

pot_direct = K_tan_direct.eval(src_s, targ) * sigma_tan;
err = norm(pot_direct - pot_ref) / norm(pot_ref);
failures = report(failures, 'direct times eval vs manual', err, 1e-13);

% =========================================================================
% 3. fmm: compare against wrapped eval
% =========================================================================

pot_eval = K_tan.eval(src_s, targ) * sigma_tan;
pot_fmm  = K_tan.fmm(eps, src_s, targ, sigma_tan);

err = norm(pot_fmm(:) - pot_eval(:)) / norm(pot_eval(:));
failures = report(failures, 'fmm vs eval', err, 1e-5);

% =========================================================================
% 4. layer_eval (off-surface): compare against smooth-quad reference
% =========================================================================

sigma_wt   = repmat(wts, 2, 1) .* reshape(sigma_tan, 2, ns);
pot_le     = K_tan.layer_eval(S, sigma_tan, targ, eps);
pot_le_ref = K_tan.eval(src_s, targ) * sigma_wt(:);

err = norm(pot_le(:) - pot_le_ref(:)) / norm(pot_le_ref(:));
failures = report(failures, 'layer_eval vs smooth quad', err, 1e-4);

% =========================================================================
% 5. layer_eval with precomp_quadrature matches layer_eval without
% =========================================================================

Q    = K_tan.getquad(S, eps, targ);
opts = struct('precomp_quadrature', Q);

pot_precomp = K_tan.layer_eval(S, sigma_tan, targ, eps, opts);

err = norm(pot_precomp(:) - pot_le(:)) / norm(pot_le(:));
failures = report(failures, 'layer_eval precomp vs plain', err, 1e-10);

% =========================================================================
% 6. On-surface layer_eval via jump relation
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
% --- chained multiply: A .* ((C .* K_cart) .* B), stok s ---
%
%   K_cart: opdims [3 3]
%   B (right, src): 3x4 constant matrix  => K_cart .* B:  opdims [3 4]
%   C (left,  tgt): 2x3 constant matrix  => C .* (...):   opdims [2 4]
%   A (left,  tgt): 5x2 constant matrix  => A .* (...):   opdims [5 4]
%
%   Reference: A*C * K_cart_mat * B * sigma_chain
% =========================================================================

fprintf('\n--- chained: A.*((C.*K).*B) ---\n');

rng(7);
B_mat = randn(3, 4);   % right: 3x4, K.opdims(2)=3
C_mat = randn(2, 3);   % left:  2x3, K.opdims(1)=3
A_mat = randn(5, 2);   % left:  5x2, prev opdims(1)=2

K_chain = A_mat .* ((C_mat .* K_cart) .* B_mat);  % opdims [5 4]

sigma_chain = randn(4*ns, 1);
wts4        = repmat(wts, 4, 1);

% Reference: (A*C) * K_cart_mat * B (applied pointwise via block structure)
AC = A_mat * C_mat;                          % 5x3
K_cart_mat3 = K_cart.eval(src_s, targ);     % (3*nt x 3*ns)
nt_c = size(targ.r, 2);

% Build block-diagonal B expansion: (3*ns x 4*ns) and AC contraction: (5*nt x 3*nt)
Bblk  = kron(speye(ns), B_mat);             % (3*ns x 4*ns)
ACblk = kron(speye(nt_c), AC);              % (5*nt x 3*nt)
pot_chain_ref = ACblk * K_cart_mat3 * Bblk * sigma_chain;

% =========================================================================
% 14. eval
% =========================================================================

pot_chain = K_chain.eval(src_s, targ) * sigma_chain;
err = norm(pot_chain(:) - pot_chain_ref(:)) / norm(pot_chain_ref(:));
failures = report(failures, 'chain eval vs manual', err, 1e-13);

% =========================================================================
% 15. layer_eval (off-surface) vs smooth-quad reference
% =========================================================================

sigma_chain_wt = wts4(:) .* sigma_chain;
pot_chain_le     = K_chain.layer_eval(S, reshape(sigma_chain, 4, ns), targ, eps);
pot_chain_le_ref = K_chain.eval(src_s, targ) * sigma_chain_wt;

err = norm(pot_chain_le(:) - pot_chain_le_ref(:)) / norm(pot_chain_le_ref(:));
failures = report(failures, 'chain layer_eval vs smooth quad', err, 1e-4);

% =========================================================================
% 16. layer_eval with precomp_quadrature
% =========================================================================

Q_chain    = K_chain.getquad(S, eps, targ);
opts_chain = struct('precomp_quadrature', Q_chain);
pot_chain_precomp = K_chain.layer_eval(S, reshape(sigma_chain, 4, ns), targ, eps, opts_chain);

err = norm(pot_chain_precomp(:) - pot_chain_le(:)) / norm(pot_chain_le(:));
failures = report(failures, 'chain layer_eval precomp vs plain', err, 1e-10);

% =========================================================================
% --- nrccie-eval kernel: partial ids + multi-column jds ---
% =========================================================================
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
% 7. eval: against manually projected Cartesian matrix
% =========================================================================

scalar_src = 4;                 % rho passes through on source side
scalar_tgt = zeros(0,1);        % no passthrough on target side (all projected)

Ps_exp_em = tangent_block(src_s, ids_em, scalar_src, 3, false);  % (4 x 3 x ns)
Pt_em     = tangent_block(targ,  jds_em, scalar_tgt, 4, true);   % (4 x 6 x nt)

nt_em = size(targ.r, 2);
K_em_mat   = K_em.eval(src_s, targ);  % (6*nt x 4*ns)
sig_em_pages  = reshape(sigma_em, 3, 1, ns);
Ps_sig_em     = reshape(pagemtimes(Ps_exp_em, sig_em_pages), 4*ns, 1);  % (4*ns x 1)
Kps_sig_em    = K_em_mat * Ps_sig_em;                                    % (6*nt x 1)
pot_em_ref    = reshape(pagemtimes(Pt_em, reshape(Kps_sig_em, 6, 1, nt_em)), 4*nt_em, 1);
pot_em_tan = K_em_tan.eval(src_s, targ) * sigma_em;

err = norm(pot_em_tan(:) - pot_em_ref(:)) / norm(pot_em_ref(:));
failures = report(failures, 'nrccie eval vs manual', err, 1e-13);

% =========================================================================
% 8. layer_eval (off-surface): against smooth-quad reference
% =========================================================================

sigma_em_wt = repmat(wts, 3, 1) .* reshape(sigma_em, 3, ns);
pot_em_le     = K_em_tan.layer_eval(S, sigma_em, targ, eps);
pot_em_le_ref = K_em_tan.eval(src_s, targ) * sigma_em_wt(:);

err = norm(pot_em_le(:) - pot_em_le_ref(:)) / norm(pot_em_le_ref(:));
failures = report(failures, 'nrccie layer_eval vs smooth', err, 1e-4);

% =========================================================================
% 9. nrccie layer_eval with precomp_quadrature
% =========================================================================

Q_em    = K_em_tan.getquad(S, eps, targ);
opts_em = struct('precomp_quadrature', Q_em);
pot_em_precomp = K_em_tan.layer_eval(S,sigma_em, targ, eps, opts_em);

% nrccie has non-trivial rsc interleave so precomp falls back to recompute;
% result should still match layer_eval to FMM accuracy.
err = norm(pot_em_precomp(:) - pot_em_le(:)) / norm(pot_em_le(:));
failures = report(failures, 'nrccie layer_eval precomp vs plain', err, 1e-4);

% =========================================================================
% --- n-dot left multiply: f(t) = n(t).', stok s -> scalar-valued ---
%
%   f_ndot(t) = reshape(t.n, 1, 3, nt)   i.e. (1 x 3 x nt)
%   K_ndot = f_ndot .* K_cart            opdims [1, 3]
%   (K_ndot * sigma)(t) = n(t).' * (K_cart * sigma)(t)
% =========================================================================

fprintf('\n--- n-dot left multiply ---\n');

f_ndot  = @(t) reshape(t.n, 1, 3, size(t.r,2));
K_ndot  = f_ndot .* K_cart;   % opdims [1, 3]

sigma_3 = randn(3*ns, 1);
wts3    = repmat(wts, 3, 1);

% =========================================================================
% 10. eval: n(t).' * K_cart(s,t) * sigma
% =========================================================================

K_cart_mat2 = K_cart.eval(src_s, targ);
nt          = size(targ.r, 2);
pot_ndot_ref = zeros(nt, 1);
for ii = 1:nt
    ni = targ.n(:,ii);
    pot_ndot_ref(ii) = ni.' * K_cart_mat2((ii-1)*3+(1:3), :) * sigma_3;
end
pot_ndot = K_ndot.eval(src_s, targ) * sigma_3;

err = norm(pot_ndot(:) - pot_ndot_ref(:)) / norm(pot_ndot_ref(:));
failures = report(failures, 'n-dot eval vs manual', err, 1e-13);

% =========================================================================
% 11. fmm: against eval
% =========================================================================

pot_ndot_fmm = K_ndot.fmm(eps, src_s, targ, sigma_3);
err = norm(pot_ndot_fmm(:) - pot_ndot(:)) / norm(pot_ndot(:));
failures = report(failures, 'n-dot fmm vs eval', err, 1e-5);

% =========================================================================
% 12. layer_eval (off-surface): against smooth-quad reference
% =========================================================================

sigma_3_wt  = wts3(:) .* sigma_3;
pot_ndot_le     = K_ndot.layer_eval(S, sigma_3, targ, eps);
pot_ndot_le_ref = K_ndot.eval(src_s, targ) * sigma_3_wt;

err = norm(pot_ndot_le(:) - pot_ndot_le_ref(:)) / norm(pot_ndot_le_ref(:));
failures = report(failures, 'n-dot layer_eval vs smooth quad', err, 1e-4);

% =========================================================================
% 13. layer_eval with precomp_quadrature
% =========================================================================

Q_ndot    = K_ndot.getquad(S, eps, targ);
opts_ndot = struct('precomp_quadrature', Q_ndot);
pot_ndot_precomp = K_ndot.layer_eval(S, sigma_3, targ, eps, opts_ndot);

err = norm(pot_ndot_precomp(:) - pot_ndot_le(:)) / norm(pot_ndot_le(:));
failures = report(failures, 'n-dot layer_eval precomp vs plain', err, 1e-10);

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

function F = local_tangent_pages(src, vds, scalar_rows, opdim, use_pinv)
%LOCAL_TANGENT_PAGES  Per-point tangent basis pages; thin wrapper around tangent_block.
F = tangent_block(src, vds, scalar_rows, opdim, use_pinv);
end


function failures = report(failures, label, err, tol)
    if err > tol, s = 'FAIL'; else, s = 'PASS'; end
    fprintf('  %-36s  err=%.2e  [%s]\n', label, err, s);
    if strcmp(s, 'FAIL')
        failures{end+1} = sprintf('%s: err=%.2e > tol=%.2e', label, err, tol);
    end
end
