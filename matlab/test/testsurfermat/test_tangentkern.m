% Test kernel3d.tangent_kern and kernel3d/times matrix-fn multiply.
% Checks stok s (full projection), chained multiply, nrccie-eval (partial
% ids/jds), and n-dot left multiply across eval, fmm, and layer_eval.

run ../../startup.m

rng(42);
eps = 1e-6;
zk  = 1.3;

S      = slicesurfer(geometries.sphere(1, 4, [0;0;0], 8, 1), 1);
ns     = S.npts;
wts    = S.wts(:).';
src_s  = struct('r', S.r(:,:), 'n', S.n(:,:), 'du', S.du(:,:), 'dv', S.dv(:,:));

targ    = [];
targ.r  = randn(3, 200);  targ.r  = targ.r  ./ vecnorm(targ.r) * 3.0;
targ.n  = randn(3, 200);  targ.n  = targ.n  ./ vecnorm(targ.n);
targ.du = randn(3, 200);  targ.dv = randn(3, 200);

%% Now run the tests

failures = {};
failures = test_stok_projection(failures, S, ns, wts, src_s, targ, eps);
failures = test_chained_multiply(failures, S, ns, wts, src_s, targ, eps);
failures = test_nrccie_projection(failures, S, ns, wts, src_s, targ, eps, zk);
failures = test_ndot_multiply(failures, S, ns, wts, src_s, targ, eps);

nfail = numel(failures);
if nfail > 0
    for k = 1:nfail, fprintf('  FAIL: %s\n', failures{k}); end
    assert(false, 'test_tangentkern: %d failure(s)', nfail);
end


function failures = test_stok_projection(failures, S, ns, wts, src_s, targ, eps)
% stok s, full 3->2 tangent projection on both src and targ:
% eval (tangent_kern and direct .* construction), fmm, layer_eval, on-surface PV.

K_cart = kernel3d('stok', 's');
ids = [1;2;3];  jds = [1;2;3];
K_tan = kernel3d.tangent_kern(K_cart, ids, jds);

f_right      = @(s) local_tangent_pages(s, ids, zeros(0,1), 2, false);
f_left       = @(t) local_tangent_pages(t, jds, zeros(0,1), 2, true);
K_tan_direct = f_left .* K_cart .* f_right;

sigma_tan = randn(2*ns, 1);

Ps_exp  = tangent_block(src_s, ids, zeros(0,1), 2, false);
Pt      = tangent_block(targ,  jds, zeros(0,1), 2, true);
nt      = size(targ.r, 2);
K_mat   = K_cart.eval(src_s, targ);
Ps_sig  = reshape(pagemtimes(Ps_exp, reshape(sigma_tan, 2, 1, ns)), 3*ns, 1);
pot_ref = reshape(pagemtimes(Pt, reshape(K_mat*Ps_sig, 3, 1, nt)), 2*nt, 1);

failures = report(failures, 'stok eval vs manual',        norm(K_tan.eval(src_s,targ)*sigma_tan        - pot_ref)/norm(pot_ref), 1e-13);
failures = report(failures, 'stok direct times vs manual', norm(K_tan_direct.eval(src_s,targ)*sigma_tan - pot_ref)/norm(pot_ref), 1e-13);
failures = report(failures, 'stok fmm vs eval',           norm(K_tan.fmm(eps,src_s,targ,sigma_tan) - K_tan.eval(src_s,targ)*sigma_tan)/norm(K_tan.eval(src_s,targ)*sigma_tan), 1e-5);

sigma_wt = repmat(wts, 2, 1) .* reshape(sigma_tan, 2, ns);
pot_le   = K_tan.layer_eval(S, sigma_tan, targ, eps);
pot_le_ref = K_tan.eval(src_s, targ) * sigma_wt(:);
failures = report(failures, 'stok layer_eval vs smooth quad', norm(pot_le(:)-pot_le_ref(:))/norm(pot_le_ref(:)), 1e-4);

Q = K_tan.getquad(S, eps, targ);
pot_precomp = K_tan.layer_eval(S, sigma_tan, targ, eps, struct('precomp_quadrature', Q));
failures = report(failures, 'stok layer_eval precomp vs plain', norm(pot_precomp(:)-pot_le(:))/norm(pot_le(:)), 1e-10);

% on-surface: stok s has no jump, so PV value = avg of one-sided limits
h = 1e-4;  tol_jump = 0.05;
[~, ic] = min(vecnorm(S.r - mean(S.r,2)));
n0 = S.n(:,ic)/norm(S.n(:,ic));  r0 = S.r(:,ic);  du0 = S.du(:,ic);  dv0 = S.dv(:,ic);
te.r = r0+h*n0; te.n = n0; te.du = du0; te.dv = dv0;
ti.r = r0-h*n0; ti.n = n0; ti.du = du0; ti.dv = dv0;
sigma_tan2 = reshape(sigma_tan, 2, ns);
u_avg      = (K_tan.layer_eval(S,sigma_tan2,te,eps) + K_tan.layer_eval(S,sigma_tan2,ti,eps)) / 2;
pot_pv_full = surfermatapply(S, K_tan, sigma_tan, eps);
od = K_tan.opdims(1);
u_pv = pot_pv_full((ic-1)*od+1 : ic*od);
failures = report(failures, 'stok on-surface PV vs avg', norm(u_pv(:)-u_avg(:))/(norm(u_avg(:))+1e-30), tol_jump);

end


function failures = test_chained_multiply(failures, S, ns, wts, src_s, targ, eps)
% Chained multiply A.*((C.*K).*B): eval and layer_eval match manual reference.

K_cart = kernel3d('stok', 's');
rng(7);
B_mat = randn(3, 4);
C_mat = randn(2, 3);
A_mat = randn(5, 2);
K_chain = A_mat .* ((C_mat .* K_cart) .* B_mat);

sigma_chain = randn(4*ns, 1);
wts4        = repmat(wts, 4, 1);
nt          = size(targ.r, 2);
AC          = A_mat * C_mat;
Bblk        = kron(speye(ns), B_mat);
ACblk       = kron(speye(nt), AC);
pot_ref     = ACblk * K_cart.eval(src_s, targ) * Bblk * sigma_chain;

failures = report(failures, 'chain eval vs manual',         norm(K_chain.eval(src_s,targ)*sigma_chain - pot_ref)/norm(pot_ref), 1e-13);

pot_le     = K_chain.layer_eval(S, reshape(sigma_chain, 4, ns), targ, eps);
pot_le_ref = K_chain.eval(src_s, targ) * (wts4(:) .* sigma_chain);
failures = report(failures, 'chain layer_eval vs smooth quad', norm(pot_le(:)-pot_le_ref(:))/norm(pot_le_ref(:)), 1e-4);

Q           = K_chain.getquad(S, eps, targ);
pot_precomp = K_chain.layer_eval(S, reshape(sigma_chain, 4, ns), targ, eps, struct('precomp_quadrature', Q));
failures = report(failures, 'chain layer_eval precomp vs plain', norm(pot_precomp(:)-pot_le(:))/norm(pot_le(:)), 1e-10);

end


function failures = test_nrccie_projection(failures, S, ns, wts, src_s, targ, eps, zk)
% nrccie-eval, partial src ids + multi-column tgt jds: eval and layer_eval.

K_em     = kernel3d('em', 'nrccie-eval', zk);
ids_em   = [1;2;3];
jds_em   = [[1;2;3], [4;5;6]];
K_em_tan = kernel3d.tangent_kern(K_em, ids_em, jds_em);

sigma_em = randn(3*ns, 1);

Ps_exp_em = tangent_block(src_s, ids_em, 4,          3, false);
Pt_em     = tangent_block(targ,  jds_em, zeros(0,1), 4, true);
nt        = size(targ.r, 2);
Ps_sig_em = reshape(pagemtimes(Ps_exp_em, reshape(sigma_em, 3, 1, ns)), 4*ns, 1);
pot_em_ref = reshape(pagemtimes(Pt_em, reshape(K_em.eval(src_s,targ)*Ps_sig_em, 6, 1, nt)), 4*nt, 1);

failures = report(failures, 'nrccie eval vs manual', norm(K_em_tan.eval(src_s,targ)*sigma_em - pot_em_ref)/norm(pot_em_ref), 1e-13);

sigma_em_wt = repmat(wts, 3, 1) .* reshape(sigma_em, 3, ns);
pot_em_le   = K_em_tan.layer_eval(S, sigma_em, targ, eps);
failures = report(failures, 'nrccie layer_eval vs smooth', norm(pot_em_le(:) - K_em_tan.eval(src_s,targ)*sigma_em_wt(:))/norm(K_em_tan.eval(src_s,targ)*sigma_em_wt(:)), 1e-4);

Q_em        = K_em_tan.getquad(S, eps, targ);
pot_em_pre  = K_em_tan.layer_eval(S, sigma_em, targ, eps, struct('precomp_quadrature', Q_em));
failures = report(failures, 'nrccie layer_eval precomp vs plain', norm(pot_em_pre(:)-pot_em_le(:))/norm(pot_em_le(:)), 1e-4);

end


function failures = test_ndot_multiply(failures, S, ns, wts, src_s, targ, eps)
% n-dot left multiply f(t)=n(t).': eval, fmm, layer_eval match manual reference.

K_cart = kernel3d('stok', 's');
f_ndot = @(t) reshape(t.n, 1, 3, size(t.r,2));
K_ndot = f_ndot .* K_cart;

sigma_3 = randn(3*ns, 1);
wts3    = repmat(wts, 3, 1);
nt      = size(targ.r, 2);
K_mat   = K_cart.eval(src_s, targ);
pot_ref = zeros(nt, 1);
for ii = 1:nt
    pot_ref(ii) = targ.n(:,ii).' * K_mat((ii-1)*3+(1:3), :) * sigma_3;
end

failures = report(failures, 'n-dot eval vs manual',       norm(K_ndot.eval(src_s,targ)*sigma_3 - pot_ref)/norm(pot_ref), 1e-13);
failures = report(failures, 'n-dot fmm vs eval',          norm(K_ndot.fmm(eps,src_s,targ,sigma_3) - K_ndot.eval(src_s,targ)*sigma_3)/norm(K_ndot.eval(src_s,targ)*sigma_3), 1e-5);

pot_le     = K_ndot.layer_eval(S, sigma_3, targ, eps);
pot_le_ref = K_ndot.eval(src_s, targ) * (wts3(:) .* sigma_3);
failures = report(failures, 'n-dot layer_eval vs smooth quad', norm(pot_le(:)-pot_le_ref(:))/norm(pot_le_ref(:)), 1e-4);

Q           = K_ndot.getquad(S, eps, targ);
pot_precomp = K_ndot.layer_eval(S, sigma_3, targ, eps, struct('precomp_quadrature', Q));
failures = report(failures, 'n-dot layer_eval precomp vs plain', norm(pot_precomp(:)-pot_le(:))/norm(pot_le(:)), 1e-10);

end


function F = local_tangent_pages(src, vds, scalar_rows, opdim, use_pinv)
F = tangent_block(src, vds, scalar_rows, opdim, use_pinv);
end

function failures = report(failures, label, err, tol)
if err > tol
    failures{end+1} = sprintf('%s: err=%.2e > tol=%.2e', label, err, tol);
end
fprintf('  %-36s  err=%.2e  [%s]\n', label, err, ternary(err>tol,'FAIL','PASS'));
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
