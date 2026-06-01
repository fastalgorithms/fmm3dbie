% Verify kernel3d eval, fmm, arithmetic operators, and zeros
% for all implemented kernels. Failures accumulate; one assert at the end.

run ../../startup.m

rng(42);

zk  = 0.1 + 0.2i;
eps = 1e-6;
ns  = 500;
nt  = 300;

src = [];
src.r  = randn(3,ns);  src.r  = src.r  ./ vecnorm(src.r) .* rand(1,ns);
src.n  = randn(3,ns);  src.n  = src.n  ./ vecnorm(src.n);
src.du = randn(3,ns);  src.dv = randn(3,ns);

targ = [];
targ.r  = randn(3,nt);  targ.r  = targ.r ./ vecnorm(targ.r) .* (2 + rand(1,nt));
targ.n  = randn(3,nt);  targ.n  = targ.n ./ vecnorm(targ.n);
targ.du = randn(3,nt);  targ.dv = randn(3,nt);

S = slicesurfer(geometries.sphere(1, 4, [0;0;0], 8, 1), 1);

targ_off.r  = randn(3,40);  targ_off.r  = targ_off.r  ./ vecnorm(targ_off.r)  * 5.0;
targ_off.n  = randn(3,40);  targ_off.n  = targ_off.n  ./ vecnorm(targ_off.n);
targ_off.du = randn(3,40);  targ_off.du = targ_off.du ./ vecnorm(targ_off.du);
targ_off.dv = randn(3,40);  targ_off.dv = targ_off.dv ./ vecnorm(targ_off.dv);

%% Now run the tests

failures = {};
failures = test_per_kernel(failures, zk, eps, ns, nt, src, targ, S, targ_off);
failures = test_transmission(failures, zk, eps, ns, nt, src, targ, S, targ_off);
failures = test_stok_sp_identity(failures, src, targ);
failures = test_combined_layer(failures, zk, ns, src, targ);
failures = test_arithmetic(failures, zk, eps, src, targ, S, targ_off);
failures = test_zeros(failures, eps, ns, nt, src, targ);

nfail = numel(failures);
if nfail > 0
    for k = 1:nfail, fprintf('  FAIL: %s\n', failures{k}); end
    assert(false, 'test_kernel3d: %d failure(s)', nfail);
end


% =========================================================================

function failures = test_per_kernel(failures, zk, eps, ~, ~, src, targ, S, targ_off)
% FD-prime and FMM checks for all basic kernels.

sigma_r  = randn(S.npts, 1);
sigma_s  = randn(S.npts, 1) + 1i*randn(S.npts, 1);
sigma_stok = randn(3, S.npts);
sigma_em3  = randn(3, S.npts) + 1i*randn(3, S.npts);
sigma_em4  = randn(4, S.npts) + 1i*randn(4, S.npts);

sigma_r_pt    = randn(size(src.r,2), 1);
sigma_pt      = randn(size(src.r,2), 1) + 1i*randn(size(src.r,2), 1);
sigma_stok_pt = randn(3*size(src.r,2), 1);
sigma_em3_pt  = randn(3*size(src.r,2), 1) + 1i*randn(3*size(src.r,2), 1);
sigma_em4_pt  = randn(4*size(src.r,2), 1) + 1i*randn(4*size(src.r,2), 1);

tol_fd = 1e-3;  tol_fmm = 1e-5;  h_fd = 1e-5;
coefs_lap  = [0.7; 1.3];
coefs_h    = [1i*zk; 1.0];
coefs_ct   = [1i*zk; 1.0];
coefs_stok = [0.7; 1.3];
zk_em = 2.5 + 0.1i;  alpha_em = 0.5;

kernels = { ...
  kernel3d('l','s'),              'lap s',       sigma_r,    sigma_r_pt,   [];                          ...
  kernel3d('l','d'),              'lap d',       sigma_r,    sigma_r_pt,   [];                          ...
  kernel3d('l','sp'),             'lap sp',      sigma_r,    sigma_r_pt,   kernel3d('l','s');            ...
  kernel3d('l','dp'),             'lap dp',      sigma_r,    sigma_r_pt,   kernel3d('l','d');            ...
  kernel3d('l','c',coefs_lap),    'lap c',       sigma_r,    sigma_r_pt,   [];                          ...
  kernel3d('l','cp',coefs_lap),   'lap cp',      sigma_r,    sigma_r_pt,   kernel3d('l','c',coefs_lap);  ...
  kernel3d('h','s',zk),           'helm s',      sigma_s,    sigma_pt,     [];                          ...
  kernel3d('h','d',zk),           'helm d',      sigma_s,    sigma_pt,     [];                          ...
  kernel3d('h','sp',zk),          'helm sp',     sigma_s,    sigma_pt,     kernel3d('h','s',zk);         ...
  kernel3d('h','dp',zk),          'helm dp',     sigma_s,    sigma_pt,     kernel3d('h','d',zk);         ...
  kernel3d('h','c',zk,coefs_h),   'helm c',      sigma_s,    sigma_pt,     [];                          ...
  kernel3d('h','cp',zk,coefs_h),  'helm cp',     sigma_s,    sigma_pt,     kernel3d('h','c',zk,coefs_h); ...
  kernel3d('h','s2trans',zk),     'helm s2trans',sigma_s,    sigma_s,      kernel3d('h','s',zk);         ...
  kernel3d('h','d2trans',zk),     'helm d2trans',sigma_s,    sigma_s,      kernel3d('h','d',zk);         ...
  kernel3d('h','c2trans',zk,coefs_ct),'helm c2trans',sigma_s,sigma_s,     kernel3d('h','c',zk,coefs_ct); ...
  kernel3d('stok','s'),           'stok s',      sigma_stok, sigma_stok_pt, [];                         ...
  kernel3d('stok','d'),           'stok d',      sigma_stok, sigma_stok_pt, [];                         ...
  kernel3d('stok','sp'),          'stok sp',     sigma_stok, sigma_stok_pt, [];                         ...
  kernel3d('stok','c',coefs_stok),'stok c',      sigma_stok, sigma_stok_pt, [];                         ...
  kernel3d('em','nrccie-bc',zk_em,alpha_em),  'em nrccie-bc',  sigma_em3, sigma_em3_pt, []; ...
  kernel3d('em','nrccie-eval',zk_em),         'em nrccie-eval', sigma_em4, sigma_em4_pt, []; ...
};

for k = 1:size(kernels, 1)
    K      = kernels{k,1};
    label  = kernels{k,2};
    sig_s  = kernels{k,3};
    sig_pt = kernels{k,4};
    K_aux  = kernels{k,5};

    is_x2trans = ismember(K.type, {'s2trans','d2trans','c2trans'});
    if is_x2trans
        switch K.type
            case 's2trans',  K_bot = kernel3d('h','sp',zk);
            case 'd2trans',  K_bot = kernel3d('h','dp',zk);
            case 'c2trans',  K_bot = kernel3d('h','cp',zk,coefs_ct);
        end
        failures = check_x2trans(failures, K, K_aux, K_bot, src, targ, label);
    else
        if ~isempty(K_aux)
            failures = check_fd_targ(failures, K, K_aux, src, targ, h_fd, tol_fd, label);
        end
        if ~isempty(K.fmm) && ~isempty(sig_pt)
            failures = check_fmm(failures, K, src, targ, sig_pt, eps, tol_fmm, label);
        end
    end
end

end


function failures = test_transmission(failures, zk, eps, ~, ~, src, targ, S, targ_off)
% Scalar diff kernels, trans_rep, trans_sys, trans_sys_diff, x2trans_diff.

zk1        = 0.3 + 0.15i;
diff_coefs = [1.2; 0.8];
coefs_ct   = [1i*zk; 1.0];
coe_22     = [0.7+0.1i, 0.3; 0.5, 1.2i];
coe_222    = cat(3, coe_22, 0.9*coe_22);
c_coefs4   = [coefs_ct; diff_coefs];
zks        = [zk; zk1];
sigma_2d   = randn(2*S.npts, 1) + 1i*randn(2*S.npts, 1);

% scalar diff kernels
for kd = { {'s_diff','s',{}}; {'d_diff','d',{}}; {'sp_diff','sp',{}}; {'dp_diff','dp',{}} }
    kd = kd{1};
    Kdiff = kernel3d('h', kd{1}, zks, diff_coefs);
    Kref  = diff_coefs(1).*kernel3d('h',kd{2},zk,kd{3}{:}) - diff_coefs(2).*kernel3d('h',kd{2},zk1,kd{3}{:});
    failures = report(failures, ['helm ' kd{1}], 'eval vs explicit diff', ...
        norm(Kdiff.eval(src,targ)-Kref.eval(src,targ),'fro')/norm(Kref.eval(src,targ),'fro'), 1e-13);
end
for kd = { {'c_diff','c'}; {'cp_diff','cp'} }
    kd = kd{1};
    Kdiff = kernel3d('h', kd{1}, zks, c_coefs4);
    Kref  = diff_coefs(1).*kernel3d('h',kd{2},zk,coefs_ct) - diff_coefs(2).*kernel3d('h',kd{2},zk1,coefs_ct);
    failures = report(failures, ['helm ' kd{1}], 'eval vs explicit diff', ...
        norm(Kdiff.eval(src,targ)-Kref.eval(src,targ),'fro')/norm(Kref.eval(src,targ),'fro'), 1e-13);
end

% trans_rep
Ktr  = kernel3d('h','trans_rep',zk,[0.7;1.3]);
M_tr = Ktr.eval(src,targ);
M_D  = 0.7*kernel3d('h','d',zk).eval(src,targ);
M_S  = 1.3*kernel3d('h','s',zk).eval(src,targ);
failures = report(failures, 'helm trans_rep', 'eval blocks vs D, S', ...
    max(norm(M_tr(:,1:2:end)-M_D,'fro')/norm(M_D,'fro'), norm(M_tr(:,2:2:end)-M_S,'fro')/norm(M_S,'fro')), 1e-13);

% trans_sys
Kts = kernel3d('h','trans_sys',zk,coe_22);
M_ts = Kts.eval(src,targ);
M_D  = kernel3d('h','d', zk).eval(src,targ);  M_S  = kernel3d('h','s', zk).eval(src,targ);
M_Dp = kernel3d('h','dp',zk).eval(src,targ);  M_Sp = kernel3d('h','sp',zk).eval(src,targ);
e = max([norm(M_ts(1:2:end,1:2:end)-coe_22(1,1)*M_D, 'fro')/norm(coe_22(1,1)*M_D,'fro'), ...
         norm(M_ts(1:2:end,2:2:end)-coe_22(1,2)*M_S, 'fro')/norm(coe_22(1,2)*M_S,'fro'), ...
         norm(M_ts(2:2:end,1:2:end)-coe_22(2,1)*M_Dp,'fro')/norm(coe_22(2,1)*M_Dp,'fro'), ...
         norm(M_ts(2:2:end,2:2:end)-coe_22(2,2)*M_Sp,'fro')/norm(coe_22(2,2)*M_Sp,'fro')]);
failures = report(failures, 'helm trans_sys', 'eval 2x2 blocks', e, 1e-13);

% trans_sys_diff
Ktsd = kernel3d('h','trans_sys_diff',zks,coe_222);
M_tsd = Ktsd.eval(src,targ);
c0 = coe_222(:,:,1);  c1 = coe_222(:,:,2);
M_D0 = kernel3d('h','d', zk ).eval(src,targ);  M_D1 = kernel3d('h','d', zk1).eval(src,targ);
M_S0 = kernel3d('h','s', zk ).eval(src,targ);  M_S1 = kernel3d('h','s', zk1).eval(src,targ);
M_Dp0= kernel3d('h','dp',zk ).eval(src,targ);  M_Dp1= kernel3d('h','dp',zk1).eval(src,targ);
M_Sp0= kernel3d('h','sp',zk ).eval(src,targ);  M_Sp1= kernel3d('h','sp',zk1).eval(src,targ);
f = max([norm(M_tsd(1:2:end,1:2:end)-(c0(1,1)*M_D0-c1(1,1)*M_D1),'fro')/norm(c0(1,1)*M_D0-c1(1,1)*M_D1,'fro'), ...
         norm(M_tsd(1:2:end,2:2:end)-(c0(1,2)*M_S0-c1(1,2)*M_S1),'fro')/norm(c0(1,2)*M_S0-c1(1,2)*M_S1,'fro'), ...
         norm(M_tsd(2:2:end,1:2:end)-(c0(2,1)*M_Dp0-c1(2,1)*M_Dp1),'fro')/norm(c0(2,1)*M_Dp0-c1(2,1)*M_Dp1,'fro'), ...
         norm(M_tsd(2:2:end,2:2:end)-(c0(2,2)*M_Sp0-c1(2,2)*M_Sp1),'fro')/norm(c0(2,2)*M_Sp0-c1(2,2)*M_Sp1,'fro')]);
failures = report(failures, 'helm trans_sys_diff', 'eval 2x2 blocks', f, 1e-13);

% x2trans_diff
for kd = { {'s2trans_diff','s_diff','sp_diff'}; {'d2trans_diff','d_diff','dp_diff'} }
    kd = kd{1};
    failures = check_x2trans_diff(failures, ...
        kernel3d('h',kd{1},zks,diff_coefs), kernel3d('h',kd{2},zks,diff_coefs), ...
        kernel3d('h',kd{3},zks,diff_coefs), src, targ, ['helm ' kd{1}]);
end
failures = check_x2trans_diff(failures, ...
    kernel3d('h','c2trans_diff',zks,c_coefs4), kernel3d('h','c_diff',zks,c_coefs4), ...
    kernel3d('h','cp_diff',zks,c_coefs4), src, targ, 'helm c2trans_diff');

end


function failures = test_stok_sp_identity(failures, src, targ)
% Stokes S': eval == D(targ,src)^T.

K_sp = kernel3d('stok','sp');
K_d  = kernel3d('stok','d');
err  = norm(K_sp.eval(src,targ) - K_d.eval(targ,src).', 'fro') / norm(K_sp.eval(src,targ),'fro');
failures = report(failures, 'stok sp', 'eval == D(targ,src)^T', err, 1e-12);

end


function failures = test_combined_layer(failures, zk, ns, src, targ)
% C.eval == coefs(1)*S.eval + coefs(2)*D.eval for Laplace, Helmholtz, Stokes.

coefs_lap  = [0.7; 1.3];
coefs_h    = [1i*zk; 1.0];
coefs_stok = [0.7; 1.3];

sigma_r = randn(ns, 1);
sigma   = randn(ns, 1) + 1i*randn(ns, 1);

combined = { ...
  kernel3d('l','c',coefs_lap),     kernel3d('l','s'),     kernel3d('l','d'),     coefs_lap,  sigma_r,       'lap c';  ...
  kernel3d('h','c',zk,coefs_h),    kernel3d('h','s',zk),  kernel3d('h','d',zk),  coefs_h,    sigma,         'helm c'; ...
  kernel3d('stok','c',coefs_stok), kernel3d('stok','s'),   kernel3d('stok','d'),  coefs_stok, randn(3*ns,1), 'stok c'; ...
};
for k = 1:size(combined,1)
    failures = check_combined(failures, combined{k,:}, src, targ);
end

end


function failures = test_arithmetic(failures, zk, eps, src, targ, ~, ~)
% Composite kernel expression tests eval and fmm.

tol_fmm = 1e-5;
coefs_h  = [1i*zk; 1.0];
ns       = size(src.r, 2);
sigma_r  = randn(ns, 1);

K_hs = kernel3d('h','s',zk);  K_hd = kernel3d('h','d',zk);
K_ls = kernel3d('l','s');      K_hc = kernel3d('h','c',zk,coefs_h);
c1 = 2.0 + 0.5i;  c2 = 1.5;

Kexpr = ((-K_hs).*c1 + K_hd - (K_ls+K_hc)/3) / c2;

pot_expr = Kexpr.eval(src,targ)*sigma_r;
pot_ref  = ((-c1)*K_hs.eval(src,targ)*sigma_r + K_hd.eval(src,targ)*sigma_r ...
           - (K_ls.eval(src,targ)*sigma_r + K_hc.eval(src,targ)*sigma_r)/3) / c2;
failures = report(failures, 'arithmetic', 'composite eval',       norm(pot_expr-pot_ref)/norm(pot_ref), 1e-13);

sigma_r_pt = randn(size(src.r,2), 1);
pot_fmm  = Kexpr.fmm(eps,src,targ,sigma_r_pt);
pot_fref = ((-c1)*K_hs.fmm(eps,src,targ,sigma_r_pt) + K_hd.fmm(eps,src,targ,sigma_r_pt) ...
           - (K_ls.fmm(eps,src,targ,sigma_r_pt) + K_hc.fmm(eps,src,targ,sigma_r_pt))/3) / c2;
failures = report(failures, 'arithmetic', 'composite fmm',        norm(pot_fmm-pot_fref)/norm(pot_fref), tol_fmm);

end


function failures = test_zeros(failures, eps, ns, nt, src, targ)
% kernel3d.zeros: correct opdims, eval size, fmm size.

for od = {[1 1], [3 3], [2 3]}
    opdims = od{1};
    Kz = kernel3d.zeros(opdims);
    assert(isequal(Kz.opdims, opdims) && Kz.iszero, 'zeros: wrong opdims or iszero');
    assert(isequal(size(Kz.eval(src,targ)),         [opdims(1)*nt, opdims(2)*ns]));
    assert(isequal(size(Kz.fmm(eps,src,targ,zeros(opdims(2)*ns,1))), [opdims(1)*nt, 1]));
end

end


% =========================================================================
% Helpers
% =========================================================================

function failures = check_fmm(failures, K, src, targ, sigma, eps, tol, label)
% Compare K.fmm against K.eval.
pot_fmm = K.fmm(eps, src, targ, sigma);
pot_ref = K.eval(src, targ) * sigma(:);
failures = report(failures, label, 'fmm vs eval', norm(pot_fmm(:)-pot_ref(:))/norm(pot_ref(:)), tol);
end

function failures = check_fd_targ(failures, K_prime, K_base, src, targ, h, tol, label)
% Verify prime kernel eval against centred FD in n_targ.
targ_p = targ;  targ_p.r = targ.r + h*targ.n;
targ_m = targ;  targ_m.r = targ.r - h*targ.n;
mat_fd = (K_base.eval(src,targ_p) - K_base.eval(src,targ_m)) / (2*h);
failures = report(failures, label, 'prime vs FD(targ)', norm(K_prime.eval(src,targ)-mat_fd,'fro')/norm(mat_fd,'fro'), tol);
end

function failures = check_combined(failures, Kc, Ks, Kd, coefs, sigma, label, src, targ)
% Verify C.eval == coefs(1)*S.eval + coefs(2)*D.eval.
pot_c   = Kc.eval(src,targ)*sigma;
pot_ref = coefs(1)*Ks.eval(src,targ)*sigma + coefs(2)*Kd.eval(src,targ)*sigma;
failures = report(failures, label, 'C = a*S + b*D', norm(pot_c-pot_ref)/norm(pot_ref), 1e-13);
end

function failures = check_x2trans(failures, K2t, K_top, K_bot, src, targ, label)
% Verify stacked x2trans: odd rows -> K_top, even rows -> K_bot.
mat = K2t.eval(src,targ);
e = max(norm(mat(1:2:end,:)-K_top.eval(src,targ),'fro')/norm(K_top.eval(src,targ),'fro'), ...
        norm(mat(2:2:end,:)-K_bot.eval(src,targ),'fro')/norm(K_bot.eval(src,targ),'fro'));
failures = report(failures, label, 'eval blocks', e, 1e-14);
end

function failures = check_x2trans_diff(failures, K2t, K_top, K_bot, src, targ, label)
% Like check_x2trans but K_top/K_bot are already difference kernels.
mat = K2t.eval(src,targ);
e = max(norm(mat(1:2:end,:)-K_top.eval(src,targ),'fro')/norm(K_top.eval(src,targ),'fro'), ...
        norm(mat(2:2:end,:)-K_bot.eval(src,targ),'fro')/norm(K_bot.eval(src,targ),'fro'));
failures = report(failures, label, 'eval blocks', e, 1e-13);
end

function failures = report(failures, label, desc, err, tol)
if err > tol
    failures{end+1} = sprintf('%s %s err=%.2e', label, desc, err);
end
fprintf('  %-24s  %-22s  err=%.2e  [%s]\n', label, desc, err, ternary(err>tol,'FAIL','PASS'));
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
