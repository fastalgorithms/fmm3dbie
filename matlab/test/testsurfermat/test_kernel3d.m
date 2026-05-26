%TEST_KERNEL3D  Verify kernel3d eval and layer_eval agree for all
% implemented kernels, and that arithmetic operators (+, -, .*, ./) work.
%
% For prime kernels, eval is checked against a finite-difference of the
% unprimed kernel in the target normal direction.
% Combined kernels are checked against coefs(1)*S + coefs(2)*D.
% FMM is checked against eval for Laplace, Helmholtz, and Stokes.
% The arithmetic section uses a single composite kernel expression.
% Transmission kernels:
%   trans_sys / trans_sys_diff: eval checked against 4 scalar sub-kernels
%     assembled into the 2x2 block.  trans_sys also gets a layer_eval check.
%   trans_rep: eval checked against D+(ik)*S applied to split densities.
%   x2trans / x2trans_diff: eval rows checked against the two scalar halves.
%
% Failures are accumulated and reported at the end.  The script asserts
% only once, at the very end, so every kernel is always exercised.

run ../../startup.m

rng(42);

zk  = 0.1 + 0.2i;
eps = 1e-6;
ns  = 500;
nt  = 300;

% Sources inside unit ball, targets in annulus [2,3] — no singularities.
src.r = randn(3,ns);  src.r = src.r ./ vecnorm(src.r) .* rand(1,ns);
src.n = randn(3,ns);  src.n = src.n ./ vecnorm(src.n);
src.du = randn(3,ns); src.dv = randn(3,ns);  % tangent vectors for EM kernels

targ.r = randn(3,nt);  targ.r = targ.r ./ vecnorm(targ.r) .* (2 + rand(1,nt));
targ.n = randn(3,nt);  targ.n = targ.n ./ vecnorm(targ.n);
targ.du = randn(3,nt); targ.dv = randn(3,nt);

sigma_r = randn(ns,1);                   % real density, for Laplace/Stokes
sigma   = randn(ns,1) + 1i*randn(ns,1); % complex density, for Helmholtz

% Surface for layer_eval tests.  Order 8 keeps the smooth-quadrature
% reference below tol_leval for both Helmholtz and the hypersingular dp kernel.
S = geometries.sphere(1, 4, [0;0;0], 8, 1);
S = slicesurfer(S, 1);

sigma_sr   = randn(S.npts, 1);                        % real surface density
sigma_s    = randn(S.npts, 1) + 1i*randn(S.npts, 1); % complex surface density
sigma_stok = randn(3, S.npts);                        % Stokes: (3, npts)
sigma_em3  = randn(3, S.npts) + 1i*randn(3, S.npts); % em nrccie-bc: (3, npts) [j_ru, j_rv, rho] (orthonormal frame)
sigma_em4  = randn(4, S.npts) + 1i*randn(4, S.npts); % em nrccie-eval: (4, npts) [Jx, Jy, Jz, rho]

% Off-surface targets at radius 5.
targ_off.r  = randn(3,40);
targ_off.r  = targ_off.r ./ vecnorm(targ_off.r) * 5.0;
targ_off.n  = randn(3,40);  targ_off.n  = targ_off.n  ./ vecnorm(targ_off.n);
targ_off.du = randn(3,40);  targ_off.du = targ_off.du ./ vecnorm(targ_off.du);
targ_off.dv = randn(3,40);  targ_off.dv = targ_off.dv ./ vecnorm(targ_off.dv);

tol_leval = 1e-4;  % layer_eval vs smooth-quadrature reference
tol_fd    = 1e-3;  % prime eval vs finite difference
tol_fmm   = 1e-5;  % fmm vs eval  (matches eps=1e-6 with margin)

h_fd = 1e-5;       % finite-difference step

% Accumulate failure messages — assert once at the very end.
failures = {};

% =========================================================================
% Kernel table
% =========================================================================
%
% Columns:
%   1  kernel object
%   2  label string
%   3  surface density   (for check_layer)
%   4  point density     (for check_fmm)
%   5  base kernel       ([] = not a prime; non-empty -> run check_fd_targ)
%
% Stokes point densities are (3,ns); all others are column vectors.
% Transmission x2trans kernels use check_x2trans instead of check_layer/fmm/fd.
%
coefs_lap  = [0.7; 1.3];
coefs_h    = [1i*zk; 1.0];
coefs_ct   = [1i*zk; 1.0];
coefs_stok = [0.7; 1.3];
zk1        = 0.3 + 0.15i;   % second wavenumber for diff kernels
zk_em      = 2.5 + 0.1i;
alpha_em   = 0.5;

sigma_stok_pt = randn(3, ns);   % Stokes point density for FMM check
sigma_em3_pt  = randn(3, ns) + 1i*randn(3,ns);
sigma_em4_pt  = randn(4, ns) + 1i*randn(4, ns);

% Transmission 2-density densities (interleaved [rho; sigma] columns).
sigma_2d   = randn(2*S.npts, 1) + 1i*randn(2*S.npts, 1);  % for trans_sys layer_eval
sigma_2d_pt = randn(2*ns, 1)   + 1i*randn(2*ns, 1);        % point version

kernels = { ...
% --- Laplace ---
  kernel3d('l','s'),              'lap s',       sigma_sr,   sigma_r,         [];                           ...
  kernel3d('l','d'),              'lap d',       sigma_sr,   sigma_r,         [];                           ...
  kernel3d('l','sp'),             'lap sp',      sigma_sr,   sigma_r,         kernel3d('l','s');             ...
  kernel3d('l','dp'),             'lap dp',      sigma_sr,   sigma_r,         kernel3d('l','d');             ...
  kernel3d('l','c',coefs_lap),    'lap c',       sigma_sr,   sigma_r,         [];                           ...
  kernel3d('l','cp',coefs_lap),   'lap cp',      sigma_sr,   sigma_r,         kernel3d('l','c',coefs_lap);   ...
% --- Helmholtz ---
  kernel3d('h','s',zk),           'helm s',      sigma_s,    sigma,           [];                           ...
  kernel3d('h','d',zk),           'helm d',      sigma_s,    sigma,           [];                           ...
  kernel3d('h','sp',zk),          'helm sp',     sigma_s,    sigma,           kernel3d('h','s',zk);          ...
  kernel3d('h','dp',zk),          'helm dp',     sigma_s,    sigma,           kernel3d('h','d',zk);          ...
  kernel3d('h','c',zk,coefs_h),   'helm c',      sigma_s,    sigma,           [];                           ...
  kernel3d('h','cp',zk,coefs_h),  'helm cp',     sigma_s,    sigma,           kernel3d('h','c',zk,coefs_h);  ...
% --- Transmission x2trans ---
  kernel3d('h','s2trans',zk),     'helm s2trans',sigma_s,    sigma_s,         kernel3d('h','s',zk);          ...
  kernel3d('h','d2trans',zk),     'helm d2trans',sigma_s,    sigma_s,         kernel3d('h','d',zk);          ...
  kernel3d('h','c2trans',zk,coefs_ct),'helm c2trans',sigma_s,sigma_s,        kernel3d('h','c',zk,coefs_ct);  ...
% --- Stokes ---
  kernel3d('stok','s'),           'stok s',      sigma_stok, sigma_stok_pt,   [];                           ...
  kernel3d('stok','d'),           'stok d',      sigma_stok, sigma_stok_pt,   [];                           ...
  kernel3d('stok','sp'),          'stok sp',     sigma_stok, sigma_stok_pt,   [];                           ...
  kernel3d('stok','c',coefs_stok),'stok c',      sigma_stok, sigma_stok_pt,   [];                           ...
% --- Maxwell ---
  kernel3d('em','nrccie-bc',zk_em,alpha_em),   'em nrccie-bc',   sigma_em3, sigma_em3_pt, [];  ...
  kernel3d('em','nrccie-eval',zk_em),           'em nrccie-eval', sigma_em4, sigma_em4_pt, [];  ...
};

% =========================================================================
% Per-kernel checks
% =========================================================================

fprintf('=== Per-kernel checks ===\n');

for k = 1:size(kernels, 1)
    K      = kernels{k,1};
    label  = kernels{k,2};
    sig_s  = kernels{k,3};
    sig_pt = kernels{k,4};
    K_aux  = kernels{k,5};  % base kernel (FD check) or top-half kernel (x2trans)

    is_x2trans = ismember(K.type, {'s2trans','d2trans','c2trans'});

    if is_x2trans
        K_top = K_aux;
        switch K.type
            case 's2trans',  K_bot = kernel3d('h','sp',zk);
            case 'd2trans',  K_bot = kernel3d('h','dp',zk);
            case 'c2trans',  K_bot = kernel3d('h','cp',zk,coefs_ct);
        end
        failures = check_x2trans(failures, K, K_top, K_bot, src, targ, label);
    else
        failures = check_layer(failures, K, S, sig_s, targ_off, eps, tol_leval, label);

        if ~isempty(K_aux)   % prime kernel: FD check
            failures = check_fd_targ(failures, K, K_aux, src, targ, h_fd, tol_fd, label);
        end

        if ~isempty(K.fmm) && ~isempty(sig_pt)
            failures = check_fmm(failures, K, src, targ, sig_pt, eps, tol_fmm, label);
        end
    end
end

% =========================================================================
% Transmission kernel checks
% =========================================================================

fprintf('\n=== Transmission kernel checks ===\n');

zks        = [zk; zk1];
diff_coefs = [1.2; 0.8];          % a0, a1 for scalar diff checks
coe_22     = [0.7+0.1i, 0.3; ...  % 2x2 coefs for trans_sys
              0.5,       1.2i];
coe_222    = cat(3, coe_22, 0.9*coe_22);  % 2x2x2 coefs for trans_sys_diff

% --- scalar diff kernels: each against explicit scalar subtraction ---
for kd = { {'s_diff',  's',  {}};  {'d_diff',  'd',  {}};  ...
            {'sp_diff', 'sp', {}}; {'dp_diff', 'dp', {}} }
    kd = kd{1};
    Kdiff = kernel3d('h', kd{1}, zks, diff_coefs);
    K0    = diff_coefs(1) .* kernel3d('h', kd{2}, zk,  kd{3}{:});
    K1    = diff_coefs(2) .* kernel3d('h', kd{2}, zk1, kd{3}{:});
    Kref  = K0 - K1;
    M     = Kdiff.eval(src, targ);
    Mref  = Kref.eval(src, targ);
    failures = report(failures, ['helm ' kd{1}], 'eval vs explicit diff', ...
                      norm(M - Mref,'fro') / norm(Mref,'fro'), 1e-13);
end
% c_diff and cp_diff (need combined coefs)
c_coefs4 = [coefs_ct; diff_coefs];
for kd = { {'c_diff', 'c'}; {'cp_diff', 'cp'} }
    kd = kd{1};
    Kdiff = kernel3d('h', kd{1}, zks, c_coefs4);
    K0    = diff_coefs(1) .* kernel3d('h', kd{2}, zk,  coefs_ct);
    K1    = diff_coefs(2) .* kernel3d('h', kd{2}, zk1, coefs_ct);
    Kref  = K0 - K1;
    M     = Kdiff.eval(src, targ);
    Mref  = Kref.eval(src, targ);
    failures = report(failures, ['helm ' kd{1}], 'eval vs explicit diff', ...
                      norm(M - Mref,'fro') / norm(Mref,'fro'), 1e-13);
end

% --- trans_rep: check blocks match coefs(1)*D and coefs(2)*S ---
Ktr      = kernel3d('h', 'trans_rep', zk, [0.7; 1.3]);
M_tr     = Ktr.eval(src, targ);            % (nt x 2*ns)
M_D      = 0.7 * kernel3d('h','d',zk).eval(src,targ);
M_S      = 1.3 * kernel3d('h','s',zk).eval(src,targ);
err_D    = norm(M_tr(:,1:2:end) - M_D,'fro') / norm(M_D,'fro');
err_S    = norm(M_tr(:,2:2:end) - M_S,'fro') / norm(M_S,'fro');
failures = report(failures, 'helm trans_rep', 'eval blocks vs D, S', ...
                  max(err_D, err_S), 1e-13);

% --- trans_sys: 2x2 interleaved eval vs coefs.*[D,S;D',S'] blocks ---
Kts      = kernel3d('h', 'trans_sys', zk, coe_22);
M_ts     = Kts.eval(src, targ);            % (2*nt x 2*ns)
M_D      = kernel3d('h','d', zk).eval(src,targ);
M_S      = kernel3d('h','s', zk).eval(src,targ);
M_Dp     = kernel3d('h','dp',zk).eval(src,targ);
M_Sp     = kernel3d('h','sp',zk).eval(src,targ);
e11 = norm(M_ts(1:2:end,1:2:end) - coe_22(1,1)*M_D, 'fro') / norm(coe_22(1,1)*M_D,'fro');
e12 = norm(M_ts(1:2:end,2:2:end) - coe_22(1,2)*M_S, 'fro') / norm(coe_22(1,2)*M_S,'fro');
e21 = norm(M_ts(2:2:end,1:2:end) - coe_22(2,1)*M_Dp,'fro') / norm(coe_22(2,1)*M_Dp,'fro');
e22 = norm(M_ts(2:2:end,2:2:end) - coe_22(2,2)*M_Sp,'fro') / norm(coe_22(2,2)*M_Sp,'fro');
failures = report(failures, 'helm trans_sys', 'eval 2x2 blocks', max([e11 e12 e21 e22]), 1e-13);

% --- trans_sys: layer_eval ---
sigma_ts = reshape(sigma_2d, 2, S.npts);
failures = check_layer(failures, Kts, S, sigma_ts, targ_off, eps, tol_leval, 'helm trans_sys');

% --- trans_sys_diff: 2x2 eval vs explicit (c0.*K0 - c1.*K1) for each block ---
Ktsd     = kernel3d('h', 'trans_sys_diff', zks, coe_222);
M_tsd    = Ktsd.eval(src, targ);           % (2*nt x 2*ns)
c0 = coe_222(:,:,1);  c1 = coe_222(:,:,2);
M_D0 = kernel3d('h','d', zk ).eval(src,targ);  M_D1 = kernel3d('h','d', zk1).eval(src,targ);
M_S0 = kernel3d('h','s', zk ).eval(src,targ);  M_S1 = kernel3d('h','s', zk1).eval(src,targ);
M_Dp0= kernel3d('h','dp',zk ).eval(src,targ);  M_Dp1= kernel3d('h','dp',zk1).eval(src,targ);
M_Sp0= kernel3d('h','sp',zk ).eval(src,targ);  M_Sp1= kernel3d('h','sp',zk1).eval(src,targ);
f11 = norm(M_tsd(1:2:end,1:2:end) - (c0(1,1)*M_D0 -c1(1,1)*M_D1),'fro') / norm(c0(1,1)*M_D0-c1(1,1)*M_D1,'fro');
f12 = norm(M_tsd(1:2:end,2:2:end) - (c0(1,2)*M_S0 -c1(1,2)*M_S1),'fro') / norm(c0(1,2)*M_S0-c1(1,2)*M_S1,'fro');
f21 = norm(M_tsd(2:2:end,1:2:end) - (c0(2,1)*M_Dp0-c1(2,1)*M_Dp1),'fro') / norm(c0(2,1)*M_Dp0-c1(2,1)*M_Dp1,'fro');
f22 = norm(M_tsd(2:2:end,2:2:end) - (c0(2,2)*M_Sp0-c1(2,2)*M_Sp1),'fro') / norm(c0(2,2)*M_Sp0-c1(2,2)*M_Sp1,'fro');
failures = report(failures, 'helm trans_sys_diff', 'eval 2x2 blocks', max([f11 f12 f21 f22]), 1e-13);

% --- x2trans_diff kernels: rows match scalar diff halves ---
for kd = { {'s2trans_diff', 's_diff', 'sp_diff'}; ...
            {'d2trans_diff', 'd_diff', 'dp_diff'} }
    kd  = kd{1};
    K2t = kernel3d('h', kd{1}, zks, diff_coefs);
    Ktop = kernel3d('h', kd{2}, zks, diff_coefs);
    Kbot = kernel3d('h', kd{3}, zks, diff_coefs);
    failures = check_x2trans_diff(failures, K2t, Ktop, Kbot, src, targ, ['helm ' kd{1}]);
end
% c2trans_diff
Kc2d  = kernel3d('h', 'c2trans_diff',  zks, c_coefs4);
Kctop = kernel3d('h', 'c_diff',  zks, c_coefs4);
Kcbot = kernel3d('h', 'cp_diff', zks, c_coefs4);
failures = check_x2trans_diff(failures, Kc2d, Kctop, Kcbot, src, targ, 'helm c2trans_diff');

% =========================================================================
% Stokes S' kernel: eval == negative transpose of DLP (src<->targ swapped)
% =========================================================================
%
% Physical kernels: t(S)_{ij}(x,y) = T_{ijk}(x,y) n_k(x) = D_{ji}(y,x)
%
% Both kern.m DLP and S' use the double-negative convention
% (they return the negative of their respective physical kernels), so:
%   kern_SP(src,targ) = -t(S)(targ,src) = -D^T(src,targ) = kern_D^T(targ,src)
%
% In full-matrix notation (block 3x3 layout, (3*nt) x (3*ns)):
%   SP_mat  ==  -D_mat(targ,src)^T
%
fprintf('\n=== Stokes S'' vs negative-transpose of DLP ===\n');
K_sp = kernel3d('stok','sp');
K_d  = kernel3d('stok','d');
% Swap roles: D evaluated with sources=targ, targets=src
% (targ needs .n for S'; src needs .n for D — both already populated)
SP_mat   = K_sp.eval(src, targ);                 % (3*nt, 3*ns)
D_swap   = K_d.eval(targ, src);                  % (3*ns, 3*nt)  [src<->targ]
err_spdlp = norm(SP_mat + D_swap.', 'fro') / norm(SP_mat, 'fro');
failures = report(failures, 'stok sp', 'eval == -D(targ,src)^T', err_spdlp, 1e-12);

% =========================================================================
% Combined-layer checks: C.eval == coefs(1)*S.eval + coefs(2)*D.eval
% =========================================================================

fprintf('\n=== Combined-layer checks ===\n');

% Table: {Kc, Ks, Kd, coefs, sigma, label}
combined = { ...
  kernel3d('l','c',coefs_lap),     kernel3d('l','s'),        kernel3d('l','d'),         coefs_lap,  sigma_r,        'lap c';    ...
  kernel3d('h','c',zk,coefs_h),    kernel3d('h','s',zk),     kernel3d('h','d',zk),      coefs_h,    sigma,          'helm c';   ...
  kernel3d('stok','c',coefs_stok), kernel3d('stok','s'),      kernel3d('stok','d'),      coefs_stok, randn(3*ns,1),  'stok c';   ...
};

for k = 1:size(combined, 1)
    failures = check_combined(failures, combined{k,:}, src, targ);
end

% =========================================================================
% Arithmetic operators — single composite expression
% =========================================================================

fprintf('\n=== Arithmetic operators ===\n');

% Kexpr = ( (-helm_s).*(2+0.5i) + helm_d - (lap_s + helm_c)/3 ) / 1.5
% exercises: uminus, times, plus, minus, mrdivide.
K_hs  = kernel3d('h','s',zk);
K_hd  = kernel3d('h','d',zk);
K_ls  = kernel3d('l','s');
K_hc  = kernel3d('h','c',zk,coefs_h);
c1    = 2.0 + 0.5i;
c2    = 1.5;

Kexpr    = ((-K_hs) .* c1 + K_hd - (K_ls + K_hc)/3) / c2;
pot_expr = Kexpr.eval(src, targ) * sigma_r;
pot_ref  = ((-c1)*K_hs.eval(src,targ)*sigma_r + K_hd.eval(src,targ)*sigma_r ...
           - (K_ls.eval(src,targ)*sigma_r + K_hc.eval(src,targ)*sigma_r)/3) / c2;
err_arith = norm(pot_expr - pot_ref) / norm(pot_ref);
status    = pass_fail(err_arith, 1e-13);
if strcmp(status, 'FAIL')
    failures{end+1} = sprintf('arithmetic composite eval error %.2e', err_arith);
end
fprintf('  %-24s  err=%.2e  [%s]\n', 'composite expression', err_arith, status);

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
    assert(false, 'test_kernel3d: %d failure(s) — see output above.', nfail);
end

% =========================================================================
% Helpers
% =========================================================================

function s = pass_fail(err, tol)
    if err > tol, s = 'FAIL'; else, s = 'PASS'; end
end

function failures = check_layer(failures, K, S, sigma_s, targ_off, eps, tol, label)
%CHECK_LAYER  Compare layer_eval against smooth-quadrature reference.
%   sigma_s is (m, npts) for vector kernels or (npts,1) for scalar.
    if isempty(K.layer_eval)
        fprintf('  %-24s  layer_eval=none  (skip)\n', label);
        return
    end
    try
        pot_leval = K.layer_eval(S, sigma_s, targ_off, eps);
    catch ME
        failures{end+1} = sprintf('%s: layer_eval threw: %s', label, ME.message);
        fprintf('  %-24s  layer_eval ERROR: %s\n', label, ME.message);
        return
    end

    src_s.r  = S.r(:,:);  src_s.n  = S.n(:,:);
    src_s.du = S.du(:,:); src_s.dv = S.dv(:,:);
    wts = S.wts(:).';
    m   = K.opdims(2);
    if m > 1
        sig_flat = repmat(wts, m, 1) .* sigma_s(:,:);
        pot_ref  = K.eval(src_s, targ_off) * sig_flat(:);
    else
        pot_ref = K.eval(src_s, targ_off) * (wts(:) .* sigma_s(:));
    end

    err    = norm(pot_leval(:) - pot_ref(:)) / norm(pot_ref(:));
    status = pass_fail(err, tol);
    if strcmp(status, 'FAIL')
        failures{end+1} = sprintf('%s: layer_eval vs eval err=%.2e > tol=%.2e', label, err, tol);
    end
    fprintf('  %-24s  layer_eval vs eval  err=%.2e  [%s]\n', label, err, status);
end

function failures = check_fmm(failures, K, src, targ, sigma, eps, tol, label)
%CHECK_FMM  Compare K.fmm against K.eval for point-to-point evaluation.
%   sigma may be (m,ns) for vector kernels or (ns,1) for scalar; the
%   reference is K.eval(src,targ)*sigma(:).
    pot_fmm = K.fmm(eps, src, targ, sigma);
    pot_ref = K.eval(src, targ) * sigma(:);
    err     = norm(pot_fmm(:) - pot_ref(:)) / norm(pot_ref(:));
    status  = pass_fail(err, tol);
    if strcmp(status, 'FAIL')
        failures{end+1} = sprintf('%s: fmm vs eval err=%.2e > tol=%.2e', label, err, tol);
    end
    fprintf('  %-24s  fmm vs eval         err=%.2e  [%s]\n', label, err, status);
end

function failures = check_fd_targ(failures, K_prime, K_base, src, targ, h, tol, label)
%CHECK_FD_TARG  Verify prime kernel eval against centred FD in n_targ.
    mat_prime = K_prime.eval(src, targ);
    targ_p = targ;  targ_p.r = targ.r + h*targ.n;
    targ_m = targ;  targ_m.r = targ.r - h*targ.n;
    mat_fd = (K_base.eval(src, targ_p) - K_base.eval(src, targ_m)) / (2*h);
    err    = norm(mat_prime - mat_fd, 'fro') / norm(mat_fd, 'fro');
    status = pass_fail(err, tol);
    if strcmp(status, 'FAIL')
        failures{end+1} = sprintf('%s: prime vs FD(targ) err=%.2e > tol=%.2e', label, err, tol);
    end
    fprintf('  %-24s  prime vs FD(targ)   err=%.2e  [%s]\n', label, err, status);
end

function failures = check_combined(failures, Kc, Ks, Kd, coefs, sigma, label, src, targ)
%CHECK_COMBINED  Verify C.eval == coefs(1)*S.eval + coefs(2)*D.eval.
    pot_c   = Kc.eval(src, targ) * sigma;
    pot_ref = coefs(1)*Ks.eval(src,targ)*sigma + coefs(2)*Kd.eval(src,targ)*sigma;
    err     = norm(pot_c - pot_ref) / norm(pot_ref);
    status  = pass_fail(err, 1e-13);
    if strcmp(status, 'FAIL')
        failures{end+1} = sprintf('%s: C != a*S + b*D, err=%.2e', label, err);
    end
    fprintf('  %-24s  C = a*S + b*D       err=%.2e  [%s]\n', label, err, status);
end

function failures = check_x2trans(failures, K2t, K_top, K_bot, src, targ, label)
%CHECK_X2TRANS  Verify a [2x1] stacked transmission kernel.
%   Odd rows of K2t.eval -> K_top.eval,  even rows -> K_bot.eval.
    mat     = K2t.eval(src, targ);
    mat_top = K_top.eval(src, targ);
    mat_bot = K_bot.eval(src, targ);
    err_top = norm(mat(1:2:end,:) - mat_top, 'fro') / norm(mat_top, 'fro');
    err_bot = norm(mat(2:2:end,:) - mat_bot, 'fro') / norm(mat_bot, 'fro');
    err     = max(err_top, err_bot);
    status  = pass_fail(err, 1e-14);
    if strcmp(status, 'FAIL')
        failures{end+1} = sprintf('%s: block mismatch top=%.2e bot=%.2e', label, err_top, err_bot);
    end
    fprintf('  %-24s  eval blocks         err=%.2e  [%s]\n', label, err, status);
end

function failures = check_x2trans_diff(failures, K2t, K_top, K_bot, src, targ, label)
%CHECK_X2TRANS_DIFF  Like check_x2trans but K_top/K_bot are already the
%   expected difference kernels (kernel3d objects with eval handles).
    mat     = K2t.eval(src, targ);
    mat_top = K_top.eval(src, targ);
    mat_bot = K_bot.eval(src, targ);
    err_top = norm(mat(1:2:end,:) - mat_top, 'fro') / norm(mat_top, 'fro');
    err_bot = norm(mat(2:2:end,:) - mat_bot, 'fro') / norm(mat_bot, 'fro');
    err     = max(err_top, err_bot);
    status  = pass_fail(err, 1e-13);
    if strcmp(status, 'FAIL')
        failures{end+1} = sprintf('%s: block mismatch top=%.2e bot=%.2e', label, err_top, err_bot);
    end
    fprintf('  %-24s  eval blocks         err=%.2e  [%s]\n', label, err, status);
end

function failures = report(failures, label, desc, err, tol)
%REPORT  Print one result line and accumulate failures.
    status = pass_fail(err, tol);
    if strcmp(status, 'FAIL')
        failures{end+1} = sprintf('%s %s err=%.2e', label, desc, err);
    end
    fprintf('  %-24s  %-22s  err=%.2e  [%s]\n', label, desc, err, status);
end
