%TEST_KERNEL3D  Verify kernel3d eval, fmm, and layer_eval agree for all
% implemented kernels, and that arithmetic operators (+, -, .*, ./) work.
% For prime kernels, eval is also checked against a finite difference of
% the unprimed kernel.
%
% Failures are accumulated and reported at the end.  The script asserts
% only once, at the very end, so every kernel is always exercised.

run ../../startup.m

rng(42);

zk   = 3.0 + 0.2i;
eps  = 1e-6;
ns   = 500;
nt   = 300;

% Sources inside unit ball, targets in annulus [2,3] — no singularities.
src.r = randn(3,ns);  src.r = src.r ./ vecnorm(src.r) .* rand(1,ns);
src.n = randn(3,ns);  src.n = src.n ./ vecnorm(src.n);

targ.r = randn(3,nt);  targ.r = targ.r ./ vecnorm(targ.r) .* (2 + rand(1,nt));
targ.n = randn(3,nt);  targ.n = targ.n ./ vecnorm(targ.n);

sigma_r = randn(ns,1);                     % real, for Laplace
sigma   = randn(ns,1) + 1i*randn(ns,1);   % complex, for Helmholtz

% Surface for layer_eval tests. Use high order so that the smooth quadrature
% sum K.eval*(wts.*sigma) is an accurate reference for layer_eval.
S = geometries.sphere(1, 4, [0;0;0], 6, 1);
sigma_sr = randn(S.npts,1);                        % real, for Laplace
sigma_s  = randn(S.npts,1) + 1i*randn(S.npts,1);  % complex, for Helmholtz

% Targets at radius 5 -- far enough that smooth quadrature converges well.
targ_off.r = randn(3,40);
targ_off.r = targ_off.r ./ vecnorm(targ_off.r) * 5.0;
targ_off.n = randn(3,40);  targ_off.n = targ_off.n ./ vecnorm(targ_off.n);

tol_fmm   = 1e-4;   % eval vs fmm
tol_leval = 1e-4;   % eval vs layer_eval  (surface quadrature limited)
tol_fd    = 1e-3;   % prime eval vs finite difference

h_fd = 1e-5;        % finite-difference step

% Accumulate failure messages — assert once at the very end.
failures = {};

% -------------------------------------------------------------------------
% Helper: check eval vs fmm
% -------------------------------------------------------------------------
function [failures] = check_fmm(failures, K, src, targ, sigma, eps, tol, label)
    mat      = K.eval(src, targ);
    pot_eval = mat * sigma;
    if isempty(K.fmm)
        fprintf('  %-32s  fmm=none  (skip)\n', label);
        return
    end
    pot_fmm  = K.fmm(eps, src, targ, sigma);
    err = norm(pot_fmm - pot_eval) / norm(pot_eval);
    status = 'PASS';
    if err > tol
        status = 'FAIL';
        failures{end+1} = sprintf('%s: eval vs fmm err=%.2e > tol=%.2e', label, err, tol);
    end
    fprintf('  %-32s  eval vs fmm       err=%.2e  [%s]\n', label, err, status);
end

% -------------------------------------------------------------------------
% Helper: check eval vs layer_eval (off-surface targets)
% -------------------------------------------------------------------------
function [failures] = check_layer(failures, K, S, sigma_s, targ_off, eps, tol, label)
    if isempty(K.layer_eval)
        fprintf('  %-32s  layer_eval=none  (skip)\n', label);
        return
    end
    pot_leval = K.layer_eval(S, sigma_s, targ_off, eps);
    src_s.r = S.r(:,:);  src_s.n = S.n(:,:);
    wts = S.wts(:).';
    if ~isempty(K.fmm)
        pot_ref = K.fmm(eps, src_s, targ_off, wts(:) .* sigma_s(:));
        refname = 'fmm';
    else
        pot_ref = K.eval(src_s, targ_off) * (wts(:) .* sigma_s(:));
        refname = 'smooth quad';
    end
    err = norm(pot_leval(:) - pot_ref(:)) / norm(pot_ref(:));
    status = 'PASS';
    if err > tol
        status = 'FAIL';
        failures{end+1} = sprintf('%s: layer_eval vs %s err=%.2e > tol=%.2e', ...
                                  label, refname, err, tol);
    end
    fprintf('  %-32s  layer_eval vs %-11s err=%.2e  [%s]\n', label, refname, err, status);
end

% -------------------------------------------------------------------------
% Helper: check prime kernel eval vs finite difference of unprimed in n_targ
% -------------------------------------------------------------------------
function [failures] = check_fd_targ(failures, K_prime, K_base, src, targ, h, tol, label)
    mat_prime = K_prime.eval(src, targ);
    targ_p = targ;  targ_p.r = targ.r + h*targ.n;
    targ_m = targ;  targ_m.r = targ.r - h*targ.n;
    mat_fd = (K_base.eval(src, targ_p) - K_base.eval(src, targ_m)) / (2*h);
    err = norm(mat_prime - mat_fd, 'fro') / norm(mat_fd, 'fro');
    status = 'PASS';
    if err > tol
        status = 'FAIL';
        failures{end+1} = sprintf('%s: prime vs FD(targ) err=%.2e > tol=%.2e', label, err, tol);
    end
    fprintf('  %-32s  prime vs FD(targ)  err=%.2e  [%s]\n', label, err, status);
end

% Helper: check prime kernel eval vs finite difference of unprimed in n_src
function [failures] = check_fd_src(failures, K_prime, K_base, src, targ, h, tol, label)
    mat_prime = K_prime.eval(src, targ);
    src_p = src;  src_p.r = src.r + h*src.n;
    src_m = src;  src_m.r = src.r - h*src.n;
    mat_fd = (K_base.eval(src_p, targ) - K_base.eval(src_m, targ)) / (2*h);
    err = norm(mat_prime - mat_fd, 'fro') / norm(mat_fd, 'fro');
    status = 'PASS';
    if err > tol
        status = 'FAIL';
        failures{end+1} = sprintf('%s: prime vs FD(src) err=%.2e > tol=%.2e', label, err, tol);
    end
    fprintf('  %-32s  prime vs FD(src)   err=%.2e  [%s]\n', label, err, status);
end

% =========================================================================
fprintf('=== Laplace kernels ===\n');

coefs_lap = [0.7; 1.3];

kernels_lap = { ...
    kernel3d('l','s'),           'lap s',    false; ...
    kernel3d('l','d'),           'lap d',    false; ...
    kernel3d('l','sp'),          'lap sp',   true;  ...   % sp = d/dn_targ S
    kernel3d('l','dp'),          'lap dp',   true;  ...   % dp = d/dn_targ D
    kernel3d('l','c',coefs_lap), 'lap c',    false; ...
    kernel3d('l','cp',coefs_lap),'lap cp',   true;  ...
};

for k = 1:size(kernels_lap,1)
    K       = kernels_lap{k,1};
    label   = kernels_lap{k,2};
    isprime = kernels_lap{k,3};

    failures = check_fmm(failures, K, src, targ, sigma_r, eps, tol_fmm, label);
    failures = check_layer(failures, K, S, sigma_sr, targ_off, eps, tol_leval, label);

    if isprime
        switch label
            case 'lap sp'   % d/dn_targ of S
                failures = check_fd_targ(failures, K, kernel3d('l','s'), src, targ, h_fd, tol_fd, label);
            case 'lap dp'   % d/dn_targ of D
                failures = check_fd_targ(failures, K, kernel3d('l','d'), src, targ, h_fd, tol_fd, label);
            case 'lap cp'   % d/dn_targ of C
                failures = check_fd_targ(failures, K, kernel3d('l','c',coefs_lap), src, targ, h_fd, tol_fd, label);
        end
    end
end

% =========================================================================
fprintf('\n=== Helmholtz kernels ===\n');

coefs_h = [1i*zk; 1.0];

kernels_helm = { ...
    kernel3d('h','s',zk),          'helm s',       false; ...
    kernel3d('h','d',zk),          'helm d',       false; ...
    kernel3d('h','sp',zk),         'helm sp',      true;  ...
    kernel3d('h','dp',zk),         'helm dprime',  true;  ...
    kernel3d('h','c',zk,coefs_h),  'helm c',       false; ...
    kernel3d('h','cp',zk,coefs_h), 'helm cp',      true;  ...
};

for k = 1:size(kernels_helm,1)
    K       = kernels_helm{k,1};
    label   = kernels_helm{k,2};
    isprime = kernels_helm{k,3};

    failures = check_fmm(failures, K, src, targ, sigma, eps, tol_fmm, label);
    failures = check_layer(failures, K, S, sigma_s, targ_off, eps, tol_leval, label);

    if isprime
        switch label
            case 'helm sp'      % d/dn_targ of S
                failures = check_fd_targ(failures, K, kernel3d('h','s',zk), src, targ, h_fd, tol_fd, label);
            case 'helm dprime'  % d/dn_targ of D
                failures = check_fd_targ(failures, K, kernel3d('h','d',zk), src, targ, h_fd, tol_fd, label);
            case 'helm cp'      % d/dn_targ of C
                failures = check_fd_targ(failures, K, kernel3d('h','c',zk,coefs_h), src, targ, h_fd, tol_fd, label);
        end
    end
end

% =========================================================================
fprintf('\n=== Transmission kernels ===\n');
% Tests for trans_rep, s2trans, d2trans, c2trans.
% All are built via interleave so we check eval dimensions and consistency
% with the individual scalar sub-kernels.

coefs_ct = [1i*zk; 1.0];   % combined-layer coefficients for c2trans

% --- s2trans: [S; S'] applied to a single-component density ---
Ks2t = kernel3d('h', 's2trans', zk);
mat_s2t = Ks2t.eval(src, targ);   % (2*nt, ns)
mat_Ks  = kernel3d('h','s', zk).eval(src, targ);
mat_Ksp = kernel3d('h','sp',zk).eval(src, targ);
err_s2t_row1 = norm(mat_s2t(1:2:end,:) - mat_Ks,  'fro') / norm(mat_Ks,  'fro');
err_s2t_row2 = norm(mat_s2t(2:2:end,:) - mat_Ksp, 'fro') / norm(mat_Ksp, 'fro');
status = 'PASS';
if err_s2t_row1 > 1e-14 || err_s2t_row2 > 1e-14
    status = 'FAIL';
    failures{end+1} = sprintf('s2trans eval block mismatch: row1=%.2e row2=%.2e', ...
                               err_s2t_row1, err_s2t_row2);
end
fprintf('  %-32s  eval shape/blocks     [%s]\n', 'helm s2trans', status);

% --- d2trans: [D; D'] ---
Kd2t = kernel3d('h', 'd2trans', zk);
mat_d2t = Kd2t.eval(src, targ);
mat_Kd  = kernel3d('h','d', zk).eval(src, targ);
mat_Kdp = kernel3d('h','dp',zk).eval(src, targ);
err_d2t_row1 = norm(mat_d2t(1:2:end,:) - mat_Kd,  'fro') / norm(mat_Kd,  'fro');
err_d2t_row2 = norm(mat_d2t(2:2:end,:) - mat_Kdp, 'fro') / norm(mat_Kdp, 'fro');
status = 'PASS';
if err_d2t_row1 > 1e-14 || err_d2t_row2 > 1e-14
    status = 'FAIL';
    failures{end+1} = sprintf('d2trans eval block mismatch: row1=%.2e row2=%.2e', ...
                               err_d2t_row1, err_d2t_row2);
end
fprintf('  %-32s  eval shape/blocks     [%s]\n', 'helm d2trans', status);

% --- c2trans: [C; C'] ---
Kc2t   = kernel3d('h', 'c2trans', zk, coefs_ct);
mat_c2t = Kc2t.eval(src, targ);
mat_Kc  = kernel3d('h','c', zk, coefs_ct).eval(src, targ);
mat_Kcp = kernel3d('h','cp',zk, coefs_ct).eval(src, targ);
err_c2t_row1 = norm(mat_c2t(1:2:end,:) - mat_Kc,  'fro') / norm(mat_Kc,  'fro');
err_c2t_row2 = norm(mat_c2t(2:2:end,:) - mat_Kcp, 'fro') / norm(mat_Kcp, 'fro');
status = 'PASS';
if err_c2t_row1 > 1e-14 || err_c2t_row2 > 1e-14
    status = 'FAIL';
    failures{end+1} = sprintf('c2trans eval block mismatch: row1=%.2e row2=%.2e', ...
                               err_c2t_row1, err_c2t_row2);
end
fprintf('  %-32s  eval shape/blocks     [%s]\n', 'helm c2trans', status);

% --- trans_rep: [D, (i*zk)*S] -- field evaluation kernel ---
Ktr = kernel3d('h', 'trans_rep', zk);
assert(isequal(Ktr.opdims, [1, 2]), 'trans_rep opdims should be [1 2]');
% Apply to a 2-component density [rho; sigma]: should equal D[rho] + (i*zk)*S[sigma]
sigma_2 = [sigma_r(1:ns); sigma_r(1:ns)];   % 2*ns vector (rho then sigma)
pot_tr  = Ktr.eval(src, targ) * sigma_2;
pot_ref = kernel3d('h','d',zk).eval(src,targ)*sigma_r(1:ns) ...
        + (1i*zk)*kernel3d('h','s',zk).eval(src,targ)*sigma_r(1:ns);
err_tr  = norm(pot_tr - pot_ref) / norm(pot_ref);
status = 'PASS';
if err_tr > 1e-13
    status = 'FAIL';
    failures{end+1} = sprintf('trans_rep eval err=%.2e', err_tr);
end
fprintf('  %-32s  eval vs D+(ik)S       err=%.2e  [%s]\n', 'helm trans_rep', err_tr, status);

% =========================================================================
fprintf('\n=== Arithmetic operators ===\n');

Ks = kernel3d('h','s',zk);
Kd = kernel3d('h','d',zk);

% --- uminus: -K ---
Kneg     = -Ks;
pot_pos  = Ks.eval(src,targ)*sigma;
pot_neg  = Kneg.eval(src,targ)*sigma;
err = norm(pot_neg + pot_pos) / norm(pot_pos);
status = 'PASS';
if err > 1e-14
    status = 'FAIL';
    failures{end+1} = sprintf('uminus eval error %.2e', err);
end
if ~isempty(Kneg.fmm)
    err_fmm = norm(Kneg.fmm(eps,src,targ,sigma) + Ks.fmm(eps,src,targ,sigma)) / norm(pot_pos);
    if err_fmm > tol_fmm
        status = 'FAIL';
        failures{end+1} = sprintf('uminus fmm error %.2e', err_fmm);
    end
end
fprintf('  %-32s  [%s]\n', 'uminus', status);

% --- times: c.*K ---
c      = 3.0 + 0.5i;
Kscale = c .* Ks;
err = norm(Kscale.eval(src,targ)*sigma - c*pot_pos) / norm(pot_pos);
status = 'PASS';
if err > 1e-14
    status = 'FAIL';
    failures{end+1} = sprintf('times eval error %.2e', err);
end
if ~isempty(Kscale.fmm)
    err_fmm = norm(Kscale.fmm(eps,src,targ,sigma) - c*Ks.fmm(eps,src,targ,sigma)) / norm(pot_pos);
    if err_fmm > tol_fmm
        status = 'FAIL';
        failures{end+1} = sprintf('times fmm error %.2e', err_fmm);
    end
end
fprintf('  %-32s  [%s]\n', 'times (c.*K)', status);

% --- mrdivide: K/c ---
Kdiv = Ks / c;
err  = norm(Kdiv.eval(src,targ)*sigma - pot_pos/c) / norm(pot_pos);
status = 'PASS';
if err > 1e-14
    status = 'FAIL';
    failures{end+1} = sprintf('mrdivide eval error %.2e', err);
end
fprintf('  %-32s  [%s]\n', 'mrdivide (K/c)', status);

% --- plus: Ks + Kd ---
Ksum  = Ks + Kd;
pot_s = Ks.eval(src,targ)*sigma;
pot_d = Kd.eval(src,targ)*sigma;
err   = norm(Ksum.eval(src,targ)*sigma - (pot_s + pot_d)) / norm(pot_s + pot_d);
status = 'PASS';
if err > 1e-14
    status = 'FAIL';
    failures{end+1} = sprintf('plus eval error %.2e', err);
end
if ~isempty(Ksum.fmm)
    pot_sum_fmm = Ksum.fmm(eps,src,targ,sigma);
    pot_ref_fmm = Ks.fmm(eps,src,targ,sigma) + Kd.fmm(eps,src,targ,sigma);
    err_fmm = norm(pot_sum_fmm - pot_ref_fmm) / norm(pot_ref_fmm);
    if err_fmm > tol_fmm
        status = 'FAIL';
        failures{end+1} = sprintf('plus fmm error %.2e', err_fmm);
    end
end
fprintf('  %-32s  [%s]\n', 'plus (Ks + Kd)', status);

% --- minus: Ks - Kd ---
Kdiff = Ks - Kd;
err   = norm(Kdiff.eval(src,targ)*sigma - (pot_s - pot_d)) / norm(pot_s - pot_d);
status = 'PASS';
if err > 1e-14
    status = 'FAIL';
    failures{end+1} = sprintf('minus eval error %.2e', err);
end
fprintf('  %-32s  [%s]\n', 'minus (Ks - Kd)', status);

% --- combined expression: (lap_s - 2*helm_c) / 3 ---
Kl    = kernel3d('l','s');
Kh    = kernel3d('h','c',zk,coefs_h);
Kexpr = (Kl - 2.*Kh) / 3;
pot_expr = Kexpr.eval(src,targ)*sigma_r;
pot_ref  = (Kl.eval(src,targ)*sigma_r - 2*Kh.eval(src,targ)*sigma_r) / 3;
err = norm(pot_expr - pot_ref) / norm(pot_ref);
status = 'PASS';
if err > 1e-14
    status = 'FAIL';
    failures{end+1} = sprintf('combined expression eval error %.2e', err);
end
fprintf('  %-32s  [%s]\n', '(lap_s - 2*helm_c)/3', status);

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
