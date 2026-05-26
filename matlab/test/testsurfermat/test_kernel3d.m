%TEST_KERNEL3D  Verify kernel3d eval and layer_eval agree for all
% implemented kernels, and that arithmetic operators (+, -, .*, ./) work.
%
% For prime kernels, eval is checked against a finite-difference of the
% unprimed kernel in the target normal direction.
% Combined kernels are checked against coefs(1)*S + coefs(2)*D.
% FMM is checked against eval for Laplace, Helmholtz, and Stokes.
% The arithmetic section uses a single composite kernel expression.
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
coefs_lap = [0.7; 1.3];
coefs_h   = [1i*zk; 1.0];
coefs_ct  = [1i*zk; 1.0];
coefs_stok = [0.7; 1.3];
zk_em     = 2.5 + 0.1i;
alpha_em  = 0.5;

sigma_stok_pt = randn(3, ns);   % Stokes point density for FMM check
sigma_em3_pt  = randn(3, ns) + 1i*randn(3,ns); % em nrccie-bc: (3, ns) [j_ru, j_rv, rho] (orthonormal frame)
sigma_em4_pt  = randn(4, ns) + 1i*randn(4, ns); % em nrccie-eval: (4, npts) [Jx, Jy, Jz, rho]

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
% --- Transmission x2trans (check_x2trans; surface sigma unused but must be present) ---
  kernel3d('h','s2trans',zk),     'helm s2trans',sigma_s,    sigma_s,              kernel3d('h','s',zk);          ...
  kernel3d('h','d2trans',zk),     'helm d2trans',sigma_s,    sigma_s,              kernel3d('h','d',zk);          ...
  kernel3d('h','c2trans',zk,coefs_ct),'helm c2trans',sigma_s,sigma_s,             kernel3d('h','c',zk,coefs_ct);  ...
% --- Stokes ---
  kernel3d('stok','s'),           'stok s',      sigma_stok, sigma_stok_pt,   [];                           ...
  kernel3d('stok','d'),           'stok d',      sigma_stok, sigma_stok_pt,   [];                           ...
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
    K_aux  = kernels{k,5};  % base kernel (FD check) or top kernel (x2trans)

    % x2trans kernels store their bottom-half kernel in K_aux (col 5).
    % The top-half kernel is stored one row up, in col 5 of that row — but
    % it is simpler to just pass both halves directly in the table (K_aux
    % holds the top kernel; the bottom is derived from the type string).
    is_x2trans = ismember(K.type, {'s2trans','d2trans','c2trans'});

    if is_x2trans
        K_top = K_aux;   % top-half kernel stored in col 5
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

% trans_rep: [1x2] interleaved kernel  u = D[rho] + (ik)*S[sigma].
% Density must be interleaved: [rho(1); sigma(1); rho(2); sigma(2); ...]
Ktr = kernel3d('h', 'trans_rep', zk);
rho_rep   = sigma_r(1:ns);
sig_rep   = 2*sigma_r(1:ns);
sigma_int = reshape([rho_rep.'; sig_rep.'], 2*ns, 1);   % interleaved layout
pot_tr    = Ktr.eval(src, targ) * sigma_int;
pot_ref   = kernel3d('h','d',zk).eval(src,targ)*rho_rep ...
          + (1i*zk)*kernel3d('h','s',zk).eval(src,targ)*sig_rep;
err_tr    = norm(pot_tr - pot_ref) / norm(pot_ref);
status    = pass_fail(err_tr, 1e-13);
if strcmp(status, 'FAIL')
    failures{end+1} = sprintf('helm trans_rep: eval vs D+(ik)S  err=%.2e', err_tr);
end
fprintf('  %-24s  eval vs D+(ik)S  err=%.2e  [%s]\n', 'helm trans_rep', err_tr, status);

% =========================================================================
% Maxwell eval checks
% =========================================================================
%
% nrccie-bc  eval check: K.eval(src,targ)*sigma == manual assembly from
%   sub-kernels NxCurlS, NxS, NxDelS, NdotS, DeldotS, S, S'.
%
% nrccie-eval eval check: K.eval(src,targ)*sigma == ik*S[J] - gradS[rho]
%   for E, and curlS[J] for H, computed directly from helm3d.green.

fprintf('\n=== Maxwell eval checks ===\n');

% --- shared setup ---
% src already has .r, .n, .du, .dv set at the top of the script

% Density for nrccie-bc: [j_u, j_v, rho]  each ns x 1, complex
ju  = randn(ns,1) + 1i*randn(ns,1);
jv  = randn(ns,1) + 1i*randn(ns,1);
rho = randn(ns,1) + 1i*randn(ns,1);
% Interleaved layout to match reshape((3,ns,3,ns),...) column-major convention
% used in em3d.kern (same as Stokes): row/col index = icomp + m*(ipoint-1).
sigma_bc = reshape([ju, jv, rho].', 3*ns, 1);   % (3*ns, 1), interleaved: [ju(1);jv(1);rho(1);ju(2);...]

% ---- nrccie-bc eval ----
K_bc = kernel3d('em', 'nrccie-bc', zk_em, alpha_em);

pot_bc = K_bc.eval(src, targ) * sigma_bc;   % (3*nt, 1)

% Reference: same computation as em3d.kern 'nrccie-bc', inlined.
% Follows get_nrccie_inteq_comps_from_potgrad in em_nrccie_pec.f90.
[Gfun, gradG] = helm3d.green(zk_em, src.r, targ.r);
gx = gradG(:,:,1);  gy = gradG(:,:,2);  gz = gradG(:,:,3);

nx  = repmat(targ.n(1,:).',  1, ns);  ny  = repmat(targ.n(2,:).', 1, ns);  nz  = repmat(targ.n(3,:).', 1, ns);
dxut = repmat(targ.du(1,:).', 1, ns); dyut = repmat(targ.du(2,:).', 1, ns); dzut = repmat(targ.du(3,:).', 1, ns);
dxvt = repmat(targ.dv(1,:).', 1, ns); dyvt = repmat(targ.dv(2,:).', 1, ns); dzvt = repmat(targ.dv(3,:).', 1, ns);
dxus = repmat(src.du(1,:), nt, 1);    dyus = repmat(src.du(2,:), nt, 1);    dzus = repmat(src.du(3,:), nt, 1);
dxvs = repmat(src.dv(1,:), nt, 1);    dyvs = repmat(src.dv(2,:), nt, 1);    dzvs = repmat(src.dv(3,:), nt, 1);
nxs  = repmat(src.n(1,:),  nt, 1);    nys  = repmat(src.n(2,:),  nt, 1);    nzs  = repmat(src.n(3,:),  nt, 1);

% j_u column
curl_x_u = gy.*dzus - gz.*dyus;  curl_y_u = gz.*dxus - gx.*dzus;  curl_z_u = gx.*dyus - gy.*dxus;
zv_x_u = ny.*curl_z_u - nz.*curl_y_u;  zv_y_u = nz.*curl_x_u - nx.*curl_z_u;  zv_z_u = nx.*curl_y_u - ny.*curl_x_u;
ndE_u = 1i*zk_em*Gfun.*(nx.*dxus + ny.*dyus + nz.*dzus);
v3x_u = alpha_em*(ndE_u.*nx - 1i*zk_em*Gfun.*dxus) - zv_x_u;
v3y_u = alpha_em*(ndE_u.*ny - 1i*zk_em*Gfun.*dyus) - zv_y_u;
v3z_u = alpha_em*(ndE_u.*nz - 1i*zk_em*Gfun.*dzus) - zv_z_u;
K_uu_r = dxut.*v3x_u + dyut.*v3y_u + dzut.*v3z_u;
K_vu_r = dxvt.*v3x_u + dyvt.*v3y_u + dzvt.*v3z_u;
K_ru_r = -1i*zk_em*Gfun.*(dxus.*nx+dyus.*ny+dzus.*nz) + alpha_em*(gx.*dxus+gy.*dyus+gz.*dzus);

% j_v column
curl_x_v = gy.*dzvs - gz.*dyvs;  curl_y_v = gz.*dxvs - gx.*dzvs;  curl_z_v = gx.*dyvs - gy.*dxvs;
zv_x_v = ny.*curl_z_v - nz.*curl_y_v;  zv_y_v = nz.*curl_x_v - nx.*curl_z_v;  zv_z_v = nx.*curl_y_v - ny.*curl_x_v;
ndE_v = 1i*zk_em*Gfun.*(nx.*dxvs + ny.*dyvs + nz.*dzvs);
v3x_v = alpha_em*(ndE_v.*nx - 1i*zk_em*Gfun.*dxvs) - zv_x_v;
v3y_v = alpha_em*(ndE_v.*ny - 1i*zk_em*Gfun.*dyvs) - zv_y_v;
v3z_v = alpha_em*(ndE_v.*nz - 1i*zk_em*Gfun.*dzvs) - zv_z_v;
K_uv_r = dxut.*v3x_v + dyut.*v3y_v + dzut.*v3z_v;
K_vv_r = dxvt.*v3x_v + dyvt.*v3y_v + dzvt.*v3z_v;
K_rv_r = -1i*zk_em*Gfun.*(dxvs.*nx+dyvs.*ny+dzvs.*nz) + alpha_em*(gx.*dxvs+gy.*dyvs+gz.*dzvs);

% rho column
ndE_r = -(gx.*nx + gy.*ny + gz.*nz);
v3x_r = alpha_em*(ndE_r.*nx + gx);  v3y_r = alpha_em*(ndE_r.*ny + gy);  v3z_r = alpha_em*(ndE_r.*nz + gz);
K_ur_r = dxut.*v3x_r + dyut.*v3y_r + dzut.*v3z_r;
K_vr_r = dxvt.*v3x_r + dyvt.*v3y_r + dzvt.*v3z_r;
Sp_ref = -(nxs.*gx + nys.*gy + nzs.*gz);
K_rr_r = Sp_ref - alpha_em*1i*zk_em*Gfun;

% Compute each output component per target (nt x 1 vectors)
out_u   = K_uu_r*ju + K_uv_r*jv + K_ur_r*rho;
out_v   = K_vu_r*ju + K_vv_r*jv + K_vr_r*rho;
out_rho = K_ru_r*ju + K_rv_r*jv + K_rr_r*rho;
% Interleave: [u(1); v(1); rho(1); u(2); v(2); rho(2); ...]
pot_bc_ref = reshape([out_u, out_v, out_rho].', 3*nt, 1);

err_bc = norm(pot_bc - pot_bc_ref) / norm(pot_bc_ref);
status = pass_fail(err_bc, 1e-10);
if strcmp(status, 'FAIL')
    failures{end+1} = sprintf('em nrccie-bc: eval vs sub-kernels err=%.2e', err_bc);
end
fprintf('  %-24s  eval vs sub-kernels  err=%.2e  [%s]\n', 'em nrccie-bc', err_bc, status);

% ---- nrccie-eval eval ----
% Density: [Jx, Jy, Jz, rho_c]  each ns columns, complex
Jx  = randn(ns,1) + 1i*randn(ns,1);
Jy  = randn(ns,1) + 1i*randn(ns,1);
Jz  = randn(ns,1) + 1i*randn(ns,1);
rhoc = randn(ns,1) + 1i*randn(ns,1);
sigma_ev = reshape([Jx, Jy, Jz, rhoc].', 4*ns, 1);   % (4*ns,1), interleaved: [Jx(1);Jy(1);Jz(1);rho(1);...]

K_ev = kernel3d('em', 'nrccie-eval', zk_em);

pot_ev = K_ev.eval(src, targ) * sigma_ev;  % (6*nt, 1)

% Reference: E = ik*S_k[J] - gradS_k[rho],  H = curl_t(S_k[J]) = gradG x J
% Gfun/gradG/gx/gy/gz already computed above for zk_em, same src/targ.
ikG = 1i*zk_em * Gfun;   % (nt, ns)

Ex_ref = ikG*Jx - gx*rhoc;
Ey_ref = ikG*Jy - gy*rhoc;
Ez_ref = ikG*Jz - gz*rhoc;
% H = gradG x J  (curl of scalar G times vector J = gradG x J)
%   H_x = gy*Jz - gz*Jy
%   H_y = gz*Jx - gx*Jz
%   H_z = gx*Jy - gy*Jx
Hx_ref = gy*Jz - gz*Jy;
Hy_ref = gz*Jx - gx*Jz;
Hz_ref = gx*Jy - gy*Jx;

% Row layout is interleaved: icomp + 6*(itarg-1), matching reshape((6,nt,4,ns),...).
% Interleave to: [Ex(1);Ey(1);Ez(1);Hx(1);Hy(1);Hz(1); Ex(2);...]
pot_ev_ref = reshape([Ex_ref, Ey_ref, Ez_ref, Hx_ref, Hy_ref, Hz_ref].', 6*nt, 1);

err_ev = norm(pot_ev - pot_ev_ref) / norm(pot_ev_ref);
status = pass_fail(err_ev, 1e-10);
if strcmp(status, 'FAIL')
    failures{end+1} = sprintf('em nrccie-eval: eval vs direct err=%.2e', err_ev);
end
fprintf('  %-24s  eval vs direct       err=%.2e  [%s]\n', 'em nrccie-eval', err_ev, status);

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
