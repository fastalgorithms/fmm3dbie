%TEST_JUMP_RELATIONS   Verify jump relations and on-surface limits.
%
%  For each kernel we pick a single point near the patch centre, displace it
%  by ±h along the outward normal, and check:
%
%  (1) Jump check:  u^+ - u^-  matches the known discontinuity.
%
%  (2) One-sided limit (usematlab = 0 and 1):
%        u^+ = surfermatapply(at that point) + jump/2
%      where surfermatapply gives the principal-value on-surface value.
%
%  Accuracy of the smooth quadrature is tested in test_kernel3d.m; here
%  we only need the near-field corrections to produce the correct O(h)-
%  accurate behaviour.  A single patch keeps the test fast.
%
%  Jump relations tested:
%    lap  d   :  [[u]] =  +sigma
%    lap  sp  :  [[u]] =  -sigma
%    lap  c   :  [[u]] =  +coefs(2)*sigma
%    lap  cp  :  [[u]] =  -coefs(1)*sigma
%    helm d   :  [[u]] =  +sigma
%    helm sp  :  [[u]] =  -sigma
%    helm c   :  [[u]] =  +coefs(2)*sigma
%    helm cp  :  [[u]] =  -coefs(1)*sigma
%    stok d   :  [[u]] =  +sigma  (3-vector)
%    stok sp  :  [[u]] =  -sigma  (3-vector)
%    stok c   :  [[u]] =  +coefs(2)*sigma
%    em nrccie-bc  :  [[u]] = -sigma  (3 frame components)
%    em nrccie-eval:
%      n x [[E]] = 0,  n x [[H]] = J,  n . [[H]] = 0,  n . [[E]] = rho

run ../../startup.m

rng(42);

% -----------------------------------------------------------------------
% Geometry: one patch of a sphere, order 6
% -----------------------------------------------------------------------
S_full = geometries.sphere(1, 4, [0;0;0], 6, 1);
S      = slicesurfer(S_full, 1);

eps = 1e-5;
tol = 0.05;    % allow O(h) error
h   = 1e-4;    % displacement for exterior / interior targets

% -----------------------------------------------------------------------
% Smooth densities (on the full patch)
% -----------------------------------------------------------------------
sigma_r = cos(S.r(1,:).' + 0.5*S.r(2,:).') .* exp(-S.r(3,:).'.^2);
sigma_r = sigma_r / norm(sigma_r);

zk      = 1.1 + 0.3i;
sigma_c = (S.r(3,:).^2 + 1i*S.r(1,:)).';
sigma_c = sigma_c / norm(sigma_c);

zk_em = 1.3;

coefs_stok = [0.7; 1.3];
sigma_stok = [cos(S.r(2,:)); sin(S.r(3,:)); cos(S.r(1,:) + S.r(2,:))];
sigma_stok = sigma_stok / norm(sigma_stok(:));

% EM nrccie-bc density: [j_ru, j_rv, rho] in orthonormal tangent frame
sigma_em3 = [cos(S.r(1,:) + S.r(2,:)); ...
             sin(S.r(2,:) + S.r(3,:)); ...
             cos(S.r(3,:))] ...
          + 1i*[sin(S.r(1,:)); cos(S.r(2,:)); sin(S.r(3,:))];
sigma_em3 = sigma_em3 / norm(sigma_em3(:));

% EM nrccie-eval density: [Jx, Jy, Jz, rho], J projected tangential
J_cart = [cos(S.r(1,:) + S.r(2,:)); ...
          sin(S.r(2,:) + S.r(3,:)); ...
          cos(S.r(3,:))] ...
       + 1i*[sin(S.r(1,:) + S.r(3,:)); ...
             cos(S.r(1,:));             ...
             sin(S.r(2,:))];
J_cart = J_cart - S.n .* sum(S.n .* J_cart, 1);
rho_em = cos(S.r(3,:)) + 1i*cos(S.r(1,:) + S.r(2,:));
sigma_em4 = [J_cart; rho_em];
sigma_em4 = sigma_em4 / norm(sigma_em4(:));

% -----------------------------------------------------------------------
% Single target point: node closest to the centroid of all patch nodes.
% This avoids boundary effects from the sorted node ordering.
% All fields are (3,1) column vectors as required by kernel evaluators.
% -----------------------------------------------------------------------
rcen = mean(S.r, 2);                         % centroid of all nodes (3,1)
[~, ic] = min(vecnorm(S.r - rcen));          % nearest node to centroid
n0  = S.n(:, ic);  n0 = n0(:) / norm(n0);
r0  = S.r(:,  ic); r0  = r0(:);
du0 = S.du(:, ic); du0 = du0(:);
dv0 = S.dv(:, ic); dv0 = dv0(:);

targ_ext.r  = r0 + h * n0;
targ_ext.n  = n0;
targ_ext.du = du0;
targ_ext.dv = dv0;

targ_int.r  = r0 - h * n0;
targ_int.n  = n0;
targ_int.du = du0;
targ_int.dv = dv0;

nfail = 0;

fprintf('=== test_jump_relations ===\n\n');

% -----------------------------------------------------------------------
% Scalar / vector kernels
% -----------------------------------------------------------------------
coefs_lap = [0.7; 1.3];
coefs_h   = [1i*zk; 1.0];

% hypersingular=true: layer_eval (usematlab=0) path does not support off-surface
% targets for hypersingular kernels; only test usematlab=1 for those.
jump_tests = {
%  label          kernel                              density      expected jump (at ic)  hypersingular
  % 'lap  d',   kernel3d('l','d'),                sigma_r,    @(s)  s(ic),               false;
  % 'lap  sp',  kernel3d('l','sp'),               sigma_r,    @(s) -s(ic),               false;
  'lap  dp',  kernel3d('l','dp'),               sigma_r,    @(s) 0*s,               false;
  'lap  c',   kernel3d('l','c',coefs_lap),      sigma_r,    @(s)  coefs_lap(2)*s(ic),  false;
  'lap  cp',  kernel3d('l','cp',coefs_lap),     sigma_r,    @(s) -coefs_lap(1)*s(ic),  true;
  'helm d',   kernel3d('h','d',zk),             sigma_c,    @(s)  s(ic),               false;
  'helm sp',  kernel3d('h','sp',zk),            sigma_c,    @(s) -s(ic),               false;
  'helm  dp', kernel3d('h','dp',zk),            sigma_c,    @(s) 0*s,               false;
  'helm c',   kernel3d('h','c',zk,coefs_h),    sigma_c,    @(s)  coefs_h(2)*s(ic),    false;
  'helm cp',  kernel3d('h','cp',zk,coefs_h),   sigma_c,    @(s) -coefs_h(1)*s(ic),    true;
  'stok d',   kernel3d('stok','d'),             sigma_stok, @(s)  s(:,ic),             false;
  'stok sp',  kernel3d('stok','sp'),            sigma_stok, @(s) -s(:,ic),             false;
  'stok c',   kernel3d('stok','c',coefs_stok), sigma_stok, @(s)  coefs_stok(2)*s(:,ic), false;
};

for k = 1:size(jump_tests, 1)
    label    = jump_tests{k,1};
    kern     = jump_tests{k,2};
    sigma    = jump_tests{k,3};
    exp_fn   = jump_tests{k,4};
    hypsing  = jump_tests{k,5};

    exp_jump = exp_fn(sigma);

    % (1) Jump check, usematlab = 0 and 1
    for usematlab = [0, 1]
        opts_ke = struct('usematlab', usematlab);
        u_ext = surferkerneval(S, kern, sigma(:), targ_ext, eps, opts_ke);
        u_int = surferkerneval(S, kern, sigma(:), targ_int, eps, opts_ke);
        jump     = u_ext - u_int;
        % err_jump = norm(jump(:) - exp_jump(:)) / norm(exp_jump(:));
        err_jump = norm(jump(:) - exp_jump(:)) / norm(abs(u_ext(:)));
        if err_jump < tol, st = 'PASS'; else, st = 'FAIL'; nfail = nfail+1; end
        fprintf('  %-12s  jump err (usematlab=%d) = %.2e  [%s]\n', label, usematlab, err_jump, st);
    end
    % Use usematlab=1 values for the one-sided check below
    u_ext = surferkerneval(S, kern, sigma(:), targ_ext, eps, struct('usematlab', 1));

    % (2) One-sided: u_ext = surfermatapply(at ic) + exp_jump/2
    % For hypersingular kernels, layer_eval (usematlab=0) does not support
    % off-surface targets, so only test usematlab=1.
    um_vals = [0, 1];
    % if hypsing, um_vals = 1; end
    for usematlab = um_vals
        opts_ma  = struct('usematlab', usematlab);
        pot_full = surfermatapply(S, kern, sigma(:), eps, [], [], opts_ma);
        od       = kern.opdims(1);
        pot_ic   = pot_full((ic-1)*od+1 : ic*od);
        u_ext_pred = pot_ic(:) + exp_jump(:)/2;
        err_on = norm(u_ext(:) - u_ext_pred(:)) / (norm(u_ext(:)) + 1e-30);
        if err_on < tol, st = 'PASS'; else, st = 'FAIL'; nfail = nfail+1; end
        fprintf('  %-12s  one-sided err (usematlab=%d) = %.2e  [%s]\n', ...
                label, usematlab, err_on, st);
    end
end

% -----------------------------------------------------------------------
% EM nrccie-bc:  [[u]] = -sigma  (3 frame components)
% -----------------------------------------------------------------------
fprintf('\n');
alpha_em = 0.5;
kern_em3 = kernel3d('em', 'nrccie-bc', zk_em, alpha_em);

exp_jump3 = -sigma_em3(:, ic);

for usematlab = [0, 1]
    opts_ke3 = struct('usematlab', usematlab);
    u_ext3 = surferkerneval(S, kern_em3, sigma_em3(:), targ_ext, eps, opts_ke3);
    u_int3 = surferkerneval(S, kern_em3, sigma_em3(:), targ_int, eps, opts_ke3);
    err_em3 = norm((u_ext3 - u_int3) - exp_jump3) / norm(exp_jump3);
    if err_em3 < tol, st = 'PASS'; else, st = 'FAIL'; nfail = nfail+1; end
    fprintf('  %-12s  jump err (usematlab=%d) = %.2e  [%s]\n', 'em nrccie-bc', usematlab, err_em3, st);
end
u_ext3 = surferkerneval(S, kern_em3, sigma_em3(:), targ_ext, eps, struct('usematlab', 1));

for usematlab = [0, 1]
    opts_ma  = struct('usematlab', usematlab);
    pot_full3 = surfermatapply(S, kern_em3, sigma_em3(:), eps, [], [], opts_ma);
    pot_ic3   = pot_full3((ic-1)*3+1 : ic*3);
    u_ext3_pred = pot_ic3(:) + exp_jump3(:)/2;
    err_on3 = norm(u_ext3(:) - u_ext3_pred(:)) / (norm(u_ext3(:)) + 1e-30);
    if err_on3 < tol, st = 'PASS'; else, st = 'FAIL'; nfail = nfail+1; end
    fprintf('  %-12s  one-sided err (usematlab=%d) = %.2e  [%s]\n', ...
            'em nrccie-bc', usematlab, err_on3, st);
end

fprintf('\n');
% -----------------------------------------------------------------------
% EM nrccie-eval: jump checks on E and H at the single target
% -----------------------------------------------------------------------
kern_em4 = kernel3d('em', 'nrccie-eval', zk_em);

J  = sigma_em4(1:3, ic);

for usematlab = [0, 1]
    opts_ke = struct('usematlab', usematlab);
    u_ext4 = surferkerneval(S, kern_em4, sigma_em4(:), targ_ext, eps, opts_ke);
    u_int4 = surferkerneval(S, kern_em4, sigma_em4(:), targ_int, eps, opts_ke);
    diff4 = u_ext4 - u_int4;
    dE = diff4(1:3);
    dH = diff4(4:6);

    % n x [[E]] = 0
    nxdE = cross(n0, dE);
    err_nxE = norm(nxdE) / (norm(dE) + norm(dH) + 1e-30);
    if err_nxE < tol, st = 'PASS'; else, st = 'FAIL'; nfail = nfail+1; end
    fprintf('  %-12s  n x [[E]]    err (usematlab=%d) = %.2e  [%s]\n', 'em nrccie-ev', usematlab, err_nxE, st);

    % n x [[H]] = J
    nxdH = cross(n0, dH);
    err_nxH = norm(nxdH - J(:)) / norm(J(:));
    if err_nxH < tol, st = 'PASS'; else, st = 'FAIL'; nfail = nfail+1; end
    fprintf('  %-12s  n x [[H]]-J  err (usematlab=%d) = %.2e  [%s]\n', 'em nrccie-ev', usematlab, err_nxH, st);

    % n . [[H]] = 0
    err_ndH = abs(dot(n0, dH)) / (norm(dH) + 1e-30);
    if err_ndH < tol, st = 'PASS'; else, st = 'FAIL'; nfail = nfail+1; end
    fprintf('  %-12s  n . [[H]]    err (usematlab=%d) = %.2e  [%s]\n', 'em nrccie-ev', usematlab, err_ndH, st);

    % n . [[E]] = rho
    rho0 = sigma_em4(4, ic);
    err_ndE = abs(dot(n0, dE) - rho0) / abs(rho0);
    if err_ndE < tol, st = 'PASS'; else, st = 'FAIL'; nfail = nfail+1; end
    fprintf('  %-12s  n . [[E]]-rho err (usematlab=%d) = %.2e  [%s]\n', 'em nrccie-ev', usematlab, err_ndE, st);
end
% Use usematlab=1 values for one-sided check below
u_ext4 = surferkerneval(S, kern_em4, sigma_em4(:), targ_ext, eps, struct('usematlab',1));
u_int4 = surferkerneval(S, kern_em4, sigma_em4(:), targ_int, eps, struct('usematlab',1));
diff4   = u_ext4 - u_int4;

% One-sided for nrccie-eval
exp_jump4 = diff4;   % known jump at this single point
for usematlab = [0, 1]
    opts_ma  = struct('usematlab', usematlab);
    pot_full4 = surfermatapply(S, kern_em4, sigma_em4(:), eps, [], [], opts_ma);
    pot_ic4   = pot_full4((ic-1)*6+1 : ic*6);
    u_ext4_pred = pot_ic4(:) + exp_jump4(:)/2;
    err_on4 = norm(u_ext4(:) - u_ext4_pred(:)) / (norm(u_ext4(:)) + 1e-30);
    if err_on4 < tol, st = 'PASS'; else, st = 'FAIL'; nfail = nfail+1; end
    fprintf('  %-12s  one-sided err (usematlab=%d) = %.2e  [%s]\n', ...
            'em nrccie-ev', usematlab, err_on4, st);
end

% -----------------------------------------------------------------------
fprintf('\nDone. %d FAILED.\n', nfail);
if nfail > 0
    error('test_jump_relations: %d failure(s)', nfail);
end
