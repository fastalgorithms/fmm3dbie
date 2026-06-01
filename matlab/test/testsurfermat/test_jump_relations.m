% Verify jump relations and on-surface limits.
% For each kernel: pick a point near the patch centre, displace ±h along n,
% and check (1) u^+ - u^- matches the known discontinuity, and (2) the
% one-sided limit u^+ = surfermatapply(PV) + jump/2.

run ../../startup.m

rng(42);

S   = slicesurfer(geometries.disk(), 1);
eps = 1e-5;
tol = 0.05;
h   = 1e-4;

sigma_r = cos(S.r(1,:).' + 0.5*S.r(2,:).') .* exp(-S.r(3,:).'.^2);
sigma_r = sigma_r / norm(sigma_r);

zk    = 1.1 + 0.3i;
zk_em = 1.3;

sigma_c = (S.r(3,:).^2 + 1i*S.r(1,:)).';
sigma_c = 1 + 0*sigma_c / norm(sigma_c);

sigma_stok = [cos(S.r(2,:)); sin(S.r(3,:)); cos(S.r(1,:)+S.r(2,:))];
sigma_stok = sigma_stok / norm(sigma_stok(:));

J_cart = [cos(S.r(1,:)+S.r(2,:)); sin(S.r(2,:)+S.r(3,:)); cos(S.r(3,:))] ...
       + 1i*[sin(S.r(1,:)+S.r(3,:)); cos(S.r(1,:)); sin(S.r(2,:))];
J_cart = J_cart - S.n .* sum(S.n .* J_cart, 1);

sigma_em3 = [cos(S.r(1,:)+S.r(2,:)); sin(S.r(2,:)+S.r(3,:)); cos(S.r(3,:))] ...
          + 1i*[sin(S.r(1,:)); cos(S.r(2,:)); sin(S.r(3,:))];
sigma_em3 = sigma_em3 / norm(sigma_em3(:));

rho_em    = cos(S.r(3,:)) + 1i*cos(S.r(1,:)+S.r(2,:));
sigma_em4 = [J_cart; rho_em] / norm([J_cart; rho_em]);

[~, ic] = min(vecnorm(S.r - mean(S.r,2)));
n0  = S.n(:,ic)  / norm(S.n(:,ic));
r0  = S.r(:,ic);
du0 = S.du(:,ic);
dv0 = S.dv(:,ic);

te.r = r0+h*n0; te.n = n0; te.du = du0; te.dv = dv0;
ti.r = r0-h*n0; ti.n = n0; ti.du = du0; ti.dv = dv0;

%% Now run the tests

nfail = 0;
nfail = test_scalar_vector(nfail, S, eps, tol, zk, ic, n0, te, ti, sigma_r, sigma_c, sigma_stok);
nfail = test_nrccie_bc(nfail,   S, eps, tol, zk_em, ic, n0, te, ti, sigma_em3);
nfail = test_nrccie_eval(nfail, S, eps, tol, zk_em, ic, n0, te, ti, sigma_em4);

fprintf('\nDone. %d FAILED.\n', nfail);
if nfail > 0
    error('test_jump_relations: %d failure(s)', nfail);
end


function nfail = test_scalar_vector(nfail, S, eps, tol, zk, ic, ~, te, ti, sigma_r, sigma_c, sigma_stok)
% Laplace, Helmholtz, and Stokes kernels: jump and one-sided limit.

coefs_lap  = [0.7; 1.3];
coefs_h    = [1i*zk; 1.0];
coefs_stok = [0.7; 1.3];

jump_tests = {
  'lap  d',  kernel3d('l','d'),               sigma_r,    @(s)  s(ic);
  'lap  sp', kernel3d('l','sp'),              sigma_r,    @(s) -s(ic);
  'lap  dp', kernel3d('l','dp'),              sigma_r,    @(s) 0*s;
  'lap  c',  kernel3d('l','c',coefs_lap),     sigma_r,    @(s)  coefs_lap(2)*s(ic);
  'lap  cp', kernel3d('l','cp',coefs_lap),    sigma_r,    @(s) -coefs_lap(1)*s(ic);
  'helm d',  kernel3d('h','d',zk),            sigma_c,    @(s)  s(ic);
  'helm sp', kernel3d('h','sp',zk),           sigma_c,    @(s) -s(ic);
  'helm dp', kernel3d('h','dp',zk),           sigma_c,    @(s) 0*s;
  'helm c',  kernel3d('h','c',zk,coefs_h),    sigma_c,    @(s)  coefs_h(2)*s(ic);
  'helm cp', kernel3d('h','cp',zk,coefs_h),   sigma_c,    @(s) -coefs_h(1)*s(ic);
  'stok d',  kernel3d('stok','d'),             sigma_stok, @(s)  s(:,ic);
  'stok sp', kernel3d('stok','sp'),            sigma_stok, @(s) -s(:,ic);
  'stok c',  kernel3d('stok','c',coefs_stok), sigma_stok, @(s)  coefs_stok(2)*s(:,ic);
};

for k = 1:size(jump_tests, 1)
    label    = jump_tests{k,1};
    kern     = jump_tests{k,2};
    sigma    = jump_tests{k,3};
    exp_jump = jump_tests{k,4}(sigma);

    for usematlab = [0, 1]
        u_ext = surferkerneval(S, kern, sigma(:), te, eps, struct('usematlab', usematlab));
        u_int = surferkerneval(S, kern, sigma(:), ti, eps, struct('usematlab', usematlab));
        err   = norm((u_ext-u_int) - exp_jump) / norm(abs(u_ext(:)));
        nfail = nfail + (err >= tol);
        fprintf('  %-10s  jump      (usematlab=%d) err=%.2e  [%s]\n', label, usematlab, err, pf(err,tol));
    end

    u_ext = surferkerneval(S, kern, sigma(:), te, eps, struct('usematlab', 1));
    for usematlab = [0, 1]
        u_pv = surfermatapply(S, kern, sigma(:), eps, [], [], struct('usematlab', usematlab));
        od   = kern.opdims(1);
        u_pv = u_pv((ic-1)*od+1 : ic*od);
        err  = norm(u_ext(:) - (u_pv(:)+exp_jump(:)/2)) / (norm(u_ext(:))+1e-30);
        nfail = nfail + (err >= tol);
        fprintf('  %-10s  one-sided (usematlab=%d) err=%.2e  [%s]\n', label, usematlab, err, pf(err,tol));
    end
end

end


function nfail = test_nrccie_bc(nfail, S, eps, tol, zk_em, ic, ~, te, ti, sigma_em3)
% EM nrccie-bc: [[u]] = -sigma.

kern      = kernel3d('em', 'nrccie-bc', zk_em, 0.5);
exp_jump  = -sigma_em3(:, ic);

for usematlab = [0, 1]
    u_ext = surferkerneval(S, kern, sigma_em3(:), te, eps, struct('usematlab', usematlab));
    u_int = surferkerneval(S, kern, sigma_em3(:), ti, eps, struct('usematlab', usematlab));
    err   = norm((u_ext-u_int) - exp_jump) / norm(exp_jump);
    nfail = nfail + (err >= tol);
    fprintf('  %-12s  jump      (usematlab=%d) err=%.2e  [%s]\n', 'em nrccie-bc', usematlab, err, pf(err,tol));
end

u_ext = surferkerneval(S, kern, sigma_em3(:), te, eps, struct('usematlab', 1));
for usematlab = [0, 1]
    u_pv = surfermatapply(S, kern, sigma_em3(:), eps, [], [], struct('usematlab', usematlab));
    u_pv = u_pv((ic-1)*3+1 : ic*3);
    err  = norm(u_ext(:) - (u_pv(:)+exp_jump(:)/2)) / (norm(u_ext(:))+1e-30);
    nfail = nfail + (err >= tol);
    fprintf('  %-12s  one-sided (usematlab=%d) err=%.2e  [%s]\n', 'em nrccie-bc', usematlab, err, pf(err,tol));
end

end


function nfail = test_nrccie_eval(nfail, S, eps, tol, zk_em, ic, n0, te, ti, sigma_em4)
% EM nrccie-eval: n x [[E]]=0, n x [[H]]=J, n.[[H]]=0, n.[[E]]=rho, and one-sided.

kern = kernel3d('em', 'nrccie-eval', zk_em);
J    = sigma_em4(1:3, ic);

for usematlab = [0, 1]
    u_ext = surferkerneval(S, kern, sigma_em4(:), te, eps, struct('usematlab', usematlab));
    u_int = surferkerneval(S, kern, sigma_em4(:), ti, eps, struct('usematlab', usematlab));
    dE = u_ext(1:3) - u_int(1:3);
    dH = u_ext(4:6) - u_int(4:6);
    denom = norm(dE) + norm(dH) + 1e-30;

    e1 = norm(cross(n0,dE))  / denom;
    e2 = norm(cross(n0,dH) - J(:)) / norm(J(:));
    e3 = abs(dot(n0,dH))     / (norm(dH)+1e-30);
    e4 = abs(dot(n0,dE) - sigma_em4(4,ic)) / abs(sigma_em4(4,ic));

    nfail = nfail + (e1>=tol) + (e2>=tol) + (e3>=tol) + (e4>=tol);
    fprintf('  %-12s  n x [[E]]    (usematlab=%d) err=%.2e  [%s]\n', 'em nrccie-ev', usematlab, e1, pf(e1,tol));
    fprintf('  %-12s  n x [[H]]-J  (usematlab=%d) err=%.2e  [%s]\n', 'em nrccie-ev', usematlab, e2, pf(e2,tol));
    fprintf('  %-12s  n . [[H]]    (usematlab=%d) err=%.2e  [%s]\n', 'em nrccie-ev', usematlab, e3, pf(e3,tol));
    fprintf('  %-12s  n . [[E]]-rho(usematlab=%d) err=%.2e  [%s]\n', 'em nrccie-ev', usematlab, e4, pf(e4,tol));
end

u_ext = surferkerneval(S, kern, sigma_em4(:), te, eps, struct('usematlab', 1));
u_int = surferkerneval(S, kern, sigma_em4(:), ti, eps, struct('usematlab', 1));
exp_jump4 = u_ext - u_int;

for usematlab = [0, 1]
    u_pv = surfermatapply(S, kern, sigma_em4(:), eps, [], [], struct('usematlab', usematlab));
    u_pv = u_pv((ic-1)*6+1 : ic*6);
    err  = norm(u_ext(:) - (u_pv(:)+exp_jump4(:)/2)) / (norm(u_ext(:))+1e-30);
    nfail = nfail + (err >= tol);
    fprintf('  %-12s  one-sided (usematlab=%d) err=%.2e  [%s]\n', 'em nrccie-ev', usematlab, err, pf(err,tol));
end

end


function s = pf(err, tol)
if err >= tol, s = 'FAIL'; else, s = 'PASS'; end
end
