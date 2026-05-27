%TEST_JUMP_RELATIONS   Verify jump relations via surferkerneval.
%
%  For each kernel with a known jump we evaluate the layer potential at
%  targets just outside (+h along normal) and just inside (-h), then check
%  that the difference matches the predicted discontinuity.
%
%  Accuracy of the smooth quadrature is already tested in test_kernel3d.m;
%  here we only need the near-field corrections to produce the right
%  O(h)-accurate jump.  A single patch is used so the test is fast.
%
%  Jump relations tested:
%    lap  d   :  u^+ - u^-  =  +sigma          (scalar)
%    lap  sp  :  u^+ - u^-  =  -sigma          (scalar)
%    lap  cp  :  u^+ - u^-  =  -coefs(1)*sigma  (S' part only; D' has no classical jump)
%    helm d   :  u^+ - u^-  =  +sigma          (scalar)
%    helm sp  :  u^+ - u^-  =  -sigma          (scalar)
%    helm cp  :  u^+ - u^-  =  -coefs(1)*sigma
%    stok d   :  u^+ - u^-  =  +sigma          (3-vector)
%    stok sp  :  u^+ - u^-  =  -sigma          (3-vector)
%    stok c   :  u^+ - u^-  =  coefs(2)*sigma  (D part only)
%    em nrccie-bc  :  u^+ - u^-  =  -sigma     (all 3 frame components;
%                      M_k jumps give -J_tan, S'_k gives -rho, alpha terms don't jump)
%    em nrccie-eval:
%      n x [[E]] = 0     (E = ik S_k[J] - grad S_k[rho]; S_k has no jump)
%      n x [[H]] = J     (J = tangential Cartesian current in density)
%      n . [[H]] = 0
%      n . [[E]] = rho   (normal component of E jumps by rho; from grad S_k[rho])

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
% Smooth densities (trigonometric, so the near-field quadrature is not
% fighting a rough density)
% -----------------------------------------------------------------------
sigma_r = cos(S.r(1,:).' + 0.5*S.r(2,:).') .* exp(-S.r(3,:).'.^2);
sigma_r = sigma_r / norm(sigma_r);

zk      = 1.1 + 0.3i;
sigma_c = (S.r(3,:).^2 + 1i*S.r(1,:)).';
sigma_c = sigma_c / norm(sigma_c);

zk_em = 1.3;

% Stokes: smooth 3-vector density (not necessarily tangential; the jump
% relation holds for any density, tangential or not)
coefs_stok = [0.7; 1.3];
sigma_stok = [cos(S.r(2,:)); sin(S.r(3,:)); cos(S.r(1,:) + S.r(2,:))];
sigma_stok = sigma_stok / norm(sigma_stok(:));

% EM nrccie-bc density: [j_ru, j_rv, rho] in orthonormal tangent frame (3 x npts)
sigma_em3 = [cos(S.r(1,:) + S.r(2,:)); ...
             sin(S.r(2,:) + S.r(3,:)); ...
             cos(S.r(3,:))] ...
          + 1i*[sin(S.r(1,:)); cos(S.r(2,:)); sin(S.r(3,:))];
sigma_em3 = sigma_em3 / norm(sigma_em3(:));

% EM nrccie-eval density: [Jx, Jy, Jz, rho] (4 x npts), J projected
% tangential so the jump relation n x [[H]] = J is clean
J_cart = [cos(S.r(1,:) + S.r(2,:)); ...
          sin(S.r(2,:) + S.r(3,:)); ...
          cos(S.r(3,:))] ...
       + 1i*[sin(S.r(1,:) + S.r(3,:)); ...
             cos(S.r(1,:));             ...
             sin(S.r(2,:))];
J_cart = J_cart - S.n .* sum(S.n .* J_cart, 1);   % project tangential
rho_em = cos(S.r(3,:)) + 1i*cos(S.r(1,:) + S.r(2,:));
sigma_em4 = [J_cart; rho_em];
sigma_em4 = sigma_em4 / norm(sigma_em4(:));

% -----------------------------------------------------------------------
% Exterior / interior targets
% -----------------------------------------------------------------------
nn = S.n ./ vecnorm(S.n);

targ_ext.r  = S.r + h * nn;
targ_ext.n  = nn;
targ_ext.du = S.du;  targ_ext.dv = S.dv;

targ_int.r  = S.r - h * nn;
targ_int.n  = nn;
targ_int.du = S.du;  targ_int.dv = S.dv;

nfail = 0;

fprintf('=== test_jump_relations ===\n\n');

% -----------------------------------------------------------------------
% Scalar / vector kernels
% -----------------------------------------------------------------------
coefs_lap = [0.7; 1.3];
coefs_h   = [1i*zk; 1.0];

jump_tests = {
%  label          kernel                              density      expected jump
  'lap  d',   kernel3d('l','d'),                sigma_r,    @(s)  s;
  'lap  sp',  kernel3d('l','sp'),               sigma_r,    @(s) -s;
  'lap  c',   kernel3d('l','c',coefs_lap),      sigma_r,    @(s) coefs_lap(2)*s;
  'lap  cp',  kernel3d('l','cp',coefs_lap),     sigma_r,    @(s) -coefs_lap(1)*s;
  'helm d',   kernel3d('h','d',zk),             sigma_c,    @(s)  s;
  'helm sp',  kernel3d('h','sp',zk),            sigma_c,    @(s) -s;
  'helm c',   kernel3d('h','c',zk,coefs_h),    sigma_c,    @(s) coefs_h(2)*s;
  'helm cp',  kernel3d('h','cp',zk,coefs_h),   sigma_c,    @(s) -coefs_h(1)*s;
  'stok d',   kernel3d('stok','d'),             sigma_stok, @(s)  s;
  'stok sp',  kernel3d('stok','sp'),            sigma_stok, @(s) -s;
  'stok c',   kernel3d('stok','c',coefs_stok), sigma_stok, @(s) coefs_stok(2)*s;
};

for k = 1:size(jump_tests, 1)
    label    = jump_tests{k,1};
    kern     = jump_tests{k,2};
    sigma    = jump_tests{k,3};
    exp_fn   = jump_tests{k,4};

    u_ext = surferkerneval(S, kern, sigma(:), targ_ext, eps);
    u_int = surferkerneval(S, kern, sigma(:), targ_int, eps);
    jump  = reshape(u_ext - u_int, size(sigma));

    exp_jump = exp_fn(sigma);
    err = norm(jump(:) - exp_jump(:)) / norm(exp_jump(:));

    if err < tol, status = 'PASS'; else, status = 'FAIL'; nfail = nfail+1; end
    fprintf('  %-12s  jump err = %.2e  [%s]\n', label, err, status);
end

% -----------------------------------------------------------------------
% EM nrccie-bc:  u^+ - u^- = -sigma  (all 3 orthonormal-frame components)
%   zvec3 = alpha*nxnxE - n x curl(S_k[J]):
%     n x curl(S_k[J]) has jump +J_tan/2 from each side => [[n x curl S_k]] = +J_tan
%     => [[-zvec]] = -J_tan components  (rows 1-2)
%     nxnxE uses S_k which has no jump => [[nxnxE]] = 0
%     => [[zvec3]] = -J_tan
%   pot_rho: contains S'_k[rho], jump -rho/2 each side => [[pot_rho]] = -rho
%   Overall: u^+ - u^- = -sigma
% -----------------------------------------------------------------------
fprintf('\n');
alpha_em = 0.5;
kern_em3 = kernel3d('em', 'nrccie-bc', zk_em, alpha_em);

u_ext3 = surferkerneval(S, kern_em3, sigma_em3(:), targ_ext, eps);
u_int3 = surferkerneval(S, kern_em3, sigma_em3(:), targ_int, eps);
jump3  = reshape(u_ext3 - u_int3, size(sigma_em3));

err_em3 = norm(jump3(:) + sigma_em3(:)) / norm(sigma_em3(:));
if err_em3 < tol, status = 'PASS'; else, status = 'FAIL'; nfail = nfail+1; end
fprintf('  %-12s  jump err = %.2e  [%s]\n', 'em nrccie-bc', err_em3, status);

fprintf('\n');
% -----------------------------------------------------------------------
% EM nrccie-eval:
%   Density: [Jx,Jy,Jz,rho] interleaved (4 components / point)
%   Output:  [Ex,Ey,Ez,Hx,Hy,Hz] interleaved (6 components / point)
%
%   E = ik S_k[J] - grad S_k[rho]:
%     S_k has no jump, but grad S_k[rho] has a jump in the normal direction
%     equal to -rho (same as S'_k), so [[E]] = (0,0,..., rho*n) in normal dir
%     => n x [[E]] = 0,  n . [[E]] = rho
%   H = curl S_k[J]:  jump of curl S_k is -n x J (Calderon),
%     so [[H]] = -n x J,  and n x [[H]] = n x (-n x J) = J_tan = J (since J tangential)
%     n . [[H]] = n . (-n x J) = 0 (cross product is tangential)
% -----------------------------------------------------------------------
kern_em4 = kernel3d('em', 'nrccie-eval', zk_em);

u_ext4 = surferkerneval(S, kern_em4, sigma_em4(:), targ_ext, eps);
u_int4 = surferkerneval(S, kern_em4, sigma_em4(:), targ_int, eps);
diff4  = reshape(u_ext4 - u_int4, 6, []);   % (6, npts): [Ex,Ey,Ez,Hx,Hy,Hz]

dE = diff4(1:3, :);
dH = diff4(4:6, :);
J  = sigma_em4(1:3, :);   % tangential Cartesian current

% n x [[E]] = 0
nxdE = cross(nn, dE, 1);
err_nxE = norm(nxdE(:)) / (norm(dE(:)) + norm(dH(:)) + 1e-30);
if err_nxE < tol, status = 'PASS'; else, status = 'FAIL'; nfail = nfail+1; end
fprintf('  %-12s  n x [[E]]  err = %.2e  [%s]\n', 'em nrccie-ev', err_nxE, status);

% n x [[H]] = J
nxdH = cross(nn, dH, 1);
err_nxH = norm(nxdH(:) - J(:)) / norm(J(:));
if err_nxH < tol, status = 'PASS'; else, status = 'FAIL'; nfail = nfail+1; end
fprintf('  %-12s  n x [[H]]-J err = %.2e  [%s]\n', 'em nrccie-ev', err_nxH, status);

% n . [[H]] = 0
ndotdH = sum(nn .* dH, 1);
err_ndH = norm(ndotdH(:)) / (norm(dH(:)) + 1e-30);
if err_ndH < tol, status = 'PASS'; else, status = 'FAIL'; nfail = nfail+1; end
fprintf('  %-12s  n . [[H]]  err = %.2e  [%s]\n', 'em nrccie-ev', err_ndH, status);

% n . [[E]] = rho  (rho = 4th row of sigma_em4, real Cartesian scalar)
rho = sigma_em4(4, :);   % (1, npts)
ndotdE = sum(nn .* dE, 1);
err_ndE = norm(ndotdE(:) - rho(:)) / norm(rho(:));
if err_ndE < tol, status = 'PASS'; else, status = 'FAIL'; nfail = nfail+1; end
fprintf('  %-12s  n . [[E]]-rho err = %.2e  [%s]\n', 'em nrccie-ev', err_ndE, status);

% -----------------------------------------------------------------------
fprintf('\nDone. %d FAILED.\n', nfail);
if nfail > 0
    error('test_jump_relations: %d failure(s)', nfail);
end
