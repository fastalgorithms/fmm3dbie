%TEST_KERNEL3D   Verify that kernel3d.fmm and kernel3d.eval agree for all
%   implemented Laplace and Helmholtz kernels.
%
%   For each kernel type, random sources and targets are placed in disjoint
%   regions so no singularity is encountered.  The potential is computed two
%   ways:
%     (1) dense matrix-vector product via obj.eval
%     (2) FMM via obj.fmm
%   and the relative l2 error is printed.  A warning is raised if any error
%   exceeds the tolerance.

run ../startup.m

zk = 3.0 + 0.0i;
eps   = 1e-6;   % FMM precision
ns    = 2e3;    % number of sources
nt    = 1e3;    % number of targets
tol   = 1e-5;   % error threshold for pass/fail

rng(42);

% Sources in unit ball, targets in an annulus [2,3] to avoid self-interaction
src.r = randn(3, ns);  src.r = src.r ./ vecnorm(src.r) .* rand(1,ns);
src.n = randn(3, ns);  src.n = src.n ./ vecnorm(src.n);

targ.r = randn(3, nt);  targ.r = targ.r ./ vecnorm(targ.r) .* (2 + rand(1,nt));
targ.n = randn(3, nt);  targ.n = targ.n ./ vecnorm(targ.n);

sigma = randn(ns, 1);

% -------------------------------------------------------------------------
fprintf('=== kernel3d lap3d tests ===\n\n');

types = {'s', 'd', 'sp'};

for k = 1:length(types)
    tp = types{k};
    K  = (-kernel3d('laplace', tp)+ 2*kernel3d('helmholtz', tp, zk))/2;

    % Dense eval: (nt x ns) matrix times sigma
    mat     = K.eval(src, targ);
    pot_mat = mat * sigma;

    % FMM eval
    pot_fmm = K.fmm(eps, src, targ, sigma);

    err = norm(pot_fmm - pot_mat) / norm(pot_mat);
    status = 'PASS';
    if err > tol, status = 'FAIL'; end
    fprintf('lap3d  type=''%s''   rel err = %.2e   [%s]\n', tp, err, status);
end

% Combined layer
coefs = [0.9; -0.5];
coefs_h = [1i*zk; 1.0];
K = kernel3d('laplace', 'c', coefs)+kernel3d('helmholtz', 'c',zk, coefs_h);

mat     = K.eval(src, targ);
pot_mat = mat * sigma;
pot_fmm = K.fmm(eps, src, targ, sigma);

err = norm(pot_fmm - pot_mat) / norm(pot_mat);
status = 'PASS';
if err > tol, status = 'FAIL'; end
fprintf('lap3d  type=''c'' coefs=[%.1f %.1f]   rel err = %.2e   [%s]\n', ...
    coefs(1), coefs(2), err, status);

% -------------------------------------------------------------------------
fprintf('\n=== kernel3d helm3d cprime tests ===\n\n');

zk2 = 1.3 + 0.2i;
coefs_cp = [1.0 + 0.5i; 2.0 - 0.3i];

kern_c  = kernel3d('h', 'c',      zk2, coefs_cp);   % alpha*S + beta*D
kern_cp = kernel3d('h', 'cprime', zk2, coefs_cp);   % alpha*S' + beta*D'

% --- Test 1: kern.eval agrees with finite-difference d/dn_x of kern_c ---
fprintf('cprime eval vs finite-difference of combined ...\n');

h = 1e-5;
mat_cp  = kern_cp.eval(src, targ);
% build targ +/- h*n
targ_p = targ; targ_p.r = targ.r + h*targ.n;
targ_m = targ; targ_m.r = targ.r - h*targ.n;
mat_fd = (kern_c.eval(src, targ_p) - kern_c.eval(src, targ_m)) / (2*h);

err_fd = norm(mat_cp - mat_fd, 'fro') / norm(mat_fd, 'fro');

status = 'PASS'; if err_fd > 1e-4, status = 'FAIL'; end
fprintf('  rel err (eval vs FD): %.2e   [%s]\n', err_fd, status);

if ~isempty(kern_cp.fmm)
sigma = randn(size(targ.r,2),1);
pot_fmm = kern_cp.fmm(eps, src, targ, sigma);
err_fmm = norm(mat_cp*sigma - pot_fmm, 'fro') / norm(mat_fd*sigma, 'fro');
status = 'PASS'; if err_fmm > 1e-4, status = 'FAIL'; end
fprintf('  rel err (eval vs fmm): %.2e   [%s]\n', err_fmm, status);
end
% --- Test 2: layer_eval of cprime agrees with finite-diff of layer_eval c ---
fprintf('cprime layer_eval vs finite-difference of combined layer_eval ...\n');

S = geometries.sphere(1, 2, [0;0;0], 4, 1);
eps_lp = 1e-7;
sigma_s = randn(S.npts, 1);

% off-surface targets, well away from surface
targ_off.r = randn(3, 20); targ_off.r = targ_off.r ./ vecnorm(targ_off.r) * 2.5;
targ_off.n = randn(3, 20); targ_off.n = targ_off.n ./ vecnorm(targ_off.n);

pot_c_p = kern_c.layer_eval(S, sigma_s, struct('r', targ_off.r + h*targ_off.n), eps_lp);
pot_c_m = kern_c.layer_eval(S, sigma_s, struct('r', targ_off.r - h*targ_off.n), eps_lp);
pot_fd   = (pot_c_p - pot_c_m) / (2*h);

pot_cp = kern_cp.layer_eval(S, sigma_s, targ_off, eps_lp);

err_lp = norm(pot_cp - pot_fd) / norm(pot_fd);
status = 'PASS'; if err_lp > 1e-4, status = 'FAIL'; end
fprintf('  rel err (layer_eval vs FD): %.2e   [%s]\n', err_lp, status);

fprintf('\nDone.\n');
