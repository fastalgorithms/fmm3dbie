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
    K  = (-kernel3d('laplace', tp)+ 2*kernel3d('helmholtz', zk, tp))/2;

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
K = kernel3d('laplace', 'c', coefs);

mat     = K.eval(src, targ);
pot_mat = mat * sigma;
pot_fmm = K.fmm(eps, src, targ, sigma);

err = norm(pot_fmm - pot_mat) / norm(pot_mat);
status = 'PASS';
if err > tol, status = 'FAIL'; end
fprintf('lap3d  type=''c'' coefs=[%.1f %.1f]   rel err = %.2e   [%s]\n', ...
    coefs(1), coefs(2), err, status);

% -------------------------------------------------------------------------
fprintf('\n=== kernel3d helm3d tests ===\n\n');

zk = 3.0 + 0.0i;

types_h = {'s', 'd', 'sp'};

for k = 1:length(types_h)
    tp = types_h{k};
    K  = kernel3d('helmholtz', zk, tp);

    mat     = K.eval(src, targ);
    pot_mat = mat * sigma;

    pot_fmm = K.fmm(eps, src, targ, sigma);

    err = norm(pot_fmm - pot_mat) / norm(pot_mat);
    status = 'PASS';
    if err > tol, status = 'FAIL'; end
    fprintf('helm3d  zk=%.1f  type=''%s''   rel err = %.2e   [%s]\n', ...
        real(zk), tp, err, status);
end

% Combined layer
coefs_h = [1i*zk; 1.0];
K = kernel3d('helmholtz', zk, 'c', coefs_h);

mat     = K.eval(src, targ);
pot_mat = mat * sigma;
pot_fmm = K.fmm(eps, src, targ, sigma);

err = norm(pot_fmm - pot_mat) / norm(pot_mat);
status = 'PASS';
if err > tol, status = 'FAIL'; end
fprintf('helm3d  zk=%.1f  type=''c'' coefs=[%.1fi %.1f]   rel err = %.2e   [%s]\n', ...
    real(zk), imag(coefs_h(1)), real(coefs_h(2)), err, status);

fprintf('\nDone.\n');
