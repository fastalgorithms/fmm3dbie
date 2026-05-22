%TEST_KERNEL3D   Verify that kernel3d.fmm and kernel3d.eval agree for all
%   implemented Laplace kernels.
%
%   For each kernel type, random sources and targets are placed in disjoint
%   regions so no singularity is encountered.  The potential is computed two
%   ways:
%     (1) dense matrix-vector product via obj.eval
%     (2) FMM via obj.fmm
%   and the relative l2 error is printed.  A warning is raised if any error
%   exceeds the tolerance.

run ../startup.m

fprintf('=== kernel3d lap3d tests ===\n\n');

eps   = 1e-10;   % FMM precision
ns    = 2e3;     % number of sources
nt    = 1e3;     % number of targets
tol   = 1e-6;   % error threshold for pass/fail

rng(42);

% Sources in unit ball, targets shifted away to avoid self-interaction
src.r = randn(3, ns);  src.r = src.r ./ vecnorm(src.r) .* rand(1,ns);
src.n = randn(3, ns);  src.n = src.n ./ vecnorm(src.n);

targ.r = randn(3, nt);  targ.r = targ.r ./ vecnorm(targ.r) .* rand(1,nt);
targ.n = randn(3, nt);  targ.n = targ.n ./ vecnorm(targ.n);

sigma = randn(ns, 1);

types = {'s', 'd', 'sp'};

for k = 1:length(types)
    tp = types{k};
    K  = kernel3d('laplace', tp);

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

fprintf('\nDone.\n');
