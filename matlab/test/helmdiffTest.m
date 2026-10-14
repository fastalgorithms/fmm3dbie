% helm3d.helmdiff_green and the helm3d _diff kernels built on it with
% kernel3d(g, type): a0*K_{k0} - a1*K_{k1} as one kernel, with the
% singular parts cancelled analytically.

% run ../startup.m

rng(42);

test_green();
test_scalar_kernels();
test_block_kernels();
test_fmm();
test_fmm_apply();
test_getquad();


function test_green()
% limits at small r, and the plain difference away from 0

k = [1.1 + 0.3i; 2.3];

% a0 = a1: G -> i(k0-k1)/(4 pi), G'' bounded, G'/r ~ -(k0^2-k1^2)/(8 pi r)
for r = [1e-4 1e-8 1e-10]
    [v, ~, h] = helm3d.helmdiff_green(k, [1; 1], [0;0;0], [r;0;0]);
    assert(abs(v - 1i*(k(1)-k(2))/(4*pi)) < 10*r, 'helmdiff_green: G limit');
    assert(abs(h(1,1,2,2)*r + (k(1)^2-k(2)^2)/(8*pi)) < 10*r, ...
           'helmdiff_green: G''/r limit');
    assert(abs(h(1,1,1,1)) < 1, 'helmdiff_green: G'''' not bounded');
end

s = randn(3, 4);  t = randn(3, 5) + [1.5;0;0];  a = [0.7; 1.2];
[v, g, h] = helm3d.helmdiff_green(k, a, s, t);
[v0, g0, h0] = helm3d.green(k(1), s, t);  [v1, g1, h1] = helm3d.green(k(2), s, t);
err = max([norm(v(:) - (a(1)*v0(:) - a(2)*v1(:)))/norm(v(:)), ...
           norm(g(:) - (a(1)*g0(:) - a(2)*g1(:)))/norm(g(:)), ...
           norm(h(:) - (a(1)*h0(:) - a(2)*h1(:)))/norm(h(:))]);
assert(err < 1e-13, 'helmdiff_green vs difference: %.2e', err);

end


function test_scalar_kernels()
% one radgreen kernel each; kernel_order = ipv-1; ipv as helm3d, or 0
% when a0 = a1; eval equal to the explicit difference

k = [1.1 + 0.3i; 2.3];  comb = [1i*k(1); 1];
[src, targ] = points();
ipv_helm = struct('s', 0, 'd', 1, 'sp', 1, 'dp', 2, 'c', 1, 'cp', 2);
for type = {'s', 'd', 'sp', 'dp', 'c', 'cp'}
    type = type{1};
    hc = {};  if any(strcmp(type, {'c', 'cp'})), hc = {comb}; end
    for a = [[1; 1], [0.7; 1.2]]
        K = kernel3d('h', [type '_diff'], k, [hc{:}; a]);
        Kh = kernel3d.helm3d(type, k(1), hc{:});
        assert(strcmp(K.name, 'radgreen'), 'helm %s_diff: not radgreen', type);
        ipv = ipv_helm.(type);
        if a(1) == a(2), ipv = 0; end
        assert(K.params.ipv == ipv, 'helm %s_diff a = %s: ipv %d', type, ...
               mat2str(a.'), K.params.ipv);
        assert(K.kernel_order == ipv - 1, 'helm %s_diff: kernel_order', type);

        B = a(1)*Kh.eval(src, targ) - a(2)*kernel3d.helm3d(type, k(2), hc{:}).eval(src, targ);
        err = norm(K.eval(src, targ) - B, 'fro')/norm(B, 'fro');
        assert(err < 1e-13, 'helm %s_diff eval: %.2e', type, err);
    end
end

end


function test_block_kernels()
% s2trans_diff, d2trans_diff, c2trans_diff and uniform trans_sys_diff are
% single radgreen kernels; kernel_order = max(ipv)-1; blocks match the
% scalar kernels

k = [1.1 + 0.3i; 2.3];  a = [1.2; 0.8];  comb = [1i*k(1); 1];
[src, targ] = points();
cases = {'s2trans', 's', 'sp', {}; 'd2trans', 'd', 'dp', {}; ...
         'c2trans', 'c', 'cp', {comb}};
for c = 1:size(cases, 1)
    [type, t1, t2, hc] = cases{c,:};
    K = kernel3d('h', [type '_diff'], k, [hc{:}; a]);
    assert(strcmp(K.name, 'radgreen') && isequal(K.opdims, [2 1]), ...
           'helm %s_diff: not one 2x1 kernel', type);
    assert(K.kernel_order == max(K.params.ipv) - 1, ...
           'helm %s_diff: kernel_order', type);
    A = K.eval(src, targ);
    B1 = kernel3d('h', [t1 '_diff'], k, [hc{:}; a]).eval(src, targ);
    B2 = kernel3d('h', [t2 '_diff'], k, [hc{:}; a]).eval(src, targ);
    err = max(norm(A(1:2:end,:) - B1, 'fro')/norm(B1, 'fro'), ...
              norm(A(2:2:end,:) - B2, 'fro')/norm(B2, 'fro'));
    assert(err < 1e-13, 'helm %s_diff blocks: %.2e', type, err);
end

C = [0.7+0.1i, 0.3; 0.5, 1.2i];
for c1 = {0.9*C, [0.2 1.1; 0.4 0.6]}               % uniform ratio, or not
    K = kernel3d('h', 'trans_sys_diff', k, cat(3, C, c1{1}));
    assert(isequal(K.opdims, [2 2]), 'helm trans_sys_diff: opdims');
    if isfield(K.params, 'ipv')                      % one radgreen kernel
        assert(K.kernel_order == max(K.params.ipv) - 1, ...
               'helm trans_sys_diff: kernel_order');
    end
    A = K.eval(src, targ);
    blk = {'d', 's'; 'dp', 'sp'};
    for i = 1:2, for j = 1:2
        B = C(i,j)*kernel3d.helm3d(blk{i,j}, k(1)).eval(src, targ) ...
            - c1{1}(i,j)*kernel3d.helm3d(blk{i,j}, k(2)).eval(src, targ);
        err = norm(A(i:2:end, j:2:end) - B, 'fro')/norm(B, 'fro');
        assert(err < 1e-12, 'helm trans_sys_diff block (%d,%d): %.2e', i, j, err);
    end, end
end
K = kernel3d('h', 'trans_diff', k, cat(3, C, 0.9*C));
assert(strcmp(K.name, 'radgreen'), 'helm trans_diff: uniform ratio not one kernel');

end


function test_fmm()
% the difference of helmholtz fmms against eval, for every type

k = [1.1 + 0.3i; 2.3];  a = [1.2; 0.8];  comb = [1i*k(1); 1];
C = [0.7+0.1i, 0.3; 0.5, 1.2i];
[src, targ] = points();
targ.r = targ.r + [3;0;0];
cases = {'s', {}; 'd', {}; 'sp', {}; 'dp', {}; 'c', {comb}; 'cp', {comb}; ...
         's2trans', {}; 'd2trans', {}; 'c2trans', {comb}};
for c = 1:size(cases, 1)
    [type, hc] = cases{c,:};
    check_fmm(kernel3d('h', [type '_diff'], k, [hc{:}; a]), src, targ, type);
end
check_fmm(kernel3d('h', 'trans_sys_diff', k, cat(3, C, 0.9*C)), src, targ, ...
          'trans_sys (one kernel)');
check_fmm(kernel3d('h', 'trans_sys_diff', k, cat(3, C, [0.2 1.1; 0.4 0.6])), ...
          src, targ, 'trans_sys (per block)');

end


function check_fmm(K, src, targ, label)
assert(~isempty(K.fmm), 'helm %s_diff: no fmm', label);
ns = size(src.r, 2);
sig = complex(randn(K.opdims(2)*ns, 1), randn(K.opdims(2)*ns, 1));
A = K.eval(src, targ)*sig;
err = norm(K.fmm(1e-10, src, targ, sig) - A)/norm(A);
assert(err < 1e-8, 'helm %s_diff fmm: %.2e', label, err);
end

function test_fmm_apply()
% the fmm with the radgreen near corrections in surfermatapply, against
% the dense matrix

k = [1.1 + 0.3i; 2.3];  a = [1; 1];  comb = [1i*k(1); 1];  eps = 1e-8;
S = geometries.sphere(1, 2, [0;0;0], 4, 1);
[~, idx] = sort(vecnorm(S.cms - S.cms(:,1)));
S = slicesurfer(S, idx(1:3));                    % three neighbouring patches
for type = {'s', 'sp', 'dp', 'c', 's2trans', 'd2trans'}
    hc = {};  if strcmp(type{1}, 'c'), hc = {comb}; end
    K = kernel3d('h', [type{1} '_diff'], k, [hc{:}; a]);
    sig = complex(randn(K.opdims(2)*S.npts, 1), randn(K.opdims(2)*S.npts, 1));
    A = surfermat(S, K, eps)*sig;
    B = surfermatapply(S, K, sig, eps);
    err = norm(A - B)/norm(A);
    assert(err < 1e-6, 'helm %s_diff surfermatapply fmm: %.2e', type{1}, err);
end

end


function test_getquad()
% near quadrature against a smooth integral on a curved patch, including
% dp_diff whose G'' is bounded while G'/r ~ 1/r

k = [1.1 + 0.3i; 2.3];
for a = [[1; 1], [0.7; 1.2]]
    for type = {'s', 'd', 'sp', 'dp'}
        K = kernel3d('h', [type{1} '_diff'], k, a);
        err = kernel3d.check_quad(K, 1e-7);
        assert(err < 1e-6, 'helm %s_diff check: %.2e', type{1}, err);
    end
end
err = kernel3d.check_quad(kernel3d('h', 's2trans_diff', k, [1; 1]), 1e-7);
assert(err < 1e-6, 'helm s2trans_diff check: %.2e', err);

end


function [src, targ] = points()
src.r  = randn(3, 30);
src.n  = randn(3, 30);  src.n  = src.n./vecnorm(src.n);
targ.r = randn(3, 20) + [2;0;0];
targ.n = randn(3, 20);  targ.n = targ.n./vecnorm(targ.n);
end
