% Test the near corrections to the oversampled smooth rule:
% smooth_sparse_quad, getquad opts.corrections (Fortran subtraction),
% propagation of getquad_sub through kernel algebra, and
% kernel3d.getquad_corrections.

%% Now run the tests

test_smooth_sparse_quad();
test_fortran_corrections();
test_getquad_sub_flags();
test_getquad_corrections();
test_kernel_order();


function test_smooth_sparse_quad()
% smooth_sparse_quad agrees with the dense oversampled smooth rule
% restricted to the near field.

S   = geometries.sphere(1, 1, [0;0;0], 4, 1);
tol = 1e-13;

kscal = kernel3d.helm3d('c', 1.3, [1; 0.5]);
kblk  = kernel3d(@(s,t) kron(lap3d.kern(s,t,'s'), [1 2i; 3 4]));
toff  = struct('r', 1.05*S.r(:,1:60), 'n', S.n(:,1:60));
nvec  = 5 + mod(1:S.npatches, 3).';

check_ssq(kscal, S, S, 7, tol, 'scalar, on surface');
check_ssq(kblk,  S, S, 7, tol, '2x2, on surface');
check_ssq(kblk,  S, toff, 6, tol, '2x2, off surface');
check_ssq(kscal, S, toff, NaN, tol, 'scalar, no oversampling');
check_ssq(kscal, S, S, nvec, tol, 'scalar, per patch nover');

end


function check_ssq(K, S, targ, nover, tol, label)

m = K.opdims(1); n = K.opdims(2);
rsc = getnear(S, targ);
[wsmth, iquad] = smooth_sparse_quad(K, targ, S, rsc.row_ptr, rsc.col_ind, nover);
assert(isequal(iquad(:), rsc.iquad(:)), 'ssq %s: iquad differs from getnear', label);
assert(isequal(size(wsmth), [m*n, rsc.nquad]), 'ssq %s: wsmth size', label);
ri = kernel3d.rsc_interleave_full(m, n);
A  = conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, wsmth, ri);

% dense reference
if any(isnan(nover))
    Sov = S; X = speye(S.npts);
else
    [Sov, X] = oversample(S, nover);
end
w   = repmat(Sov.wts(:).', n, 1);
Ad  = (K.eval_mask(Sov, targ) .* w(:).') * kron(X, eye(n));
msk = conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, ones(m*n, rsc.nquad), ri) ~= 0;
Ad  = Ad .* msk;

err = full(max(abs(A(:) - Ad(:)))) / max(abs(Ad(:)));
assert(err < tol, 'ssq %s: error %g', label, err);

end


function test_fortran_corrections()
% getquad with opts.corrections equals getquad minus smooth_sparse_quad.

S     = geometries.sphere(1, 2, [0;0;0], 4, 1);
eps   = 1e-6;
nover = 7;
tol   = 1e-11;
toff  = struct('r', 1.05*S.r(:,1:60), 'n', S.n(:,1:60));

kerns = {kernel3d.lap3d('s'), kernel3d.lap3d('c', [1; 0.7]), ...
         kernel3d.lap3d('cp', [0.6; 1.3]), kernel3d.helm3d('d', 1.3), ...
         kernel3d.helm3d('cp', 1.3, [1; 0]), kernel3d.belpde('klb'), ...
         kernel3d.belpde('rhb', 2.1)};
targs = {S, toff};
for k = 1:numel(kerns)
    K = kerns{k};
    for it = 1:numel(targs)
        t = targs{it};
        if strcmp(K.name, 'belpde') && it == 2, continue; end
        A    = K.getquad(S, eps, t);
        rsc  = getnear(S, t);
        ws   = smooth_sparse_quad(K, t, S, rsc.row_ptr, rsc.col_ind, nover);
        Aref = A - conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, ws);
        Asub = K.getquad(S, eps, t, struct('corrections', true, 'nover', nover));
        err  = full(max(abs(Asub(:) - Aref(:)))) / max(abs(Aref(:)));
        assert(err < tol, 'corrections: %s %s, error %g', K.name, K.type, err);
    end
end

% a different oversampling order on each patch
nvec = 5 + mod(1:S.npatches, 3).';
K    = kernel3d.helm3d('c', 1.3, [1; 0.5]);
A    = K.getquad(S, eps, S);
rsc  = getnear(S, S);
ws   = smooth_sparse_quad(K, S, S, rsc.row_ptr, rsc.col_ind, nvec);
Aref = A - conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, ws);
Asub = K.getquad(S, eps, S, struct('corrections', true, 'nover', nvec));
err  = full(max(abs(Asub(:) - Aref(:)))) / max(abs(Aref(:)));
assert(err < tol, 'corrections: per patch nover, error %g', err);

% getnearquad_mex directly, rsc output
Q  = getnearquad_mex(S, 'lap3d_d', eps);
Qc = getnearquad_mex(S, 'lap3d_d', eps, [], [], struct('corrections', true, 'nover', nover));
assert(Qc.corrections && isequal(Qc.nover, nover*ones(S.npatches,1)), ...
    'corrections: Q fields');
assert(isequal(Q.row_ptr, Qc.row_ptr) && isequal(Q.iquad, Qc.iquad), ...
    'corrections: near field changed');

end


function test_getquad_sub_flags()
% getquad_sub is true only when getquad can subtract the smooth rule itself.

L = kernel3d.lap3d('d');
H = kernel3d.helm3d('s', 1.1);
F = kernel3d(@(s,t) lap3d.kern(s,t,'s'));

B = kernel3d.belpde('klb');
assert(L.getquad_sub && H.getquad_sub && B.getquad_sub, 'flags: base kernels');
assert(~F.getquad_sub, 'flags: function handle kernel');
K = kernel3d.stok3d('s');    assert(~K.getquad_sub, 'flags: stokes');
K = kernel3d.eye(1);         assert(K.getquad_sub, 'flags: eye');
K = kernel3d.zeros([1 1]);   assert(K.getquad_sub, 'flags: zeros');
% the smooth rule drops coincident nodes, so ones has a nonzero correction
K = kernel3d.ones(1,1);      assert(~K.getquad_sub, 'flags: ones');

K = L + 0.5*H;               assert(K.getquad_sub, 'flags: plus');
K = L - kernel3d.eye(0.5);   assert(K.getquad_sub, 'flags: minus eye');
K = L + F;                   assert(~K.getquad_sub, 'flags: plus with handle kernel');
K = -L;                      assert(K.getquad_sub, 'flags: uminus');
K = L/2;                     assert(K.getquad_sub, 'flags: mrdivide');
K = 2.*L;                    assert(K.getquad_sub, 'flags: times');
K = [1; 2] * L;              assert(K.getquad_sub, 'flags: left multiply');
K = L * [1 2];               assert(~K.getquad_sub, 'flags: right multiply');
K = kernel3d.interleave([L, H; H, L]); assert(K.getquad_sub, 'flags: interleave');
O = kernel3d.ones(1,1);
K = kernel3d.interleave([L, O; O, L]); assert(~K.getquad_sub, 'flags: interleave mixed');
K = kernel3d.helm3d('s_diff', [1.1; 2.0]); assert(K.getquad_sub, 'flags: helm3d diff');

end


function test_getquad_corrections()
% kernel3d.getquad_corrections gives the same answer through the Fortran
% subtraction and through the smooth_sparse_quad fallback.

S     = geometries.sphere(1, 2, [0;0;0], 4, 1);
eps   = 1e-6;
nover = 7;
tol   = 1e-11;

kerns = {kernel3d.lap3d('d') + 0.5*kernel3d.lap3d('s'), ...
         kernel3d.helm3d('c', 1.3, [1; -1i]) - kernel3d.eye(0.5), ...
         [1; 2] * kernel3d.lap3d('sp'), ...
         kernel3d.interleave([kernel3d.lap3d('s'), kernel3d.lap3d('d'); ...
                              kernel3d.lap3d('sp'), kernel3d.lap3d('s')])};
for k = 1:numel(kerns)
    K = kerns{k};
    assert(K.getquad_sub, 'getquad_corrections: case %d should use Fortran', k);
    A1 = kernel3d.getquad_corrections(K, S, eps, S, nover);
    K.getquad_sub = false;
    A2 = kernel3d.getquad_corrections(K, S, eps, S, nover);
    err = full(max(abs(A1(:) - A2(:)))) / max(abs(A2(:)));
    assert(err < tol, 'getquad_corrections: case %d, error %g', k, err);
end

% ones falls back to smooth_sparse_quad: the corrections are the weights
% of the source nodes that coincide with targets
K  = kernel3d.ones(1,1);
A  = kernel3d.getquad_corrections(K, S, eps, S, nover);
rsc = getnear(S, S);
ws = smooth_sparse_quad(K, S, S, rsc.row_ptr, rsc.col_ind, nover);
Aref = K.getquad(S, eps, S) - conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, ws);
assert(full(max(abs(A(:) - Aref(:)))) < tol, 'getquad_corrections: ones');

% kernels without built-in corrections, against getquad minus
% smooth_sparse_quad of the whole kernel
kerns = {kernel3d.lap3d('d') + kernel3d.ones(1,1), ...
         kernel3d.lap3d('s') * [1 2]};
for k = 1:numel(kerns)
    K = kerns{k};
    assert(~K.getquad_sub, 'getquad_corrections: fallback case %d flag', k);
    A  = kernel3d.getquad_corrections(K, S, eps, S, nover);
    ws = smooth_sparse_quad(K, S, S, rsc.row_ptr, rsc.col_ind, nover);
    Aref = K.getquad(S, eps, S) - conv_rsc_to_spmat(S, rsc.row_ptr, ...
        rsc.col_ind, ws, kernel3d.rsc_interleave_full(K.opdims(1), K.opdims(2)));
    err = full(max(abs(A(:) - Aref(:)))) / max(abs(Aref(:)));
    assert(err < tol, 'getquad_corrections: fallback case %d, error %g', k, err);
end

% NaN nover: no oversampling
K  = kernel3d.lap3d('d');
A1 = kernel3d.getquad_corrections(K, S, eps, S, NaN);
K.getquad_sub = false;
A2 = kernel3d.getquad_corrections(K, S, eps, S, NaN);
assert(full(max(abs(A1(:) - A2(:)))) < tol*max(abs(A2(:))), ...
    'getquad_corrections: NaN nover');

% surfermat with corrections uses the same path
K  = kernel3d.lap3d('d') + 0.5*kernel3d.lap3d('s');
B1 = surfermat({S}, K, eps, struct('corrections', true, 'nonsmoothonly', true));
K.getquad_sub = false;
B2 = surfermat({S}, K, eps, struct('corrections', true, 'nonsmoothonly', true));
assert(full(max(abs(B1(:) - B2(:)))) < tol*max(abs(B2(:))), ...
    'getquad_corrections: surfermat');

% surferkernevalmat with corrections, off-surface targets
t  = struct('r', 1.05*S.r(:,1:60));
K  = kernel3d.lap3d('d') + 0.5*kernel3d.lap3d('s');
B1 = surferkernevalmat({S}, K, t, eps, struct('corrections', true, 'nonsmoothonly', true));
K.getquad_sub = false;
B2 = surferkernevalmat({S}, K, t, eps, struct('corrections', true, 'nonsmoothonly', true));
assert(full(max(abs(B1(:) - B2(:)))) < tol*max(abs(B2(:))), ...
    'getquad_corrections: surferkernevalmat');

end


function test_kernel_order()
% Function handle kernels default to kernel_order -1, and +/- keep the
% larger order.

F = kernel3d(@(s,t) lap3d.kern(s,t,'s'));
assert(F.kernel_order == -1, 'kernel_order: default');

S = kernel3d.lap3d('s');
D = kernel3d.lap3d('d');
K = S + D;  assert(K.kernel_order == 0, 'kernel_order: S + D');
K = D + S;  assert(K.kernel_order == 0, 'kernel_order: D + S');
K = S - D;  assert(K.kernel_order == 0, 'kernel_order: S - D');

K = kernel3d.interleave([D, F; F, D]);
assert(K.kernel_order == 0, 'kernel_order: interleave with handle kernel');

end

