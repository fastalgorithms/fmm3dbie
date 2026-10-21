% Test the near corrections to the oversampled smooth rule:
% smooth_sparse_quad, propagation of getquad_sub through kernel algebra,
% kernel3d.getquad_corrections, and the Fortran subtraction for radpp.

%% Now run the tests

test_smooth_sparse_quad();
test_getquad_sub_flags();
test_getquad_corrections();
test_kernel_order();
test_radpp_corrections();


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


function test_getquad_sub_flags()
% getquad_sub is true only when getquad can subtract the smooth rule itself.

L = kernel3d.lap3d('d');
F = kernel3d(@(s,t) lap3d.kern(s,t,'s'));
E = kernel3d.eye(1);
Z = kernel3d.zeros([1 1]);

assert(~L.getquad_sub && ~F.getquad_sub, 'flags: default false');
assert(E.getquad_sub, 'flags: eye');
assert(Z.getquad_sub, 'flags: zeros');
% the smooth rule drops coincident nodes, so ones has a nonzero correction
K = kernel3d.ones(1,1);      assert(~K.getquad_sub, 'flags: ones');

K = E + Z;                   assert(K.getquad_sub, 'flags: plus');
K = E - Z;                   assert(K.getquad_sub, 'flags: minus');
K = L - kernel3d.eye(0.5);   assert(~K.getquad_sub, 'flags: minus mixed');
K = -E;                      assert(K.getquad_sub, 'flags: uminus');
K = E/2;                     assert(K.getquad_sub, 'flags: mrdivide');
K = 2.*E;                    assert(K.getquad_sub, 'flags: times');
K = [1; 2] * E;              assert(K.getquad_sub, 'flags: left multiply');
K = E * [1 2];               assert(~K.getquad_sub, 'flags: right multiply');
K = kernel3d.interleave([E, Z; Z, E]); assert(K.getquad_sub, 'flags: interleave');
K = kernel3d.interleave([E, L; L, E]); assert(~K.getquad_sub, 'flags: interleave mixed');

end


function test_getquad_corrections()
% kernel3d.getquad_corrections against getquad minus smooth_sparse_quad.

S     = geometries.sphere(1, 2, [0;0;0], 4, 1);
eps   = 1e-6;
nover = 7;
tol   = 1e-11;
rsc   = getnear(S, S);
toff  = struct('r', 1.05*S.r(:,1:60), 'n', S.n(:,1:60));

kerns = {kernel3d.lap3d('d') + 0.5*kernel3d.lap3d('s'), ...
         kernel3d.helm3d('c', 1.3, [1; -1i]) - kernel3d.eye(0.5), ...
         [1; 2] * kernel3d.lap3d('sp'), ...
         kernel3d.lap3d('d') + kernel3d.ones(1,1), ...
         kernel3d.lap3d('s') * [1 2], ...
         kernel3d.ones(1,1)};
for k = 1:numel(kerns)
    K = kerns{k};
    assert(~K.getquad_sub, 'getquad_corrections: case %d flag', k);
    A  = kernel3d.getquad_corrections(K, S, eps, S, nover);
    ws = smooth_sparse_quad(K, S, S, rsc.row_ptr, rsc.col_ind, nover);
    Aref = K.getquad(S, eps, S) - conv_rsc_to_spmat(S, rsc.row_ptr, ...
        rsc.col_ind, ws, kernel3d.rsc_interleave_full(K.opdims(1), K.opdims(2)));
    err = full(max(abs(A(:) - Aref(:)))) / max(abs(Aref(:)));
    assert(err < tol, 'getquad_corrections: case %d, error %g', k, err);
end

% kernels whose getquad already gives the corrections (eye, zeros)
K = kernel3d.eye(0.5);
A = kernel3d.getquad_corrections(K, S, eps, S, nover);
B = K.getquad(S, eps, S);
assert(isequal(size(A), size(B)) && full(max(abs(A(:) - B(:)))) == 0, ...
    'getquad_corrections: eye');

% NaN nover: no oversampling
K  = kernel3d.lap3d('d');
A  = kernel3d.getquad_corrections(K, S, eps, toff, NaN);
rt = getnear(S, toff);
ws = smooth_sparse_quad(K, toff, S, rt.row_ptr, rt.col_ind, NaN);
Aref = K.getquad(S, eps, toff) - conv_rsc_to_spmat(S, rt.row_ptr, rt.col_ind, ws);
err = full(max(abs(A(:) - Aref(:)))) / max(abs(Aref(:)));
assert(err < tol, 'getquad_corrections: NaN nover, error %g', err);

% surfermat and surferkernevalmat with corrections run through the same path
K  = kernel3d.lap3d('d') + 0.5*kernel3d.lap3d('s');
B  = surfermat({S}, K, eps, struct('corrections', true, 'nonsmoothonly', true));
assert(all(isfinite(nonzeros(B))), 'getquad_corrections: surfermat');
B  = surferkernevalmat({S}, K, toff, eps, struct('corrections', true, 'nonsmoothonly', true));
assert(all(isfinite(nonzeros(B))), 'getquad_corrections: surferkernevalmat');

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


function test_radpp_corrections()
% radpp getquad with opts.corrections (getnearquad_radpp_sub) equals
% getquad minus smooth_sparse_quad, and getquad_corrections agrees with
% its fallback.

zk    = 1.1;
eps   = 1e-7;
nover = 7;
tol   = 1e-11;
S = geometries.sphere(1, 1, [0;0;0], 4, 1);
toff.r = randn(3, 40);  toff.r = toff.r./vecnorm(toff.r)*1.05;
toff.n = toff.r./vecnorm(toff.r);
nvec = 5 + mod(1:S.npatches, 3).';

g  = @(s,t) helm3d.green(zk, s, t);
fk = kernel3d.helm3d('s', zk);
kerns = {kernel3d.radpp(fk, [1e-9 4], 's', 'opts', struct('eps', 1e-13)), ...
         kernel3d(g, 'cp'), ...
         kernel3d(g, 's') + 2*kernel3d(g, 'd'), ...
         kernel3d(g, 'trans_sys')};
for k = 1:numel(kerns)
    K = kerns{k};
    assert(K.getquad_sub, 'radpp corrections: case %d flag', k);
    m = K.opdims(1); n = K.opdims(2);
    ri = kernel3d.rsc_interleave_full(m, n);
    for t = {S, toff}
        for nv = {nover, nvec}
            rsc  = getnear(S, t{1});
            A    = K.getquad(S, eps, t{1});
            ws   = smooth_sparse_quad(K, t{1}, S, rsc.row_ptr, rsc.col_ind, nv{1});
            Aref = A - conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, ws, ri);
            Asub = K.getquad(S, eps, t{1}, ...
                struct('corrections', true, 'nover', nv{1}));
            err  = full(max(abs(Asub(:) - Aref(:)))) / max(abs(Aref(:)));
            assert(err < tol, 'radpp corrections: case %d, error %g', k, err);
        end
    end

    A1 = kernel3d.getquad_corrections(K, S, eps, S, nover);
    K.getquad_sub = false;
    A2 = kernel3d.getquad_corrections(K, S, eps, S, nover);
    err = full(max(abs(A1(:) - A2(:)))) / max(abs(A2(:)));
    assert(err < tol, 'radpp getquad_corrections: case %d, error %g', k, err);
end

end
