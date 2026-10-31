% Test byindex.kernbyindex entries against surfermat.

run ../startup.m

rng(42);

%% Now run the tests

test_mixed_opdims();
test_stok_with_corrections();
test_additive_corrections();
test_diag_term();
test_extra_targ_field();
test_shared_nodes();
test_oversampled_sources();
test_old_signature_rejected();


function test_mixed_opdims()
% Two surfers with a 2x2 kernel array of opdims [3 3], [3 1], [1 3] and
% [1 1], including a zero block: kernbyindex matches surfermat.

eps = 1e-8;  tol = 1e-9;
opts_sm = struct('selfquad', false, 'adaptive_correction', false, 'ifoversamp', 0);

S1 = geometries.sphere(1, 2, [0;0;0], 3, 1);
S2 = geometries.sphere(1, 2, [3;0;0], 3, 1);
srfrs = [S1, S2];

kerns = mixed_kerns();

Asmth = surfermat(srfrs, kerns, eps, opts_sm);
ii = randi(size(Asmth,1), 100, 1);
jj = randi(size(Asmth,2), 100, 1);

err = norm(byindex.kernbyindex(ii, jj, srfrs, kerns) - Asmth(ii,jj), 'fro') ...
    / norm(Asmth(ii,jj), 'fro');
assert(err < tol, 'mixed opdims: %.2e', err);

end


function test_stok_with_corrections()
% Two-surfer 2x2 Stokes kernels: kernbyindex with nonsmoothonly corrections
% matches the full surfermat reference.

eps = 1e-8;  tol = 1e-9;

S3a = geometries.sphere(1, 2, [0;0;0],   3, 1);
S3b = geometries.sphere(1, 2, [2.5;0;0], 3, 1);
srfrs = [S3a, S3b];

kerns(2,2) = kernel3d();
kerns(1,1) = kernel3d('stok', 's');
kerns(2,1) = kernel3d('stok', 'sp');
kerns(1,2) = kernel3d('stok', 'd');
kerns(2,2) = kernel3d('stok', 's');

Afull = surfermat(srfrs, kerns, eps, struct('ifoversamp', 0));
Qns   = surfermat(srfrs, kerns, eps, struct('nonsmoothonly', 1, 'ifoversamp', 0));

ii = randi(size(Afull,1), 100, 1);
jj = randi(size(Afull,2), 100, 1);

err = norm(byindex.kernbyindex(ii, jj, srfrs, kerns, Qns) - Afull(ii,jj), 'fro') ...
    / norm(Afull(ii,jj), 'fro');
assert(err < tol, 'stok with corrections: %.2e', err);

end


function test_additive_corrections()
% corrections=1 output from surfermat, added with opts.replace_quadcorr = false.

eps = 1e-8;  tol = 1e-9;

S    = geometries.ellipsoid([1,1,1], [3,3,3], [], 4);
kern = kernel3d('l', 'd');

Afull = surfermat(S, kern, eps, struct('ifoversamp', 0));
Qc    = surfermat(S, kern, eps, struct('corrections', 1, 'ifoversamp', 0));

ii = randi(size(Afull,1), 80, 1);
jj = randi(size(Afull,2), 80, 1);

Atest = byindex.kernbyindex(ii, jj, S, kern, Qc, struct('replace_quadcorr', false));
err = norm(Atest - Afull(ii,jj), 'fro') / norm(Afull(ii,jj), 'fro');
assert(err < tol, 'additive corrections: %.2e', err);

end


function test_diag_term()
% Identity term with opdims [3 3] comes from Qsparse in both modes,
% including nonsmoothonly output without quadrature.

eps = 1e-8;  tol = 1e-9;

S    = geometries.sphere(1, 2, [0;0;0], 3, 1);
kern = kernel3d('stok', 'd') + kernel3d.eye(-0.5*eye(3));

ii = randi(3*S.npts, 60, 1);
jj = [randi(3*S.npts, 60, 1); ii];
ii = [ii; ii];                     % include diagonal entries

opts_sm = struct('selfquad', false, 'adaptive_correction', false, 'ifoversamp', 0);
Asmth = surfermat(S, kern, eps, opts_sm);
opts_sm.nonsmoothonly = true;
Qd = surfermat(S, kern, eps, opts_sm);
assert(nnz(Qd) == 3*S.npts, 'expected only the identity term in Qd');
err = norm(byindex.kernbyindex(ii, jj, S, kern, Qd) - Asmth(ii,jj), 'fro') ...
    / norm(Asmth(ii,jj), 'fro');
assert(err < tol, 'diag term, identity only: %.2e', err);

Afull = surfermat(S, kern, eps, struct('ifoversamp', 0));
Qns   = surfermat(S, kern, eps, struct('nonsmoothonly', 1, 'ifoversamp', 0));
Qc    = surfermat(S, kern, eps, struct('corrections', 1, 'ifoversamp', 0));

err = norm(byindex.kernbyindex(ii, jj, S, kern, Qns) - Afull(ii,jj), 'fro') ...
    / norm(Afull(ii,jj), 'fro');
assert(err < tol, 'diag term, replacement: %.2e', err);

Atest = byindex.kernbyindex(ii, jj, S, kern, Qc, struct('replace_quadcorr', false));
err = norm(Atest - Afull(ii,jj), 'fro') / norm(Afull(ii,jj), 'fro');
assert(err < tol, 'diag term, additive: %.2e', err);

end


function test_extra_targ_field()
% Kernel declaring a target field beyond the geometry (belpde reads
% mean_curv, which the surfer stores as an npts x 1 column). Its identity
% term comes from Qsparse.

eps = 1e-8;  tol = 1e-9;
opts_sm = struct('selfquad', false, 'adaptive_correction', false, 'ifoversamp', 0);

S1 = geometries.ellipsoid([1,1,1.2], [2,2,2], [],      4);
S2 = geometries.ellipsoid([1,1,1],   [2,2,2], [4;0;0], 4);
srfrs = [S1, S2];

kern = kernel3d('bel', 'rlb');

Asmth = surfermat(srfrs, kern, eps, opts_sm);
opts_sm.nonsmoothonly = true;
Qd = surfermat(srfrs, kern, eps, opts_sm);

ii = randi(size(Asmth,1), 80, 1);
jj = [randi(size(Asmth,2), 80, 1); ii];
ii = [ii; ii];

err = norm(byindex.kernbyindex(ii, jj, srfrs, kern, Qd) - Asmth(ii,jj), 'fro') ...
    / norm(Asmth(ii,jj), 'fro');
assert(err < tol, 'extra targ field: %.2e', err);

end


function test_shared_nodes()
% Two copies of the same surfer: nodes shared between surfers give zero
% entries rather than Inf/NaN.

S    = geometries.sphere(1, 2, [0;0;0], 3, 1);
kern = kernel3d('l', 'd');

ii = randi(S.npts, 60, 1);
jj = [randi(S.npts, 60, 1); ii];

A11 = byindex.kernbyindex(ii, jj, [S, S], kern);
A12 = byindex.kernbyindex(ii, jj + S.npts, [S, S], kern);
assert(all(isfinite(A12(:))), 'shared nodes should give finite entries');
assert(norm(A12 - A11, 'fro') <= 1e-14*norm(A11, 'fro'), ...
    'shared nodes should match the self block');

end


function test_oversampled_sources()
% Oversampled on-surface matrix built from kernbyindexeval with the
% oversampled surfer as sources and the original surfer as targets,
% applying the interpolation matrix outside. opdims = [3 3].

tol = 1e-12;

S      = geometries.sphere(1, 2, [0;0;0], 3, 1);
kern   = kernel3d('stok', 's');
opdims = kern.opdims;

[Sover, xinterp] = S.oversample((S.norders(1) + 2) * ones(S.npatches, 1));
xinterp = kron(xinterp, eye(opdims(2)));
wtsover = repmat(Sover.wts(:).', opdims(2), 1);

Aover_ref = kern.eval(Sover, S) .* wtsover(:).';
Aover_ref(~isfinite(Aover_ref)) = 0;
Aref = Aover_ref * xinterp;

ii = randi(3*S.npts, 60, 1);
jj = randi(3*S.npts, 60, 1);

Aover = byindex.kernbyindexeval(ii, 1:3*Sover.npts, Sover, kern, S);
err = norm(Aover * xinterp(:, jj) - Aref(ii,jj), 'fro') / norm(Aref(ii,jj), 'fro');
assert(err < tol, 'oversampled sources: %.2e', err);

end


function test_old_signature_rejected()
% The old (..., eps, objover) arguments must be rejected rather than
% being taken as Qsparse.

S    = geometries.sphere(1, 1, [0;0;0], 3, 1);
kern = kernel3d('l', 's');

caught = false;
try
    byindex.kernbyindex(1:3, 1:3, S, kern, 1e-10, {NaN*zeros(S.npatches,1)});
catch
    caught = true;
end
assert(caught, 'eps in place of Qsparse should raise an error');

end


function kerns = mixed_kerns()
% Helper: 2x2 kernel array with opdims [3 3], [3 1]; [1 3], [1 1]
kerns(2,2) = kernel3d();
kerns(1,1) = kernel3d('stok', 's');
kerns(1,2) = grad(kernel3d('l', 's'));
kerns(2,1) = kernel3d.zeros([1 3]);
kerns(2,2) = kernel3d('l', 's');
end
