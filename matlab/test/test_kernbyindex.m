% Test byindex.kernbyindex entries against surfermat.

run ../startup.m

rng(42);

%% Now run the tests

test_two_surfer_lap();
test_stok_with_corrections();
test_additive_corrections();
test_single_surfer_lap();
test_extra_targ_field();
test_oversampled_sources();
test_old_signature_rejected();


function test_two_surfer_lap()
% Two-surfer 2x2 Laplace kernel array: kernbyindex matches surfermat.

eps = 1e-10;  tol = 1e-9;
opts_sm = struct('selfquad', false, 'adaptive_correction', false, 'ifoversamp', 0);

S1a = geometries.ellipsoid([1,1,1], [3,3,3], [],      6);
S1b = geometries.ellipsoid([1,1,1], [3,3,3], [4;0;0], 6);
srfrs = [S1a, S1b];

kerns(2,2) = kernel3d();
kerns(1,1) = kernel3d('l', 's');
kerns(2,1) = kernel3d('l', 'sp');
kerns(1,2) = kernel3d('l', 'd');
kerns(2,2) = kernel3d('l', 'dp');

Asmth = surfermat(srfrs, kerns, eps, opts_sm);
ii = randi(size(Asmth,1), 80, 1);
jj = randi(size(Asmth,2), 80, 1);

err = norm(byindex.kernbyindex(ii, jj, srfrs, kerns) - Asmth(ii,jj), 'fro') ...
    / (norm(Asmth(ii,jj), 'fro') + eps);
assert(err < tol, 'two-surfer lap: %.2e', err);

end


function test_stok_with_corrections()
% Two-surfer 2x2 Stokes kernels: kernbyindex with nonsmoothonly corrections
% matches the full surfermat reference.

eps = 1e-10;  tol = 1e-9;

S3a = geometries.ellipsoid([1,1,1], [1,1,1], [],        5);
S3b = geometries.ellipsoid([1,1,1], [1,1,1], [2.5;0;0], 5);
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
    / (norm(Afull(ii,jj), 'fro') + eps);
assert(err < tol, 'stok with corrections: %.2e', err);

end


function test_additive_corrections()
% corrections=1 output from surfermat, added with opts.replace_quadcorr = false.

eps = 1e-10;  tol = 1e-9;

S    = geometries.ellipsoid([1,1,1], [3,3,3], [], 5);
kern = kernel3d('l', 'd');

Afull = surfermat(S, kern, eps, struct('ifoversamp', 0));
Qc    = surfermat(S, kern, eps, struct('corrections', 1, 'ifoversamp', 0));

ii = randi(size(Afull,1), 80, 1);
jj = randi(size(Afull,2), 80, 1);

Atest = byindex.kernbyindex(ii, jj, S, kern, Qc, struct('replace_quadcorr', false));
err = norm(Atest - Afull(ii,jj), 'fro') / (norm(Afull(ii,jj), 'fro') + eps);
assert(err < tol, 'additive corrections: %.2e', err);

end


function test_single_surfer_lap()
% Single surfer, Laplace SLP: kernbyindex matches surfermat.

eps = 1e-10;  tol = 1e-9;
opts_sm = struct('selfquad', false, 'adaptive_correction', false, 'ifoversamp', 0);

S    = geometries.ellipsoid([1,1,1], [3,3,3], [], 5);
kern = kernel3d('l', 's');
od   = kern.opdims;

Asmth = surfermat(S, kern, eps, opts_sm);
ii = randi(S.npts * od(1), 80, 1);
jj = randi(S.npts * od(2), 80, 1);

err = norm(byindex.kernbyindex(ii, jj, S, kern) - Asmth(ii,jj), 'fro') ...
    / (norm(Asmth(ii,jj), 'fro') + eps);
assert(err < tol, 'single-surfer lap: %.2e', err);

end


function test_extra_targ_field()
% Kernel declaring a target field beyond the geometry (belpde reads
% mean_curv, which the surfer stores as an npts x 1 column).

eps = 1e-10;  tol = 1e-9;
opts_sm = struct('selfquad', false, 'adaptive_correction', false, 'ifoversamp', 0);

S1 = geometries.ellipsoid([1,1,1.2], [2,2,2], [],      5);
S2 = geometries.ellipsoid([1,1,1],   [2,2,2], [4;0;0], 5);
srfrs = [S1, S2];

kern = kernel3d('bel', 'rlb');
kern.diag = [];   % compare smooth rule only; identity term is tested in test_kernel3d_eye

Asmth = surfermat(srfrs, kern, eps, opts_sm);
ii = randi(size(Asmth,1), 80, 1);
jj = randi(size(Asmth,2), 80, 1);

err = norm(byindex.kernbyindex(ii, jj, srfrs, kern) - Asmth(ii,jj), 'fro') ...
    / (norm(Asmth(ii,jj), 'fro') + eps);
assert(err < tol, 'extra targ field: %.2e', err);

end


function test_oversampled_sources()
% Oversampled on-surface matrix built from kernbyindexeval with the
% oversampled surfer as sources and the original surfer as targets,
% applying the interpolation matrix outside.

eps = 1e-10;  tol = 1e-9;

S      = geometries.ellipsoid([1,1,1], [3,3,3], [], 4);
kern   = kernel3d('l', 's');
opdims = kern.opdims;

nover = (S.norders(1) + 2) * ones(S.npatches, 1);
[Sover, xinterp] = S.oversample(nover);
xinterp = kron(xinterp, eye(opdims(2)));
wtsover = repmat(Sover.wts(:).', opdims(2), 1);
wtsover = wtsover(:).';

Aover_ref = kern.eval(Sover, S) .* wtsover;
src_norm  = max(vecnorm(Sover.r), 1);
for q = 1:size(Sover.r, 2)
    diff_norms = vecnorm(S.r - Sover.r(:,q));
    self_mask  = diff_norms < 1e-14 * src_norm(q);
    if any(self_mask)
        row_off = (1:opdims(1)) + opdims(1)*(find(self_mask).' - 1);
        col_off = (1:opdims(2)) + opdims(2)*(q-1);
        Aover_ref(row_off(:), col_off) = 0;
    end
end
Aref = Aover_ref * xinterp;

ii = randi(S.npts, 60, 1);
jj = randi(S.npts, 60, 1);

Aover = byindex.kernbyindexeval(ii, 1:Sover.npts*opdims(2), Sover, kern, S);
err = norm(Aover * xinterp(:, jj) - Aref(ii,jj), 'fro') / (norm(Aref(ii,jj), 'fro') + eps);
assert(err < tol, 'oversampled sources: %.2e', err);

end


function test_old_signature_rejected()
% The old (..., eps, objover) arguments must be rejected rather than
% being taken as Qsparse.

S    = geometries.ellipsoid([1,1,1], [3,3,3], [], 4);
kern = kernel3d('l', 's');

caught = false;
try
    byindex.kernbyindex(1:3, 1:3, S, kern, 1e-10, {NaN*zeros(S.npatches,1)});
catch
    caught = true;
end
assert(caught, 'eps in place of Qsparse should raise an error');

end
