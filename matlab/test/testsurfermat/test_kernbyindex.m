% Test byindex.kernbyindex entries against surfermat.

run ../../startup.m

rng(42);

%% Now run the tests

test_two_surfer_lap();
test_two_surfer_lap_ordervec_mode();
test_recompute_on_the_fly();
test_explicit_novers();
test_stok_with_corrections();
test_single_surfer_stok();
test_extra_targ_field();
test_wrong_novers();


function test_two_surfer_lap()
% Two-surfer 2x2 Laplace kernel array: kernbyindex matches surfermat.

eps = 1e-10;  tol = 1e-9;
opts_sm = struct('selfquad', false, 'adaptive_correction', false);

S1a = geometries.ellipsoid([1,1,1], [3,3,3], [],      6);
S1b = geometries.ellipsoid([1,1,1], [3,3,3], [4;0;0], 6);
srfrs = [S1a, S1b];

kerns(2,2) = kernel3d();
kerns(1,1) = kernel3d('l', 's');
kerns(2,1) = kernel3d('l', 'sp');
kerns(1,2) = kernel3d('l', 'd');
kerns(2,2) = kernel3d('l', 'dp');

[Asmth, novers] = surfermat(srfrs, kerns, eps, opts_sm);
ii = randi(size(Asmth,1), 80, 1);
jj = randi(size(Asmth,2), 80, 1);

err = norm(byindex.kernbyindex(ii, jj, srfrs, kerns, eps, novers) - Asmth(ii,jj), 'fro') ...
    / (norm(Asmth(ii,jj), 'fro') + eps);
assert(err < tol, 'two-surfer lap: %.2e', err);

end


function test_two_surfer_lap_ordervec_mode()
% Same two-surfer setup as test_two_surfer_lap, but forcing surfermat to
% return the numeric norders (opts.ifreturnovers=0).

eps = 1e-10;  tol = 1e-9;
opts_sm = struct('selfquad', false, 'adaptive_correction', false, ...
    'ifreturnovers', 0);

S1a = geometries.ellipsoid([1,1,1], [3,3,3], [],      6);
S1b = geometries.ellipsoid([1,1,1], [3,3,3], [4;0;0], 6);
srfrs = [S1a, S1b];

kerns(2,2) = kernel3d();
kerns(1,1) = kernel3d('l', 's');
kerns(2,1) = kernel3d('l', 'sp');
kerns(1,2) = kernel3d('l', 'd');
kerns(2,2) = kernel3d('l', 'dp');

[Asmth, novers] = surfermat(srfrs, kerns, eps, opts_sm);

% Confirm we actually got the order-vector form.
assert(iscell(novers) && isnumeric(novers{1,1}), ...
    'expected order-vector form of objover');
has_oversamp = false;
for k = 1:numel(novers)
    if ~any(isnan(novers{k})) && ~isempty(novers{k})
        has_oversamp = has_oversamp || any(novers{k} ~= 0);
    end
end

ii = randi(size(Asmth,1), 80, 1);
jj = randi(size(Asmth,2), 80, 1);

err = norm(byindex.kernbyindex(ii, jj, srfrs, kerns, eps, novers) - Asmth(ii,jj), 'fro') ...
    / (norm(Asmth(ii,jj), 'fro') + eps);
assert(err < tol, 'two-surfer lap (order-vector mode): %.2e', err);

end


function test_recompute_on_the_fly()
% Passing objover=[] (or omitting it) makes kernbyindex recompute
% oversampling orders on the fly.

eps = 1e-10;  tol = 1e-9;
opts_sm = struct('selfquad', false, 'adaptive_correction', false);

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

% objover omitted entirely -> defaults to [] -> on-the-fly recomputation
err = norm(byindex.kernbyindex(ii, jj, srfrs, kerns, eps) - Asmth(ii,jj), 'fro') ...
    / (norm(Asmth(ii,jj), 'fro') + eps);
assert(err < tol, 'recompute on the fly (objover omitted): %.2e', err);

% objover explicitly [] -> same on-the-fly recomputation
err2 = norm(byindex.kernbyindex(ii, jj, srfrs, kerns, eps, []) - Asmth(ii,jj), 'fro') ...
    / (norm(Asmth(ii,jj), 'fro') + eps);
assert(err2 < tol, 'recompute on the fly (objover=[]): %.2e', err2);

end


function test_explicit_novers()
% Explicit novers vector: kernbyindex matches manually oversampled smooth matrix.

eps = 1e-10;  tol = 1e-9;

S      = geometries.ellipsoid([1,1,1], [3,3,3], [], 4);
kern   = kernel3d('l', 's');
opdims = kern.opdims;

nover = (S.norders(1) + 2) * ones(S.npatches, 1);
[Sover, xinterp] = S.oversample(nover);
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
Aref = Aover_ref * kron(xinterp, eye(opdims(2)));

ii = randi(S.npts, 60, 1);
jj = randi(S.npts, 60, 1);

err = norm(byindex.kernbyindex(ii, jj, S, kern, eps, {nover}) - Aref(ii,jj), 'fro') ...
    / (norm(Aref(ii,jj), 'fro') + eps);
assert(err < tol, 'explicit novers: %.2e', err);

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

[Afull, novers] = surfermat(srfrs, kerns, eps);
[Qns,   ~]      = surfermat(srfrs, kerns, eps, struct('nonsmoothonly', 1));

ii = randi(size(Afull,1), 100, 1);
jj = randi(size(Afull,2), 100, 1);

err = norm(byindex.kernbyindex(ii, jj, srfrs, kerns, eps, novers, Qns) - Afull(ii,jj), 'fro') ...
    / (norm(Afull(ii,jj), 'fro') + eps);
assert(err < tol, 'stok with corrections: %.2e', err);

end


function test_single_surfer_stok()
% Single surfer, Laplace SLP: kernbyindex matches surfermat.

eps = 1e-10;  tol = 1e-9;
opts_sm = struct('selfquad', false, 'adaptive_correction', false);

S    = geometries.ellipsoid([1,1,1], [3,3,3], [], 5);
kern = kernel3d('l', 's');
od   = kern.opdims;

[Asmth, novers] = surfermat(S, kern, eps, opts_sm);
ii = randi(S.npts * od(1), 80, 1);
jj = randi(S.npts * od(2), 80, 1);

err = norm(byindex.kernbyindex(ii, jj, S, kern, eps, novers) - Asmth(ii,jj), 'fro') ...
    / (norm(Asmth(ii,jj), 'fro') + eps);
assert(err < tol, 'single-surfer stok: %.2e', err);

end


function test_extra_targ_field()
% Kernel declaring a target field beyond the geometry (belpde reads
% mean_curv, which the surfer stores as an npts x 1 column).

eps = 1e-10;  tol = 1e-9;
opts_sm = struct('selfquad', false, 'adaptive_correction', false);

S1 = geometries.ellipsoid([1,1,1.2], [2,2,2], [],      5);
S2 = geometries.ellipsoid([1,1,1],   [2,2,2], [4;0;0], 5);
srfrs = [S1, S2];

kern = kernel3d('bel', 'rlb');

[Asmth, novers] = surfermat(srfrs, kern, eps, opts_sm);
ii = randi(size(Asmth,1), 80, 1);
jj = randi(size(Asmth,2), 80, 1);

err = norm(byindex.kernbyindex(ii, jj, srfrs, kern, eps, novers) - Asmth(ii,jj), 'fro') ...
    / (norm(Asmth(ii,jj), 'fro') + eps);
assert(err < tol, 'extra targ field (objover): %.2e', err);

[Asmth2, novers2] = surfermat(srfrs, kern, eps, ...
    struct('selfquad', false, 'adaptive_correction', false, 'ifreturnovers', 0));
err2 = norm(byindex.kernbyindex(ii, jj, srfrs, kern, eps, novers2) - Asmth2(ii,jj), 'fro') ...
    / (norm(Asmth2(ii,jj), 'fro') + eps);
assert(err2 < tol, 'extra targ field (order-vector mode): %.2e', err2);

end

function test_wrong_novers()
% Wrong novers should give an error above tolerance (negative test).
% Reuses the same small geometry as test_explicit_novers.

eps = 1e-10;  tol = 1e-9;

S    = geometries.ellipsoid([1,1,1], [3,3,3], [], 4);
kern = kernel3d('l', 's');
od   = kern.opdims;

nover_good  = (S.norders(1) + 2) * ones(S.npatches, 1);
nover_wrong = S.norders(1)       * ones(S.npatches, 1);

[Sover, xinterp] = S.oversample(nover_good);
wtsover = repmat(Sover.wts(:).', od(2), 1);  wtsover = wtsover(:).';
Aover_ref = kern.eval(Sover, S) .* wtsover;
Aref = Aover_ref * kron(xinterp, eye(od(2)));

ii = randi(S.npts, 60, 1);
jj = randi(S.npts, 60, 1);

err = norm(byindex.kernbyindex(ii, jj, S, kern, eps, {nover_wrong}) - Aref(ii,jj), 'fro') ...
    / (norm(Aref(ii,jj), 'fro') + eps);
assert(err > tol, 'wrong novers should give error > tol but got %.2e', err);

end
