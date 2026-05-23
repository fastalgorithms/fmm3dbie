%TEST_KERNBYINDEX  Verify flam.kernbyindex entries match surfermat entries.
%
% Scenario 1: single surfer, no novers passed (uses get_overs_orders internally)
% Scenario 2: two-surfer array, same
% Scenario 3: explicit novers vector, compared against manually built reference

run ../startup.m

tol = 1e-10;
rng(42);

fprintf('=== test_kernbyindex ===\n\n');

opts_sm = struct('selfquad', false, 'adaptive_correction', false);

%% Scenario 1: single surfer
fprintf('Scenario 1: single surfer, lap combined ...\n');

S = geometries.ellipsoid([1,1,1], [3,3,3], [], 6);
kern = kernel3d('l', 'c', [1,1]);
opdims = kern.opdims;

Asmth = surfermat(S, kern, opts_sm);

n = S.npts * opdims(1);
ii = randi(n, 50, 1);
jj = randi(n, 50, 1);

Acheck = flam.kernbyindex(ii, jj, S, kern, opdims);

err1 = norm(Acheck - Asmth(ii, jj), 'fro') / (norm(Asmth(ii, jj), 'fro') + eps);
fprintf('  relative error: %.2e\n', err1);
assert(err1 < tol, 'FAIL: %.2e', err1);
fprintf('  PASS\n\n');

%% Scenario 2: two-surfer array
fprintf('Scenario 2: two-surfer array, lap single-layer ...\n');

S2 = geometries.ellipsoid([1,1,1], [3,3,3], [4;0;0], 6);
srfrs = [S, S2];
kern_sl = kernel3d('l', 's');
opdims_sl = kern_sl.opdims;

Asmth2 = surfermat(srfrs, kern_sl, opts_sm);
n2 = size(Asmth2, 1);
m2 = size(Asmth2, 2);
ii2 = randi(n2, 80, 1);
jj2 = randi(m2, 80, 1);

Acheck2 = flam.kernbyindex(ii2, jj2, srfrs, kern_sl, opdims_sl);

err2 = norm(Acheck2 - Asmth2(ii2, jj2), 'fro') / (norm(Asmth2(ii2, jj2), 'fro') + eps);
fprintf('  relative error: %.2e\n', err2);
assert(err2 < tol, 'FAIL: %.2e', err2);
fprintf('  PASS\n\n');

%% Scenario 3: explicit novers vector (one order per surfer)
fprintf('Scenario 3: explicit novers, oversampled smooth rule ...\n');

S3 = geometries.ellipsoid([1,1,1], [3,3,3], [], 4);
kern3 = kernel3d('l', 's');
opdims3 = kern3.opdims;

% novers as vector of length nsurfers=1
nover3 = S3.norders(1) + 2;

% Build reference matching what kernbyindex does:
% oversample uniformly, evaluate kern, zero self, apply xinterp
[S3over, xinterp3] = S3.oversample(nover3 * ones(S3.npatches, 1));
wtsover3 = repmat(S3over.wts(:).', opdims3(2), 1);
wtsover3 = wtsover3(:).';

Aover_ref = kern3.eval(S3over, S3) .* wtsover3;

src_norm = max(vecnorm(S3over.r), 1);
for q = 1:size(S3over.r, 2)
    diff_norms = vecnorm(S3.r - S3over.r(:,q));
    self_mask = diff_norms < 1e-14 * src_norm(q);
    if any(self_mask)
        row_off = (1:opdims3(1)) + opdims3(1)*(find(self_mask).' - 1);
        col_off = (1:opdims3(2)) + opdims3(2)*(q-1);
        Aover_ref(row_off(:), col_off) = 0;
    end
end

Aref3 = Aover_ref * kron(xinterp3, eye(opdims3(2)));

n3 = S3.npts;
ii3 = randi(n3, 60, 1);
jj3 = randi(n3, 60, 1);

Acheck3 = flam.kernbyindex(ii3, jj3, S3, kern3, opdims3, nover3);

err3 = norm(Acheck3 - Aref3(ii3, jj3), 'fro') / (norm(Aref3(ii3, jj3), 'fro') + eps);
fprintf('  relative error: %.2e\n', err3);
assert(err3 < tol, 'FAIL: %.2e', err3);
fprintf('  PASS\n\n');

%% Scenario 4: wrong novers should give wrong answer
fprintf('Scenario 4: wrong novers gives wrong answer ...\n');

% Use norders(1) (no oversampling) when the correct answer needs more
nover_wrong = S3.norders(1);
Acheck_wrong = flam.kernbyindex(ii3, jj3, S3, kern3, opdims3, nover_wrong);

err_wrong = norm(Acheck_wrong - Aref3(ii3, jj3), 'fro') / (norm(Aref3(ii3, jj3), 'fro') + eps);
fprintf('  relative error with wrong novers: %.2e\n', err_wrong);
assert(err_wrong > tol, 'FAIL: wrong novers should give error > tol but got %.2e', err_wrong);
fprintf('  PASS (correctly wrong)\n\n');

fprintf('All tests passed.\n');
