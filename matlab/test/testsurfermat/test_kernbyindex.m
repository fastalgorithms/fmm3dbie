%TEST_KERNBYINDEX  Verify byindex.kernbyindex entries match surfermat entries.
%
% Scenario 1: two-surfer array, 2x2 lap kernel array (scalar, 1x1 opdims)
% Scenario 2: explicit novers vector, compared against manually built reference
% Scenario 3: two-surfer array, 2x2 stok kernel array (vector, 3x3 opdims)
%             with nonsmoothonly corrections passed in as Qsparse
% Scenario 4: single surfer, stok s (vector-valued, 3x3 opdims)
% Scenario 5: wrong novers should give wrong answer (negative test)

run ../../startup.m

eps = 1e-10;
tol = 1e-9;
rng(42);

fprintf('=== test_kernbyindex ===\n\n');

opts_sm = struct('selfquad', false, 'adaptive_correction', false);

%% Scenario 1: two-surfer array, 2x2 lap kernel array
fprintf('Scenario 1: two-surfer array, 2x2 lap kernels ...\n');

S1a = geometries.ellipsoid([1,1,1], [3,3,3], [],      6);
S1b = geometries.ellipsoid([1,1,1], [3,3,3], [4;0;0], 6);
srfrs1 = [S1a, S1b];

kerns1(2,2) = kernel3d();
kerns1(1,1) = kernel3d('l', 's');
kerns1(2,1) = kernel3d('l', 'sp');
kerns1(1,2) = kernel3d('l', 'd');
kerns1(2,2) = kernel3d('l', 'dp');

[Asmth1, novers1] = surfermat(srfrs1, kerns1, eps, opts_sm);
n1 = size(Asmth1, 1);
m1 = size(Asmth1, 2);
ii1 = randi(n1, 80, 1);
jj1 = randi(m1, 80, 1);

Acheck1 = byindex.kernbyindex(ii1, jj1, srfrs1, kerns1, eps, novers1);

err1 = norm(Acheck1 - Asmth1(ii1,jj1), 'fro') / (norm(Asmth1(ii1,jj1), 'fro') + eps);
fprintf('  relative error: %.2e\n', err1);
assert(err1 < tol, 'FAIL: %.2e', err1);
fprintf('  PASS\n\n');

%% Scenario 2: explicit novers vector, compared against manually built reference
fprintf('Scenario 2: explicit novers, oversampled smooth rule ...\n');

S2  = geometries.ellipsoid([1,1,1], [3,3,3], [], 4);
kern2   = kernel3d('l', 's');
opdims2 = kern2.opdims;

nover2 = (S2.norders(1) + 2) * ones(S2.npatches, 1);

[S2over, xinterp2] = S2.oversample(nover2);
wtsover2 = repmat(S2over.wts(:).', opdims2(2), 1);
wtsover2 = wtsover2(:).';

Aover_ref = kern2.eval(S2over, S2) .* wtsover2;

src_norm = max(vecnorm(S2over.r), 1);
for q = 1:size(S2over.r, 2)
    diff_norms = vecnorm(S2.r - S2over.r(:,q));
    self_mask  = diff_norms < 1e-14 * src_norm(q);
    if any(self_mask)
        row_off = (1:opdims2(1)) + opdims2(1)*(find(self_mask).' - 1);
        col_off = (1:opdims2(2)) + opdims2(2)*(q-1);
        Aover_ref(row_off(:), col_off) = 0;
    end
end

Aref2 = Aover_ref * kron(xinterp2, eye(opdims2(2)));

n2  = S2.npts;
ii2 = randi(n2, 60, 1);
jj2 = randi(n2, 60, 1);

Acheck2 = byindex.kernbyindex(ii2, jj2, S2, kern2, eps, {nover2});

err2 = norm(Acheck2 - Aref2(ii2,jj2), 'fro') / (norm(Aref2(ii2,jj2), 'fro') + eps);
fprintf('  relative error: %.2e\n', err2);
assert(err2 < tol, 'FAIL: %.2e', err2);
fprintf('  PASS\n\n');

%% Scenario 3: two-surfer, 2x2 stok kernels (3x3 opdims) + nonsmoothonly Qsparse
fprintf('Scenario 3: two-surfer, 2x2 stok kernels + nonsmoothonly corrections ...\n');

% All stok kernels have opdims [3,3], so row/col opdims are consistent
% across the 2x2 array: every row block has opdims(1)=3, every col block
% has opdims(2)=3.
S3a = geometries.ellipsoid([1,1,1], [1,1,1], [],      5);
S3b = geometries.ellipsoid([1,1,1], [1,1,1], [2.5;0;0], 5);
srfrs3 = [S3a, S3b];

kerns3(2,2) = kernel3d();
kerns3(1,1) = kernel3d('stok', 's');
kerns3(2,1) = kernel3d('stok', 'sp');
kerns3(1,2) = kernel3d('stok', 'd');
kerns3(2,2) = kernel3d('stok', 's');   % reuse s for (2,2) — opdims still [3,3]

% Full surfermat (default opts: selfquad on) as reference
[Afull3,novers3] = surfermat(srfrs3, kerns3, eps);

% Nonsmoothonly sparse matrix: full near-singular quadrature values
opts_ns = struct('nonsmoothonly', 1);
[Qns3, ~] = surfermat(srfrs3, kerns3, eps, opts_ns);

n3 = size(Afull3, 1);
m3 = size(Afull3, 2);
ii3 = randi(n3, 100, 1);
jj3 = randi(m3, 100, 1);

% With Qsparse: smooth entries overwritten by Qns3 where near-singular (default)
Acorr3 = byindex.kernbyindex(ii3, jj3, srfrs3, kerns3, eps, novers3, Qns3);

err3 = norm(Acorr3 - Afull3(ii3,jj3), 'fro') / (norm(Afull3(ii3,jj3), 'fro') + eps);
fprintf('  relative error with corrections: %.2e\n', err3);
assert(err3 < tol, 'FAIL: %.2e', err3);
fprintf('  PASS\n\n');

%% Scenario 4: single surfer, lap s
fprintf('Scenario 4: single surfer, lap s ...\n');

S4   = geometries.ellipsoid([1,1,1], [3,3,3], [], 5);
kern4   = kernel3d('l', 's');
opdims4 = kern4.opdims;   % [3, 3]

[Asmth4, novers4] = surfermat(S4, kern4, eps, opts_sm);

n4  = S4.npts * opdims4(1);
m4  = S4.npts * opdims4(2);
ii4 = randi(n4, 80, 1);
jj4 = randi(m4, 80, 1);

Acheck4 = byindex.kernbyindex(ii4, jj4, S4, kern4, eps, novers4);

err4 = norm(Acheck4 - Asmth4(ii4,jj4), 'fro') / (norm(Asmth4(ii4,jj4), 'fro') + eps);
fprintf('  relative error: %.2e\n', err4);
assert(err4 < tol, 'FAIL: %.2e', err4);
fprintf('  PASS\n\n');

%% Scenario 5: wrong novers should give wrong answer (negative test)
fprintf('Scenario 5: wrong novers gives wrong answer ...\n');

nover_wrong  = S2.norders(1) * ones(S2.npatches, 1);
Acheck_wrong = byindex.kernbyindex(ii2, jj2, S2, kern2, eps, {nover_wrong});

err_wrong = norm(Acheck_wrong - Aref2(ii2,jj2), 'fro') / (norm(Aref2(ii2,jj2), 'fro') + eps);
fprintf('  relative error with wrong novers: %.2e\n', err_wrong);
assert(err_wrong > tol, 'FAIL: wrong novers should give error > tol but got %.2e', err_wrong);
fprintf('  PASS (correctly wrong)\n\n');
