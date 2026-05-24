%TEST_UNIF_NOVER   Verify that the unif_nover option in surfermat makes
%   oversampling orders uniform across target surfaces for each source j.
%
%   With two source surfaces S1, S2 and unif_nover=1, the returned novers
%   cell must satisfy novers{1,j} == novers{2,j} for j=1,2, and each
%   value must be >= the corresponding entry without unif_nover.

run ../startup.m

fprintf('=== test_unif_nover ===\n\n');

S1 = geometries.sphere(1, 4, [0;0;0], 4);
S2 = geometries.sphere(1, 4, [4;0;0], 4);
srfrs = [S1, S2];

kern = kernel3d('laplace', 'c', [1;1]);
eps  = 1e-6;

% --- without unif_nover ---
[~, novers_base] = surfermat(srfrs, kern, eps, struct('corrections', 1));

% --- with unif_nover ---
opts = struct('corrections', 1, 'unif_nover', 1);
[~, novers_unif] = surfermat(srfrs, kern, eps, opts);

nfail = 0;

% Test 1: novers{1,j} == novers{2,j} for each source surface j
fprintf('Test 1: novers constant across target surfaces for each source j\n');
ok = true;
for j = 1:2
    if ~isequal(novers_unif{1,j}, novers_unif{2,j})
        fprintf('  FAIL: novers{1,%d} ~= novers{2,%d}\n', j, j);
        ok = false;
    end
end
if ok
    fprintf('  PASS\n');
else
    nfail = nfail + 1;
end

% Test 2: uniform novers >= base novers entry-wise (never decreases)
fprintf('Test 2: unif_nover never decreases any oversampling order\n');
ok = true;
for i = 1:2
    for j = 1:2
        if any(novers_unif{i,j} < novers_base{i,j})
            fprintf('  FAIL: novers_unif{%d,%d} < novers_base{%d,%d} on some patch\n', i, j, i, j);
            ok = false;
        end
    end
end
if ok
    fprintf('  PASS\n');
else
    nfail = nfail + 1;
end

% Test 3: uniform value equals the element-wise max of the base column
fprintf('Test 3: uniform value equals element-wise max over target surfaces\n');
ok = true;
for j = 1:2
    expected = max(novers_base{1,j}, novers_base{2,j});
    for i = 1:2
        if ~isequal(novers_unif{i,j}, expected)
            fprintf('  FAIL: novers_unif{%d,%d} does not equal column max\n', i, j);
            ok = false;
        end
    end
end
if ok
    fprintf('  PASS\n');
else
    nfail = nfail + 1;
end

fprintf('\n--- %d test(s) failed ---\n', nfail);
if nfail == 0, fprintf('All tests PASSED.\n'); end
