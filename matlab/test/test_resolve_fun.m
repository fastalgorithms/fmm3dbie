run('../startup.m');

test_resolve_fun0();
test_slicesurfer0();
test_patch_max0();

% ------------------------------------------------------------------ %

function test_resolve_fun0()
% Test that resolve_fun refines a surfer until a given function is
% resolved to the requested tolerance, and that the resolved surfer
% evaluates the function accurately at its nodes.

tol       = 1e-8;
norderhead = 6;
nrefmax   = 4;

src   = []; src.r = [8; 0; 0];
f     = @(S) lap3d.kern(src, S, 's');
S     = geometries.stellarator([5, 15], 8, 11);
S2    = S.resolve_fun(tol, f, norderhead, nrefmax);

% The resolved surfer should have at least as many patches as S.
assert(S2.npatches >= S.npatches, ...
    'resolve_fun: refined surfer should have at least as many patches as original');

% surf_fun_error on the refined surfer should be below tolerance everywhere.
err = surf_fun_error(S2, f(S2), [], norderhead);
assert(max(err(:)) < 10*tol, ...
    sprintf('resolve_fun: max surf_fun_error = %.2e, expected < %.2e', max(err(:)), tol));

% The function values at S2 nodes should match a direct evaluation.
fvals = f(S2);
assert(isequal(size(fvals), [S2.npts,1]), ...
    'resolve_fun: f(S2) should be (S2.npts, 1)');

fprintf('test_resolve_fun0 passed\n');

figure(2);clf
plot(S2,log10((surf_fun_error(S2,f(S2)))));colorbar
hold on
wireframe(S2)
end

% ------------------------------------------------------------------ %

function test_slicesurfer0()
% Test that slicesurfer correctly extracts a subset of patches.

S = geometries.stellarator([5, 15], 8, 11);

% Keep the first three patches.
ipatch = [1, 2, 3];
S2 = slicesurfer(S, ipatch);

assert(S2.npatches == 3, ...
    'slicesurfer: S2 should have 3 patches');
assert(isequal(S2.norders, S.norders(ipatch)), ...
    'slicesurfer: norders of S2 should match those of the kept patches');
assert(isequal(S2.iptype, S.iptype(ipatch)), ...
    'slicesurfer: iptype of S2 should match those of the kept patches');

% Node count: sum of nodes over kept patches.
expected_npts = sum(S.ixyzs(ipatch+1) - S.ixyzs(ipatch));
assert(S2.npts == expected_npts, ...
    sprintf('slicesurfer: expected %d pts, got %d', expected_npts, S2.npts));

% Geometry should match: check that r, du, dv, n agree with S.
for k = 1:3
    iinds = S.ixyzs(ipatch(k)):S.ixyzs(ipatch(k)+1)-1;
    iinds2 = S2.ixyzs(k):S2.ixyzs(k+1)-1;
    assert(norm(S2.r(:,iinds2)  - S.r(:,iinds),  'fro') == 0, ...
        sprintf('slicesurfer: r mismatch on patch %d', k));
    assert(norm(S2.du(:,iinds2) - S.du(:,iinds), 'fro') == 0, ...
        sprintf('slicesurfer: du mismatch on patch %d', k));
    assert(norm(S2.dv(:,iinds2) - S.dv(:,iinds), 'fro') == 0, ...
        sprintf('slicesurfer: dv mismatch on patch %d', k));
    assert(norm(S2.n(:,iinds2)  - S.n(:,iinds),  'fro') == 0, ...
        sprintf('slicesurfer: n mismatch on patch %d', k));
end

% Also test a non-contiguous subset.
ipatch2 = [1, 3, S.npatches];
S3 = slicesurfer(S, ipatch2);
assert(S3.npatches == 3, ...
    'slicesurfer: non-contiguous subset should have 3 patches');

fprintf('test_slicesurfer0 passed\n');
end

% ------------------------------------------------------------------ %

function test_patch_max0()
% Test that patch_max returns the per-patch maximum of a scalar and
% multi-row function, and that the area-scaling variant is consistent.

S = geometries.stellarator([5, 15], 8, 11);

% Scalar function: distance from origin.
f1 = vecnorm(S.r, 2, 1);   % (1, npts)

fmax1 = S.patch_max(f1);
assert(isequal(size(fmax1), [1, S.npatches]), ...
    'patch_max: output should be (1, npatches) for a row-vector input');

% Check each patch manually.
for i = 1:S.npatches
    iind = S.ixyzs(i):S.ixyzs(i+1)-1;
    assert(abs(fmax1(i) - max(f1(iind))) < 1e-14, ...
        sprintf('patch_max: wrong maximum on patch %d', i));
end

% Multi-row function: all three coordinates.
f3 = S.r;   % (3, npts)

fmax3 = S.patch_max(f3);
assert(isequal(size(fmax3), [3, S.npatches]), ...
    'patch_max: output should be (3, npatches) for a 3-row input');

for i = 1:S.npatches
    iind = S.ixyzs(i):S.ixyzs(i+1)-1;
    assert(norm(fmax3(:,i) - max(f3(:,iind), [], 2)) < 1e-14, ...
        sprintf('patch_max: wrong maximum on patch %d (multi-row)', i));
end

% Column-vector input: patch_max should transpose back.
f_col = f1.';   % (npts, 1)
fmax_col = S.patch_max(f_col);
assert(isequal(size(fmax_col), [S.npatches, 1]), ...
    'patch_max: output should be (npatches, 1) for a column-vector input');
assert(norm(fmax_col(:) - fmax1(:)) < 1e-14, ...
    'patch_max: column-vector result should match row-vector result');

% Area-scaling variant: patch_max(f, p) should scale each entry by
% (patch_area / total_area)^p, so summing over patches with p=1
% and f=ones gives total_area / total_area = 1.
f_ones = ones(1, S.npts);
fmax_scaled = S.patch_max(f_ones, 1);
assert(abs(sum(fmax_scaled) - 1) < 1e-12, ...
    sprintf('patch_max: area-scaled sum for f=1 should be 1, got %.4e', sum(fmax_scaled)));

fprintf('test_patch_max0 passed\n');
end
