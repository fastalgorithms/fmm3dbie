% Tests for vertex_surfer, oversample (nover == norder), wireframe, and split_patches
% assumes pwd is the directory this script is in

addpath '../';

%% oversample with nover == norder

% Triangular patches (iptype=1)
S = geometries.sphere(1, 2, [0;0;0], 4, 1);
norder = S.norders(1);
[Sover, val2over] = oversample(S, norder);

assert(Sover.npts == S.npts, 'oversample nover==norder tria: npts should be unchanged');
rerr = norm(val2over - speye(S.npts), 'fro') / S.npts;
assert(rerr < 1e-12, 'oversample nover==norder tria: val2over should be identity');

% Quadrilateral patches (iptype=11)
S = geometries.sphere(1, 2, [0;0;0], 4, 11);
norder = S.norders(1);
[Sover, val2over] = oversample(S, norder);

assert(Sover.npts == S.npts, 'oversample nover==norder quad: npts should be unchanged');
rerr = norm(val2over - speye(S.npts), 'fro') / S.npts;
assert(rerr < 1e-12, 'oversample nover==norder quad: val2over should be identity');


%% vertex_surfer

S = geometries.sphere(1, 2, [0;0;0], 4, 1);
vert = S.end_pt_verts{1}(:,1);

% Auto-detect patches sharing vert
[Sout, islice] = vertex_surfer(S, vert);
assert(~isempty(Sout), 'vertex_surfer: should find at least one patch at vertex');

% Explicit idvertpatch covering all patches
[Sout2, ~] = vertex_surfer(S, vert, 1:S.npatches);
assert(Sout2.npatches == S.npatches, 'vertex_surfer: explicit idvertpatch should cover all patches');

% Sshift vertex-patch positions should match Sout
[Sout3, islice3, Sshift] = vertex_surfer(S, vert);
rerr = norm(Sshift.r(:, islice3) - Sout3.r, 'fro');
assert(rerr < 1e-12, 'vertex_surfer: Sshift vertex-patch positions should match Sout');

% Non-vertex patches should be S.r - vert
non_islice = setdiff(1:S.npts, islice3);
rerr = norm(Sshift.r(:, non_islice) - (S.r(:, non_islice) - vert), 'fro');
assert(rerr < 1e-12, 'vertex_surfer: Sshift non-vertex patches should equal S.r - vert');

% Tight tolerance should find no patches
[Sout_notol, ~] = vertex_surfer(S, vert, [], 1e-15);
assert(isempty(Sout_notol), 'vertex_surfer: tight tol should find no patches');


%% wireframe

S = geometries.sphere(1, 2, [0;0;0], 4, 1);
fig = figure('Visible','off');
wireframe(S);
wireframe(S, struct('nsign', -1, 'wfill', true));
close(fig);

S = geometries.sphere(1, 2, [0;0;0], 4, 11);
fig = figure('Visible','off');
wireframe(S);
close(fig);


%% split_patches

S = geometries.sphere(1, 2, [0;0;0], 4, 1);
npatches0 = S.npatches;
npols = S.ixyzs(2) - S.ixyzs(1);

% Split no patches: output should match S
Sout = split_patches(S, []);
assert(Sout.npatches == npatches0, 'split_patches: no split should leave npatches unchanged');
assert(Sout.npts == S.npts,        'split_patches: no split should leave npts unchanged');

% Split one patch: should add 3 extra patches
[Sout, val2split] = split_patches(S, 1);
assert(Sout.npatches == npatches0 + 3,    'split_patches: one split should add 3 patches');
assert(Sout.npts == S.npts + 3*npols,     'split_patches: one split should add 3*npols points');
assert(size(val2split,1) == Sout.npts,    'split_patches: val2split rows should equal objout.npts');
assert(size(val2split,2) == S.npts,       'split_patches: val2split cols should equal obj.npts');

% val2split should map S positions to Ssplit positions
[Sout, val2split] = split_patches(S, 1);
rerr = norm(Sout.r - S.r * val2split.', 'fro') / norm(Sout.r, 'fro');
assert(rerr < 1e-12, 'split_patches: Ssplit.r should equal S.r * val2split.''');

% Split all patches
Sout = split_patches(S, 1:npatches0);
assert(Sout.npatches == 4*npatches0, 'split_patches: splitting all patches should quadruple npatches');
assert(Sout.npts == 4*S.npts,        'split_patches: splitting all patches should quadruple npts');
