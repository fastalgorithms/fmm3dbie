% Tests for vertex_surfer, oversample (nover == norder), wireframe, and split_patches
% assumes pwd is the directory this script is in

run('../startup.m');

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


%% triangle endpoint order and alignment

uvs = koorn.rv_nodes(1);
r = [uvs; zeros(1,size(uvs,2))];
du = repmat([1;0;0], 1, size(uvs,2));
dv = repmat([0;1;0], 1, size(uvs,2));
n = repmat([0;0;1], 1, size(uvs,2));
Slin = surfer(1, 1, [r;du;dv;n], 1);

tri_verts = [0 1 0; 0 0 1; 0 0 0];
assert(norm(Slin.end_pt_verts{1} - tri_verts, 'fro') < 1e-12, ...
    'surfer: triangle end_pt_verts should use Fortran vertex order');

for k = 1:3
    distmin = ones(3,1);
    distmin(k) = 0;
    Svert = align_patches(Slin, distmin, zeros(2,1), 102);
    assert(norm(Svert.end_pt_verts{1}(:,3) - tri_verts(:,k)) < 1e-12, ...
        'align_patches: nearest triangle vertex should rotate to (0,1)');
end

near_edges = [1 2; 1 3; 2 3];
edge_dists = [0 0 1; 0 1 0; 1 0 0];
for k = 1:3
    Sedge = align_patches(Slin, edge_dists(k,:).', zeros(2,1), 101);
    edge_got = Sedge.end_pt_verts{1}(:,1:2);
    edge_exp = tri_verts(:,near_edges(k,:));
    err = min(norm(edge_got - edge_exp, 'fro'), ...
              norm(edge_got - fliplr(edge_exp), 'fro'));
    assert(err < 1e-12, ...
        'align_patches: near triangle edge should rotate to leading edge');
end


%% vertex_surfer

S = geometries.sphere(1, 2, [0;0;0], 6, 1);
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
wireframe(S, struct('nfac', -1, 'wfill', true));
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


%% edge alignment


S = geometries.disk([1,1,1],[],[],8);

figure(10);clf
subplot(1,2,1)
plot(S,S.uvs_targ(2,:));colorbar

vs = cell2mat(S.end_pt_verts(:).');
ds = 1-vecnorm(vs);

S_rot = align_patches(S,ds,[]);
subplot(1,2,2)
plot(S_rot,S_rot.uvs_targ(2,:));colorbar

% test alignment
% TODO test vertex-on patches
thetas = 2*pi*rand(1,100);
[~,~,uvs,dist,flag] = get_closest_pts(S_rot,[cos(thetas);sin(thetas);0*thetas]);

assert(norm(uvs(2,:))<1e-4)
assert(max(dist)<1e-4)
assert(max(flag)==0)

S = geometries.disk([1,1,1],[],[],8,11);

dstemp = ones(1,S.npatches*4);
dstemp(1:4:end) = 0;
dstemp(2:4:end) = 0;
S = align_patches(S,dstemp,[],[],[2,3,1]);

figure(11);clf
subplot(1,2,1)
plot(S,S.uvs_targ(2,:));colorbar

vs = cell2mat(S.end_pt_verts(:).');
ds = 1-vecnorm(vs);

S_rot = align_patches(S,ds,[]);
subplot(1,2,2)
plot(S_rot,S_rot.uvs_targ(2,:));colorbar

% test alignment
thetas = 2*pi*rand(1,100);
[~,~,uvs,dist,flag] = get_closest_pts(S_rot,[cos(thetas);sin(thetas);0*thetas]);
assert(norm(uvs(2,:)+1)<1e-4)
assert(max(dist)<1e-4)
assert(max(flag)==0)








