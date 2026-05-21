% Tests for surfer.load_from_msh covering every supported gmsh variant:
%   - format versions 2.x, 4.0, 4.1, both ASCII and binary (LE and BE)
%   - triangle orders 1..10 and quad orders 1..10 (full Lagrange)
%   - Q8 serendipity projected to Q9
%   - mixed-shape and mixed-order meshes
%   - parser robustness: blank lines, tabs, CRLF, ignored sections
%   - polynomial recovery (degree-p input recovered at the target nodes
%     to machine precision)
%   - large stress mesh (5000 triangles)
%   - defensive error when $MeshFormat is missing
%
% Flat patches in z=0 are exactly representable at any polynomial order,
% so geometry-based assertions can use a tight tolerance.

addpath '../'

tol = 1e-13;

%% 1. Order-1 (3-node) single triangle, vertices (0,0,0),(1,0,0),(0,1,0)
f1 = [tempname '.msh'];
fid = fopen(f1, 'w');
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');
fprintf(fid, '$Nodes\n3\n');
fprintf(fid, '1 0 0 0\n2 1 0 0\n3 0 1 0\n');
fprintf(fid, '$EndNodes\n$Elements\n1\n');
fprintf(fid, '1 2 2 1 1 1 2 3\n');
fprintf(fid, '$EndElements\n');
fclose(fid);
S1 = surfer.load_from_msh(f1);
delete(f1);

assert(S1.npatches == 1, 'order-1: wrong npatches');
assert(S1.norders(1) == 1, 'order-1: wrong norder');
assert(S1.npts == 3, 'order-1: wrong npts');
% flat in z=0
assert(all(abs(S1.r(3,:)) < tol), 'order-1: z not zero');
% du,dv constant equal to edge vectors
assert(norm(S1.du - repmat([1;0;0],1,3), 'fro') < tol, 'order-1: du wrong');
assert(norm(S1.dv - repmat([0;1;0],1,3), 'fro') < tol, 'order-1: dv wrong');
% unit normal (0,0,1) everywhere
assert(norm(S1.n - repmat([0;0;1],1,3), 'fro') < tol, 'order-1: normal wrong');
% total area = 0.5
assert(abs(sum(S1.wts) - 0.5) < tol, 'order-1: area wrong');
fprintf('order-1 single triangle  : OK\n');

%% 2. Order-2 (6-node) same flat triangle (corners + edge midpoints)
f2 = [tempname '.msh'];
fid = fopen(f2, 'w');
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');
fprintf(fid, '$Nodes\n6\n');
fprintf(fid, '1 0 0 0\n2 1 0 0\n3 0 1 0\n');
fprintf(fid, '4 0.5 0 0\n5 0.5 0.5 0\n6 0 0.5 0\n');
fprintf(fid, '$EndNodes\n$Elements\n1\n');
fprintf(fid, '1 9 2 1 1 1 2 3 4 5 6\n');
fprintf(fid, '$EndElements\n');
fclose(fid);
S2 = surfer.load_from_msh(f2);
delete(f2);

assert(S2.npatches == 1, 'order-2: wrong npatches');
assert(S2.norders(1) == 2, 'order-2: wrong norder');
assert(S2.npts == 6, 'order-2: wrong npts');
assert(all(abs(S2.r(3,:)) < tol), 'order-2: z not zero');
assert(norm(S2.du - repmat([1;0;0],1,6), 'fro') < tol, 'order-2: du wrong');
assert(norm(S2.dv - repmat([0;1;0],1,6), 'fro') < tol, 'order-2: dv wrong');
assert(norm(S2.n - repmat([0;0;1],1,6), 'fro') < tol, 'order-2: normal wrong');
assert(abs(sum(S2.wts) - 0.5) < tol, 'order-2: area wrong');
fprintf('order-2 single triangle  : OK\n');

%% 3. Two order-1 triangles forming the unit square, mixed line/triangle elems
f3 = [tempname '.msh'];
fid = fopen(f3, 'w');
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');
fprintf(fid, '$Nodes\n4\n');
fprintf(fid, '1 0 0 0\n2 1 0 0\n3 1 1 0\n4 0 1 0\n');
fprintf(fid, '$EndNodes\n$Elements\n3\n');
% an unsupported 1-line element (type 1) -- must be ignored
fprintf(fid, '1 1 2 1 1 1 2\n');
fprintf(fid, '2 2 2 1 1 1 2 3\n');
fprintf(fid, '3 2 2 1 1 1 3 4\n');
fprintf(fid, '$EndElements\n');
fclose(fid);
S3 = surfer.load_from_msh(f3);
delete(f3);

assert(S3.npatches == 2, 'square: wrong npatches (unsupported elem not ignored?)');
assert(S3.npts == 6, 'square: wrong npts');
assert(abs(sum(S3.wts) - 1.0) < tol, 'square: total area should be 1');
fprintf('square (2 tris, +ignored line element): OK\n');

%% 4. gmsh v4.1 ASCII, same flat triangle (entity-block layout, tags then xyz)
f4 = [tempname '.msh'];
fid = fopen(f4, 'w');
fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
fprintf(fid, '$Nodes\n1 3 1 3\n2 1 0 3\n1\n2\n3\n0 0 0\n1 0 0\n0 1 0\n$EndNodes\n');
fprintf(fid, '$Elements\n1 1 1 1\n2 1 2 1\n1 1 2 3\n$EndElements\n');
fclose(fid);
S4 = surfer.load_from_msh(f4);
delete(f4);

assert(S4.npatches == 1, 'v4.1: wrong npatches');
assert(S4.norders(1) == 1, 'v4.1: wrong norder');
assert(S4.npts == 3, 'v4.1: wrong npts');
assert(all(abs(S4.r(3,:)) < tol), 'v4.1: z not zero');
assert(norm(S4.n - repmat([0;0;1],1,3), 'fro') < tol, 'v4.1: normal wrong');
assert(abs(sum(S4.wts) - 0.5) < tol, 'v4.1: area wrong');
fprintf('v4.1 single triangle     : OK\n');

%% 5. gmsh v4.0 ASCII, same flat triangle (interleaved tag-xyz per node line)
f5 = [tempname '.msh'];
fid = fopen(f5, 'w');
fprintf(fid, '$MeshFormat\n4.0 0 8\n$EndMeshFormat\n');
fprintf(fid, '$Nodes\n1 3\n1 2 0 3\n1 0 0 0\n2 1 0 0\n3 0 1 0\n$EndNodes\n');
fprintf(fid, '$Elements\n1 1\n1 2 2 1\n1 1 2 3\n$EndElements\n');
fclose(fid);
S5 = surfer.load_from_msh(f5);
delete(f5);

assert(S5.npatches == 1, 'v4.0: wrong npatches');
assert(S5.norders(1) == 1, 'v4.0: wrong norder');
assert(S5.npts == 3, 'v4.0: wrong npts');
assert(all(abs(S5.r(3,:)) < tol), 'v4.0: z not zero');
assert(norm(S5.n - repmat([0;0;1],1,3), 'fro') < tol, 'v4.0: normal wrong');
assert(abs(sum(S5.wts) - 0.5) < tol, 'v4.0: area wrong');
fprintf('v4.0 single triangle     : OK\n');

%% 6-8. Binary round-trip in v2.2 / v4.0 / v4.1, both little- and big-endian.
% The big-endian iterations exercise the endian-probe + byte-swap branch
% (host is typically little-endian).
for en_cell = {'ieee-le', 'ieee-be'}
en = en_cell{1};
for variant = {'2.2','4.0','4.1'}
    f = [tempname '.msh'];
    fid = fopen(f, 'wb');
    fprintf(fid, '$MeshFormat\n%s 1 8\n', variant{1});
    fwrite(fid, 1, 'int32', 0, en);                 % endian probe
    fprintf(fid, '\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n');
    switch variant{1}
      case '2.2'
        fprintf(fid, '3\n');
        for k = 1:3
            fwrite(fid, k, 'int32', 0, en);
            switch k
              case 1, p = [0 0 0]; case 2, p = [1 0 0]; case 3, p = [0 1 0];
            end
            fwrite(fid, p, 'double', 0, en);
        end
      case '4.0'
        fwrite(fid, [1 3], 'uint64', 0, en);        % nBlocks, nNodes
        fwrite(fid, [1 2 0], 'int32', 0, en);       % tagEnt dimEnt typeNode
        fwrite(fid, 3, 'uint64', 0, en);            % m
        for k = 1:3
            fwrite(fid, k, 'int32', 0, en);
            switch k
              case 1, p = [0 0 0]; case 2, p = [1 0 0]; case 3, p = [0 1 0];
            end
            fwrite(fid, p, 'double', 0, en);
        end
      case '4.1'
        fwrite(fid, [1 3 1 3], 'uint64', 0, en);    % nBlocks nNodes minTag maxTag
        fwrite(fid, [2 1 0], 'int32', 0, en);       % entityDim entityTag parametric
        fwrite(fid, 3, 'uint64', 0, en);            % m
        fwrite(fid, [1 2 3], 'uint64', 0, en);      % tags
        fwrite(fid, [0 0 0  1 0 0  0 1 0], 'double', 0, en);  % coords
    end
    fprintf(fid, '\n$EndNodes\n$Elements\n');
    switch variant{1}
      case '2.2'
        fprintf(fid, '1\n');
        fwrite(fid, [2 1 2], 'int32', 0, en);       % elemType numOfType numTags
        fwrite(fid, [1 1 1 1 2 3], 'int32', 0, en); % elemTag tag1 tag2 n1 n2 n3
      case '4.0'
        fwrite(fid, [1 1], 'uint64', 0, en);        % nBlocks nElems
        fwrite(fid, [1 2 2], 'int32', 0, en);       % tagEnt dimEnt typeEle
        fwrite(fid, 1, 'uint64', 0, en);            % m
        fwrite(fid, [1 1 2 3], 'int32', 0, en);     % elemTag n1 n2 n3
      case '4.1'
        fwrite(fid, [1 1 1 1], 'uint64', 0, en);    % nBlocks nElems minTag maxTag
        fwrite(fid, [2 1 2], 'int32', 0, en);       % entityDim entityTag elementType
        fwrite(fid, 1, 'uint64', 0, en);            % m
        fwrite(fid, [1 1 2 3], 'uint64', 0, en);    % elemTag n1 n2 n3
    end
    fprintf(fid, '\n$EndElements\n');
    fclose(fid);

    S = surfer.load_from_msh(f);
    delete(f);
    msgsuf = sprintf('v%s %s bin', variant{1}, en);
    assert(S.npatches == 1,                              [msgsuf ': wrong npatches']);
    assert(S.norders(1) == 1,                            [msgsuf ': wrong norder']);
    assert(S.npts == 3,                                  [msgsuf ': wrong npts']);
    assert(all(abs(S.r(3,:)) < tol),                     [msgsuf ': z not zero']);
    assert(norm(S.n - repmat([0;0;1],1,3),'fro') < tol,  [msgsuf ': normal wrong']);
    assert(abs(sum(S.wts) - 0.5) < tol,                  [msgsuf ': area wrong']);
    fprintf('v%s binary (%s) single triangle: OK\n', variant{1}, en);
end
end

%% 9. Quad patches Q4/Q9/Q8/Q16/Q25, each as a flat unit square in z=0.
%    Map each gmsh reference uv in [-1,1]^2 to (u+1)/2 in [0,1]^2.
%    Default opts.iptype = 11 (Gauss-Legendre tensor product).
cases = { ...
    {3 , 1, 4 , 'Q4 (type 3)'},  {10, 2, 9 , 'Q9 (type 10)'}, ...
    {16, 2, 8 , 'Q8 serendipity (type 16, projected)'}, ...
    {36, 3, 16, 'Q16 (type 36)'}, {37, 4, 25, 'Q25 (type 37)'}, ...
};
for c = 1:numel(cases)
    [qtype, p, nn, desc] = deal(cases{c}{:});
    uv_full = i_build_gmsh_quad_uvs_local(p);
    uv_use  = uv_full(:, 1:nn);                    % Q8 keeps only outer 8 of 9
    xyz     = [(uv_use + 1)/2; zeros(1, nn)];      % map to [0,1]^2 x {0}

    f = [tempname '.msh'];
    fid = fopen(f, 'w');
    fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n%d\n', nn);
    for k = 1:nn, fprintf(fid, '%d %.15g %.15g %.15g\n', k, xyz(:,k)); end
    fprintf(fid, '$EndNodes\n$Elements\n1\n');
    fprintf(fid, '1 %d 2 1 1', qtype);
    fprintf(fid, ' %d', 1:nn);
    fprintf(fid, '\n$EndElements\n');
    fclose(fid);

    S = surfer.load_from_msh(f);
    delete(f);

    nq = (p+1)^2;
    assert(S.npatches == 1,                              [desc ': npatches']);
    assert(S.iptype(1) == 11,                            [desc ': default iptype']);
    assert(S.norders(1) == p,                            [desc ': norder']);
    assert(S.npts == nq,                                 [desc ': npts']);
    assert(all(abs(S.r(3,:)) < tol),                     [desc ': z']);
    assert(norm(S.n - repmat([0;0;1],1,nq),'fro') < tol, [desc ': normal']);
    assert(abs(sum(S.wts) - 1.0) < tol,                  [desc ': area']);
    fprintf('%-44s: OK\n', desc);
end

%% 10. opts.iptype override (Chebyshev quads) on a Q9 patch.
qtype = 10;
p = 2;
uv_full = i_build_gmsh_quad_uvs_local(p);
xyz = [(uv_full + 1)/2; zeros(1, 9)];
f = [tempname '.msh'];
fid = fopen(f, 'w');
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n9\n');
for k = 1:9, fprintf(fid, '%d %.15g %.15g %.15g\n', k, xyz(:,k)); end
fprintf(fid, '$EndNodes\n$Elements\n1\n1 %d 2 1 1', qtype);
fprintf(fid, ' %d', 1:9);
fprintf(fid, '\n$EndElements\n');
fclose(fid);
S = surfer.load_from_msh(f, struct('iptype', 12));
delete(f);
assert(S.iptype(1) == 12, 'opts.iptype=12 override not honored');
assert(abs(sum(S.wts) - 1.0) < tol, 'opts.iptype=12: area should be 1');
fprintf('opts.iptype=12 override                     : OK\n');


%% 11. Higher-order tri (p=5..10) and quad (p=5..10) on flat unit patches.
%     Tri patch:  vertices (0,0,0),(1,0,0),(0,1,0) -- standard simplex, area 0.5.
%     Quad patch: maps [-1,1]^2 reference to [0,1]^2 in z=0 -- area 1.
%     Tests exercise the shell algorithm at every order for both shapes.
gmsh_tri_type = containers.Map( ...
    {1,2,3,4,5,6,7,8,9,10}, {2,9,21,23,25,42,43,44,45,46});
gmsh_quad_type = containers.Map( ...
    {1,2,3,4,5,6,7,8,9,10}, {3,10,36,37,38,47,48,49,50,51});
for p = 5:10
    % ---- triangle ----
    uvt = i_build_gmsh_tri_uvs_local(p);
    nt  = size(uvt, 2);
    xyz = [uvt; zeros(1, nt)];           % flat tri on standard simplex, z=0
    f = [tempname '.msh'];
    fid = fopen(f, 'w');
    fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n', nt);
    for k = 1:nt, fprintf(fid, '%d %.15g %.15g %.15g\n', k, xyz(:,k)); end
    fprintf(fid, '$EndNodes\n$Elements\n1\n1 %d 2 1 1', gmsh_tri_type(p));
    fprintf(fid, ' %d', 1:nt); fprintf(fid, '\n$EndElements\n');
    fclose(fid);
    S = surfer.load_from_msh(f);
    delete(f);
    nrv = (p+1)*(p+2)/2;
    assert(S.npatches == 1 && S.norders(1) == p && S.npts == nrv, ...
           sprintf('tri p=%d: shape/order/npts', p));
    assert(all(abs(S.r(3,:)) < tol),                       sprintf('tri p=%d: z', p));
    assert(norm(S.n - repmat([0;0;1],1,nrv),'fro') < tol,  sprintf('tri p=%d: normal', p));
    assert(abs(sum(S.wts) - 0.5) < tol,                    sprintf('tri p=%d: area', p));
    fprintf('tri  p=%-2d  (type %2d, %d nodes): OK\n', p, gmsh_tri_type(p), nt);

    % ---- quad ----
    uvq = i_build_gmsh_quad_uvs_local(p);
    nq  = size(uvq, 2);
    xyz = [(uvq + 1)/2; zeros(1, nq)];   % map ref [-1,1]^2 -> [0,1]^2, z=0
    f = [tempname '.msh'];
    fid = fopen(f, 'w');
    fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n', nq);
    for k = 1:nq, fprintf(fid, '%d %.15g %.15g %.15g\n', k, xyz(:,k)); end
    fprintf(fid, '$EndNodes\n$Elements\n1\n1 %d 2 1 1', gmsh_quad_type(p));
    fprintf(fid, ' %d', 1:nq); fprintf(fid, '\n$EndElements\n');
    fclose(fid);
    S = surfer.load_from_msh(f);
    delete(f);
    npt = (p+1)^2;
    assert(S.npatches == 1 && S.norders(1) == p && S.npts == npt, ...
           sprintf('quad p=%d: shape/order/npts', p));
    assert(all(abs(S.r(3,:)) < tol),                       sprintf('quad p=%d: z', p));
    assert(norm(S.n - repmat([0;0;1],1,npt),'fro') < tol,  sprintf('quad p=%d: normal', p));
    assert(abs(sum(S.wts) - 1.0) < tol,                    sprintf('quad p=%d: area', p));
    fprintf('quad p=%-2d  (type %2d, %d nodes): OK\n', p, gmsh_quad_type(p), nq);
end

%% 12. Mixed-shape file: 1 triangle + 1 quad (with a shared edge).
%     Triangle vertices (0,0,0),(1,0,0),(1,1,0); quad (1,0,0),(2,0,0),(2,1,0),(1,1,0).
f = [tempname '.msh'];
fid = fopen(f, 'w');
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');
fprintf(fid, '$Nodes\n5\n');
fprintf(fid, '1 0 0 0\n2 1 0 0\n3 1 1 0\n4 2 0 0\n5 2 1 0\n');
fprintf(fid, '$EndNodes\n$Elements\n2\n');
fprintf(fid, '1 2 2 1 1 1 2 3\n');               % tri3
fprintf(fid, '2 3 2 1 1 2 4 5 3\n');             % Q4
fprintf(fid, '$EndElements\n');
fclose(fid);
S = surfer.load_from_msh(f);
delete(f);
assert(S.npatches == 2,                            'mixed: npatches');
assert(isequal(sort(S.iptype(:).'), [1 11]),       'mixed: iptype vector');
assert(isequal(sort(S.norders(:).'), [1 1]),       'mixed: norders');
assert(abs(sum(S.wts) - (0.5 + 1.0)) < tol,        'mixed: total area = 1.5');
fprintf('mixed-shape (1 tri + 1 quad)                : OK\n');

%% 13. Mixed-order triangles in one file: a p=1 tri sharing an edge with a p=2 tri.
%     Verifies the per-patch (shape, order) dispatch and per-patch norders vector.
f = [tempname '.msh'];
fid = fopen(f, 'w');
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');
% Nodes: 1,2,3 are p=1 tri vertices; 4,5,6 share edge 2-3 of the first tri and add
% a vertex to form a p=2 tri (nodes 7,8,9 are the edge midpoints in gmsh order).
fprintf(fid, '$Nodes\n9\n');
fprintf(fid, '1 0 0 0\n2 1 0 0\n3 0 1 0\n');         % vertices of tri-A (p=1)
fprintf(fid, '4 1 0 0\n5 1 1 0\n6 0 1 0\n');         % corners of tri-B (p=2)
fprintf(fid, '7 1 0.5 0\n8 0.5 1 0\n9 0.5 0.5 0\n'); % edge mids of tri-B (in gmsh order)
fprintf(fid, '$EndNodes\n$Elements\n2\n');
fprintf(fid, '1 2 2 1 1 1 2 3\n');                   % tri3, p=1
fprintf(fid, '2 9 2 1 1 4 5 6 7 8 9\n');             % tri6, p=2
fprintf(fid, '$EndElements\n');
fclose(fid);
S = surfer.load_from_msh(f);
delete(f);
assert(S.npatches == 2,                  'mixed-order: npatches');
assert(isequal(S.norders(:).', [1 2]),   'mixed-order: per-patch norders');
assert(S.npts == 3 + 6,                  'mixed-order: total points');
fprintf('mixed-order (p=1 tri + p=2 tri)             : OK\n');

%% 14. Blank lines, tabs, varying whitespace in v2 ASCII.
f = [tempname '.msh'];
fid = fopen(f, 'w');
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');
fprintf(fid, '\n\n');                            % blank lines between sections
fprintf(fid, '$Nodes\n   3       \n');           % trailing spaces on count line
fprintf(fid, '1 0 0 0\n\n');                     % blank line in nodes
fprintf(fid, '2\t1\t0\t0\n');                    % tabs
fprintf(fid, '3   0   1   0\n');                 % multiple spaces
fprintf(fid, '$EndNodes\n\n');
fprintf(fid, '$Elements\n1\n');
fprintf(fid, '1   2   2   1   1   1   2   3\n');
fprintf(fid, '$EndElements\n');
fclose(fid);
S = surfer.load_from_msh(f);
delete(f);
assert(S.npatches == 1 && abs(sum(S.wts) - 0.5) < tol, ...
       'whitespace-tolerance test failed');
fprintf('blank lines / tabs / extra spaces           : OK\n');

%% 15. CRLF line endings (Windows-style).
f = [tempname '.msh'];
fid = fopen(f, 'wb');
crlf = sprintf('\r\n');
fwrite(fid, ['$MeshFormat'   crlf '2.2 0 8' crlf '$EndMeshFormat' crlf]);
fwrite(fid, ['$Nodes' crlf '3' crlf '1 0 0 0' crlf '2 1 0 0' crlf '3 0 1 0' crlf]);
fwrite(fid, ['$EndNodes' crlf '$Elements' crlf '1' crlf '1 2 2 1 1 1 2 3' crlf '$EndElements' crlf]);
fclose(fid);
S = surfer.load_from_msh(f);
delete(f);
assert(S.npatches == 1 && abs(sum(S.wts) - 0.5) < tol, ...
       'CRLF line ending test failed');
fprintf('CRLF line endings                           : OK\n');

%% 16. Extra unknown sections ($PhysicalNames, $Entities) before $Nodes.
f = [tempname '.msh'];
fid = fopen(f, 'w');
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');
fprintf(fid, '$PhysicalNames\n2\n');
fprintf(fid, '2 1 "outer surface"\n2 2 "inner surface"\n');
fprintf(fid, '$EndPhysicalNames\n');
fprintf(fid, '$Nodes\n3\n1 0 0 0\n2 1 0 0\n3 0 1 0\n');
fprintf(fid, '$EndNodes\n$Elements\n1\n1 2 2 1 1 1 2 3\n$EndElements\n');
fclose(fid);
S = surfer.load_from_msh(f);
delete(f);
assert(S.npatches == 1 && abs(sum(S.wts) - 0.5) < tol, ...
       'extra-section ignore test failed');
fprintf('extra unknown sections ignored              : OK\n');

%% 17. Polynomial recovery: a degree-p polynomial in (u,v) is captured
%     exactly at the surfer's target nodes for a tri at p=2 and a quad at p=2.
% Triangle at p=2 (Koornwinder spans total degree <= 2)
p = 2;
uv_g = i_build_gmsh_tri_uvs_local(p);
poly = @(u, v) u.^2 + 2*u.*v + 0.5*v.^2 + 3;     % total degree 2
xyz_g = [uv_g(1,:); uv_g(2,:); poly(uv_g(1,:), uv_g(2,:))];
f = [tempname '.msh'];
fid = fopen(f, 'w');
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n6\n');
for k = 1:6, fprintf(fid, '%d %.16g %.16g %.16g\n', k, xyz_g(:,k)); end
fprintf(fid, '$EndNodes\n$Elements\n1\n1 9 2 1 1 1 2 3 4 5 6\n$EndElements\n');
fclose(fid);
S = surfer.load_from_msh(f);
delete(f);
% At each RV node, the loaded patch's z should equal poly(uv_RV).
z_expected = poly(S.uvs_targ(1,:), S.uvs_targ(2,:));
err = max(abs(S.r(3,:) - z_expected));
assert(err < 1e-12, sprintf('tri polynomial recovery: err=%g', err));
fprintf('tri p=2 polynomial recovery (z = u^2+2uv+v^2/2+3): err=%.2g\n', err);

% Quad at p=2 (tensor Legendre spans bidegree (2,2))
p = 2;
uv_g = i_build_gmsh_quad_uvs_local(p);
poly = @(u, v) u.^2 .* v.^2 + 2*u.*v + 3;        % bidegree (2,2)
xyz_g = [uv_g(1,:); uv_g(2,:); poly(uv_g(1,:), uv_g(2,:))];
f = [tempname '.msh'];
fid = fopen(f, 'w');
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n9\n');
for k = 1:9, fprintf(fid, '%d %.16g %.16g %.16g\n', k, xyz_g(:,k)); end
fprintf(fid, '$EndNodes\n$Elements\n1\n1 10 2 1 1 1 2 3 4 5 6 7 8 9\n$EndElements\n');
fclose(fid);
S = surfer.load_from_msh(f);
delete(f);
% Surfer reference uvs on [-1,1]^2 — surface r(:,k) = (u_k, v_k, poly(u_k, v_k))
err = max(abs(S.r(3,:) - poly(S.uvs_targ(1,:), S.uvs_targ(2,:))));
assert(err < 1e-12, sprintf('quad polynomial recovery: err=%g', err));
fprintf('quad p=2 polynomial recovery (z = u^2 v^2+2uv+3): err=%.2g\n', err);

%% 18. v4.1 ASCII with two entity blocks (each one triangle), verify both load.
f = [tempname '.msh'];
fid = fopen(f, 'w');
fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
fprintf(fid, '$Nodes\n2 6 1 6\n');                    % 2 blocks, 6 nodes, tags 1..6
fprintf(fid, '2 1 0 3\n1\n2\n3\n0 0 0\n1 0 0\n0 1 0\n');    % block 1: nodes 1,2,3
fprintf(fid, '2 2 0 3\n4\n5\n6\n2 0 0\n3 0 0\n2 1 0\n');    % block 2: nodes 4,5,6
fprintf(fid, '$EndNodes\n');
fprintf(fid, '$Elements\n2 2 1 2\n');                 % 2 blocks, 2 elements, tags 1..2
fprintf(fid, '2 1 2 1\n1 1 2 3\n');                   % block 1: 1 tri
fprintf(fid, '2 2 2 1\n2 4 5 6\n');                   % block 2: 1 tri
fprintf(fid, '$EndElements\n');
fclose(fid);
S = surfer.load_from_msh(f);
delete(f);
assert(S.npatches == 2,                'v4 two-block: npatches');
assert(abs(sum(S.wts) - 1.0) < tol,    'v4 two-block: total area = 1.0 (two 0.5 tris)');
fprintf('v4.1 two entity blocks                      : OK\n');

%% 19. Stress: a 50x50 grid of order-1 triangles (5000 patches) — round-tripped.
nside = 50;
[xv, yv] = meshgrid(linspace(0,1,nside+1));
xs = xv(:); ys = yv(:); zs = zeros(numel(xs), 1);
ntri = 2 * nside * nside;
ix = reshape(1:(nside+1)*(nside+1), nside+1, nside+1);
tris = zeros(ntri, 3);
ki = 0;
for jj = 1:nside
    for ii = 1:nside
        a = ix(ii,   jj  ); b = ix(ii+1, jj  );
        c = ix(ii+1, jj+1); d = ix(ii,   jj+1);
        ki = ki + 1; tris(ki, :) = [a, b, d];
        ki = ki + 1; tris(ki, :) = [b, c, d];
    end
end
f = [tempname '.msh'];
fid = fopen(f, 'w');
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n', numel(xs));
fprintf(fid, '%d %.16g %.16g %.16g\n', [1:numel(xs); xs.'; ys.'; zs.']);
fprintf(fid, '$EndNodes\n$Elements\n%d\n', ntri);
fprintf(fid, '%d 2 2 1 1 %d %d %d\n', [1:ntri; tris.']);
fprintf(fid, '$EndElements\n');
fclose(fid);
tic; S = surfer.load_from_msh(f); tload = toc;
delete(f);
assert(S.npatches == ntri,                         'stress: npatches');
assert(all(S.norders == 1),                        'stress: orders');
assert(abs(sum(S.wts) - 1.0) < 1e-12,              'stress: area');
fprintf('stress: 50x50 grid (5000 tris)  load=%.2fs  area=%.10g  : OK\n', tload, sum(S.wts));

%% 20. Missing $MeshFormat -> defensive error.
f = [tempname '.msh'];
fid = fopen(f, 'w');
fprintf(fid, '$Nodes\n3\n1 0 0 0\n2 1 0 0\n3 0 1 0\n$EndNodes\n');
fprintf(fid, '$Elements\n1\n1 2 2 1 1 1 2 3\n$EndElements\n');
fclose(fid);
try
    surfer.load_from_msh(f);
    assert(false, 'expected an error');
catch e
    assert(contains(e.identifier, 'load_from_msh:format') ...
        || contains(e.message, '$MeshFormat'), ...
        ['unexpected error: ' e.message]);
end
delete(f);
fprintf('missing $MeshFormat raises clear error      : OK\n');


function uv = i_build_gmsh_quad_uvs_local(p)
% Mirrors i_gmsh_quad_uvs inside load_from_msh.m. Independent copy so the
% test cross-checks the reader's reference-uv table by construction.
    n = p + 1;
    s = linspace(-1, 1, n);
    uv = zeros(2, n*n);
    idx = 0; lo = 1; hi = n;
    while lo <= hi
        if lo == hi
            idx = idx + 1; uv(:, idx) = [s(lo); s(lo)];
            break
        end
        corners = [s(lo), s(hi), s(hi), s(lo); ...
                   s(lo), s(lo), s(hi), s(hi)];
        uv(:, idx+1:idx+4) = corners; idx = idx + 4;
        ne = hi - lo - 1;
        if ne > 0
            uv(:, idx+1:idx+ne) = [s(lo+1:hi-1); repmat(s(lo), 1, ne)]; idx = idx + ne;
            uv(:, idx+1:idx+ne) = [repmat(s(hi), 1, ne); s(lo+1:hi-1)]; idx = idx + ne;
            uv(:, idx+1:idx+ne) = [s(hi-1:-1:lo+1); repmat(s(hi), 1, ne)]; idx = idx + ne;
            uv(:, idx+1:idx+ne) = [repmat(s(lo), 1, ne); s(hi-1:-1:lo+1)]; idx = idx + ne;
        end
        lo = lo + 1; hi = hi - 1;
    end
end

function uv = i_build_gmsh_tri_uvs_local(p)
% Mirrors i_gmsh_tri_uvs inside load_from_msh.m. Independent copy so the
% test cross-checks the reader's reference-uv table by construction.
    uv = zeros(2, (p+1)*(p+2)/2);
    idx = 0;
    for shell = 0:floor(p/3)
        ps = p - 3*shell;
        A = [shell/p; shell/p];
        B = [(p-2*shell)/p; shell/p];
        C = [shell/p; (p-2*shell)/p];
        if ps == 0
            idx = idx + 1; uv(:, idx) = (A + B + C) / 3;
            break
        end
        uv(:, idx+1:idx+3) = [A, B, C]; idx = idx + 3;
        ne = ps - 1;
        if ne > 0
            t = (1:ne) / ps;
            uv(:, idx+1:idx+ne) = A + (B - A) * t; idx = idx + ne;
            uv(:, idx+1:idx+ne) = B + (C - B) * t; idx = idx + ne;
            uv(:, idx+1:idx+ne) = C + (A - C) * t; idx = idx + ne;
        end
    end
end
