function obj = load_from_msh(fname, opts)
%LOAD_FROM_MSH  Create a surfer object from a Gmsh .msh file.
%
% S = surfer.load_from_msh(fname) reads a Gmsh mesh (format version 2.x
%  or 4.x, ASCII or binary, little- or big-endian) and builds a surfer
%  object with one patch per supported triangle or quadrilateral element.
%
% S = surfer.load_from_msh(fname, opts) accepts:
%    opts.iptype = 11 (Gauss-Legendre tensor product, default)
%                or 12 (Chebyshev tensor product) for quad patches.
%  Triangles always use Rokhlin-Vioreanu nodes (iptype = 1).
%
% Supported gmsh element type codes:
%    triangles  : 2 (tri3),  9 (tri6), 21 (tri10), 23 (tri15)
%    quads      : 3 (Q4),   10 (Q9),  16 (Q8 serendipity),
%                36 (Q16), 37 (Q25)
%  Q8 is projected to Q9 by reconstructing the centre node:
%    centre = -1/4 sum(corners) + 1/2 sum(edge_midpoints).
%  Other element types (lines, tetrahedra, points, ...) are ignored.
%  Mixed-shape and mixed-order meshes are supported.
%
% Algorithm: xyz at each element's gmsh reference nodes is converted to
% coefficients in the natural polynomial basis on the reference element
% (Koornwinder for triangles, tensor-product Legendre/Chebyshev for
% quads), then resampled at the surfer's target nodes. Tangents du,dv
% are obtained from the same expansion by spectral differentiation; the
% unit normal is normalize(du x dv).
%
% See also: surfer.load_from_file, surfer

    if nargin < 2 || isempty(opts), opts = struct(); end
    iptype_q = 11;
    if isfield(opts, 'iptype'), iptype_q = opts.iptype; end
    if iptype_q ~= 11 && iptype_q ~= 12
        error('surfer:load_from_msh:iptype', ...
              'opts.iptype must be 11 (Gauss-Legendre) or 12 (Chebyshev)');
    end

    [nodes, eshapes, eorders, enodes] = i_parse_msh(fname);

    npatches = numel(eorders);
    if npatches == 0
        error('surfer:load_from_msh:no_supported_elements', ...
              'No supported elements (tri 2/9/21/23 or quad 3/10/16/36/37) found in %s', fname);
    end

    % Precompute per-(shape,order) resampling operators (small, npols x npols).
    Mg = struct('T', {cell(4,1)}, 'Q', {cell(4,1)});
    Pt = struct('T', {cell(4,1)}, 'Q', {cell(4,1)});
    Du = struct('T', {cell(4,1)}, 'Q', {cell(4,1)});
    Dv = struct('T', {cell(4,1)}, 'Q', {cell(4,1)});
    for s = unique(eshapes)
        for p = unique(eorders(eshapes == s))
            if s == 'T', uv_gmsh = i_gmsh_tri_uvs(p);
            else,        uv_gmsh = i_gmsh_quad_uvs(p);
            end
            uv_tgt    = i_target_nodes(s, iptype_q, p);
            Mg.(s){p} = i_vals2coefs(s, iptype_q, p, uv_gmsh);
            [Pt.(s){p}, Du.(s){p}, Dv.(s){p}] = i_ders(s, iptype_q, p, uv_tgt);
        end
    end

    % Per-patch point counts and ranges in srcvals.
    nps = zeros(1, npatches);
    for k = 1:npatches
        if eshapes(k) == 'T', nps(k) = (eorders(k)+1)*(eorders(k)+2)/2;
        else,                 nps(k) = (eorders(k)+1)^2;
        end
    end
    offs    = [0, cumsum(nps)];
    srcvals = zeros(12, offs(end));
    iptypes = ones(npatches, 1);

    for k = 1:npatches
        s = eshapes(k);  p = eorders(k);
        gxyz = nodes(:, enodes{k});                       % 3 x ngmsh
        if s == 'Q' && p == 2 && size(gxyz,2) == 8
            % Q8 serendipity -> Q9: centre = -1/4 sum(corners) + 1/2 sum(edge_mids).
            gxyz(:, 9) = -0.25*sum(gxyz(:,1:4),2) + 0.5*sum(gxyz(:,5:8),2);
        end
        coefs = gxyz * Mg.(s){p}.';                       % 3 x np basis coefs
        rng   = offs(k)+1 : offs(k+1);
        du    = coefs * Du.(s){p};
        dv    = coefs * Dv.(s){p};
        nrm   = cross(du, dv, 1);
        nrm   = nrm ./ repmat(sqrt(sum(nrm.*nrm, 1)), [3, 1]);
        srcvals(1:3,   rng) = coefs * Pt.(s){p};
        srcvals(4:6,   rng) = du;
        srcvals(7:9,   rng) = dv;
        srcvals(10:12, rng) = nrm;
        if s == 'Q', iptypes(k) = iptype_q; end
    end

    obj = surfer(npatches, eorders(:), srcvals, iptypes);
end

% ------------------------------------------------------------------------
function [nodes, eshapes, eorders, enodes] = i_parse_msh(fname)
% Read $Nodes and $Elements from a gmsh ASCII or binary file (v2.x or v4.x).
%   nodes   : 3 x maxNodeTag (indexed by gmsh node tag)
%   eshapes : 1 x ne char vector, 'T' (triangle) or 'Q' (quad), per element
%   eorders : 1 x ne polynomial orders, per element
%   enodes  : 1 x ne cell of gmsh node-tag vectors per element

    fid = fopen(fname, 'rb');                   % binary mode; fgetl still works
    if fid < 0
        error('surfer:load_from_msh:fopen', 'cannot open %s', fname);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    nodes   = zeros(3, 0);
    eshapes = char(zeros(1, 0));
    eorders = zeros(1, 0);
    enodes  = cell(1, 0);
    major   = 0;
    is_v40  = false;
    is_bin  = false;
    en      = 'ieee-le';

    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        line = strtrim(line);
        switch line
            case '$MeshFormat'
                fmt = strsplit(strtrim(fgetl(fid)));
                ver = sscanf(fmt{1}, '%g');
                is_bin = (numel(fmt) >= 2) && strcmp(fmt{2}, '1');
                major = floor(ver);
                is_v40 = (major == 4) && (ver < 4.05);
                if major ~= 2 && major ~= 4
                    error('surfer:load_from_msh:version', ...
                          'unsupported .msh major version %g (need 2.x or 4.x)', ver);
                end
                if is_bin
                    % Endian probe: int32 written as 1 in native byte order
                    probe = fread(fid, 1, 'int32=>int32', 0, 'ieee-le');
                    if probe ~= 1
                        fseek(fid, -4, 'cof');
                        probe = fread(fid, 1, 'int32=>int32', 0, 'ieee-be');
                        if probe ~= 1
                            error('surfer:load_from_msh:endian', ...
                                  'binary endianness probe failed');
                        end
                        en = 'ieee-be';
                    end
                end
            case '$Nodes'
                if major == 0
                    error('surfer:load_from_msh:format', ...
                          '$MeshFormat must precede $Nodes');
                end
                if     is_bin && major == 2,  nodes = i_read_nodes_v2_bin (fid, en);
                elseif is_bin && is_v40,      nodes = i_read_nodes_v40_bin(fid, en);
                elseif is_bin,                nodes = i_read_nodes_v41_bin(fid, en);
                elseif           major == 2,  nodes = i_read_nodes_v2 (fid);
                elseif           is_v40,      nodes = i_read_nodes_v40(fid);
                else,                         nodes = i_read_nodes_v41(fid);
                end
            case '$Elements'
                if major == 0
                    error('surfer:load_from_msh:format', ...
                          '$MeshFormat must precede $Elements');
                end
                if     is_bin && major == 2,  [eshapes, eorders, enodes] = i_read_elems_v2_bin(fid, en);
                elseif is_bin,                [eshapes, eorders, enodes] = i_read_elems_v4_bin(fid, en, is_v40);
                elseif           major == 2,  [eshapes, eorders, enodes] = i_read_elems_v2(fid);
                else,                         [eshapes, eorders, enodes] = i_read_elems_v4(fid, is_v40);
                end
        end
    end
end

% ------------------------------------------------------------------------
function nodes = i_read_nodes_v2(fid)
% v2 layout: numnodes, then numnodes lines each "tag x y z".
    numnodes = sscanf(fgetl(fid), '%d', 1);
    C = textscan(fid, '%d %f %f %f', numnodes);
    ids = double(C{1});
    nodes = zeros(3, max(ids));
    nodes(:, ids) = [C{2}, C{3}, C{4}].';
end

% ------------------------------------------------------------------------
function [eshapes, eorders, enodes] = i_read_elems_v2(fid)
% v2 layout: numelem, then per element "id type ntags tag1..tagN node1..nodeM".
    numelem = sscanf(fgetl(fid), '%d', 1);
    eshapes = char(zeros(1, numelem));
    eorders = zeros(1, numelem);
    enodes  = cell(1, numelem);
    ke = 0;
    for i = 1:numelem
        row = sscanf(fgetl(fid), '%d').';
        [s, ord, np] = i_elem_type(row(2));
        if s == ' ', continue; end
        ntags = row(3);
        ke = ke + 1;
        eshapes(ke) = s;
        eorders(ke) = ord;
        enodes{ke}  = row((4+ntags):(3+ntags+np));
    end
    eshapes = eshapes(1:ke);
    eorders = eorders(1:ke);
    enodes  = enodes(1:ke);
end

% ------------------------------------------------------------------------
function nodes = i_read_nodes_v41(fid)
% v4.1 layout:
%   numEntityBlocks numNodes minNodeTag maxNodeTag
%   per block: entityDim entityTag parametric numNodesInBlock
%              tag_1 ... tag_M  (then) x_1 y_1 z_1 ... x_M y_M z_M
    hdr = fscanf(fid, '%d', 4);
    nBlocks = hdr(1); maxTag = hdr(4);
    nodes = zeros(3, max(maxTag, 1));
    for b = 1:nBlocks
        bhdr = fscanf(fid, '%d', 4);
        if bhdr(3) ~= 0
            error('surfer:load_from_msh:parametric', ...
                  'parametric nodes not supported');
        end
        m = bhdr(4);
        ids = fscanf(fid, '%d', m);
        xyz = reshape(fscanf(fid, '%g', 3*m), 3, m);
        nodes(:, ids) = xyz;
    end
end

% ------------------------------------------------------------------------
function nodes = i_read_nodes_v40(fid)
% v4.0 layout (differs from v4.1):
%   numEntityBlocks numNodes                                 (2 ints)
%   per block: entityTag entityDim parametric numNodesInBlock
%              per node: tag x y z   (interleaved one per line)
    hdr = fscanf(fid, '%d', 2);
    nBlocks = hdr(1); numNodes = hdr(2);
    tags = zeros(1, numNodes);
    xyz  = zeros(3, numNodes);
    idx = 0;
    for b = 1:nBlocks
        bhdr = fscanf(fid, '%d', 4);
        if bhdr(3) ~= 0
            error('surfer:load_from_msh:parametric', ...
                  'parametric nodes not supported');
        end
        m = bhdr(4);
        data = reshape(fscanf(fid, '%g', 4*m), 4, m);
        rng = (idx+1):(idx+m);
        tags(rng)  = data(1, :);
        xyz(:, rng) = data(2:4, :);
        idx = idx + m;
    end
    nodes = zeros(3, max(tags));
    nodes(:, tags) = xyz;
end

% ------------------------------------------------------------------------
function [eshapes, eorders, enodes] = i_read_elems_v4(fid, is_v40)
% v4.0/v4.1 element layout (differ only in section header width).
    if is_v40
        hdr = fscanf(fid, '%d', 2);
    else
        hdr = fscanf(fid, '%d', 4);
    end
    nBlocks = hdr(1); nElems = hdr(2);
    eshapes = char(zeros(1, nElems));
    eorders = zeros(1, nElems);
    enodes  = cell(1, nElems);
    ke = 0;
    for b = 1:nBlocks
        bhdr = fscanf(fid, '%d', 4);
        nb = bhdr(4);
        [s, ord, np] = i_elem_type(bhdr(3));
        if s == ' '
            fgetl(fid); % finish current line started by last fscanf
            for j = 1:nb, fgetl(fid); end
            continue
        end
        for j = 1:nb
            row = fscanf(fid, '%d', 1+np);
            ke = ke + 1;
            eshapes(ke) = s;
            eorders(ke) = ord;
            enodes{ke}  = row(2:end).';
        end
    end
    eshapes = eshapes(1:ke);
    eorders = eorders(1:ke);
    enodes  = enodes(1:ke);
end

% ------------------------------------------------------------------------
function [s, ord, np] = i_elem_type(etype)
% Map gmsh element type code to (shape char, polynomial order, #nodes).
% Full-Lagrange triangles and quads of orders 1..10. Returns (' ',0,0)
% for unsupported types ("incomplete" / serendipity variants other than
% Q8 are not handled).
    switch etype
        case 2,  s = 'T'; ord = 1;  np = 3;
        case 9,  s = 'T'; ord = 2;  np = 6;
        case 21, s = 'T'; ord = 3;  np = 10;
        case 23, s = 'T'; ord = 4;  np = 15;
        case 25, s = 'T'; ord = 5;  np = 21;
        case 42, s = 'T'; ord = 6;  np = 28;
        case 43, s = 'T'; ord = 7;  np = 36;
        case 44, s = 'T'; ord = 8;  np = 45;
        case 45, s = 'T'; ord = 9;  np = 55;
        case 46, s = 'T'; ord = 10; np = 66;
        case 3,  s = 'Q'; ord = 1;  np = 4;     % Q4
        case 10, s = 'Q'; ord = 2;  np = 9;     % Q9
        case 16, s = 'Q'; ord = 2;  np = 8;     % Q8 serendipity
        case 36, s = 'Q'; ord = 3;  np = 16;
        case 37, s = 'Q'; ord = 4;  np = 25;
        case 38, s = 'Q'; ord = 5;  np = 36;
        case 47, s = 'Q'; ord = 6;  np = 49;
        case 48, s = 'Q'; ord = 7;  np = 64;
        case 49, s = 'Q'; ord = 8;  np = 81;
        case 50, s = 'Q'; ord = 9;  np = 100;
        case 51, s = 'Q'; ord = 10; np = 121;
        otherwise, s = ' '; ord = 0; np = 0;
    end
end

% ------------------------------------------------------------------------
function uv = i_gmsh_tri_uvs(p)
% Reference (u,v) of gmsh's order-p triangular nodes on the standard
% simplex (0,0)/(1,0)/(0,1), shell-by-shell from the outside in: each
% shell adds 3 corners + (p_s-1) interior nodes along each of 3 edges
% (CCW), then recurses into the inner triangle of order p_s = p - 3s.
% When p_s = 0 the recursion ends with a single centroid. The p=1..4
% cases agree with load_uv_gmsh2_all in src/surface_routs/in_gmsh2.f90.
    if p < 1 || p > 10
        error('surfer:load_from_msh:order', ...
              'unsupported gmsh triangle order %d (must be 1..10)', p);
    end
    uv = zeros(2, (p+1)*(p+2)/2);
    idx = 0;
    for s = 0:floor(p/3)
        ps = p - 3*s;
        A = [s/p; s/p];
        B = [(p-2*s)/p; s/p];
        C = [s/p; (p-2*s)/p];
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

function uv = i_gmsh_quad_uvs(p)
% Gmsh order-p quadrilateral reference nodes on [-1,1]^2, in gmsh's
% shell traversal: outer shell first (4 corners CCW from (-1,-1), then
% p-1 interior nodes along each of the 4 edges, CCW), then recursively
% the inner shell on the (p-2)x(p-2) sub-square, etc. For odd inner
% sizes the recursion ends at a single centre point.
    if p < 1 || p > 10
        error('surfer:load_from_msh:order', ...
              'unsupported gmsh quad order %d (must be 1..10)', p);
    end
    n = p + 1;
    s = linspace(-1, 1, n);
    uv = zeros(2, n*n);
    idx = 0;
    lo = 1; hi = n;
    while lo <= hi
        if lo == hi
            idx = idx + 1; uv(:, idx) = [s(lo); s(lo)];
            break
        end
        % four corners CCW from (s(lo), s(lo))
        corners = [s(lo), s(hi), s(hi), s(lo); ...
                   s(lo), s(lo), s(hi), s(hi)];
        uv(:, idx+1:idx+4) = corners; idx = idx + 4;
        ne = hi - lo - 1;  % # interior nodes per edge of this shell
        if ne > 0
            % bottom edge L->R
            uv(:, idx+1:idx+ne) = [s(lo+1:hi-1); repmat(s(lo), 1, ne)]; idx = idx + ne;
            % right edge B->T
            uv(:, idx+1:idx+ne) = [repmat(s(hi), 1, ne); s(lo+1:hi-1)]; idx = idx + ne;
            % top edge R->L
            uv(:, idx+1:idx+ne) = [s(hi-1:-1:lo+1); repmat(s(hi), 1, ne)]; idx = idx + ne;
            % left edge T->B
            uv(:, idx+1:idx+ne) = [repmat(s(lo), 1, ne); s(hi-1:-1:lo+1)]; idx = idx + ne;
        end
        lo = lo + 1; hi = hi - 1;
    end
end

% ------------------------------------------------------------------------
% Basis dispatch: triangles use Koornwinder on the std simplex; quads use
% tensor-product Legendre (iptype=11) or Chebyshev (iptype=12) on [-1,1]^2.

function uv = i_target_nodes(shape, iptype_q, p)
    if shape == 'T'
        uv = koorn.rv_nodes(p);
    elseif iptype_q == 11
        uv = polytens.lege.nodes(p);
    else
        uv = polytens.cheb.nodes(p);
    end
end

function M = i_vals2coefs(shape, iptype_q, p, uvs)
    if shape == 'T'
        M = koorn.vals2coefs(p, uvs);
    elseif iptype_q == 11
        M = polytens.lege.vals2coefs(p, uvs);
    else
        M = polytens.cheb.vals2coefs(p, uvs);
    end
end

function [P, Du, Dv] = i_ders(shape, iptype_q, p, uvs)
    if shape == 'T'
        [P, Du, Dv] = koorn.ders(p, uvs);
    elseif iptype_q == 11
        [P, Du, Dv] = polytens.lege.ders(p, uvs);
    else
        [P, Du, Dv] = polytens.cheb.ders(p, uvs);
    end
end

% ========================================================================
% Binary readers (v2.x / v4.0 / v4.1).  The int32=1 probe in $MeshFormat
% sets `en`; fread's machinefmt then handles byte order for single-type
% reads. For the int32+double[3] interleaved records used by v2 and v4.0
% nodes we read raw bytes and byte-reverse blocks before typecast.
% ========================================================================

function [tags, xyz] = i_read_records_tag_xyz(fid, n, en)
% n records of [int32 tag, double[3] xyz] = 28 bytes each.
    raw = fread(fid, 28*n, 'uint8=>uint8');
    raw = reshape(raw, 28, n);
    tagb = raw(1:4, :);
    xyzb = reshape(raw(5:28, :), 8, 3*n);
    if i_must_swap(en)
        tagb = flipud(tagb);
        xyzb = flipud(xyzb);
    end
    tags = double(typecast(tagb(:), 'int32'));
    xyz  = reshape(typecast(xyzb(:), 'double'), 3, n);
end

function tf = i_must_swap(en)
    [~, ~, host] = computer;            % 'L' little, 'B' big
    tf = (host == 'L' && strcmp(en, 'ieee-be')) || ...
         (host == 'B' && strcmp(en, 'ieee-le'));
end

% ------------------------------------------------------------------------
function nodes = i_read_nodes_v2_bin(fid, en)
% v2 binary nodes: ASCII numnodes line, then n records of (int32 + 3 doubles).
    numnodes = sscanf(fgetl(fid), '%d', 1);
    [tags, xyz] = i_read_records_tag_xyz(fid, numnodes, en);
    nodes = zeros(3, max(tags));
    nodes(:, tags) = xyz;
end

% ------------------------------------------------------------------------
function [eshapes, eorders, enodes] = i_read_elems_v2_bin(fid, en)
% v2 binary elements: ASCII numelem line, then chunks of
%   3*int32 header [elemType numOfType ntagsPerElem]
%   numOfType records of int32[1 + ntagsPerElem + nodesPerType]
    numelem = sscanf(fgetl(fid), '%d', 1);
    eshapes = char(zeros(1, numelem));
    eorders = zeros(1, numelem);
    enodes  = cell(1, numelem);
    ke    = 0;
    seen  = 0;
    while seen < numelem
        hdr = fread(fid, 3, 'int32=>double', 0, en);
        etype = hdr(1); nt = hdr(2); ntags = hdr(3);
        [sh, ord, np] = i_elem_type(etype);
        if sh == ' '
            % Unsupported type: still consume its bytes so we stay aligned.
            nptt = i_nodes_per_etype(etype);
            if nptt < 0
                error('surfer:load_from_msh:elem_type', ...
                      'cannot skip unsupported binary v2 elem type %d', etype);
            end
            fread(fid, nt * (1 + ntags + nptt), 'int32=>int32', 0, en);
            seen = seen + nt;
            continue
        end
        rec = fread(fid, nt * (1 + ntags + np), 'int32=>double', 0, en);
        rec = reshape(rec, 1 + ntags + np, nt);
        nrows = rec(2+ntags:end, :);    % node indices per element
        for j = 1:nt
            ke = ke + 1;
            eshapes(ke) = sh;
            eorders(ke) = ord;
            enodes{ke}  = nrows(:, j).';
        end
        seen = seen + nt;
    end
    eshapes = eshapes(1:ke);
    eorders = eorders(1:ke);
    enodes  = enodes(1:ke);
end

function n = i_nodes_per_etype(t)
% Nodes-per-element for gmsh types we may need to skip past in v2 binary
% blocks. -1 if unknown. Only types likely to appear in surface meshes.
    switch t
        case 1,  n = 2;     % line2
        case 4,  n = 4;     % tet4
        case 5,  n = 8;     % hex8
        case 6,  n = 6;     % prism6
        case 7,  n = 5;     % pyramid5
        case 8,  n = 3;     % line3
        case 11, n = 10;    % tet10
        case 12, n = 27;    % hex27
        case 15, n = 1;     % point
        case 17, n = 20;    % hex20
        otherwise
            % triangles / quads currently supported in i_elem_type
            [~, ~, n] = i_elem_type(t);
            if n == 0, n = -1; end
    end
end

% ------------------------------------------------------------------------
function nodes = i_read_nodes_v40_bin(fid, en)
% v4.0 binary nodes: 2 c_ulong header, per block (3 int32 + 1 c_ulong) +
% per-node interleaved (int32 + 3 doubles). v4.0 lacks maxTag in the
% section header, so `nodes` grows dynamically as larger tags appear.
    hdr = fread(fid, 2, 'uint64=>double', 0, en);
    nBlocks = hdr(1);
    nodes = zeros(3, 0);
    for b = 1:nBlocks
        bh = fread(fid, 3, 'int32=>double', 0, en);
        if bh(3) ~= 0
            error('surfer:load_from_msh:parametric', ...
                  'parametric nodes not supported');
        end
        m = fread(fid, 1, 'uint64=>double', 0, en);
        [tags, xyz] = i_read_records_tag_xyz(fid, m, en);
        mt = max(tags);
        if mt > size(nodes, 2), nodes(:, mt) = 0; end    % zero-extend
        nodes(:, tags) = xyz;
    end
end

% ------------------------------------------------------------------------
function nodes = i_read_nodes_v41_bin(fid, en)
% v4.1 binary nodes: 4 size_t header, per block (3 int32 + 1 size_t)
% + size_t[m] tags + double[3*m] coords.
    hdr = fread(fid, 4, 'uint64=>double', 0, en);
    nBlocks = hdr(1); maxTag = hdr(4);
    nodes = zeros(3, max(maxTag, 1));
    for b = 1:nBlocks
        bh = fread(fid, 3, 'int32=>double', 0, en);
        if bh(3) ~= 0
            error('surfer:load_from_msh:parametric', ...
                  'parametric nodes not supported');
        end
        m = fread(fid, 1, 'uint64=>double', 0, en);
        tags = fread(fid, m, 'uint64=>double', 0, en);
        xyz  = reshape(fread(fid, 3*m, 'double=>double', 0, en), 3, m);
        nodes(:, tags) = xyz;
    end
end

% ------------------------------------------------------------------------
function [eshapes, eorders, enodes] = i_read_elems_v4_bin(fid, en, is_v40)
% v4 binary elements: header (2 or 4 size_t), per block (3 int32 + 1
% size_t), then per-element records as int32 (v4.0) or size_t (v4.1)
% of width 1 + nodesPerType.
    if is_v40
        nHdr = 2;
        elemPrec = 'int32=>double';
    else
        nHdr = 4;
        elemPrec = 'uint64=>double';
    end
    hdr = fread(fid, nHdr, 'uint64=>double', 0, en);
    nBlocks = hdr(1); nElems = hdr(2);
    eshapes = char(zeros(1, nElems));
    eorders = zeros(1, nElems);
    enodes  = cell(1, nElems);
    ke = 0;
    for b = 1:nBlocks
        bh = fread(fid, 3, 'int32=>double', 0, en);
        m  = fread(fid, 1, 'uint64=>double', 0, en);
        [sh, ord, np] = i_elem_type(bh(3));
        if sh == ' '
            nptt = i_nodes_per_etype(bh(3));
            if nptt < 0
                error('surfer:load_from_msh:elem_type', ...
                      'cannot skip unsupported binary v4 elem type %d', bh(3));
            end
            fread(fid, m * (1+nptt), elemPrec, 0, en);
            continue
        end
        rec = reshape(fread(fid, m * (1+np), elemPrec, 0, en), 1+np, m);
        for j = 1:m
            ke = ke + 1;
            eshapes(ke) = sh;
            eorders(ke) = ord;
            enodes{ke}  = rec(2:end, j).';
        end
    end
    eshapes = eshapes(1:ke);
    eorders = eorders(1:ke);
    enodes  = enodes(1:ke);
end
