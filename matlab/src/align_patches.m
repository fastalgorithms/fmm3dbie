function S = align_patches(S, distmin, normals, align_type, align_id)
% ALIGN_PATCHES  Rotate patch parametric coordinates so that a specified
%   edge or vertex aligns with the singular direction.
%
% Syntax:
%   S = align_patches(S, distmin, normals)
%   S = align_patches(S, distmin, normals, align_type)
%   S = align_patches(S, distmin, normals, align_type, align_id)
%
% Input arguments:
%   S          - surfer object (iptype 1, 11, or 12; may be a mixed mesh)
%   distmin    - vector of length sum(nv(i)) over all patches, where
%                nv(i) = 3 for a triangular patch and nv(i) = 4 for a quad
%                patch.  The entries for patch i are the distances from
%                its vertices to the nearest boundary feature.
%
%                Vertex ordering within each patch:
%                  iptype=1  (triangle): 3 Rokhlin-Vioreanu vertices are
%                  ordered as
%                       (0,0), (1,0), (0,1).
%                  iptype=11/12 (quad):  4 corners are ordered as
%                       (-1,-1), (+1,-1), (+1,+1), (-1,+1)
%
%   normals    - (2, npatches) or (3, npatches) a representative boundary
%                pseudo normal for each patch. Currently unused.
%   align_type - (1, npatches) integer array:
%                  101  edge alignment (triangles): rotate so that the
%                       edge opposite the farthest-from-boundary vertex
%                       becomes the leading edge. [default for iptype=1]
%                  102  vertex alignment (triangles): rotate so that the
%                       vertex closest to the boundary is the image of
%                       (0,1)
%                  201  edge alignment (quads): rotate so that the edge 
%                       containing the closest-to-boundary corner becomes
%                       the leading edge. [default for iptype=11/12]
%                Scalar input is broadcast to all patches.
%   align_id   - length-3 vector [offset101, offset102, offset201]. Offest
%                the automatic irot by this value (mod 3 for triangles, mod
%                4 for quads). Default: [0, 0, 0].
%
% Output:
%   S          - aligned surfer object

% Threshold for counting "small" distances: a vertex distance is considered
% near-zero if it is less than THRESH * max(distmin_i).
%   nsmall == 1  =>  vertex alignment (one corner touching the feature)
%   nsmall == 2  =>  edge alignment   (one edge touching the feature)
%   nsmall >= 3  =>  no rotation needed (patch is on the feature)
THRESH = 1e-7;

npatches = S.npatches;

% --- defaults ---
auto_align_type = (nargin < 4 || isempty(align_type));
if ~auto_align_type
    % broadcast scalar
    if isscalar(align_type)
        align_type = align_type * ones(1, npatches);
    end
end
if nargin < 5 || isempty(align_id)
    align_id = [0, 0, 0];
end

% --- build per-patch start indices into the flat distmin vector ---
nv_per_patch = 3 * ones(npatches, 1);
nv_per_patch(S.iptype == 11 | S.iptype == 12) = 4;
distmin_ptr = [1; cumsum(nv_per_patch) + 1];   % length npatches+1

distmin = distmin(:);   % ensure column vector

% --- allocate output arrays (same size as input) ---
npts = S.npts;
rnew  = zeros(3, npts);
dunew = zeros(3, npts);
dvnew = zeros(3, npts);
nnew  = zeros(3, npts);

% --- per-patch loop ---
for i = 1:npatches
    iinds = S.ixyzs(i):S.ixyzs(i+1)-1;

    % extract this patch's vertex distances from the flat vector
    dptr = distmin_ptr(i):distmin_ptr(i+1)-1;
    distmin_i = distmin(dptr);

    % count how many vertex distances are negligibly small relative to the
    % largest distance on this patch
    finite_dists = distmin_i(isfinite(distmin_i));
    if isempty(finite_dists) || max(finite_dists) == 0
        nsmall = 0;
    else
        nsmall = sum(finite_dists < THRESH * max(finite_dists));
    end

    % determine align_type for this patch if not supplied by caller
    iptype_i = S.iptype(i);
    is_quad  = (iptype_i == 11 || iptype_i == 12);

    if auto_align_type
        if nsmall == 1
            % exactly one vertex on the feature: vertex alignment
            align_type_i = 102 * ~is_quad + 201 * is_quad;
        elseif nsmall == 2
            % exactly one edge on the feature: edge alignment
            align_type_i = 101 * ~is_quad + 201 * is_quad;
        else
            align_type_i = 0;   % patch interior or fully on feature, skip
        end
    else
        align_type_i = align_type(i);
    end

    do_rotate = (nsmall > 0) && (nsmall < 3) && (align_type_i > 0);

    if ~do_rotate
        rnew(:,iinds)  = S.r(:,iinds);
        dunew(:,iinds) = S.du(:,iinds);
        dvnew(:,iinds) = S.dv(:,iinds);
        nnew(:,iinds)  = S.n(:,iinds);
        continue
    end

    % --- determine irot ---
    if align_type_i == 101
        % Triangle edge alignment. Vertices: 1:(0,0), 2:(1,0), 3:(0,1).
        % Rotate so the near edge (opposite the far vertex) is leading.
        far_vertex = find(distmin_i > 1e-8);
        if isempty(far_vertex)
            irot = 0;
        else
            tri_vertex_to_irot = [2, 1, 0];
            irot = tri_vertex_to_irot(far_vertex(1));
        end
        irot = mod(irot + align_id(1), 3);

    elseif align_type_i == 102
        % Triangle vertex alignment. Vertices: 1:(0,0), 2:(1,0), 3:(0,1).
        % Rotate the nearest vertex to the (0,1) corner.
        [~, k] = min(distmin_i);
        tri_vertex_to_irot = [2, 1, 0];
        irot = tri_vertex_to_irot(k);
        irot = mod(irot + align_id(2), 3);

    elseif align_type_i == 201
        % Quad alignment. 
        % Rotate (CCW, (u,v)->(-v,u)) so near edge lands at bottom {1,2}.
        small = find(distmin_i < THRESH * max(finite_dists));
        if numel(small) >= 2
            s = sort(small(1:2))';
            if     isequal(s, [1 2]),  irot = 0;
            elseif isequal(s, [1 4]),  irot = 1;
            elseif isequal(s, [3 4]),  irot = 2;
            elseif isequal(s, [2 3]),  irot = 3;
            else  % diagonal pair, fall back to nearest vertex
                [~, k] = min(distmin_i);
                corner_to_irot = [0, 3, 2, 1];
                irot = corner_to_irot(k);
            end
        else
            corner_to_irot = [0, 3, 2, 1];
            irot = corner_to_irot(small(1));
        end
        irot = mod(irot + align_id(3), 4);

    else
        irot = 0;
    end

    % --- apply rotation or copy ---
    if irot > 0
        uvs_use = S.uvs_targ(:, iinds);
        [rnew(:,iinds), dunew(:,iinds), dvnew(:,iinds)] = ...
            rotate_patch(S.srccoefs{i}, S.norders(i), iptype_i, irot, uvs_use);

        rtmp  = cross(dunew(:,iinds), dvnew(:,iinds));
        vnorm = repmat(vecnorm(rtmp, 2), [3, 1]);
        nnew(:,iinds) = rtmp ./ vnorm;
    else
        rnew(:,iinds)  = S.r(:,iinds);
        dunew(:,iinds) = S.du(:,iinds);
        dvnew(:,iinds) = S.dv(:,iinds);
        nnew(:,iinds)  = S.n(:,iinds);
    end
end

% --- reconstruct surfer with same iptype and norders ---
srcvals = [rnew; dunew; dvnew; nnew];
S = surfer(npatches, S.norders, srcvals, S.iptype);

end


function [r,du,dv] = rotate_patch(srccoefs, norder, iptype, irot, uvstmp)
% ROTATE_PATCH  Re-evaluate patch geometry after a parametric rotation.
%
% For triangular patches (iptype=1), irot in {1,2} applies the two
% non-trivial cyclic permutations of the Rokhlin-Vioreanu triangle.
%
% For quad patches (iptype=11/12), irot in {1,2,3} applies 90-, 180-,
% 270-degree CCW rotations of [-1,1]^2, i.e. (u,v) -> (-v,u)^irot.

if iptype == 1
    % --- triangular patch ---
    if irot == 1
        uvs = [uvstmp(2,:); 1-uvstmp(1,:)-uvstmp(2,:)];
        pols = koorn.pols(norder, uvs).';
        r  = srccoefs(1:3,:)*pols';
        dv = srccoefs(4:6,:)*pols' - srccoefs(7:9,:)*pols';
        du = -srccoefs(7:9,:)*pols';
    elseif irot == 2
        uvs = [1-uvstmp(1,:)-uvstmp(2,:); uvstmp(1,:)];
        pols = koorn.pols(norder, uvs).';
        r  = srccoefs(1:3,:)*pols';
        du = -srccoefs(4:6,:)*pols' + srccoefs(7:9,:)*pols';
        dv = -srccoefs(4:6,:)*pols';
    else
        uvs = uvstmp;
        pols = koorn.pols(norder, uvs).';
        r  = srccoefs(1:3,:)*pols';
        du = srccoefs(4:6,:)*pols';
        dv = srccoefs(7:9,:)*pols';
    end

else
    % --- quad patch (iptype 11 or 12): rotations of [-1,1]^2 ---
    %   irot=1: u_new=-v, v_new= u  => du_new=-dv_old, dv_new= du_old
    %   irot=2: u_new=-u, v_new=-v  => du_new=-du_old, dv_new=-dv_old
    %   irot=3: u_new= v, v_new=-u  => du_new= dv_old, dv_new=-du_old

    % apply coordinate transform to the sample points
    uvs = uvstmp;
    for k = 1:irot
        uvs = [-uvs(2,:); uvs(1,:)];
    end

    if iptype == 11
        pols = polytens.lege.pols(norder, uvs).';
    else
        pols = polytens.cheb.pols(norder, uvs).';
    end

    r   = srccoefs(1:3,:)*pols';
    du_orig = srccoefs(4:6,:)*pols';
    dv_orig = srccoefs(7:9,:)*pols';

    % rotate the tangent vectors by the same irot quarter-turns
    du = du_orig;
    dv = dv_orig;
    for k = 1:irot
        du_tmp = -dv;
        dv     =  du;
        du     =  du_tmp;
    end
end

end
