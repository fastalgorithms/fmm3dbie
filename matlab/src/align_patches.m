function S = align_patches(S, distmin, normals, align_type, align_id)
% ALIGN_PATCHES  Rotate patch parametric coordinates so that a specified
%   edge or vertex aligns with the singular direction, without changing
%   iptype or expansion orders.
%
% Syntax:
%   S = align_patches(S, distmin, normals)
%   S = align_patches(S, distmin, normals, align_type)
%   S = align_patches(S, distmin, normals, align_type, align_id)
%
% Input arguments:
%   S          - surfer object
%   distmin    - (3, npatches) array of distances from each patch vertex
%                to the nearest boundary feature. Patches for which all
%                entries of distmin(:,i) are non-finite are left unchanged.
%   normals    - (2, npatches) or (3, npatches) boundary-normal direction
%                array used when align_type == 102. If (2, npatches), a
%                third row of zeros is appended before projecting.
%                May be [] if align_type is always 101.
%   align_type - (1, npatches) integer array:
%                  101  edge alignment: rotate so that the edge opposite
%                       the farthest-from-boundary vertex becomes the
%                       leading edge (same logic as type 101 in
%                       convert_to_edge). [default]
%                  102  vertex alignment: rotate so that the vertex closest
%                       to the boundary feature becomes vertex 0 of the
%                       parametric triangle (same logic as type 102 in
%                       convert_to_edge).
%                Scalar input is broadcast to all patches.
%   align_id   - length-2 vector [offset101, offset102]. After the
%                automatic irot is computed, offset101 is added (mod 3)
%                for type-101 patches and offset102 for type-102 patches.
%                Default: [0, 0].
%
% Output:
%   S          - surfer object with the same iptype and norders as input,
%                but with node positions, derivatives, and normals
%                re-evaluated after the parametric rotation.

npatches = S.npatches;

% --- defaults ---
if nargin < 4 || isempty(align_type)
    align_type = 101 * ones(1, npatches);
end
if nargin < 5 || isempty(align_id)
    align_id = [0, 0];
end

% broadcast scalars
if isscalar(align_type)
    align_type = align_type * ones(1, npatches);
end

% --- allocate output arrays (same size as input) ---
npts = S.npts;
rnew  = zeros(3, npts);
dunew = zeros(3, npts);
dvnew = zeros(3, npts);
nnew  = zeros(3, npts);

% --- per-patch loop ---
for i = 1:npatches
    iinds = S.ixyzs(i):S.ixyzs(i+1)-1;

    % Only rotate if distmin is finite and positive for this patch
    do_rotate = any(isfinite(distmin(:,i)) & distmin(:,i) > 0);

    % Only triangular (iptype==1) patches support parametric rotation
    if ~do_rotate || S.iptype(i) ~= 1
        rnew(:,iinds)  = S.r(:,iinds);
        dunew(:,iinds) = S.du(:,iinds);
        dvnew(:,iinds) = S.dv(:,iinds);
        nnew(:,iinds)  = S.n(:,iinds);
        continue
    end

    % --- determine irot ---
    if align_type(i) == 101
        % Edge alignment: rotate so the edge opposite the farthest vertex
        % from the boundary becomes the leading edge.
        distmin_i = distmin(:, i);
        irot = mod(find(distmin_i > 1e-8) + 1, 3);
        if isempty(irot)
            irot = 0;
        else
            irot = irot(1);
        end
        irot = mod(irot + align_id(1), 3);

    elseif align_type(i) == 102
        % Vertex alignment: find the vertex closest to the boundary.
        distmin_i = distmin(:, i);
        [~, k] = min(distmin_i);
        irot = mod(k+1, 3);
        irot = mod(irot + align_id(2), 3);

    else
        irot = 0;
    end

    % --- apply rotation or copy ---
    if irot > 0
        uvs_use = S.uvs_targ(:, iinds);
        [rnew(:,iinds), dunew(:,iinds), dvnew(:,iinds)] = ...
            rotate_patch(S.srccoefs{i}, S.norders(i), irot, uvs_use);

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


function [r,du,dv,uvs] = rotate_patch(srccoefs,norder,irot,uvstmp)

if irot == 1
    uvs = [uvstmp(2,:);1-uvstmp(1,:)-uvstmp(2,:)];

    pols = koorn.pols(norder,uvs).';

    r = srccoefs(1:3,:)*pols';
    dv = srccoefs(4:6,:)*pols'-srccoefs(7:9,:)*pols';
    du = -srccoefs(7:9,:)*pols';
elseif irot == -1 || irot ==2
    uvs = [1-uvstmp(1,:)-uvstmp(2,:);uvstmp(1,:)];

    pols = koorn.pols(norder,uvs).';

    r = srccoefs(1:3,:)*pols';
    du = -srccoefs(4:6,:)*pols'+srccoefs(7:9,:)*pols';
    dv = -srccoefs(4:6,:)*pols';
else
    uvs = uvstmp;
    pols = koorn.pols(norder,uvs).';

    r = srccoefs(1:3,:)*pols';
    du = srccoefs(4:6,:)*pols';
    dv = srccoefs(7:9,:)*pols';
end

end
