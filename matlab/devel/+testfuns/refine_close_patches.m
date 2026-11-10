function [S1refine,S2refine,info] = refine_close_patches(S1,S2,prefactor,maxlevels,max_npts)
%REFINE_CLOSE_PATCHES Refine patches whose opposing centroid is nearby.
%
% Refine patch i when
%
%   nearest opposing-centroid distance < prefactor*sqrt(patch area).
%

if nargin < 3 || isempty(prefactor)
    error('prefactor must be supplied.')
end

if nargin < 4 || isempty(maxlevels)
    maxlevels = 20;
end

if ~isscalar(prefactor) || ~isfinite(prefactor) || prefactor <= 0
    error('prefactor must be a positive finite scalar.')
end

if ~isscalar(maxlevels) || ~isfinite(maxlevels) || ...
        maxlevels < 1 || maxlevels ~= floor(maxlevels)
    error('maxlevels must be a positive integer.')
end

if nargin < 5, max_npts = Inf; end

S1refine = S1;
S2refine = S2;
nlevels = 0;

while true
    np1 = double(S1refine.npatches);
    np2 = double(S2refine.npatches);

    A1 = accumarray(double(S1refine.patch_id(:)), ...
        S1refine.wts(:),[np1,1]);

    A2 = accumarray(double(S2refine.patch_id(:)), ...
        S2refine.wts(:),[np2,1]);

    h1 = sqrt(A1);
    h2 = sqrt(A2);

    % Bound temporary distance matrices independently of the total patch count.
    d1 = inf(np1,1);
    d2 = inf(np2,1);
    for first = 1:256:np1
        idx = first:min(first+255,np1);
        dx = S1refine.cms(1,idx).' - S2refine.cms(1,:);
        dy = S1refine.cms(2,idx).' - S2refine.cms(2,:);
        dz = S1refine.cms(3,idx).' - S2refine.cms(3,:);
        D = sqrt(dx.^2 + dy.^2 + dz.^2);
        d1(idx) = min(D,[],2);
        d2 = min(d2,min(D,[],1).');
    end

    refinelist1 = find(d1 < prefactor*h1);
    refinelist2 = find(d2 < prefactor*h2);

    if isempty(refinelist1) && isempty(refinelist2)
        converged = true;
        break
    end

    if nlevels >= maxlevels
        converged = false;
        warning('refine_close_patches:maxlevels', ...
            'Maximum refinement level reached.')
        break
    end

    points1 = diff(S1refine.ixyzs(:));
    points2 = diff(S2refine.ixyzs(:));
    next_npts = S1refine.npts+S2refine.npts + ...
        3*(sum(points1(refinelist1))+sum(points2(refinelist2)));
    if next_npts > max_npts
        error('refine_close_patches:PointLimit', ...
            'Refinement needs %d points (limit %d).',next_npts,max_npts);
    end

    if ~isempty(refinelist1)
        S1refine = split_patches(S1refine,refinelist1);
    end

    if ~isempty(refinelist2)
        S2refine = split_patches(S2refine,refinelist2);
    end

    nlevels = nlevels + 1;
end

info.prefactor = prefactor;
info.nlevels = nlevels;
info.converged = converged;
info.min_centroid_distance = min(d1);
info.npatches1 = S1refine.npatches;
info.npatches2 = S2refine.npatches;

end


function [S1refine,S2refine,info] = refine_close_patches33(S1,S2,theta,maxlevels)
%REFINE_CLOSE_PATCHES Refine two nearby surfaces until their patches are small.
%
%   [S1REFINE,S2REFINE,INFO] = REFINE_CLOSE_PATCHES(S1,S2)
%   repeatedly splits patches until
%
%       sqrt(patch area) < theta * nearest opposing-patch centroid distance
%
%   on both surfaces. The default is theta = 0.5.
%
%   [...] = REFINE_CLOSE_PATCHES(S1,S2,THETA,MAXLEVELS) also sets a
%   safety limit on the number of simultaneous refinement rounds.
%
%   All patch marking is vectorized; no for-loops are used.

if nargin < 3 || isempty(theta)
    theta = 0.5;
end

if nargin < 4 || isempty(maxlevels)
    maxlevels = 20;
end

if ~isscalar(theta) || ~isfinite(theta) || theta <= 0
    error('theta must be a positive finite scalar.')
end

if ~isscalar(maxlevels) || maxlevels < 1 || maxlevels ~= floor(maxlevels)
    error('maxlevels must be a positive integer.')
end

S1refine = S1;
S2refine = S2;
nlevels = 0;

while true
    % Patch centroids as npatches-by-3 arrays.
    cms1 = S1refine.cms.';
    cms2 = S2refine.cms.';

    % Pairwise distances between centroids on opposing surfaces.
    dist2 = sum(cms1.^2,2) ...
          + sum(cms2.^2,2).' ...
          - 2*(cms1*cms2.');
    patch_dist = sqrt(max(dist2,0));

    d1 = min(patch_dist,[],2);
    d2 = min(patch_dist,[],1).';

    % Characteristic patch sizes h = sqrt(area).
    A1 = full(sparse(S1refine.patch_id(:),1,S1refine.wts(:), ...
        S1refine.npatches,1));
    A2 = full(sparse(S2refine.patch_id(:),1,S2refine.wts(:), ...
        S2refine.npatches,1));

    h1 = sqrt(A1);
    h2 = sqrt(A2);

    % Refine patches that fail h < theta*d.
    refinelist1 = find(h1 >= theta*d1);
    refinelist2 = find(h2 >= theta*d2);

    if isempty(refinelist1) && isempty(refinelist2)
        break
    end

    if nlevels >= maxlevels
        error(['Refinement did not satisfy the size-distance criterion ', ...
            'within %d levels.'],maxlevels)
    end

    if ~isempty(refinelist1)
        S1refine = split_patches(S1refine,refinelist1);
    end

    if ~isempty(refinelist2)
        S2refine = split_patches(S2refine,refinelist2);
    end

    nlevels = nlevels + 1;
end

info.theta = theta;
info.nlevels = nlevels;
info.npatches1 = S1refine.npatches;
info.npatches2 = S2refine.npatches;
info.max_ratio1 = max(h1./d1);
info.max_ratio2 = max(h2./d2);

end
