function fmax = patch_max(obj, f, p)
% PATCH_MAX  Maximum of a function over each patch.
%
% fmax = patch_max(obj, f) returns the maximum value of f over each patch.
% fmax = patch_max(obj, f, p) additionally scales each patch value by
%   (patch_area / total_area)^p, where patch areas are computed by summing
%   the quadrature weights (obj.wts) over each patch.
%
% Input:
%   f    - function values at obj's nodes, size [m, obj.npts] or [obj.npts, 1]
%   p    - (optional) area scaling exponent
%
% Output:
%   fmax - maximum value of f over each patch, size [m, obj.npatches]
%
% See also SURF_FUN_ERROR

[m1, m2] = size(f);

if m2 == obj.npts
    fuse = f;
elseif m1 == obj.npts
    fuse = f.';
else
    error('patch_max: f must have obj.npts columns or rows');
end

m = size(fuse, 1);
fmax = zeros(m, obj.npatches);

if nargin >= 3
    wts = obj.wts(:);
    total_area = sum(wts);
    scale = zeros(1, obj.npatches);
    for i = 1:obj.npatches
        iind = obj.ixyzs(i):obj.ixyzs(i+1)-1;
        scale(i) = (sum(wts(iind)) / total_area)^p;
    end
else
    scale = ones(1, obj.npatches);
end

for i = 1:obj.npatches
    iind = obj.ixyzs(i):obj.ixyzs(i+1)-1;
    fmax(:, i) = max(fuse(:, iind), [], 2) * scale(i);
end

if m1 == obj.npts
    fmax = fmax.';
end

end
