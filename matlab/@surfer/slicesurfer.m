function obj2 = slicesurfer(obj, ipatchkeep)
% SLICESURFER Extract a subset of patches from a surfer object.
%
% S2 = slicesurfer(S, ipatchkeep) returns a new surfer object S2 containing
%  only the patches of S whose indices are listed in ipatchkeep.
%
% Input:
%   S           - surfer object
%   ipatchkeep  - integer vector of patch indices to retain (1-based)
%
% Output:
%   S2          - surfer object with npatches = numel(ipatchkeep)
%
% Example:
%   S2 = slicesurfer(S, 1:3);   % keep first three patches
%
% See also MERGE, SPLIT_PATCHES

ipatchkeep = ipatchkeep(:).';

npatches2 = length(ipatchkeep);

npts_S2 = 0;
for i = ipatchkeep
    npts_S2 = npts_S2 + (obj.ixyzs(i+1) - obj.ixyzs(i));
end

srcvals_S = zeros(12, npts_S2);
col = 1;
for i = ipatchkeep
    iinds = obj.ixyzs(i):obj.ixyzs(i+1)-1;
    nipatch     = length(iinds);
    srcvals_S(:, col:col+nipatch-1) = [obj.r(:,iinds); obj.du(:,iinds); obj.dv(:,iinds); obj.n(:,iinds)];
    col = col + nipatch;
end

obj2 = surfer(npatches2, obj.norders(ipatchkeep), srcvals_S, obj.iptype(ipatchkeep));
end