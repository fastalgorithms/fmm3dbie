function val = eval_mask(obj,src,targ,L)
% KERNEL3D.EVAL_MASK call obj.eval and zero out self interactions
% Optional L is the largest box edge of the full FMM source and target set.
% If omitted or empty, use the bounding box of src and targ.
val = obj.eval(src,targ);

% As in FMM3D: ignore pairs with |t-s| <= thresh, where thresh is 2^-51
% times the largest edge of the bounding box of the sources and targets
src_r  = src.r(:,:);
targ_r = targ.r(:,:);
targ_r(end+1:size(src_r,1),:) = 0;   % e.g. 2d targets of 3d sources
if nargin < 4 || isempty(L)
    r = [src_r, targ_r];
    L = max(max(r, [], 2) - min(r, [], 2));
end
thresh = 2^-51 * L;
if thresh == 0, thresh = realmin; end

% Candidates within thresh in every coordinate, then |t-s| <= thresh
[tf, loc] = ismembertol(src_r.', targ_r.', thresh, 'ByRows', true, ...
    'OutputAllIndices', true, 'DataScale', 1);
if ~any(tf), return; end
q = find(tf);
t = vertcat(loc{q});
s = repelem(q, cellfun(@numel, loc(q)));
keep = vecnorm(targ_r(:,t) - src_r(:,s)) <= thresh;
t = t(keep);
s = s(keep);

% Zero every component of each coincident (target, source) pair
m = obj.opdims(1);
n = obj.opdims(2);
[a, b] = ndgrid(1:m, 1:n);
rows = a(:) + m*(t(:).' - 1);
cols = b(:) + n*(s(:).' - 1);
val(rows + size(val,1)*(cols - 1)) = 0;

end
