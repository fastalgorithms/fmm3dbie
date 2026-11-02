function val = eval_mask(obj,src,targ)
% KERNEL3D.EVAL_MASK call obj.eval and zero out self interactions
val = obj.eval(src,targ);

% self interactions: |t - s| <= 2^-51 * bounding box size, as in FMM3D
src_r  = src.r(:,:);
targ_r = targ.r(:,:);
d = max(size(src_r,1), size(targ_r,1));
src_r(end+1:d,:)  = 0;
targ_r(end+1:d,:) = 0;
allr   = [src_r, targ_r];
thresh = 2^-51 * max(max(allr,[],2) - min(allr,[],2));
for q = 1:size(src_r, 2)
    diff_norms = vecnorm(targ_r - src_r(:,q));
    self_mask  = diff_norms <= thresh;
    if any(self_mask)
        row_off = (1:obj.opdims(1)) + obj.opdims(1)*(find(self_mask).' - 1);
        col_off = (1:obj.opdims(2)) + obj.opdims(2)*(q-1);
        val(row_off(:), col_off) = 0;
    end
end
            
end