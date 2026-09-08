function val = eval_mask(obj,src,targ)
% KERNEL3D.EVAL_MASK call obj.eval and zero out self interactions
val = obj.eval(src,targ);

src_norm = max(vecnorm(src.r(:,:)), 1);
targ_r   = targ.r(:,:);
for q = 1:size(src.r(:,:), 2)
    diff_norms = vecnorm(targ_r - src.r(:,q));
    self_mask  = diff_norms < 1e-14 * src_norm(q);
    if any(self_mask)
        row_off = (1:obj.opdims(1)) + obj.opdims(1)*(find(self_mask).' - 1);
        col_off = (1:obj.opdims(2)) + obj.opdims(2)*(q-1);
        val(row_off(:), col_off) = 0;
    end
end
            
end