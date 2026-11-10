function [err_lp_abs,err_lp_rel,fun_lp] = comp_err_upsamp_ref( ...
    S,Sref,fun,fun_ref,lptype)

% interp_lp_error maps data from its first surface onto its second.
[err_lp_abs,fun_int2S] = testfuns.interp_lp_error( ...
    S,Sref,fun,fun_ref,lptype);

if lptype == 0
    err_ref_lp = max(abs(fun_ref),[],2);
    fun_lp = max(abs(fun),[],2);
else
    err_ref_lp = ( ...
        sum(abs(fun_ref).^lptype.*Sref.wts(:).',2)).^(1/lptype);
    fun_lp = sum(abs(fun).^lptype.*S.wts(:).',2).^(1/lptype);    
end

err_lp_rel = zeros(size(err_lp_abs));
nonzero_ref = err_ref_lp ~= 0;
err_lp_rel(nonzero_ref) = ...
    err_lp_abs(nonzero_ref)./err_ref_lp(nonzero_ref);
err_lp_rel(~nonzero_ref & err_lp_abs ~= 0) = Inf;

end