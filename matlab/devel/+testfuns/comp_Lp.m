
function [ffun_lp] = comp_Lp(S,ffun,lptype)

if lptype == 0
    ffun_lp = max(abs(ffun),[],1);
else
    ffun_lp = sum(abs(ffun).^lptype.*S.wts(:),1).^(1/lptype);
end

end
