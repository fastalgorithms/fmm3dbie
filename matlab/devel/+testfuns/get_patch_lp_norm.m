function [lp_norm] = get_patch_lp_norm(S,fun,patch_ids,lptype)


if nargin < 3 || isempty(patch_ids)
    patch_ids = 1:S.npatches;
end

if nargin < 4 || isempty(lptype)
    lptype = 1;   % 0 means infinity norm
end

[~,~,~,ixyzs,~,wts] = extract_arrays(S);
nd = size(fun,1);
lp_norm = zeros(nd,S.npatches);

for ic = 1:numel(patch_ids)
    ip = patch_ids(ic);
    idxpat = ixyzs(ip):(ixyzs(ip+1)-1);

    if lptype == 0
        lp_norm(:,ip) = max(abs(fun(:,idxpat)),[],2);
    elseif lptype == 1
        lp_norm(:,ip) = sum(abs(fun(:,idxpat)).*wts(idxpat).',2);
    elseif lptype == 2
        lp_norm(:,ip) = sqrt(sum(abs(fun(:,idxpat)).^2.*wts(idxpat).',2));
    else
        error('Unsupported lptype')
    end
end

return

end