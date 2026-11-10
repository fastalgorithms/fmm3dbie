function [refine_list,errq_list] = find_patch_refine_list(S,errlp,funlp,tol,patch_ids,lptype,eta)

if nargin < 5 || isempty(patch_ids)
    patch_ids = 1:S.npatches;
end

if nargin < 6 || isempty(lptype)
    lptype = 1;   % 0 means infinity norm
end

if nargin < 7 || isempty(eta)
    eta = 0;
end

patcharea = full(sparse(S.patch_id(:),1,S.wts(:),S.npatches,1));
hpatch = sqrt(patcharea);

if lptype == 0
    surfarealp = 1;
    patcharealp = ones(size(patcharea));
else
    surfarealp = (S.area)^(1/lptype);
    patcharealp = patcharea.^(1/lptype);
end

nd = numel(funlp);

refine = false(S.npatches,1);
errq = zeros(S.npatches,1);

for ic = 1:numel(patch_ids)
    ip = patch_ids(ic);
    for idim = 1:nd
        rhs = tol * funlp(idim) * patcharealp(ip)/ surfarealp;
        lhs = errlp(idim,ip) * hpatch(ip)^eta;

        if lhs > rhs
            refine(ip) = true;
            if (lhs/rhs) > errq(ip)
                errq(ip) = lhs/rhs;
            end
            % break
        end
    end
end

refine_list = find(refine);
errq_list = errq(refine_list);

end

