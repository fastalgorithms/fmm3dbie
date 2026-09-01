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
        rhs = tol * funlp(idim) * patcharealp(ip) / surfarealp;
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

function list = find_patch_refine_list3(S,errs,funlp,tol,patch_ids,lptype,eta)
% FIND_REFINE_PATCHES_LIST Generates a list of patches failing an error check criterion
%
% Syntax:
%   [list] = find_refine_patches_list(S, errp)
%   [list,nlist] = find_refine_patches_list(S, errp, patch_ids)
%
%   Add description
%
% Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * errp: (nd,npts) array, pointwise errors of the nd functions
%    * funlp: (nd) global Lp norm of the nd functions: ||f||_p
%    * tol: tolerance for refinement criterion, such that a patch Pj
%           is flagged for refinement if errp(:,j) < tol*funlp(:)*(|Pj|/|P|)^(1/p)
%           where |Pj| is the area of patch Pj and P=U_j Pj
%    * (optional) patch_ids: indices of patches to check
%    * (optional) eta: scaling factor of error    
%
% Output arguments:
%    - list(npatches): integer
%        indices of patches to refine
%    - nlist: integer
%        number of entries in list, number of patches to refine
%
if nargin < 5 || isempty(patch_ids)
    patch_ids = 1:S.npatches;   % default: all patches
end

if nargin < 6 || isempty(lptype)
    lptype = 0;   % default: 0, p=infinity norm
end

if nargin < 7 || isempty(eta)
    eta = 0;   % default: h^eta = 1
end
%
%   approximate surface area as sum of patchwise area
%
     [~,~,~,ixyzs,~,wts] = extract_arrays(S);

     ncheck = numel(patch_ids);
     list = zeros(ncheck,1);
     nlist = 0;
     nd = numel(funlp);

     errlp = zeros(ncheck,1);     
     patcharealp = zeros(ncheck,1);
     hpweta = zeros(ncheck,1);
     
     if lptype==0
         for ipat = 1:ncheck
             ip = patch_ids(ipat);
             idxpat = ixyzs(ip):(ixyzs(ip+1)-1);
             for idim=1:nd
                 errlp(idim,ipat) = max(abs(errs(idim,idxpat)));
             end
             patcharealp(ipat) = 1;
             hpweta(ipat) = sqrt(sum(wts(idxpat)))^eta;
         end         
     end
     
     if lptype==1
         for ipat = 1:ncheck
             ip = patch_ids(ipat);
             idxpat = ixyzs(ipat):(ixyzs(ipat+1)-1);
             for idim=1:nd
                 errlp(idim,ipat) = sum(abs(errs(idim,idxpat)).*wts(idxpat));
             end
             patcharealp(ipat) = sum(wts(idxpat));
             hpweta(ipat) = sqrt(sum(wts(idxpat)))^eta;
         end         
     end

     if lptype==2
         for ipat = 1:ncheck
             ip = patch_ids(ipat);
             idxpat = ixyzs(ipat):(ixyzs(ipat+1)-1);
             for idim=1:nd
                 errlp(idim,ipat) = sqrt(sum(abs(errs(idim,idxpat)).^2.*wts(idxpat)));
             end
             patcharealp(ipat) = sqrt(sum(wts(idxpat)));
             hpweta(ipat) = sqrt(sum(wts(idxpat)))^eta;
         end         
     end

     surfarealp = (S.area)^(1/lptype);


     
     %
     %     check for each patch if refinement criterion is satisfied. If not
     %     satisfied, label the patch for refinement.
     %
     
     ict = 0;
      for idim=1:nd
         rsc = funlp(idim) * tol / surfarealp  ;
         for ic=1:ncheck
             ip = patch_ids(ic);
             errcheck = errlp(idim,ip) * hpweta;
            if(errp(idim,ip) > patcharealp(ip)*rsc)
               ict = ict + 1;
               list(ict) = ip;
            end
         end
      end

      list = list(1:ict);

      return
      
end

