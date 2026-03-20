function [list] = find_patch_refine_list(S,errlp,funlp,tol,patch_ids,lptype,eta)

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

for ic = 1:numel(patch_ids)
    ip = patch_ids(ic);
    for idim = 1:nd
        rhs = tol * funlp(idim) * patcharealp(ip) / surfarealp;
        lhs = errlp(idim,ip) * hpatch(ip)^eta;

        if lhs > rhs
            refine(ip) = true;
            break
        end
    end
end

list = find(refine);


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

%===========================================================================

function [list, errchecks] = find_refine_patches_list_ttt(S,fun,tol,lp)
% FIND_REFINE_PATCHES_LIST  find patches to refine based on tolerance criterion
%
% [list, nlist] = find_refine_patches_list(obj,fun,p,tol,patch_ids)
% [errps] = surf_fun_error(obj,fun,p)
%  
%       This routine returns a list of patches failing a tolerance
%       based criterion. Let ej be the approximation of the error in
%       L^p norm over the patch Pj. The patch Pj is flagged for refinement
%       if |e|_j^p < tol^p |fun|^p |Pj|/|P| where TOL is the given tolerance,
%       |Pj| is the area of Pj, and |P| is the sum of all |Pj|, and fun is the
%       function where are estimating the error of.
%
%       Inputs:
%       obj     --- the surfer object
%       fun     --- (nfuns,npts) array containing function values
%                   at the points in obj
%       tol     --- tolerance
%       lp       --- optional: the Lp-norm to use in error criterion
%                   (default is infinity)
%

%
% Output arguments:
%    - list(npatches): integer
%        indices of patches to refine
%    - nlist: integer
%        number of entries in list, number of patches to refine
%
if nargin < 4 || isempty(p)
    lp = 1:S.npatches;   % default: all patches
end
%
%   approximate surface area as sum of patchwise area
%


tail_errs = S.surf_fun_error(fun,lp);

     surfarea = S.area;
     patcharea = full(sparse(S.patch_id, 1, S.wts));

     %     reflist = funlp.*tol_dens/S.area;

     ncheck = length(patch_ids);
     list = zeros(ncheck,1);
     errchecks = zeros(nd,ncheck);
     nlist = 0;
     nd = numel(funlp,1);
 
     %
     %     check for each patch if refinement criterion is satisfied. If not
     %     satisfied, label the patch for refinement.
     %
     ict = 0;
      for idim=1:nd
         rsc = funlp(idim) * tol/surfarea;
         for ic=1:ncheck
             ip = patch_ids(ic);
             errchecks(idim,ic) = patcharea(ip)*rsc;
            if(errp(idim,ip) > patcharea(ip)*rsc)
               ict = ict + 1;
               list(ict) = ip;
            end
         end
      end

      nlist = ict;
      list = list(1:nlist);

      return
      
end



