function [coefs] = vals2coefs(obj,vals)
% VALS2COEFS  convert values to coeffs on given surfer object
%
% [coefs] = vals2coefs(S,vals) returns a concatenated array of poly coeffs of
%   the vector of vals on the nodes of S, a surfer object.
%
% Inputs: S is a surfer object
%         vals is (n, S.npts) or (S.npts, n) array (n=1,2... is the number of scalar funcs)
% Outputs:
%         coeffs is array of coeffs, in the native patch ordering, same size of vals
%
% *** to check this doc
    npatches = obj.npatches;
    npts = obj.npts;
    
    vals_in = vals;
    if size(vals,1) == npts
    elseif size(vals,2) == npts
       vals = vals.';
    else
        error("VALS2COEFS:error in size");
        coefs = [];
        return
    end
    
    ikrn = find(obj.iptype == 1);
    ileg = find(obj.iptype == 11);
    iche = find(obj.iptype == 12);

%%%%%%%%%%%%%%
%       .   .   .   first the koornwinder

    nords = obj.norders(ikrn);


%%%%%%%%%%%%%%

    ntmp = zeros(npatches,2);
    ntmp(:,1) = obj.norders;
    ntmp(:,2) = obj.iptype;
    
    
    coefs = zeros(size(vals));
    
    [ntmp_uni,~,intmp] = unique(ntmp,'rows');
    nuni = size(ntmp_uni,1);
    umats = cell(nuni,1);
    for i=1:nuni
        norder = ntmp_uni(i,1);
        iptype = ntmp_uni(i,2);
        if(iptype == 1)
            rnodes = koorn.rv_nodes(norder);
            umats{i} = koorn.vals2coefs(norder,rnodes);       
        elseif(iptype == 11)
            rnodes = polytens.lege.nodes(norder);
            umats{i} = polytens.lege.vals2coefs(norder,rnodes);
        elseif(iptype == 12)
            rnodes = polytens.cheb.nodes(norder);
            umats{i} = polytens.cheb.vals2coefs(norder,rnodes);
        end
        npp  = size(umats{i},2);
        inds = (intmp == i);
        indv = obj.ixyzs(inds);
        [I,V] = meshgrid(indv,1:npp);
        itot = I(:)+V(:)-1;
        cftmp = umats{i}*reshape(vals(itot,:),npp,[]);
        coefs(itot,:) = reshape(cftmp,[],size(vals,2));
    end
    
    coefs = reshape(coefs, size(vals_in));
    
end
