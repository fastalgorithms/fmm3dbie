function [cfs] = vals2coefs(obj,vec)
%   
%   .   .   .   computes the expansion coefficients of a vector "vec" 
%               defined on a surfer object in the appropriate basis on
%               each patch (currently only Koornwinder). It is assumed 
%               that the entries of vec are ordered so contiguous points 
%               correspond to the same patch, and the ordering of the 
%               patches is the same as the ordering in the surfer object.
%
    if (obj.npts ~= size(vec,1))
        disp("The size of the input is incompatible with the surfer object");
        cfs = [];
        return
    end

    if (sum(obj.iptype ~= 1) ~= 0)
        disp("The surfer object has an unsupported patch type for this feature");
        sum(obj.iptype == 1)
        cfs = [];
    end

    umats = {};
    vmats = {};
    rnodes= {};
    norders = unique(obj.norders);

    for iord = 1:numel(norders)
        norder = norders(iord);
        rnode = koorn.rv_nodes(norder);
        umats{norder} = koorn.coefs2vals(norder,rnode);
        vmats{norder} = koorn.vals2coefs(norder,rnode);
        rnodes{norder} = rnode;
    end

    cfs = zeros(size(vec));
    iind = 0;
    
    for ipatch = 1:obj.npatches
        norder = obj.norders(ipatch);
        vmat = vmats{norder};
        npts = size(vmat,2);
        vals = vec(iind + (1:npts),:);
        coef = vmat*vals;
        cfs(iind+(1:npts),:) = coef;
        iind = iind + npts;
    end
        
    
end
