function [objout,varargout] = oversample(obj,nover)
%
%  Oversamples the given geometry given a current discretization
%  and a patchwise oversampling parameter, or a common oversampling
%  order for all patches
%
    npatches = obj.npatches;
    npts = obj.npts;
    if(length(nover)==1)
        noveruse = nover*ones(npatches,1);
    elseif(length(nover)==npatches)
        noveruse = nover;
    else
        fprintf('Incorrect size of npatches, returning same obj\n');
        objout = obj;
        varargout{1} = speye(npts);
    end
    
    ntmp = zeros(npatches,2);
    ntmp(:,1) = obj.norders;
    ntmp(:,2) = noveruse;
    
    npts_per_patch = (noveruse+1).*(noveruse+2)/2;
    npts_per_patch = npts_per_patch(:);
    npts_per_patch = [1; npts_per_patch];
    ixyzso = cumsum(npts_per_patch);
    nptso = ixyzso(end)-1;
    srcover = zeros(12,nptso);
    
    
    [ntmp_uni,~,intmp] = unique(ntmp,'rows');
    nuni = length(ntmp_uni);
    
    vmats = cell(nuni,1);
    umats = cell(nuni,1);
    ispl = cell(nuni,1);
    for i=1:nuni
        
        norder = ntmp_uni(i,1);
        nover = ntmp_uni(i,2);
        rnodes = koorn.rv_nodes(norder);
        rnodes_over = koorn.rv_nodes(nover);
        
        vmats{i} = koorn.coefs2vals(norder,rnodes_over);
        umats{i} = koorn.vals2coefs(norder,rnodes);  
        n1 = mod(norder,3);
        n2 = mod(nover,3);
        if(~n1 && ~n2) 
            i1 = 1;
            i2 = 1;
            if(norder == 15), i1 = 76; end
            if(norder == 18), i1 = 106; end
            if(nover == 15), i2 = 76; end
            if(nover == 18), i2 = 106; end
            
            ispl{i} = [i1;i2];
        end
    end
       
    
    
    for i=1:npatches
        istart = ixyzso(i);
        iend = ixyzso(i+1)-1;
        iind = istart:iend;
        srcover(1:9,iind) = obj.srccoefs{i}*vmats{intmp(i)}';
        ru = srcover(4:6,iind);
        rv = srcover(7:9,iind);    
        rtmp = cross(ru,rv);
        vnorm = repmat(vecnorm(rtmp,2),[3,1]);
        srcover(10:12,iind) = rtmp./vnorm;
        if(~isempty(ispl{intmp(i)}))
            srcover(:,ixyzso(i)+i2-1) = obj.srcvals{i}(:,i1);
        end
    end
    
    objout = surfer(npatches,noveruse,srcover);
    
    if(nargout>1)
        xinterpmat = cell(nuni,1);
        isysmat_all = cell(nuni,1);
        jsysmat_all = cell(nuni,1);
        vsysmat_all = cell(nuni,1);
        for i=1:nuni
            xinterpmat{i} = vmats{i}*umats{i};
            if(~isempty(ispl{i}))
                xinterpmat{i}(i2,:) = 0;
                xinterpmat{i}(i2,i1) = 1;
            end   
            [isysmat_all{i},jsysmat_all{i},vsysmat_all{i}] = find(xinterpmat{i});
        end
        isysmat = [];
        jsysmat = [];
        vsysmat = [];
        for i=1:npatches
            isysmat = [isysmat;isysmat_all{intmp(i)}+ixyzso(i)-1];
            jsysmat = [jsysmat;jsysmat_all{intmp(i)}+obj.ixyzs(i)-1];
            vsysmat = [vsysmat;vsysmat_all{intmp(i)}];            
        end
        varargout{1} = sparse(isysmat,jsysmat,vsysmat,nptso,obj.npts);  
    end
    
    
    
end
