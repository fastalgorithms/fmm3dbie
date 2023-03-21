function [objout,varargout] = affine_transf(obj,mat,shift)
%
%  performs an affine transformation of the given object
%
    objout = obj;

 	npatches = obj.npatches;
    norders  = obj.norders;
    
    srcvals  = obj.srcvals;
    srccoefs = obj.srccoefs;
    weights  = obj.weights;
    
    if (nargin < 2)
        return
    elseif (nargin < 3)
        shift = [0;0;0];
    end
    
    if (nargin == 2 && ~isequal(size(mat),[3,3]))
        fprintf('Incompatible scaling matrix, returning same object \n');
        return
    end
    if (nargin == 3 && ~isequal(size(shift),[3,1]))
        fprintf('Incompatible shift vector, returning same object \n');
        return
    end
    
    if(length(norders)==1)
        
        obj.norders = ones(npatches,1)*norders;
        
        rwts = cell(1);
        umats = cell(1);
        rnodes = koorn.rv_nodes(norders);
        rwts{1} = koorn.rv_weights(norders);
        umats{1} = koorn.vals2coefs(norders,rnodes);
        
        iuse = ones(npatches,1);
        
        
        
    elseif(length(norders)==npatches)
        obj.norders = norders;
        [norders_uni,~,iuse] = unique(norders);
        nuni = length(norders_uni);
        rwts = cell(nuni,1);
        umats = cell(nuni,1);
        for i=1:nuni
            rnodes = koorn.rv_nodes(norders_uni(i));
            rwts{i} = koorn.rv_weights(norders_uni(i));
            umats{i} = koorn.vals2coefs(norders_uni(i),rnodes);
        end
    else
        fprintf('Incompatible size of norders, returning\n');
        return;
    end
    
    for i=1:npatches
        sptch = [srcvals{i}];
        r = sptch(1:3,:);
        du= sptch(4:6,:);
        dv= sptch(7:9,:);
        n = sptch(10:12,:);
        
        r = mat*r+shift;
        du= mat*du;
        dv= mat*dv;
        n = mat*n;
        da = vecnorm(cross(du,dv),2).*rwts{iuse(i)};
        weights{i} = da(:);
        
        sptch = [r;du;dv;n];
        srccoefs{i} = sptch(1:9,:)*umats{iuse(i)}';
        srcvals{i}  = sptch;
    end    
    
    objout.weights = weights;
    objout.srcvals = srcvals;
    objout.srccoefs = srccoefs;
    sv = [objout.srcvals{:}];
    objout.r  = sv(1:3,:);
    objout.du = sv(4:6,:);
    objout.dv = sv(7:9,:);
    objout.n  = sv(10:12,:);
    
end
