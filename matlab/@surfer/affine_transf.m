function [objout,varargout] = affine_transf(obj,mat,shift)
%
%  performs an affine transformation of the given object
%
    objout = obj;

 	npatches = obj.npatches;
    norders  = obj.norders;
    
    srcvals  = obj.srcvals;
    iptype   = obj.iptype;
    
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
    

    
    for i=1:npatches
        sptch = [srcvals{i}];
        r = sptch(1:3,:);
        du= sptch(4:6,:);
        dv= sptch(7:9,:);
        n = sptch(10:12,:);
        
        r = mat*r+shift;
        du= mat*du;
        dv= mat*dv;
        n = cross(du,dv);
        nnorm = vecnorm(n,2,1);
        n = n./nnorm;
        sptch = [r;du;dv;n];
        srcvals{i}  = sptch;
    end    
    srcvals = [srcvals{:}];
    objout = surfer(npatches,norders,srcvals,iptype);
    
end
