function [objout,varargout] = translate(obj,r)

% translates the surfer
% Input: obj    - a surfer object
%        r(3,n) - a matrix of 'shifts'
% Output: a single surfer object containing the translated copies

    sz = size(r);
    if (sz(1) ~=3) 
        fprintf('incompatible sizes');
        objout = surfer.empty;
        return
    end


    
    Sarray = surfer.empty;
    for ii=1:sz(2)
        objtmp = affine_transf(obj,eye(3),r(:,ii));
        Sarray(ii)= objtmp;
    end

    objout = merge(Sarray);
    
end
