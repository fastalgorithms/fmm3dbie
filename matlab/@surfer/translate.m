function [objout,varargout] = translate(obj, r)
% TRANSLATE translate a surfer object by given vector(s)
%
% S2 = translate(S,v) where v is a 3x1 vector returns S translated by v.
%  If v is 3xn with n>1, duplicates the object S at each translation
%  given by each column of v.
%
% Input: obj    - a surfer object
%        r(3,n) - a matrix of 'shifts' (one column vector per shift)
% Output: a single surfer object containing the translated copies

    sz = size(r);
    if (sz(1) ~=3) 
        error('incompatible sizes\n');
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
