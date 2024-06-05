function [objout,varargout] = scale(obj,sf)

% translates the surfer
% Input: obj    - a surfer object
%        sf(3)  - scaling parameters
% Output: a single surfer object containing the scaled surfer object

    nsf = numel(sf);
    if (nsf ~=3) 
        fprintf('incompatible sizes');
        objout = surfer.empty;
        return
    end

    objout = affine_transf(obj,diag(sf),[0;0;0]);

end