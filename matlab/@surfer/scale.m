function [objout,varargout] = scale(obj,sf)
% SCALE  Create a rescaled copy of a surfer object.
%
% Snew = scale(S,sf) creates a new surfer object sf times bigger than S
%  (all geometry is similarly scaled). If sf is a 3-vector, the scale
%  factors apply in x,y,z axes respectively.

  nsf = numel(sf);
  if numel(sf) == 1
    sf = sf*ones(3,1);
  end
  if (length(sf) ~=3) 
    fprintf('incompatible sizes');
    objout = surfer.empty;
    return
  end
  
  objout = affine_transf(obj,diag(sf),[0;0;0]);
end
