function a = area(obj)
% AREA Surface area of surfer object (sum of quadrature weights).
%
% a = area(S)
    a = sum(cat(1,obj.weights{:}));
end
