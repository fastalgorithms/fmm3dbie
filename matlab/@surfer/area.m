function a = area(obj)
% AREA Surface area of surface object (sum of weights).
    a = sum(cat(1,obj.weights{:}));
end
