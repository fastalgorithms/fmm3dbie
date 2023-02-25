function a = area(obj)
    a = sum(cat(1,obj.weights{:}));
end