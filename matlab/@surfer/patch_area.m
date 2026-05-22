function as = patch_area(obj)
% PATCH_AREA  Area of each patch.
%
% as = patch_area(S)
as = zeros(1,obj.npatches);
for i = 1:obj.npatches
    as(i) = sum(obj.weights{i});
end
end
