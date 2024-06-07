function [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(obj)
% EXTRACT_ARRAYS  extract and concatenate panel arrays from surfer object.
%
% [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(obj)
    srcvals = [obj.srcvals{:}];
    srccoefs = [obj.srccoefs{:}];
    ixyzs = obj.ixyzs;
    iptype = obj.iptype;
    norders = obj.norders;
    wts = cat(1,obj.weights{:});
end