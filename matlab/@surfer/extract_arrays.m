function [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(obj)
    srcvals = horzcat(obj.srcvals{:});
    srccoefs = horzcat(obj.srccoefs{:});
    ixyzs = obj.ixyzs;
    iptype = obj.iptype;
    norders = obj.norders;
    wts = vertcat(obj.weights{:});
end