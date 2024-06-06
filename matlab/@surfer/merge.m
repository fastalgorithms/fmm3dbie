function [objout] = merge(Sarray)
% MERGE  Combine array of surfer objects into a single such object.
  
    npout  = 0;
    noout  = [];
    srout  = [];
    ipout  = [];
    for ii=1:numel(Sarray)
        S = Sarray(ii);
 	    npatches = S.npatches;
        norders  = S.norders;        
        srcvals  = S.srcvals;
        iptype   = S.iptype;
        srcvals = [srcvals{:}];
        npout = npout + npatches;
        noout = [noout;norders];
        srout = [srout,srcvals];
        ipout = [ipout;iptype];
    end

    objout = surfer(npout,noout,srout,ipout);
end