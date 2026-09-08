function ids = surferids(srfobj, isrfs)
%surferids return the indices corresponding to points in srfobj(isrfs)
if iscell(srfobj)
    npts = arrayfun(@(i)size(srfobj{i}.r(:,:),2),1:length(srfobj));
else
    npts = arrayfun(@(i)size(srfobj(i).r(:,:),2),1:length(srfobj));
end
% get starting index for each edge
idstart = [1,cumsum(npts)+1];

ids = [];
for i = 1:length(isrfs)
    ids = [ids, idstart(isrfs(i)):(idstart(isrfs(i)+1)-1)];
end
end