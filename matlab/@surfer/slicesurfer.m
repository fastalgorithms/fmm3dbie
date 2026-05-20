function obj2 = slicesurfer(obj, ipatchkeep)

ipatchkeep = ipatchkeep(:).';

npatches2 = length(ipatchkeep);

npts_S2 = 0;
for i = ipatchkeep
    npts_S2 = npts_S2 + (obj.S.ixyzs(i+1) - obj.S.ixyzs(i));
end

srcvals_S = zeros(12, npts_S2);
col = 1;
for i = ipatchkeep
    iinds = obj.S.ixyzs(i):obj.S.ixyzs(i+1)-1;
    nipatch     = length(iinds);
    srcvals_S(:, col:col+nipatch-1) = [obj.S.r(:,iinds); obj.S.du(:,iinds); obj.S.dv(:,iinds); obj.S.n(:,iinds)];
    col = col + nipatch;
end

obj2 = surfer(npatches2, obj.S.norders(ipatchkeep), srcvals_S, obj.S.iptype(ipatchkeep));
end