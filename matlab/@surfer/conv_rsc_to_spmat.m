function spmat = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear)
    [~,ndim] = size(wnear);
    ixyzs = S.ixyzs(:);
    npatches = S.npatches;
    npols = ixyzs(2:end)-ixyzs(1:end-1);


    npts_col = npols(col_ind);

    [nt,~] = size(row_ptr);
    ntarg = nt-1;
    nrep = zeros(ntarg,1);
    istarts = row_ptr(1:end-1);
    iends = row_ptr(2:end)-1;
    icol_ind = zeros(sum(npts_col),1);
    isrcinds = cell(npatches,1);
    for i=1:npatches
        isrcinds{i} = ixyzs(i):ixyzs(i+1)-1;
    end

    istart = 1;
    for i=1:ntarg
        nrep(i) = sum(npts_col(istarts(i):iends(i))); 
        iinds = horzcat(isrcinds{col_ind(istarts(i):iends(i))});
        nelem = length(iinds);
        icol_ind(istart:istart+nelem-1) = iinds;
        istart = istart+nelem;
    end
    irow_ind = repelem((1:ntarg)',nrep);
    
    if(ndim > 1)
        spmat = cell(ndim,1);
        for i=1:ndim 
            spmat{i} = sparse(irow_ind,icol_ind,wnear(:,i),ntarg,npts);
        end
    else
        spmat = sparse(irow_ind,icol_ind,wnear,ntarg,npts);
    end    
end