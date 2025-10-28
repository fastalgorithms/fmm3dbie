function Asmth = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,norders,lbat)
    if nargin < 7
        lbat = 1e3;
    end
    if isa(kern,'kernel')
        kern = kern.eval;
    end
    ntarg = length(row_ptr)-1;
    
    S_over = oversample(S,norders);
    
    [irow_ind_ov,icol_ind_ov,npts_col_ov] = unpack_rsc(row_ptr,col_ind,S_over);
    nnz_ov = length(irow_ind_ov);
    vals = zeros(1,nnz_ov);
    
    nbat = ceil(nnz_ov/lbat);

    for k = 1:nbat
        ks = (lbat*(k-1)+1):min(lbat*k,nnz_ov);
        rsrc = S_over.r(:,icol_ind_ov(ks));
        rtarg = targs(:,irow_ind_ov(ks));
        rnear = rtarg - rsrc;
        vals(ks) = kern(struct('r',[0;0;0]),struct('r',rnear)).*S_over.wts(icol_ind_ov(ks));
    end

    [irow_ind,icol_ind,npts_col] = unpack_rsc(row_ptr,col_ind,S);
    nnz = length(irow_ind);

    wnear = zeros(1,nnz);
    [xinterpmat,intmp] = get_oversamps(S,noveruse);


    strt = 0; strt_ov = 0;
    for i = 1:ntarg
        js = row_ptr(i):(row_ptr(i+1)-1);
        for j = 1:length(js)
            jdloc = strt + (1:npts_col(j));
            jdloc_ov = strt_ov + (1:npts_col_ov(j));

            xinterp = xinterpmat{intmp(col_ind(js(j)))};
            wnear(jdloc) = (xinterp*vals(jdloc_ov).').';
            
            strt = strt + npts_col(j);
            strt_ov = strt_ov + npts_col_v(j);
        end
    end

end



function [irow_ind,icol_ind,npts_col] = unpack_rsc(row_ptr,col_ind,S)

    ixyzs = S.ixyzs;
    npols = ixyzs(2:end)-ixyzs(1:end-1);
    npts_col = npols(col_ind);

    [nt,~] = size(row_ptr);
    ntarg = nt-1;
    nrep = zeros(ntarg,1);
    istarts = row_ptr(1:end-1);
    iends = row_ptr(2:end)-1;
    icol_ind = zeros(sum(npts_col),1);
    isrcinds = cell(S_over.npatches,1);
    for i=1:S_over.npatches
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
   
end

function [xinterpmat,intmp] = get_oversamps(S,noveruse)


    xinterpmat = cell(nuni,1);

    ntmp = zeros(npatches,3);
    ntmp(:,1) = S.norders;
    ntmp(:,2) = noveruse;
    ntmp(:,3) = S.iptype;
    
    [ntmp_uni,~,intmp] = unique(ntmp,'rows');
    for i=1:nuni        
        norder = ntmp_uni(i,1);
        nover = ntmp_uni(i,2);
        iptype = ntmp_uni(i,3);
        if(iptype == 1)
            rnodes = koorn.rv_nodes(norder);
            rnodes_over = koorn.rv_nodes(nover);

            vmat = koorn.coefs2vals(norder,rnodes_over);
            umat = koorn.vals2coefs(norder,rnodes);  
            
        elseif(iptype == 11)
            rnodes = polytens.lege_nodes(norder);
            rnodes_over = polytens.lege_nodes(nover);
            vmat = polytens.lege_coefs2vals(norder,rnodes_over);
            umat = polytens.lege_vals2coefs(norder,rnodes);
        elseif(iptype == 12)
            rnodes = polytens.cheb_nodes(norder);
            rnodes_over = polytens.cheb_nodes(nover);
            vmat = polytens.cheb_coefs2vals(norder,rnodes_over);
            umat = polytens.cheb_vals2coefs(norder,rnodes);
        end
        
        xinterpmat{i} = vmat*umat;
    end

end