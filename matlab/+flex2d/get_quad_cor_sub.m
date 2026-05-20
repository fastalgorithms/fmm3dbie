function xmat = get_quad_cor_sub(S, gkern, eps, zk, type, nover, targinfo)

    if nargin < 6
        nover = 0;
    end

    if nargin < 7
        targinfo = [];
        targinfo.r = S.r;
        targinfo.patch_id = S.patch_id;
        targinfo.uvs_targ = S.uvs_targ;
    end
    targs = targinfo.r;
    [ndtarg,ntarg] = size(targs);
    
    if isfield(targinfo,'patch_id')
      patch_id = targinfo.patch_id;
    else
      patch_id = zeros(ntarg,1);
    end

    if isprop(targinfo,'uvs_targ')
      uvs_targ = targinfo.uvs_targ;
    else
      uvs_targ = zeros(2,ntarg);
    end


    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    npatches = S.npatches;
    npts = S.npts;

    iptype_avg = floor(sum(iptype)/(npatches+0.0d0));
    norder_avg = floor(sum(norders)/(npatches+0.0d0));
    [~, rfac0] = get_rfacs(norder_avg,iptype_avg);
    [rsc] = getnear(S,targinfo);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    nnz = length(rsc.col_ind);
    nquad = rsc.iquad(nnz+1)-1;

    iquadtype = 1;


    wnear = flex2d.getnearquad_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targs,patch_id, uvs_targ, eps,...
          iquadtype,nnz, rsc.row_ptr,rsc.col_ind,rsc.iquad,rfac0, ...
          zk,nquad,'type');

    npols = ixyzs(2:end)-ixyzs(1:end-1);
    npts_col = npols( rsc.col_ind);

    [nt,~] = size( rsc.row_ptr);
    ntarg = nt-1;
    nrep = zeros(ntarg,1);
    istarts =  rsc.row_ptr(1:end-1);
    iends =  rsc.row_ptr(2:end)-1;
    icol_ind = zeros(sum(npts_col),1);
    isrcinds = cell(npatches,1);
    for i=1:npatches
        isrcinds{i} = ixyzs(i):ixyzs(i+1)-1;
    end

    istart = 1;
    for i=1:ntarg
        nrep(i) = sum(npts_col(istarts(i):iends(i))); 
        iinds = horzcat(isrcinds{ rsc.col_ind(istarts(i):iends(i))});
        nelem = length(iinds);
        icol_ind(istart:istart+nelem-1) = iinds;
        istart = istart+nelem;
    end
    irow_ind = repelem((1:ntarg)',nrep);    
    
    xmat = sparse(irow_ind,icol_ind, wnear, ntarg, S.npts);

    Asmth_over = smooth_sparse_quad(gkern,targs,S, rsc.row_ptr, rsc.col_ind,nover); 

    xmat = xmat - Asmth_over;

end
%
%
