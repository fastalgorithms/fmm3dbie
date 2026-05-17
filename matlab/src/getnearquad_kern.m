function [xmat,norderup] = getnearquad_kern(obj, kern, eps, getnearquad, targinfo, norderup)

    try
        kernuse = kern.eval;
    catch
        kernuse = kern;
    end

    if nargin < 6
        norderup = 2;
    end

    if nargin < 5
        targinfo = [];
        targinfo.r = obj.r;
        targinfo.patch_id = obj.patch_id;
        targinfo.uvs_targ = obj.uvs_targ;
    end

    [ndtarg,ntarg] = size(targinfo.r);
    
    if isfield(targinfo,'patch_id') || isprop(targinfo,'patch_id')
      patch_id = targinfo.patch_id;
    else
      patch_id = zeros(ntarg,1);
    end
    if isempty(patch_id), patch_id = zeros(ntarg,1); end

    if isfield(targinfo,'uvs_targ') || isprop(targinfo,'uvs_targ')
      uvs_targ = targinfo.uvs_targ;
    else
      uvs_targ = zeros(2,ntarg);
    end
    if isempty(uvs_targ), uvs_targ = zeros(2,ntarg); end

    loc_field = fieldnames(targinfo)';
    loc_field_use = {};
    for i = 1:length(loc_field)
        if mod(size(targinfo.(loc_field{i}),2), size(targinfo.r(:,:),2)) == 0 && ~isempty(targinfo.(loc_field{i}))
        loc_field_use{length(loc_field_use)+1} = loc_field{i};
        end
    end

    targdim = [];
    for field = loc_field_use
        targdim.(field{1}) = targinfo.(field{1})(:,1);
    end
    tmp = kernuse(struct('r',[0;0;0]),targdim);
    opdims = size(tmp);

    [srcvals,srccoefs,norders,ixyzs,iptype,~] = extract_arrays(obj);
    npatches = obj.npatches;
    npts = obj.npts;

    rsc = getnear(obj, targinfo);

    nnz = length(rsc.col_ind);
    nquad = rsc.iquad(end) - 1;
    wnear = zeros(prod(opdims), nquad);

    iquadtype = 1;

    if nnz > 0
        wnear = getnearquad(npatches, norders, ixyzs, ...
            iptype, npts, srccoefs, srcvals, targinfo, patch_id, uvs_targ, eps, ...
            iquadtype, nnz, rsc.row_ptr, rsc.col_ind, rsc.iquad, nquad);
        if size(wnear, 1) == nquad
            wnear = wnear.';
        end
    end

    xmat = unpack_wnear(obj, rsc, wnear, opdims);

    Asmth_over = smooth_sparse_quad(kernuse, targinfo, obj, rsc.row_ptr, rsc.col_ind, obj.norders(1)+norderup);

    xmat = xmat - Asmth_over;


end


function xmat = unpack_wnear(obj,rsc,wnear,opdims)
    npols = obj.ixyzs(2:end)-obj.ixyzs(1:end-1);
    npts_col = npols(rsc.col_ind);

    [nt,~] = size(rsc.row_ptr);
    ntarg = nt-1;
    nrep = zeros(ntarg,1);
    istarts =  rsc.row_ptr(1:end-1);
    iends =  rsc.row_ptr(2:end)-1;
    icol_ind = zeros(sum(npts_col),1);
    isrcinds = cell(obj.npatches,1);
    for i=1:obj.npatches
        isrcinds{i} = obj.ixyzs(i):obj.ixyzs(i+1)-1;
    end

    istart = 1;
    for i=1:ntarg
        nrep(i) = sum(npts_col(istarts(i):iends(i))); 
        iinds = horzcat(isrcinds{rsc.col_ind(istarts(i):iends(i))});
        nelem = length(iinds);
        icol_ind(istart:istart+nelem-1) = iinds;
        istart = istart+nelem;
    end
    irow_ind = repelem((1:ntarg)',nrep);
    irow_indopdim = (1:opdims(1)).' + opdims(1)*(irow_ind(:).'-1);
    icol_indopdim = 0*(1:opdims(1)).' + (icol_ind.');
    irow_indopdim = irow_indopdim(:);
    icol_indopdim = icol_indopdim(:);
    irow_indopdim = 0*(1:opdims(2)).' + (irow_indopdim.');
    icol_indopdim = (1:opdims(2)).' + opdims(2)*(icol_indopdim.'-1);
    irow_indopdim = irow_indopdim(:);
    icol_indopdim = icol_indopdim(:);
   xmat = sparse(irow_indopdim,icol_indopdim, wnear(:), opdims(1)*ntarg, opdims(2)*obj.npts);


end


