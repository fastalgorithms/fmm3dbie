function xmat = get_quad_cor_sub(S, type, zpars, eps, gs_kern, ivpp)

    if strcmp(type,'gs')
        iker = 1;
    elseif strcmp(type,'gphi')
        iker = 2;
    elseif strcmp(type,'bilapgs')
        iker = 3;
    elseif strcmp(type,'bilapgphi')
        iker = 4;
    elseif strcmp(type,'s3dgphi')
        iker = 5;
    elseif strcmp(type,'s3d')
        iker = 6;
    else
        error('kernel name not recognized') 
    end 

    if nargin < 6
        ivpp = 1;
    end

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    npatches = S.npatches;
    npp1 = npatches+1;
    npts = S.npts;
    ntp1 = npts+1;
    n9 = 9;
    n12 = 12;
    ndtarg = 3;

    % [patch_id, uvs_targ] = get_patch_id_uvs(S);

    %this might need fixing

    iptype_avg = floor(sum(iptype)/(npatches+0.0d0));
    norder_avg = floor(sum(norders)/(npatches+0.0d0));

    [rfac, rfac0] = get_rfacs(norder_avg,iptype_avg);
    % prin2('rfac = *',rfac,1)
    % prin2('rfac0 = *',rfac0,1)

        % allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

    [cms, rads] = get_centroid_rads(npatches,norders,ixyzs,iptype,npts, ... 
       srccoefs);

    rad_near = rads*rfac;

    % 
    % find near quadrature correction interactions
    % 

    targs = srcvals(1:3,:);
    ntarg = npts;

    % prinf('entering find near mem',0,0)
    nnz = findnearmem(cms, npatches, rad_near, ndtarg, targs, npts);
    nnzp1 = nnz+1;
    % prinf('nnz = *',nnz,1);

    [row_ptr,col_ind] = findnear(cms,npatches,rad_near,ndtarg,targs,npts,nnz);

    iquad = get_iquad_rsc(npatches,ixyzs,npts,ntarg,nnz,row_ptr,col_ind);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wts = get_qwts_f(npatches,norders,ixyzs,iptype,ntarg,srcvals);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    nquad = iquad(nnz+1)-1;

    iquadtype = 1;

    wnear = cell(length(iker),1);
    for i = 1:length(iker)
    if iker < 6
        A = zeros(nquad,1,'like',1i);
        if ivpp
            maxdist = 1.5*sqrt((max(S.r(1,:)) - min(S.r(1,:))).^2 + (max(S.r(2,:)) - min(S.r(2,:))).^2);
            mex_id_ = 'getnearquad_flex_all_vpp(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i double[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, maxdist, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, 1, nquad);
        else
            mex_id_ = 'getnearquad_flex_all(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, nquad);
        end
        wnear{i} = A;
    elseif iker == 6
        opts.rep = 'eval';
        opts.format = 'rsc';
        wstruc = lap3d.neumann.get_quadrature_correction(S,eps,S,opts);
        wnear{i} = wstruc.wnear;
        row_ptr = wstruc.row_ptr;
        iquad = wstruc.iquad;
        col_ind = wstruc.col_ind;
        nnz = length(iquad)-1;
        nquad = iquad(nnz+1)-1;
    end

    end

    if length(iker) == 1
        wnear = wnear{1};
    end

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

   nnz = length(irow_ind);
    lbat = 1e3;
    nbat = ceil(nnz/lbat);
    
    for k = 1:nbat
        ks = (lbat*(k-1)+1):min(lbat*k,nnz);
        rsrc = S.r(:,icol_ind(ks));
        rtarg = S.r(:,irow_ind(ks));
        rnear = rtarg - rsrc;
        wnear(ks) = wnear(ks)-gs_kern(struct('r',[0;0;0]),struct('r',rnear)).*S.wts(icol_ind(ks));
    end
   xmat = sparse(irow_ind,icol_ind, wnear, ntarg, S.npts);

end
%
%
%%
