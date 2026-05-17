function [xmat1,xmat2] = get_quad_cor_v2b(S, targs, kern, eps, zpars, uv_bndry, patch_id, norderup)

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npatches = S.npatches;
    npp1 = npatches+1;

    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;


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

    % prinf('entering find near mem',0,0)
    nnz = findnearmem(cms, npatches, rad_near, ndtarg, targs, ntarg);
    nnzp1 = nnz+1;
    % prinf('nnz = *',nnz,1);

    [row_ptr,col_ind] = findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,nnz);

    iquad = get_iquad_rsc(npatches,ixyzs,npts,ntarg,nnz,row_ptr,col_ind);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wts = get_qwts_f(npatches,norders,ixyzs,iptype,ntarg,srcvals);

    if nargin < 6
        ipatch_id = zeros(ntarg,1);
        uvs_targ = zeros(2,ntarg);
    else
        if ~isempty(patch_id)
            ipatch_id = patch_id; 
            uvs_targ = uv_bndry;
        else
            ipatch_id = zeros(ntarg,1);
            uvs_targ = zeros(2,ntarg);
        end
    end

    targinfo = [];
    targinfo.r = targs(1:3,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    nquad = iquad(nnz+1)-1;

    iquadtype = 1;

    wnear = zeros(2,nquad,'like',1i);

    mex_id_ = 'getnearquad_capillary_bcs(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[xx])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 6, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 2, nquad);
    xmat1 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(1,:).');
    xmat2 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(2,:).');

    if nargin < 8
        norderup = 4;
    end
    nover = S.norders(1) + norderup; 
    S_over = oversample(S,nover);
    Asmth = cell(1,1);

    rnodes = koorn.rv_nodes(S.norders(1));
    rnodes_over = koorn.rv_nodes(nover);

    vmat = koorn.coefs2vals(S.norders(1),rnodes_over);
    umat = koorn.vals2coefs(S.norders(1),rnodes);  
    Ainterp = vmat*umat;

    for k = 1:2
    Asmth{k} = sparse(ntarg,npts);

    for ii = 1:ntarg
        rtarg = [];
        rtarg.r = targs(1:2,ii);
        rtarg.n = targs(10:11,ii);
        for j = row_ptr(ii):(row_ptr(ii+1)-1)
            jj = col_ind(j);
            jds_over = S_over.ixyzs(jj):(S_over.ixyzs(jj+1)-1);
            jds = S.ixyzs(jj):(S.ixyzs(jj+1)-1);
            rsrc = [];
            rsrc.r = S_over.r(:,jds_over);
            Aloc = kern(k).eval(rsrc,rtarg).*S_over.wts(jds_over).';
            Aloc = Aloc*Ainterp;
            Asmth{k}(ii,jds) = Aloc;
        end
    end

    end

    xmat1 = xmat1 - Asmth{1};
    xmat2 = xmat2 - Asmth{2};

    % wnear = cell(length(iker),1);
    % for i = 1:length(iker)
    % wnear{i} = surfwave.flex.getnearquad_flexural(npatches,norders,ixyzs, ...
    %       iptype,npts,srccoefs,srcvals,eps,iquadtype,nnz, ...
    %       row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker(i)) ;
    % 
    % end    
    % prin2('quadrature generation time=*',t2,1)
    
    % fprintf('done generating near quad \n')
    
    %
    %           .   .   .   now correct
        
    % fprintf('correcting quadrature \n')

   
end
%
%
%
