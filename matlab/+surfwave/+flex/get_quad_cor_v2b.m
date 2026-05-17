function [xmat1,xmat2,wnear] = get_quad_cor_v2b(S, targ, kern, eps, zpars, uv_bndry, patch_id)

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npatches = S.npatches;
    npp1 = npatches+1;
    
    rx = targ.r(1,:); 
    ry = targ.r(2,:); 
    
    nx = targ.n(1,:); 
    ny = targ.n(2,:);
    
    dx = targ.d(1,:);
    dy = targ.d(2,:);
    
    ds = sqrt(dx.*dx+dy.*dy);
    taux = (dx./ds); % normalization
    tauy = (dy./ds);
    
    targs = zeros(13,length(targ.r(1,:)));
    targs(1:2,:) = [rx; ry];
    targs(4:5,:) = [taux; tauy];
    targs(10:11,:) = [nx; ny];
    targs(13,:) = targ.kappa(:);

    [ndtarg,ntarg] = size(targs);
    ntp1 = ntarg+1;

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

    targinfo = [];
    targinfo.r = targs(1:3,:);

    if nargin < 6
        ipatch_id = zeros(ntarg,1);
        uvs_targ = zeros(2,ntarg);
    else
        ipatch_id = patch_id; 
        uvs_targ = uv_bndry;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    nquad = iquad(nnz+1)-1;

    iquadtype = 1;

    wnear = zeros(2,nquad,'like',1i);

    mex_id_ = 'getnearquad_flex_bcs(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i int[x], i double[xx], i double[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[x], i int[x], io dcomplex[xx])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 2, nquad);

    xmat1 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(1,:).');
    xmat2 = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear(2,:).');

    nover = S.norders(1) + 4; 
    S_over = oversample(S,nover);
    Asmth = cell(2,1);

    rnodes = koorn.rv_nodes(S.norders(1));
    rnodes_over = koorn.rv_nodes(nover);

    vmat = koorn.coefs2vals(S.norders(1),rnodes_over);
    umat = koorn.vals2coefs(S.norders(1),rnodes);  
    Ainterp = vmat*umat;

    fprintf('here\n');

    for k = 1:2
        fprintf('k=%d\n',k)
    
        I = [];
        J = [];
        v = [];
        for ii = 1:ntarg
            rtarg = [];
            rtarg.r = targ.r(:,ii);
            rtarg.n = targ.n(:,ii);
            rtarg.d = targ.d(:,ii);
            rtarg.d2 = targ.d2(:,ii);
            for j = row_ptr(ii):(row_ptr(ii+1)-1)
                jj = col_ind(j);
                jds_over = S_over.ixyzs(jj):(S_over.ixyzs(jj+1)-1);
                jds = S.ixyzs(jj):(S.ixyzs(jj+1)-1);
                rsrc = [];
                rsrc.r = S_over.r(:,jds_over);
                Aloc = kern(k).eval(rsrc,rtarg).*S_over.wts(jds_over).';
                Aloc = Aloc*Ainterp;
                I = [I ii*ones(1,length(jds))];
                J = [J jds];
                v = [v Aloc(:).'];  
            end
        end
        Asmth{k} = sparse(I,J,v,ntarg,npts);

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
%
%
%
%
%
%
