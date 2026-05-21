function wnear = getnearquad_v2b_flex(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,targ,ipatch_id, uvs_targ, eps,...
            iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad)
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    npp1 = npatches+1;
    nnzp1 = nnz+1;
    ndz = length(zpars);

    
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    wnear = zeros(2,nquad,'like',1i);

    mex_id_ = 'getnearquad_flex_bcs(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[xx])';
[wnear] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 2, nquad);


end


%
%
