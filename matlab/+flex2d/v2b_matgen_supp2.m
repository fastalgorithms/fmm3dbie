function A = v2b_matgen_supp2(S,zk,nu,targinfo,eps,ipatch_id,uvs_targ)
%
%  flex2d.v2b_matgen_supp2
%
%  Syntax
%   A = flex2d.v2b_matgen_supp2(S,zk,nu,targinfo,eps,ipatch_id,uvs_targ)
%
%  Integral representation
%     pot = (Update representation here)
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * eps: precision requested
%

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;
    npp1 = npatches+1;

    [targs] = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;

    n3 = 3;
    row_ptr = 1:npatches:(npatches*ntarg+1);
    col_ind = repmat(1:npatches,[1,ntarg]);
    row_ptr = row_ptr(:);
    col_ind = col_ind(:);
    zpuse = 1j;
    nnz = npatches*ntarg;
    nnzp1 = nnz + 1;
    iquadtype = 1;
    ntp1 = npts + 1;
    npols = ixyzs(2) - ixyzs(1);

    rfac0 = 1.25;

    if nargin < 6
    ipatch_id = zeros(ntarg,1);
    uvs_targ = zeros(2,ntarg);
    end

    iquad = 1:npols:(npts*ntarg+1);
    nquad = iquad(nnzp1)-1;


    dpars = nu;
    
    if isscalar(zk)
        if abs(zk) < 1e-6
            A = zeros(1,nquad);
            zpars = 0;
        mex_id_ = 'getnearquad_bh2d_supp2(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        
        else 

            zpars = zeros(1,2);
            zpars(1) = zk; 
            zpars(2) = zk*1i; 
            A = zeros(1,nquad,'like',1i);
        mex_id_ = 'getnearquad_flex2d_supp2(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 2, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);

        end


    elseif length(zk) == 2

        if any(abs(zk) < 1e-6)

            error('No support for Stokes problem')
            
             

        else 

            zpars = zeros(1,2);
            zpars(1) = zk(1); 
            zpars(2) = zk(2); 
            A = zeros(1,nquad,'like',1i);
        mex_id_ = 'getnearquad_flex2d_supp2(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 2, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);


        end

    end 

    
    
    
    A = reshape(A,[S.npts size(targinfo.r,2)]).';

end
%
%
%
%
function A = v2b_matgen_free2(S,zk,nu,targinfo,eps,ipatch_id,uvs_targ)
%
%  flex2d.v2b_matgen_free2
%
%  Syntax
%   A = flex2d.v2b_matgen_free2(S,zk,nu,targinfo,eps,ipatch_id,uvs_targ)
%
%  Integral representation
%     pot = (Update representation here)
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * eps: precision requested
%

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;
    npp1 = npatches+1;

    [targs] = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;

    n3 = 3;
    row_ptr = 1:npatches:(npatches*ntarg+1);
    col_ind = repmat(1:npatches,[1,ntarg]);
    row_ptr = row_ptr(:);
    col_ind = col_ind(:);
    zpuse = 1j;
    nnz = npatches*ntarg;
    nnzp1 = nnz + 1;
    iquadtype = 1;
    ntp1 = npts + 1;
    npols = ixyzs(2) - ixyzs(1);

    rfac0 = 1.25;

    if nargin < 6
    ipatch_id = zeros(ntarg,1);
    uvs_targ = zeros(2,ntarg);
    end

    iquad = 1:npols:(npts*ntarg+1);
    nquad = iquad(nnzp1)-1;


    dpars = nu;
       
    if isscalar(zk)
        if abs(zk) < 1e-6
            A = zeros(1,nquad);
            zpars = 0;
        mex_id_ = 'getnearquad_bh2d_free2(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        
        else 

            zpars = zeros(1,2);
            zpars(1) = zk; 
            zpars(2) = zk*1i; 
            A = zeros(1,nquad,'like',1i);
        mex_id_ = 'getnearquad_flex2d_free2(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 2, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);

        end


    elseif length(zk) == 2

        if any(abs(zk) < 1e-6)

            error('No support for Stokes problem')
            
             

        else 

            zpars = zeros(1,2);
            zpars(1) = zk(1); 
            zpars(2) = zk(2); 
            A = zeros(1,nquad,'like',1i);
        mex_id_ = 'getnearquad_flex2d_free2(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 2, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);


        end

    end 
    
    A = reshape(A,[S.npts size(targinfo.r,2)]).';

end
%
%
%
%
