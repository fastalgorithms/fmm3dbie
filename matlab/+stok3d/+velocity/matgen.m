function mat = matgen(S, eps, dpars)
%  stok3d.velocity.matgen
%    This subroutine generates the matrix for the combined
%    field layer potential
%  
%  Syntax
%    mat = stok3d.velocity.matgen(S, eps, dpars)
%
%  Integral representation
%     pot = \alpha S_{stok} [\sigma] + \beta D_{stok} [\sigma]
%
%  S_{stok}, D_{stok}: stokes single and double layer potential
%  
%  \alpha, beta = dpars(1:2)
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * eps: precision requested
%    * dpars: kernel parameters
%        dpars(1) - single layer strength
%        dpars(2) - double layer strength
%
%  Note that the matrix only returns the principal
%  part of the discretization, the identity terms need to be included
%  by the user
%
%  Output:
%    * mat: double(3*n, 3*n)
%        discretized matrix ordered where components of velocity
%        and density at a point are continunous in memory
%        if \sigma_{i,j} is the ith component at the jth point
%        then the ordering is \sigma_{1,1}, \sigma_{2,1}, \sigma_{3,1}
%        \sigma{1,2}, \sigma_{2,2,} \ldots \sigma_{1,n}, \sigma_{2,n}
%        \sigma_{3,n}

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;
    npp1 = npatches+1;
    n3 = 3;
    ntp1 = npts+1;
    
    nnz = S.npatches*S.npts;
    row_ptr = 1:S.npatches:nnz+1;
    col_ind = repmat(1:S.npatches, [1,npts]);
    nnzp1 = nnz+1;
    iquad = zeros(nnz+1,1);
    mex_id_ = 'get_iquad_rsc(c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c io int[x])';
[iquad] = fmm3dbie_routs(mex_id_, npatches, ixyzs, npts, nnz, row_ptr, col_ind, iquad, 1, npp1, 1, 1, ntp1, nnz, nnzp1);
     
    iquadtype = 1;
    rfac0 = 1.25;
    uvs = S.uvs_targ;
    patch_id = S.patch_id;
    nker = 6;
    nquad = iquad(nnz+1)-1;
    wnear = zeros(nker,nquad);

    mex_id_ = 'getnearquad_stok_comb_vel_eval(c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[xx], c i double[xx], c i int[x], c i int[x], c i double[xx], c i int[x], c i double[xx], c i double[x], c i double[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[x], c i int[x], c io double[xx])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, n12, npts, srcvals, patch_id, uvs, eps, dpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, n12, npts, npts, 2, npts, 1, 2, 1, 1, ntp1, nnz, nnzp1, 1, 1, nker, nquad);

    n = S.npts;
    mat = zeros(3*n,3*n);
    mat(1:3:end,1:3:end) = reshape(wnear(1,:),[n,n]).';
    mat(2:3:end,1:3:end) = reshape(wnear(2,:),[n,n]).';
    mat(1:3:end,2:3:end) = reshape(wnear(2,:),[n,n]).';
    mat(3:3:end,1:3:end) = reshape(wnear(3,:),[n,n]).';
    mat(1:3:end,3:3:end) = reshape(wnear(3,:),[n,n]).';
    mat(2:3:end,2:3:end) = reshape(wnear(4,:),[n,n]).';
    mat(3:3:end,2:3:end) = reshape(wnear(5,:),[n,n]).';
    mat(2:3:end,3:3:end) = reshape(wnear(5,:),[n,n]).';
    mat(3:3:end,3:3:end) = reshape(wnear(6,:),[n,n]).';

end  
%------------------
