function [novers,varargout] = get_oversampling_parameters(S,Q,eps)
%
%  get_oversampling_parameters
%     subroutine to estimate the oversampling paramters for a given
%     surface, and set of quadrature corrections stored in Q
%  
%  Syntax
%    Q = get_oversampling_parameters(S,Q,eps)
%
%  Input arguments
%    * S: surfer object, see README.md in matlab for details
%    * Q: quadrature correction struct, necessary components are
%            Q.targinfo, Q.rfac, Q.wavenumber, Q.row_ptr, Q.col_ind
%    * eps: tolerance
%
    
%
%  extract arrays
%
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;

    n3 = 3;

    cms = zeros(3,npatches);
    rads = zeros(npatches,1);
    mex_id_ = 'get_centroid_rads(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], io double[xx], io double[x])';
[cms, rads] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, cms, rads, 1, npatches, npatp1, npatches, 1, n9, npts, n3, npatches, npatches);

    novers = zeros(npatches,1);
    ixyzso = zeros(npatp1,1);

    targinfo = Q.targinfo;

    targs = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);

    ntp1 = ntarg + 1;
    zk = complex(Q.wavenumber);
    ikerorder = Q.kernel_order;
    row_ptr = Q.row_ptr;
    col_ind = Q.col_ind;
    nnz = length(col_ind);
    rfac = Q.rfac;


    mex_id_ = 'get_far_order(i double[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[x], i int[x], i double[xx], i int[x], i int[x], i double[xx], i int[x], i dcomplex[x], i int[x], i int[x], i int[x], i double[x], io int[x], io int[x])';
[novers, ixyzso] = fmm3dbie_routs(mex_id_, eps, npatches, norders, ixyzs, iptype, cms, rads, npts, srccoefs, ndtarg, ntarg, targs, ikerorder, zk, nnz, row_ptr, col_ind, rfac, novers, ixyzso, 1, 1, npatches, npatp1, npatches, 3, npatches, npatches, 1, n9, npts, 1, 1, ndtarg, ntarg, 1, 1, 1, ntp1, nnz, 1, npatches, npatp1);
    varargout{1} = ixyzso;

end


%
%
%
%
%
%-------------------------------------------------
