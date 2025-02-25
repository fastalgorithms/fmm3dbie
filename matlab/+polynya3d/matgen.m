function A = matgen(S, eps)
%
%  polynya3d.dirichlet.get_quadrature_correction
%
%  Syntax
%   A = polynya3d.dirichlet.get_quadrature_correction(S,eps)
%   Q = polynya3d.dirichlet.get_quadrature_correction(S,eps,targinfo)
%   Q = polynya3d.dirichlet.get_quadrature_correction(S,eps,targinfo,opts)
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
    n3 = 3;
    row_ptr = 1:npatches:(npatches*npts+1);
    col_ind = repmat(1:npatches,[1,npts]);
    row_ptr = row_ptr(:);
    col_ind = col_ind(:);
    zpuse = 1j;
    nnz = npatches*npts;
    nnzp1 = nnz + 1;
    iquadtype = 1;
    ntp1 = npts + 1;
    npols = ixyzs(2) - ixyzs(1);
    iquad = 1:npols:(npts*npts+1);
    A = complex(zeros(npts,npts));
    rfac0 = 1.25;
    nquad = npts*npts;
      mex_id_ = 'getnearquad_polynya_lps_neu(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[x], i int[x], io dcomplex[xx])';
[A] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpuse, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, 1, 1, ntp1, nnz, nnzp1, 1, 1, npts, npts);

    
end
%
%
%
