function A = matgen(S,type, zpars0, eps, ivpp)
%
%  surfwave.flex.matgen
%
%  Syntax
%   A = surfwave.flex.matgen(S,zpars,eps)
%
%  Integral representation
%     pot = (Update representation here)
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * zpars: parameters. zpars(1:3) are the roots, zpars(4:6) are the residues
%    * eps: precision requested
%

    if strcmp(type,'gs')
        iker = 1;
    elseif strcmp(type,'gphi')
        iker = 2;
    elseif strcmp(type,'bilapgs')
        iker = 3;
    elseif strcmp(type,'bilapgphi')
        iker = 4;
    else
        error('kernel name not recognized') 
    end 

    if nargin < 5
        ivpp = 1;
    end

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npp1 = npatches+1;
    row_ptr = 1:npatches:(npatches*npts+1);
    col_ind = repmat(1:npatches,[1,npts]);
    row_ptr = row_ptr(:);
    col_ind = col_ind(:);
    nnz = npatches*npts;
    nnzp1 = nnz + 1;
    iquadtype = 1;
    ntp1 = npts + 1;
    npols = ixyzs(2) - ixyzs(1);
    iquad = 1:npols:(npts*npts+1);
    A = complex(zeros(npts,npts));
    rfac0 = 1.25;
    nquad = npts*npts;
    
    if ivpp
    maxdist = 1.5*sqrt((max(S.r(1,:)) - min(S.r(1,:))).^2 + (max(S.r(2,:)) - min(S.r(2,:))).^2);
    mex_id_ = 'getnearquad_flex_all_vpp(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i double[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars0, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, maxdist, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, 1, nquad);
    else
    mex_id_ = 'getnearquad_flex_all(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars0, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, nquad);
    end

    A = reshape(A, [S.npts,S.npts]).';

end
%
%
