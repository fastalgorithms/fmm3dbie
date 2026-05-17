function A = v2b_neu(S,targinfo,eps,ipatch_id,uvs_targ)
%
%  lap2d.v2b_neu
%
%  Syntax
%   A = lap2d.v2b_neu(S,targinfo,zpars,eps)
%
%  Integral representation
%     pot = \int G_S(r,r') \sigma(r') dA(r')
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * targinfo: target info
%       targinfo.r = (3,nt) target locations
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%          patch if on-surface (optional)
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

    if nargin < 5
    ipatch_id = zeros(ntarg,1);
    uvs_targ = zeros(2,ntarg);
    end

    iquad = 1:npols:(npts*ntarg+1);
    nquad = iquad(nnzp1)-1;
    A = zeros(1,nquad);
    zpars = 0;

    mex_id_ = 'getnearquad_lap2d_gdn(i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i double[xx], i int64_t[x], i int64_t[x], i double[xx], i int64_t[x], i double[xx], i double[x], i dcomplex[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[x], i int64_t[x], io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);

    A = reshape(A,[S.npts size(targinfo.r,2)]).';
end
%
%
%
%
%
%
