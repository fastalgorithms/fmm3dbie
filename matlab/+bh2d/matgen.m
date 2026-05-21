function A = matgen(S, zk, type, targinfo, eps, ipatch_id, uvs_targ)
%
%  bh2d.matgen
%
%  Syntax
%   A = bh2d.matgen(S, zk, type, targinfo, eps)
%   A = bh2d.matgen(S, zk, type, targinfo, eps, ipatch_id, uvs_targ)
%
%  Computes near-quadrature weight matrices for the biharmonic equation
%  (zk = 0 only). The type string selects the operator:
%
%   v2v types (no targinfo needed, pass []):
%     'v2v'        volume-to-volume (G_S single-layer, self-interaction)
%
%   v2b types (targinfo = boundary target struct):
%     'dir'        volume-to-boundary Dirichlet (SLP)
%     'neu'        volume-to-boundary Neumann (SLP normal derivative)
%     'supp2'      volume-to-boundary simply-supported, part 2 (needs nu in type as 3rd token, or pass as targinfo.nu)
%     'free2'      volume-to-boundary free plate, part 2 (needs nu)
%     'clamped'    volume-to-boundary clamped (returns [npts x 2*ntarg])
%     'gdn'        normal derivative of G_S to boundary targets
%
%   b2v types (targinfo = volume target struct, sources on boundary):
%     'b2v_dir'     boundary-to-volume Dirichlet (SLP)
%     'b2v_neu'     boundary-to-volume Neumann (SLP normal derivative)
%     'b2v_clamped' boundary-to-volume clamped (returns [ntarg x 2*npts])
%
%  Input arguments:
%    * S: surfer object
%    * zk: wavenumber (must be 0)
%    * type: string selecting the operator (see above)
%    * targinfo: target struct with fields r, patch_id, uvs_targ (optional for v2v)
%    * eps: precision requested
%    * nu: Poisson ratio, required for supp2 and free2 types; pass as
%          targinfo.nu or as an extra element in type, e.g. 'supp2'
%          with targinfo.nu set.
%
%  Note: for on-surface targets, only the principal value part is returned.
%

    if abs(zk) > 1e-8
        error('bh2d.matgen: only zk = 0 (biharmonic) is supported')
    end

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~]     = size(srccoefs);
    [npatches,~] = size(norders);
    npp1 = npatches + 1;
    zpars = 0;
    rfac0 = 1.25;
    iquadtype = 1;
    npols = ixyzs(2) - ixyzs(1);

    % --- v2v self-interaction (no targinfo) ---
    if strcmp(type, 'v2v')
        ntp1 = npts + 1;
        row_ptr = 1:npatches:(npatches*npts+1);
        col_ind = repmat(1:npatches,[1,npts]);
        row_ptr = row_ptr(:);
        col_ind = col_ind(:);
        nnz     = npatches*npts;
        nnzp1   = nnz + 1;
        iquad   = 1:npols:(npts*npts+1);
        nquad   = npts*npts;
        A = zeros(npts, npts);

        mex_id_ = 'getnearquad_bh2d_gv2v(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);

        A = reshape(A, [npts, npts]).';
        return
    end

    % --- all other types: need targinfo ---
    [targs]       = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg + 1;

    row_ptr = 1:npatches:(npatches*ntarg+1);
    col_ind = repmat(1:npatches,[1,ntarg]);
    row_ptr = row_ptr(:);
    col_ind = col_ind(:);
    nnz   = npatches*ntarg;
    nnzp1 = nnz + 1;

    if nargin < 6
        ipatch_id = zeros(ntarg,1);
        uvs_targ  = zeros(2,ntarg);
    end

    iquad = 1:npols:(npts*ntarg+1);
    nquad = iquad(nnzp1)-1;

    if strcmp(type, 'dir')
        A = zeros(1, nquad);
        mex_id_ = 'getnearquad_bh2d_dir(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        A = reshape(A, [npts, ntarg]).';

    elseif strcmp(type, 'neu')
        A = complex(zeros(1, nquad));
        mex_id_ = 'getnearquad_bh2d_neu(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        A = reshape(A, [npts, ntarg]).';

    elseif strcmp(type, 'supp2')
        dpars = targinfo.nu;
        A = zeros(1, nquad);
        mex_id_ = 'getnearquad_bh2d_supp2(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        A = reshape(A, [npts, ntarg]).';

    elseif strcmp(type, 'free2')
        dpars = targinfo.nu;
        A = zeros(1, nquad);
        mex_id_ = 'getnearquad_bh2d_free2(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        A = reshape(A, [npts, ntarg]).';

    elseif strcmp(type, 'clamped')
        A = zeros(2, nquad);
        mex_id_ = 'getnearquad_bh2d_clamped(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[xx])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 2, nquad);
        A = reshape(A, [npts, 2*ntarg]).';

    elseif strcmp(type, 'gdn')
        A = zeros(1, nquad);
        mex_id_ = 'getnearquad_bh2d_gdn(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        A = reshape(A, [npts, ntarg]).';

    elseif strcmp(type, 'b2v_dir')
        A = zeros(1, nquad);
        mex_id_ = 'getnearquad_bh2d_b2v_dir(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        A = reshape(A, [npts, ntarg]).';

    elseif strcmp(type, 'b2v_neu')
        A = zeros(1, nquad);
        mex_id_ = 'getnearquad_bh2d_v2b_neu(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, nquad);
        A = reshape(A, [npts, ntarg]).';

    elseif strcmp(type, 'b2v_clamped')
        A = zeros(2, nquad);
        mex_id_ = 'getnearquad_bh2d_clamped(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[xx])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntargp1, nnz, nnzp1, 1, 1, 2, nquad);
        A = reshape(A, [2*npts, ntarg]).';

    else
        error('bh2d.matgen: unknown type ''%s''', type)
    end

end
