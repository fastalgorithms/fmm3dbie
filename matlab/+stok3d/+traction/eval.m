function p = eval(S, sigma, targinfo, eps, varargin)
%
%  stok3d.traction.eval
%    Evaluates the traction of the Stokes single layer (S' kernel)
%    at a collection of targets.
%
%  Syntax
%   pot = stok3d.traction.eval(S, sigma, targinfo, eps)
%   pot = stok3d.traction.eval(S, sigma, targinfo, eps, opts)
%
%  Operator:
%     pot = S'_{stok}[\sigma](x)
%         = n_k(x) sigma_{ik}(S[\sigma])(x)
%
%  where sigma_{ij}(u) = -p delta_{ij} + (du_i/dx_j + du_j/dx_i)
%  is the Stokes stress tensor and n(x) is the target outward normal.
%
%  Note: for targets on surface, only principal value part is returned.
%  The jump term +- sigma/2 must be added separately by the caller.
%
%  Input arguments:
%    * S: surfer object
%    * sigma: (3, ns) layer potential density
%    * targinfo: target info
%       targinfo.r = (3,nt) target locations
%       targinfo.n = (3,nt) outward unit normals at targets  [REQUIRED]
%       targinfo.patch_id (nt,) patch id, = -1 if off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv coords on patch (optional)
%    * eps: precision requested
%    * opts: options struct (optional)
%        opts.nonsmoothonly - use smooth quadrature rule (false)
%        opts.precomp_quadrature: precomputed quadrature struct in rsc format
%

    if(nargin < 5)
      opts = [];
    else
      opts = varargin{1};
    end

    isprecompq = false;
    if isfield(opts, 'precomp_quadrature')
      Q = opts.precomp_quadrature;
      isprecompq = true;
    end

    nonsmoothonly = false;
    if(isfield(opts,'nonsmoothonly'))
      nonsmoothonly = opts.nonsmoothonly;
    end

    if(isprecompq)
      if ~(strcmpi(Q.format,'rsc'))
        fprintf('Invalid precomputed quadrature format\n');
        fprintf('Ignoring quadrature corrections\n');
        opts_qcorr = [];
        opts_qcorr.type = 'double';
        Q = init_empty_quadrature_correction(targinfo,opts_qcorr);
      end
    end


% Extract arrays
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;

    ff = 'rsc';

    [targs] = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;

% Compute quadrature corrections
    if ~isprecompq
      if ~nonsmoothonly
        opts_quad = [];
        opts_quad.format = 'rsc';
        [Q] = stok3d.traction.get_quadrature_correction(S,eps,[],targinfo,opts_quad);
      else
        opts_qcorr = [];
        opts_qcorr.type = 'double';
        Q = init_empty_quadrature_correction(targinfo,opts_qcorr);
      end
    end
    nquad = Q.iquad(end)-1;
    nnz = length(Q.col_ind);
    nnzp1 = nnz+1;

    [novers] = get_oversampling_parameters(S,Q,eps);
    Sover = oversample(S,novers);


% Extract oversampled arrays

    [srcover,~,~,ixyzso,~,wover] = extract_arrays(Sover);
    nptso = Sover.npts;

% Extract quadrature arrays
    row_ptr = Q.row_ptr;
    col_ind = Q.col_ind;
    iquad = Q.iquad;
    wnear = Q.wnear;

    p = zeros(3,ntarg);

    ndd = 1;
    dpars = [0.0];
    ndz = 0;
    zpars = [];
    ndi = 1;
    ipars = [0];
    nker = 6;
    lwork = 0;
    work = [];
    ndim = 3;
% Call the layer potential evaluator
    mex_id_ = 'stok_sprime_eval_addsub(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i double[xx], c io double[xx])';
[p] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, sigma, p, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, 1, 1, ndd, 1, ndz, 1, ndi, 1, ntargp1, nnz, nnzp1, 1, 1, nker, nquad, npatches, 1, npatp1, 12, nptso, nptso, 1, lwork, 1, ndim, npts, ndim, ntarg);
end
%
%
%
%----------------------------------
%

