function Q = getnearquad_mex(S, name, eps, pars, targinfo, opts)
%GETNEARQUAD_MEX  Near-field quadrature for a registered scalar kernel.
%
%  Syntax:
%    Q = getnearquad_mex(S, name, eps)
%    Q = getnearquad_mex(S, name, eps, pars)
%    Q = getnearquad_mex(S, name, eps, pars, targinfo)
%    Q = getnearquad_mex(S, name, eps, pars, targinfo, opts)
%
%  Kernels are registered with !@kernel annotations in src/ and listed
%  by LIST_QUAD_KERNELS. Only scalar kernels are supported.
%
%  Input arguments:
%    * S: surfer object
%    * name: registry name of the kernel, e.g. 'lap3d_s', 'helm3d_c'
%    * eps: requested precision
%    * pars: (optional) struct with fields dpars, zpars, ipars
%    * targinfo: (optional) target info, default S
%    * opts: (optional) struct
%        opts.format   - 'rsc' (default), 'csc' or 'sparse'
%        opts.corrections - if true, subtract the smooth rule oversampled
%                        to opts.nover (default false)
%        opts.nover    - oversampling order, scalar or (npatches,1)
%
%  Output arguments:
%    * Q: struct
%        Q.name, Q.ifcomplex, Q.targinfo, Q.rfac, Q.nquad, Q.format,
%        Q.corrections, Q.nover
%        rsc:    Q.row_ptr, Q.col_ind, Q.iquad, Q.wnear (nquad,1)
%        csc:    Q.col_ptr, Q.row_ind, Q.iper, Q.iquad, Q.wnear
%        sparse: Q.spmat, (ntarg x S.npts)
%
%  See also LIST_QUAD_KERNELS, GETNEAR, CONV_RSC_TO_SPMAT

    if nargin < 4 || isempty(pars), pars = []; end
    if nargin < 5 || isempty(targinfo), targinfo = S; end
    if nargin < 6 || isempty(opts), opts = []; end

    if ~(ischar(name) || (isstring(name) && isscalar(name)))
        error('GETNEARQUAD_MEX: name must be a character vector.');
    end
    name = char(name);
    kname = double(name(:));
    nchar = numel(kname);

    ipv = 0;
    ifcomplex = 0;
    ier = 0;
    mex_id_ = 'kernel_registry_info(c i int64_t[x], c i int64_t[x], c io int64_t[x], c io int64_t[x], c io int64_t[x])';
[ipv, ifcomplex, ier] = fmm3dbie_routs(mex_id_, nchar, kname, ipv, ifcomplex, ier, 1, nchar, 1, 1, 1);
    if ier ~= 0
        error('GETNEARQUAD_MEX: unknown kernel ''%s'' (registered: %s).', ...
            name, strjoin(list_quad_kernels(), ', '));
    end

    % mwrap needs nonempty arrays
    dpars = 0; zpars = complex(0); ipars = 0;
    if isfield(pars,'dpars') && ~isempty(pars.dpars), dpars = pars.dpars(:); end
    if isfield(pars,'zpars') && ~isempty(pars.zpars), zpars = complex(pars.zpars(:)); end
    if isfield(pars,'ipars') && ~isempty(pars.ipars), ipars = double(pars.ipars(:)); end
    ndd = numel(dpars); ndz = numel(zpars); ndi = numel(ipars);

    ff = 'rsc';
    if isfield(opts,'format'), ff = opts.format; end
    if ~(strcmpi(ff,'rsc') || strcmpi(ff,'csc') || strcmpi(ff,'sparse'))
        fprintf('invalid quadrature format, reverting to rsc format\n');
        ff = 'rsc';
    end

    ifsub = 0;
    if isfield(opts,'corrections') && opts.corrections, ifsub = 1; end

    [srcvals,srccoefs,norders,ixyzs,iptype,~] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npp1 = npatches+1;

    novers = norders;
    if ifsub
        if ~isfield(opts,'nover') || isempty(opts.nover) || any(isnan(opts.nover))
            error('GETNEARQUAD_MEX: opts.nover is required when opts.corrections is true.');
        end
        novers = opts.nover(:);
        if numel(novers) == 1, novers = novers*ones(npatches,1); end
        if numel(novers) ~= npatches
            error('GETNEARQUAD_MEX: opts.nover must be a scalar or have npatches entries.');
        end
    end

    targs = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);

    if (isfield(targinfo,'patch_id') || isprop(targinfo,'patch_id')) && ...
            ~isempty(targinfo.patch_id)
        patch_id = targinfo.patch_id(:);
    else
        patch_id = -ones(ntarg,1);
    end
    if (isfield(targinfo,'uvs_targ') || isprop(targinfo,'uvs_targ')) && ...
            ~isempty(targinfo.uvs_targ)
        uvs_targ = targinfo.uvs_targ;
    else
        uvs_targ = zeros(2,ntarg);
    end
    if numel(patch_id) ~= ntarg || ~isequal(size(uvs_targ), [2 ntarg])
        error('GETNEARQUAD_MEX: patch_id/uvs_targ do not match the number of targets.');
    end

    rsc = getnear(S, targinfo);
    row_ptr = rsc.row_ptr; col_ind = rsc.col_ind; iquad = rsc.iquad;
    rfac0 = rsc.rfac0; nnz = rsc.nnz; nquad = rsc.nquad;
    ntp1 = ntarg+1;
    nnzp1 = nnz+1;

    ier = 0;
    if ifcomplex
        wnear = complex(zeros(nquad,1));
        mex_id_ = 'zgetnearquad_surf_generic(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io dcomplex[x], c io int64_t[x])';
[wnear, ier] = fmm3dbie_routs(mex_id_, nchar, kname, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, ifsub, novers, wnear, ier, 1, nchar, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, ndd, 1, ndz, 1, ndi, 1, ntp1, nnz, nnzp1, 1, 1, 1, npatches, nquad, 1);
    else
        wnear = zeros(nquad,1);
        mex_id_ = 'dgetnearquad_surf_generic(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[x], c io int64_t[x])';
[wnear, ier] = fmm3dbie_routs(mex_id_, nchar, kname, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, ifsub, novers, wnear, ier, 1, nchar, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, ndd, 1, ndz, 1, ndi, 1, ntp1, nnz, nnzp1, 1, 1, 1, npatches, nquad, 1);
    end
    if ier ~= 0
        error('GETNEARQUAD_MEX: kernel ''%s'' failed with ier = %d.', name, ier);
    end

    Q = [];
    Q.name = name;
    Q.ifcomplex = ifcomplex;
    Q.targinfo = targinfo;
    Q.rfac = rsc.rfac;
    Q.nquad = nquad;
    Q.format = ff;
    Q.corrections = logical(ifsub);
    Q.nover = [];
    if ifsub, Q.nover = novers; end

    if strcmpi(ff,'rsc')
        Q.iquad = iquad;
        Q.wnear = wnear;
        Q.row_ptr = row_ptr;
        Q.col_ind = col_ind;
    elseif strcmpi(ff,'csc')
        col_ptr = zeros(npatches+1,1);
        row_ind = zeros(nnz,1);
        iper = zeros(nnz,1);
        npatp1 = npatches+1;
        mex_id_ = 'rsc_to_csc(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io int64_t[x], c io int64_t[x], c io int64_t[x])';
[col_ptr, row_ind, iper] = fmm3dbie_routs(mex_id_, npatches, ntarg, nnz, row_ptr, col_ind, col_ptr, row_ind, iper, 1, 1, 1, ntp1, nnz, npatp1, nnz, nnz);
        Q.iquad = iquad;
        Q.iper = iper;
        Q.wnear = wnear;
        Q.col_ptr = col_ptr;
        Q.row_ind = row_ind;
    else
        irow = zeros(nquad,1);
        jcol = zeros(nquad,1);
        mex_id_ = 'rsc_to_ij(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io int64_t[x], c io int64_t[x])';
[irow, jcol] = fmm3dbie_routs(mex_id_, npatches, ixyzs, ntarg, nnz, row_ptr, col_ind, iquad, nquad, irow, jcol, 1, npp1, 1, 1, ntp1, nnz, nnzp1, 1, nquad, nquad);
        Q.spmat = sparse(irow, jcol, wnear, ntarg, npts);
    end

end
%
%
%
%
%-------------------------------------------------
