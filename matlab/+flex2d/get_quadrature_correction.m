function Q = get_quadrature_correction(S, type, zk, nu, eps, targinfo, opts)
%
%  flex2d.get_quadrature_correction
%    This subroutine returns the near quadrature correction for the
%    2D flexural (thin-plate) volume potential and its plate boundary
%    traces, with density supported on the flat surface S (in the
%    z = 0 plane), and targets given by targinfo, as a sparse
%    matrix/rsc format
%
%  Syntax
%   Q = flex2d.get_quadrature_correction(S,type,zk,nu,eps)
%   Q = flex2d.get_quadrature_correction(S,type,zk,nu,eps,targinfo)
%   Q = flex2d.get_quadrature_correction(S,type,zk,nu,eps,targinfo,opts)
%
%  Kernels, with G the flexural Green's function of flex2d.green:
%     type = 's':                    G
%     type = 'clamped_plate_bcs':    [G; \nabla_{n_x} G]
%     type = 'supported_plate_bcs':  [G; G_{nn} + nu G_{tt}]
%     type = 'free_plate_bcs':       [G_{nn} + nu G_{tt}; free-plate shear]
%  The two-row kernels return Q.wnear of size (2,nquad).
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * type: kernel type, see above
%    * zk: plate wavenumbers [zk1, zk2], or a scalar zk which is
%        promoted to [zk, 1i*zk]
%    * nu: Poisson ratio, required for the supported and free types
%    * eps: precision requested
%    * targinfo: target info (optional, defaults to S)
%       targinfo.r = (3,nt) target locations, only r(1:2,:) is used
%       targinfo.n = (3,nt) target normals, required for the plate bcs types
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv coordinates of target on
%          patch if on-surface (optional)
%    * opts: options struct
%        opts.format - Storage format for sparse matrices
%           'rsc' - row sparse compressed format
%           'csc' - column sparse compressed format
%           'sparse' - sparse matrix format
%        opts.quadtype - quadrature type, currently only 'ggq' supported
%        opts.rsc - precomputed near-field structure from getnear
%          (optional)
%
%       targinfo.kappa (nt,) signed curvature of the boundary at the
%          target, needed for 'free_plate_bcs'; computed from
%          targinfo.d and targinfo.d2 if absent
%

    zpars = complex(zk(:));
    if numel(zpars) == 1
      zpars = [zpars; 1i*zpars];
    end
    if numel(zpars) ~= 2
      error('FLEX2D.GET_QUADRATURE_CORRECTION: zk must be [zk1, zk2] or a scalar.');
    end

    switch lower(type)
      case {'s', 'single'}
        kernel_order = -1;
      case {'clamped_plate_bcs'}
        kernel_order = -1;
      case {'supported_plate_bcs', 'free_plate_bcs'}
        if isempty(nu)
          error('FLEX2D.GET_QUADRATURE_CORRECTION: type ''%s'' requires nu.', type);
        end
        kernel_order = -1;
      otherwise
        error('FLEX2D.GET_QUADRATURE_CORRECTION: unknown type ''%s''.', type);
    end
    dpars = 0;
    if ~isempty(nu)
      dpars = nu;
    end

    [srcvals,srccoefs,norders,ixyzs,iptype,~] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npp1 = npatches+1;

    if nargin < 6 || isempty(targinfo)
      targinfo = S;
    end

    if nargin < 7
      opts = [];
    end

    ff = 'rsc';
    if(isfield(opts,'format'))
       ff = opts.format;
    end

    if(~(strcmpi(ff,'rsc') || strcmpi(ff,'csc') || strcmpi(ff,'sparse')))
       fprintf('invalid quadrature format, reverting to rsc format\n');
       ff = 'rsc';
    end

    if isa(targinfo, 'surfer')
      targinfo = struct('r', targinfo.r, 'du', targinfo.du, ...
          'dv', targinfo.dv, 'n', targinfo.n, ...
          'patch_id', targinfo.patch_id, 'uvs_targ', targinfo.uvs_targ);
    end

    if strcmpi(type,'free_plate_bcs') && ~isfield(targinfo,'kappa')
      if ~(isfield(targinfo,'d') && isfield(targinfo,'d2'))
        error(['FLEX2D.GET_QUADRATURE_CORRECTION: ''free_plate_bcs'' ' ...
               'needs targinfo.kappa, or targinfo.d and targinfo.d2.']);
      end
      d = targinfo.d; d2 = targinfo.d2;
      targinfo.kappa = (d(1,:).*d2(2,:) - d2(1,:).*d(2,:)) ./ ...
          sqrt(d(1,:).^2 + d(2,:).^2).^3;
    end

    targs = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);

    if(isfield(targinfo,'patch_id') && ~isempty(targinfo.patch_id))
      patch_id = targinfo.patch_id;
    else
      patch_id = -1*ones(ntarg,1);
    end

    if(isfield(targinfo,'uvs_targ') && ~isempty(targinfo.uvs_targ))
      uvs_targ = targinfo.uvs_targ;
    else
      uvs_targ = zeros(2,ntarg);
    end

    if(length(patch_id)~=ntarg)
      fprintf('Incorrect size of patch id in target info struct. Aborting! \n');
    end

    [n1,n2] = size(uvs_targ);
    if(n1 ~=2 && n2 ~=ntarg)
      fprintf('Incorrect size of uvs_targ array in targinfo struct. Aborting! \n');
    end

    if isfield(opts,'rsc') && ~isempty(opts.rsc)
        rsc = opts.rsc;
    else
        rsc = getnear(S, targinfo);
    end
    row_ptr = rsc.row_ptr; col_ind = rsc.col_ind; iquad   = rsc.iquad;
    rfac    = rsc.rfac;    rfac0   = rsc.rfac0;   nnz     = rsc.nnz;
    nquad   = rsc.nquad;
    ntp1    = ntarg+1;
    nnzp1   = nnz+1;

    iquadtype = 1;
    if(isfield(opts,'quadtype'))
      if(strcmpi(opts.quadtype,'ggq'))
         iquadtype = 1;
      else
        fprintf('Unsupported quadrature type, reverting to ggq\n');
        iquadtype = 1;
      end
    end

    w1 = complex(zeros(nquad,1));
    w2 = complex(zeros(nquad,1));
    switch lower(type)
      case {'s', 'single'}
        mex_id_ = 'getnearquad_flex2d_dir(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[w1] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, w1, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 2, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
        wnear = w1;
        nker = 1;
      case {'clamped_plate_bcs'}
        mex_id_ = 'getnearquad_flex2d_dir(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[w1] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, w1, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 2, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
        mex_id_ = 'getnearquad_flex2d_neu(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[w2] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, w2, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 2, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
        wnear = [w1.'; w2.'];
        nker = 2;
      case {'supported_plate_bcs'}
        mex_id_ = 'getnearquad_flex2d_dir(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[w1] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, w1, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 2, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
        mex_id_ = 'getnearquad_flex2d_supp2(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[w2] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, w2, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 2, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
        wnear = [w1.'; w2.'];
        nker = 2;
      case {'free_plate_bcs'}
        mex_id_ = 'getnearquad_flex2d_supp2(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[w1] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, w1, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 2, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
        mex_id_ = 'getnearquad_flex2d_free2(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[w2] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, w2, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 2, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
        wnear = [w1.'; w2.'];
        nker = 2;
    end

    Q = [];
    Q.targinfo = targinfo;
    Q.ifcomplex = 1;
    Q.wavenumber = zpars;
    Q.kernel_order = kernel_order;
    Q.rfac = rfac;
    Q.nquad = nquad;
    Q.format = ff;

    if(strcmpi(ff,'rsc'))
        Q.iquad = iquad;
        Q.wnear = wnear;
        Q.row_ptr = row_ptr;
        Q.col_ind = col_ind;
    elseif(strcmpi(ff,'csc'))
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
        if nker == 1
          spmat = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear);
        else
          spmat = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear, ...
              kernel3d.rsc_interleave_full(nker, 1));
        end
        Q.spmat = spmat;
    end

end
%
%
%
%
%-------------------------------------------------
