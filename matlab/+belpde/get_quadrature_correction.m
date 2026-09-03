function Q = get_quadrature_correction(S, type, zk, eps, targinfo, opts)
%
%  belpde.get_quadrature_correction
%    This subroutine returns the near quadrature correction
%    for the parametrix or the remainder of the Laplace-Beltrami or
%    Helmholtz-Beltrami operator, with density supported on the
%    surface S, and targets given by targinfo, as a sparse matrix/rsc
%    format
%
%  Syntax
%   Q = belpde.get_quadrature_correction(S,type,zk,eps)
%   Q = belpde.get_quadrature_correction(S,type,zk,eps,targinfo)
%   Q = belpde.get_quadrature_correction(S,type,zk,eps,targinfo,opts)
%
%  Kernels
%     type = 'klb':  K(x,y) = log|x-y|^2/(4 pi)
%     type = 'rlb':  R(x,y) = \Delta_{\Gamma} K(x,y)
%     type = 'khb':  K(x,y) = H_0(zk|x-y|)/(4i)
%     type = 'rhb':  R(x,y) = (\Delta_{\Gamma} + zk^2) K(x,y)
%
%  The remainder kernels depend on the mean curvature at the target, which
%  is read from targinfo.kappa or targinfo.mean_curv (defaults to
%  S.mean_curv)
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * type: kernel type, one of 'klb', 'rlb', 'khb', 'rhb'
%    * zk: Helmholtz wavenumber, ignored for 'klb' and 'rlb'
%    * eps: precision requested
%    * targinfo: target info (optional)
%       targinfo.r = (3,nt) target locations
%       targinfo.du = u tangential derivative info
%       targinfo.dv = v tangential derivative info
%       targinfo.n = normal info
%       targinfo.kappa or targinfo.mean_curv (nt,) mean curvature at
%          target
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
%

    switch lower(type)
      case {'klb'}
        islap = true;  iktype = 2; kernel_order = -1;
      case {'rlb'}
        islap = true;  iktype = 1; kernel_order = -1;
      case {'khb'}
        islap = false; iktype = 2; kernel_order = -1;
      case {'rhb'}
        islap = false; iktype = 1; kernel_order = -1;
      otherwise
        error('BELPDE.GET_QUADRATURE_CORRECTION: unknown type ''%s''.', type);
    end

    if islap
      zk = 0;
    end

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;
    npp1 = npatches+1;
    n3 = 3;

    if nargin < 5 || isempty(targinfo)
      targinfo = [];
      targinfo.r = S.r;
      targinfo.du = S.du;
      targinfo.dv = S.dv;
      targinfo.n = S.n;
      targinfo.patch_id = S.patch_id;
      targinfo.uvs_targ = S.uvs_targ;
      targinfo.kappa = S.mean_curv;
    end

    if nargin < 6
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
          'patch_id', targinfo.patch_id, 'uvs_targ', targinfo.uvs_targ, ...
          'kappa', targinfo.mean_curv);
    end

    if isfield(targinfo,'kappa')
      kappa = targinfo.kappa;
    elseif isfield(targinfo,'mean_curv')
      kappa = targinfo.mean_curv;
      targinfo = rmfield(targinfo,'mean_curv');
    else
      error(['BELPDE.GET_QUADRATURE_CORRECTION: mean curvature at the ' ...
             'targets must be supplied in targinfo.kappa or ' ...
             'targinfo.mean_curv.']);
    end
    targinfo.kappa = reshape(kappa, 1, []);

    targs = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;

    if(isfield(targinfo,'patch_id') || isprop(targinfo,'patch_id'))
      patch_id = targinfo.patch_id;
    else
      patch_id = -1*ones(ntarg,1);
    end

    if(isfield(targinfo,'uvs_targ') || isprop(targinfo,'uvs_targ'))
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

    rsc = getnear(S, targinfo);
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

    if islap
        wnear = zeros(nquad,1);
        mex_id_ = 'getnearquad_lap_bel_log(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, iktype, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
    else
        zpars = complex(zk);
        wnear = complex(zeros(nquad,1));
        mex_id_ = 'getnearquad_helm_bel_hank(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[x])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, zpars, iktype, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
    end

    Q = [];
    Q.targinfo = targinfo;
    Q.ifcomplex = double(~islap);
    Q.wavenumber = zk;
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
        spmat = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear);
        Q.spmat = spmat;
    end

end
%
%
%
%
%-------------------------------------------------
