function Q = get_quadrature_correction(S, type, eps, targinfo, opts)
%
%  lap2d.get_quadrature_correction
%    This subroutine returns the near quadrature correction for the
%    2D Laplace volume potential, with density supported on the flat
%    surface S (in the z = 0 plane), and targets given by targinfo, as
%    a sparse matrix/rsc format
%
%  Syntax
%   Q = lap2d.get_quadrature_correction(S,type,eps)
%   Q = lap2d.get_quadrature_correction(S,type,eps,targinfo)
%   Q = lap2d.get_quadrature_correction(S,type,eps,targinfo,opts)
%
%  Kernels
%     type = 's':   G(x,y) = -log|x-y|/(2 pi)
%     type = 'sp':  \nabla_{n_x} G(x,y)
%     type = 'sgrad': \nabla_x G(x,y), wnear is (2,nquad) with the
%                     d/dx and d/dy corrections
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * type: kernel type, one of 's', 'sp', 'sgrad'
%    * eps: precision requested
%    * targinfo: target info (optional, defaults to S)
%       targinfo.r = (3,nt) target locations, only r(1:2,:) is used
%       targinfo.n = (3,nt) target normals, required for 'sp'
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

    switch lower(type)
      case {'s', 'single'}
        issp = false; isgrad = false; kernel_order = -1;
      case {'sp', 'sprime'}
        issp = true;  isgrad = false; kernel_order = -1;
      case {'sg', 'sgrad'}
        issp = true;  isgrad = true;  kernel_order = -1;
      otherwise
        error('LAP2D.GET_QUADRATURE_CORRECTION: unknown type ''%s''.', type);
    end

    [srcvals,srccoefs,norders,ixyzs,iptype,~] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npp1 = npatches+1;

    if nargin < 4 || isempty(targinfo)
      targinfo = S;
    end

    if nargin < 5
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

    zpars = complex(0);
    % the gradient is sprime with the target normal set to e1, then e2
    ndir = 1;
    if isgrad
        ndir = 2;
        if size(targs,1) < 12, targs(12,:) = 0; end
        ndtarg = size(targs,1);
    end
    wall = zeros(ndir,nquad);
    for idir = 1:ndir
    if isgrad
        targs(10:12,:) = 0;
        targs(9+idir,:) = 1;
    end
    wnear = zeros(nquad,1);
    if issp
        mex_id_ = 'getnearquad_lap2d_gdn(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
    else
        mex_id_ = 'getnearquad_lap2d_dir(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
    end
    wall(idir,:) = wnear;
    end
    if isgrad, wnear = wall; end

    Q = [];
    Q.targinfo = targinfo;
    Q.ifcomplex = 0;
    Q.wavenumber = 0;
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
        if isgrad
          spmat = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear, ...
              kernel3d.rsc_interleave_full(2, 1));
        else
          spmat = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear);
        end
        Q.spmat = spmat;
    end

end
%
%
%
%
%-------------------------------------------------

