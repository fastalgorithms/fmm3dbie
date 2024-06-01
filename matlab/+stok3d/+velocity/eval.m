function p = eval(S,sigma,eps,dpars,varargin)
%
%  stok3d.eval
%    Evaluates the stokes layer potential at a collection 
%    of targets
%
%  Syntax
%   pot = stok3d.velocity.eval(S,sigma,eps,dpars)
%   pot = stok3d.velocity.eval(S,sigma,eps,dpars,targinfo)
%   pot = stok3d.velocity.eval(S,sigma,eps,dpars,targinfo,Q)
%   pot = stok3d.velocity.eval(S,sigma,eps,dpars,targinfo,Q,opts)
%
%  Integral representation
%     pot = \alpha S_{stok} [\sigma] + \beta D_{stok} [\sigma]
%
%  S_{stok}, D_{stok}: stokes single and double layer potential
%  
%  \alpha, beta = dpars(1:2)
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * sigma: (3, ns) layer potential density
%    * eps: precision requested
%    * dpars: kernel parameters
%        dpars(1) - single layer strength
%        dpars(2) - double layer strength
%    * targinfo: target info (optional)
%       targinfo.r = (3,nt) target locations
%       targinfo.du = u tangential derivative info
%       targinfo.dv = v tangential derivative info
%       targinfo.n = normal info
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%          patch if on-surface (optional)
%    * Q: precomputed quadrature corrections struct (optional)
%           currently only supports quadrature corrections
%           computed in rsc format 
%    * opts: options struct
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%    

    if(nargin < 7) 
      opts = [];
    else
      opts = varargin{3};
    end

    isprecompq = true;
    if(nargin < 6)
       Q = [];
       isprecompq = false;
    else
       Q = varargin{2}; 
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

    nonsmoothonly = false;
    if(isfield(opts,'nonsmoothonly'))
      nonsmoothonly = opts.nonsmoothonly;
    end

% Extract arrays
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;

    if(nargin < 5)
      targinfo = [];
      targinfo.r = S.r;
      targinfo.du = S.du;
      targinfo.dv = S.dv;
      targinfo.n = S.n;
      patch_id  = zeros(npts,1);
      uvs_targ = zeros(2,npts);
      mex_id_ = 'get_patch_id_uvs(i int[x], i int[x], i int[x], i int[x], i int[x], io int[x], io double[xx])';
[patch_id, uvs_targ] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, patch_id, uvs_targ, 1, npatches, npatp1, npatches, 1, npts, 2, npts);
      targinfo.patch_id = patch_id;
      targinfo.uvs_targ = uvs_targ;
      opts = [];
    else
      targinfo = varargin{1};
    end

    ff = 'rsc';

    [targs] = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;

% Compute quadrature corrections   
    if ~isprecompq
      if ~nonsmoothonly
        opts_quad = [];
        opts_quad.format = 'rsc';
%
%  For now Q is going to be a struct with 'quad_format', 
%  'nkernels', 'pde', 'bc', 'kernel', 'ker_order',
%  and either, 'wnear', 'row_ind', 'col_ptr', or
%  with 'spmat' as a sparse matrix or a cell array of wnear/spmat
%  if nkernel is >1
%

        [Q] = stok3d.velocity.get_quadrature_correction(S,eps,dpars,targinfo,opts_quad);
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

    ndd = 2;
    ndz = 0;
    zpars = [];
    ndi = 0;
    ipars = [];
    nker = 6;
    lwork = 0;
    work = [];
    ndim_s = 3;
    idensflag = 0;
    ipotflag = 0;
    ndim_p = 3;
% Call the layer potential evaluator
    mex_id_ = 'stok_comb_vel_eval_addsub(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i double[x], i int[x], i double[x], i int[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i int[x], i int[x], i int[x], i double[xx], i double[x], i int[x], i double[x], i int[x], i int[x], i double[xx], i int[x], i int[x], io double[xx])';
[p] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, ixyzso, srcover, wover, lwork, work, idensflag, ndim_s, sigma, ipotflag, ndim_p, p, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, 1, 1, ndd, 1, ndz, 1, ndi, 1, ntargp1, nnz, nnzp1, 1, 1, nker, nquad, npatches, 1, npatp1, 12, nptso, nptso, 1, lwork, 1, 1, ndim_s, npts, 1, 1, ndim_p, ntarg);
end    
%
%
%
%----------------------------------
%
