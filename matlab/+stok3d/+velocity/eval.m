function p = eval(S,sigma,targinfo,eps,dpars,varargin)
%
%  stok3d.eval
%    Evaluates the stokes layer potential at a collection 
%    of targets
%
%  Syntax
%   pot = stok3d.velocity.eval(S,sigma,targinfo,eps,dpars)
%   pot = stok3d.velocity.eval(S,sigma,targinfo,eps,dpars,opts)
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
%    * targinfo: target info 
%       targinfo.r = (3,nt) target locations
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%          patch if on-surface (optional)
%    * eps: precision requested
%    * dpars: kernel parameters
%        dpars(1) - single layer strength
%        dpars(2) - double layer strength
%    * opts: options struct
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.precomp_quadrature: precomputed quadrature corrections struct 
%           currently only supports quadrature corrections
%           computed in rsc format 
%    

    if(nargin < 6) 
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
    mex_id_ = 'stok_comb_vel_eval_addsub(c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[xx], c i double[xx], c i int[x], c i int[x], c i double[xx], c i double[x], c i int[x], c i double[x], c i int[x], c i dcomplex[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[xx], c i int[x], c i int[x], c i int[x], c i double[xx], c i double[x], c i int[x], c i double[x], c i int[x], c i int[x], c i double[xx], c i int[x], c i int[x], c io double[xx])';
[p] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, ixyzso, srcover, wover, lwork, work, idensflag, ndim_s, sigma, ipotflag, ndim_p, p, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, 1, 1, ndd, 1, ndz, 1, ndi, 1, ntargp1, nnz, nnzp1, 1, 1, nker, nquad, npatches, 1, npatp1, 12, nptso, nptso, 1, lwork, 1, 1, ndim_s, npts, 1, 1, ndim_p, ntarg);
end    
%
%
%
%----------------------------------
%
%
