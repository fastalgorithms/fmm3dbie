function p = eval(S,zpars,sigma,eps,varargin)
%
%  helm3d.dirichlet.eval
%    Evaluates the helmholtz dirichlet layer potential at a collection 
%    of targets
%
%  Syntax
%   pot = helm3d.dirichlet.eval(S,zpars,sigma,eps)
%   pot = helm3d.dirichlet.eval(S,zpars,sigma,eps,targinfo)
%   pot = helm3d.dirichlet.eval(S,zpars,sigma,eps,targinfo,opts)
%
%  Integral representation
%     pot = \alpha S_{k} [\sigma] + \beta D_{k} [\sigma]
%
%  S_{k}, D_{k}: helmholtz single and double layer potential
%  
%  k, \alpha, beta = zpars(1:3)
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * zpars: kernel parameters
%        zpars(1) - wave number
%        zpars(2) - single layer strength
%        zpars(3) - double layer strength
%    * sigma: layer potential density
%    * eps: precision requested
%    * targinfo: target info (optional)
%       targinfo.r = (3,nt) target locations
%       targinfo.du = u tangential derivative info
%       targinfo.dv = v tangential derivative info
%       targinfo.n = normal info
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%          patch if on-surface (optional)
%    * opts: options struct
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%    
%


%
%
% Todo: Fix varargin
%
    if(nargin < 6) 
      opts = [];
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
      opts = []
    end

    ff = 'rsc';

    [targs] = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;

% Compute quadrature corrections    
    if(~nonsmoothonly) 
      opts_quad = [];
      opts_quad.format = 'rsc';
%
%  For now Q is going to be a struct with 'quad_format', 
%  'nkernels', 'pde', 'bc', 'kernel', 'ker_order',
%  and either, 'wnear', 'row_ind', 'col_ptr', or
%  with 'spmat' as a sparse matrix or a cell array of wnear/spmat
%  if nkernel is >1
%

      [Q] = helm3d.dirichlet.get_quadrature_correction(S,zpars,eps,targinfo,opts_quad);
      nnz = length(Q.col_ind);
      nquad = Q.iquad(end)-1;
    else
      opts_qcorr = [];
      opts_qcorr.type = 'complex';
      Q = init_empty_quadrature_correction(targinfo,opts_qcorr);
    end
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


    p = complex(zeros(ntarg,1));

% Call the layer potential evaluator
    mex_id_ = 'lpcomp_helm_comb_dir_addsub(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i double[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i dcomplex[x], i dcomplex[x], i int[x], i int[x], i int[x], i double[xx], i double[x], io dcomplex[x])';
[p] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, zpars, nnz, row_ptr, col_ind, iquad, nquad, wnear, sigma, novers, nptso, ixyzso, srcover, wover, p, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, 1, 3, 1, ntargp1, nnz, nnzp1, 1, nquad, npts, npatches, 1, npatp1, 12, nptso, nptso, ntarg);
end    
%
%
%
%----------------------------------
%
