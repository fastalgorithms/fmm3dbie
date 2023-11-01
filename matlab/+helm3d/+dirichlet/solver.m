function [sigma,varargout] = solver(S,zpars,rhs,eps,varargin)
%
%  helm3d.dirichlet.solver
%    Solve the helmholtz dirichlet boundary value problem
%
%  Syntax
%   sigma = helm3d.dirichlet.solver(S,zpars,sigma,eps)
%   sigma = helm3d.dirichlet.solver(S,zpars,sigma,eps,opts)
%
%  Integral representation
%     pot = \alpha S_{k} [\sigma] + \beta D_{k} [\sigma]
%
%  S_{k}, D_{k}: helmholtz single and double layer potential
%  
%  k, \alpha, beta = zpars(1:3)
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * zpars: kernel parameters
%        zpars(1) - wave number
%        zpars(2) - single layer strength
%        zpars(3) - double layer strength
%    * rhs: boundary data 
%    * eps: precision requested
%    * opts: options struct
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.eps_gmres - tolerance to which linear system is to be
%           solved (eps_gmres = eps)
%        opts.maxit - maximum number of gmres iterations (200)
%        opts.ifout - whether to solve interior problem or not (1)
%        
%
%  Output arguemnts:
%    * sigma: layer potential density
%    
%
    
    if(nargin < 5) 
      opts = [];
    end

    nonsmoothonly = false;
    if(isfield(opts,'nonsmoothonly'))
      nonsmoothonly = opts.nonsmoothonly;
    end

    eps_gmres = eps;
    if(isfield(opts,'eps_gmres'))
      eps_gmres = opts.eps_gmres;
    end

    maxit = 200;
    if(isfield(opts,'maxit'))
      maxit = opts.maxit;
    end

    ifout = 1;
    if(isfield(opts,'ifout'))
      ifout = opts.ifout;
    end

% Extract arrays
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;

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

    ff = 'rsc';

    nptsp1 = npts+1;

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


    sigma = complex(zeros(npts,1));
    niter = 0;
    errs = zeros(maxit+1,1);
    maxitp1 = maxit + 1;
    rres = 0;

% Call the layer potential evaluator
    mex_id_ = 'helm_comb_dir_solver_guru(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[x], i dcomplex[x], i int[x], i int[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i dcomplex[x], i int[x], i int[x], i int[x], i double[xx], i double[x], i double[x], io int[x], io double[x], io double[x], io dcomplex[x])';
[niter, errs, rres, sigma] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, maxit, ifout, rhs, nnz, row_ptr, col_ind, iquad, nquad, wnear, novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, sigma, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 3, 1, 1, npts, 1, nptsp1, nnz, nnzp1, 1, nquad, npatches, 1, npatp1, 12, nptso, nptso, 1, 1, maxitp1, 1, npts);

    errs = errs(1:niter);
    varargout{1} = errs;
    varargout{2} = rres;
    varargout{3} = Q;
end    
%
%
%
%
%
%%
%%   Helmholtz Neumann routines
%% 
%
%
%----------------------------------
%
