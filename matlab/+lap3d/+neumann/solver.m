function [sigma, varargout] = solver(S, rhs, eps, opts)
%
%  lap3d.neumann.solver
%    Solve the Laplace neumann boundary value problem
%
%  Syntax
%   [sigma] = lap3d.neumann.solver(S,rhs,eps)
%   [sigma] = lap3d.neumann.solver(S,rhs,eps,opts)
%
%  Integral representation
%     pot = S_{0} [\sigma] 
%
%  S_{0}: Laplace single layer potential
%  
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * rhs: boundary data 
%    * eps: precision requested
%    * opts: options struct
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.eps_gmres - tolerance to which linear system is to be
%           solved (eps_gmres = eps)
%        opts.maxit - maximum number of gmres iterations (200)
%        opts.ifout - whether to solve interior problem or not (1)
%        opts.quadrature_correction - precomputed quadrature correction ([])
%        
%
%  Output arguemnts:
%    * sigma: layer potential density
%    
%
    
    if(nargin < 4) 
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
    mex_id_ = 'get_patch_id_uvs(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io int64_t[x], c io double[xx])';
[patch_id, uvs_targ] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, patch_id, uvs_targ, 1, npatches, npatp1, npatches, 1, npts, 2, npts);
    targinfo.patch_id = patch_id;
    targinfo.uvs_targ = uvs_targ;

    ff = 'rsc';

    nptsp1 = npts+1;

% Compute quadrature corrections    
    if(~nonsmoothonly)

      if isfield(opts, 'quadrature_correction')
         Q = opts.quadrature_correction;
      else
        opts_quad = [];
        opts_quad.format = 'rsc';
        opts_quad.rep = 'bc';
%
%  For now Q is going to be a struct with 'quad_format', 
%  'nkernels', 'pde', 'bc', 'kernel', 'ker_order',
%  and either, 'wnear', 'row_ind', 'col_ptr', or
%  with 'spmat' as a sparse matrix or a cell array of wnear/spmat
%  if nkernel is >1
%

        [Q] = lap3d.neumann.get_quadrature_correction(S,eps,targinfo,opts_quad);
      end
    else
      opts_qcorr = [];
      opts_qcorr.type = 'double';
      opts_qcorr.nker = 1;
      Q = init_empty_quadrature_correction(targinfo,opts_qcorr);
    end
    nnz = length(Q.col_ind);
    nquad = Q.iquad(end)-1;
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


    sigma = zeros(npts,1);
    niter = 0;
    errs = zeros(maxit+1,1);
    maxitp1 = maxit + 1;
    rres = 0;
    nker = 1;


% Call the layer potential evaluator
    mex_id_ = 'lap_s_neu_solver_guru(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c io int64_t[x], c io double[x], c io double[x], c io double[x])';
[niter, errs, rres, sigma] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, maxit, ifout, rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, sigma, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, 1, npts, 1, nptsp1, nnz, nnzp1, 1, 1, nquad, npatches, 1, npatp1, 12, nptso, nptso, 1, 1, maxitp1, 1, npts);

    errs = errs(1:niter);
    varargout{1} = errs;
    varargout{2} = rres;
    varargout{3} = Q;
end    
%
%

%-------------------------------------------------
%
%%
%%   Helmholtz dirichlet routines
%
%
%-------------------------------------------------

