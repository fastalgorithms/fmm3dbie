function [sigma,varargout] = solver(S, rhs, eps, dpars, opts)
%
%  stok3d.solver
%    Solve the Stokes boundary value problem
%
%  Syntax
%   sigma = stok3d.velocity.solver(S,rhs,eps,dpars)
%   sigma = stok3d.velocity.solver(S,rhs,eps,dpars,opts)
%
%  Integral representation
%     pot = \alpha S_{stok} [\sigma] + \beta D_{stok} [\sigma]
%
%  S_{stok}, D_{stok}: Stokes single and double layer potential
%  
%  \alpha, beta = dpars(1:2)
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * rhs: (3, npts) boundary data 
%    * eps: precision requested
%    * dpars: kernel parameters
%        dpars(1) - single layer strength
%        dpars(2) - double layer strength
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
%    * sigma: (3,npts) layer potential density
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

      if isfield(opts, 'quadrature_correction')
         Q = opts.quadrature_correction;
      else
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
      end
    else
      opts_qcorr = [];
      opts_qcorr.type = 'double';
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


    sigma = zeros(3,npts);
    niter = 0;
    errs = zeros(maxit+1,1);
    maxitp1 = maxit + 1;
    rres = 0;
    nker = 6;
% Call the layer potential evaluator
    mex_id_ = 'stok_comb_vel_solver_guru(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[x], i double[x], i int[x], i int[x], i double[xx], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i int[x], i int[x], i int[x], i double[xx], i double[x], i double[x], io int[x], io double[x], io double[x], io double[xx])';
[niter, errs, rres, sigma] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, dpars, maxit, ifout, rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, sigma, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 2, 1, 1, 3, npts, 1, nptsp1, nnz, nnzp1, 1, 1, nker, nquad, npatches, 1, npatp1, 12, nptso, nptso, 1, 1, maxitp1, 1, 3, npts);

    errs = errs(1:niter);
    varargout{1} = errs;
    varargout{2} = rres;
    varargout{3} = Q;
end  
%
%
%
%
%----------------------------------------
%%
%%   Maxwell pec routines
%
%
%-------------------------------------------------

