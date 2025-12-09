function [densities, varargout] = solver(S, rhs, eps, zks, rep_params, opts)
%
%  helm3d.transmission.solver
%    Solve the helmholtz transmission boundary value problem
%      (\Delta + k0^2)u0 = 0   in the exterior
%      (\Delta + k1^2)u1 = 0   in the interior
%
%      alpha0 u0 - alpha1 u1 = f 
%      beta0 dudn0 - beta1 dudn1 = g
%
%  Syntax
%   [densities] = helm3d.transmission.solver(S,rhs,eps,zks,rep_params)
%   [densities] = helm3d.transmission.solver(S,rhs,eps,zks,rep_params,opts)
%
%  Integral representation
%     pot0 = 1/beta0 (ik_{0} S_{k0} [\sigma] + D_{k0} [\rho]) (exterior representation)
%     pot1 = 1/beta1 (ik_{1} S_{k1} [\sigma] + D_{k1} [\rho]) (interior representation)
%
%  S_{k}, D_{k}: helmholtz single and double layer potential
%  
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * rhs(2,S.npts): boundary data, 
%        the first component should be data for 
%          alpha0 u0 - alpha1 u1 (f);
%        the second component should be data for
%          beta0 dudn0 - beta1 dudn1 (g);
%    * eps: precision requested
%    * zks: wave numbers, must be of size (2,1) 
%    * rep_params: [alpha0; beta0; alpha1; beta1] above 
%    * opts: options struct
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.eps_gmres - tolerance to which linear system is to be
%           solved (eps_gmres = eps)
%        opts.maxit - maximum number of gmres iterations (200)
%        opts.quadrature_correction - precomputed quadrature correction ([])
%        
%
%  Output arguemnts:
%    * densities: solution to integral equation
%         densities(1,:) = \rho 
%         densities(2,:) = \sigma
%    
%
    
    if(nargin < 6) 
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
    mex_id_ = 'get_patch_id_uvs(c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c io int[x], c io double[xx])';
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
        opts_quad.rep = 'trans-bc';
%
%  For now Q is going to be a struct with 'quad_format', 
%  'nkernels', 'pde', 'bc', 'kernel', 'ker_order',
%  and either, 'wnear', 'row_ind', 'col_ptr', or
%  with 'spmat' as a sparse matrix or a cell array of wnear/spmat
%  if nkernel is >1
%

        [Q] = helm3d.transmission.get_quadrature_correction(S,eps,zks,rep_params,targinfo,opts_quad);
      end
    else
      opts_qcorr = [];
      opts_qcorr.type = 'complex';
      opts_qcorr.nker = 4;
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


    densities = complex(zeros(2,npts));

    rhs_use = rhs; 
    z = -(1j*zks(1) + 1j*zks(2))/(rep_params(1)/rep_params(2) + rep_params(3)/rep_params(4));
    rhs_use(1,:) = z*rhs(1,:);

    niter = 0;
    errs = zeros(maxit+1,1);
    maxitp1 = maxit + 1;
    rres = 0;
    nker = 4;

    zpars = complex(zeros(6,1));
    zpars(1) = zks(1);
    zpars(2) = rep_params(1);
    zpars(3) = rep_params(2);
    zpars(4) = zks(2);
    zpars(5) = rep_params(3);
    zpars(6) = rep_params(4);


% Call the layer potential evaluator
    mex_id_ = 'helm_comb_trans_solver_guru(c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int[x], c i dcomplex[xx], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i dcomplex[xx], c i int[x], c i int[x], c i int[x], c i double[xx], c i double[x], c i double[x], c io int[x], c io double[x], c io double[x], c io dcomplex[xx])';
[niter, errs, rres, densities] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, maxit, rhs_use, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, densities, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 6, 1, 2, npts, 1, nptsp1, nnz, nnzp1, 1, 1, nker, nquad, npatches, 1, npatp1, 12, nptso, nptso, 1, 1, maxitp1, 1, 2, npts);
    
    errs = errs(1:niter);
    varargout{1} = errs;
    varargout{2} = rres;
    varargout{3} = Q;
end    
%
%

%-------------------------------------------------
%
%-------------------------------------------------
%
%%
%%   Stokes routines
%
%
%-------------------------------------------------

