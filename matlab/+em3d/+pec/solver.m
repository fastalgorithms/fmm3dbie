function [densities, varargout] = solver(S, einc, hinc, eps, zk, rep_params, opts)
%
%  em3d.pec.solver
%    Solve the Maxwell pec boundary value problem
%
%  Syntax
%   [densities] = em3d.pec.solver(S, einc, hinc, eps, zk, rep_params)
%   [densities] = em3d.pec.solver(S, einc, hinc, eps, zk, rep_params, opts)
%
%  This routine will support the following representations:
%  * nrccie   (Non-resonant charge current integral equation)
%  * dpie     (Decoupled potential integral equation)
%  * mfie     (Magnetic field integral equation)
%  * aumfie   (Augmented magnetic field integral equation)
%  * aurcsie  (Augmented regularized combined source integral equation)
%  * gendeb   (Generalized Debye)
%
%  For notes on the specific representations, boundary integral equations,
%  and order of kernels returned by this routine, checkout
%  em3d.pec.Contents.m
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * einc, hinc: incident electric and magnetic fields 
%    * eps: precision requested
%    * zk: wave number
%    * rep_params: parameters for integral representation 
%                  for nrccie, it should be a scalar
%    * opts: options struct
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.eps_gmres - tolerance to which linear system is to be
%           solved (eps_gmres = eps)
%        opts.maxit - maximum number of gmres iterations (200)
%        opts.quadrature_correction - precomputed quadrature correction ([])
%        opts.rep - integral representation being used
%                         Supported representations
%
%  Output arguemnts:
%    * densities: layer potential density
%    
%
    
    if(nargin < 7) 
      opts = [];
    end

    nonsmoothonly = false;
    if(isfield(opts,'nonsmoothonly'))
      nonsmoothonly = opts.nonsmoothonly;
    end

    rep = 'nrccie';
    if isfield(opts, 'rep')
      rep = opts.rep;
    end

    if strcmpi(rep, 'nrccie')
      nker = 9;
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
        opts_quad.rep = [rep '-bc'];
%
%  For now Q is going to be a struct with 'quad_format', 
%  'nkernels', 'pde', 'bc', 'kernel', 'ker_order',
%  and either, 'wnear', 'row_ind', 'col_ptr', or
%  with 'spmat' as a sparse matrix or a cell array of wnear/spmat
%  if nkernel is >1
%

        [Q] = em3d.pec.get_quadrature_correction(S,eps,zk,rep_params,targinfo,opts_quad);
      end
    else
      opts_qcorr = [];
      opts_qcorr.type = 'complex';
      opts_qcorr.nker = nker;
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


    niter = 0;
    errs = zeros(maxit+1,1);
    maxitp1 = maxit + 1;
    rres = 0;

    if strcmpi(rep, 'nrccie')
      zjvec = complex(zeros(3,npts));
      rho = complex(zeros(1,npts));
      zpars = complex(zeros(2,1));
      zpars(1) = zk;
      zpars(2) = rep_params;
% Call the layer potential evaluator
      mex_id_ = 'em_nrccie_pec_solver_guru(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[x], i dcomplex[x], i int[x], i dcomplex[xx], i dcomplex[xx], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i dcomplex[xx], i int[x], i int[x], i int[x], i double[xx], i double[x], i double[x], io int[x], io double[x], io double[x], io dcomplex[xx], io dcomplex[x])';
[niter, errs, rres, zjvec, rho] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, maxit, einc, hinc, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, zjvec, rho, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 2, 1, 3, npts, 3, npts, 1, nptsp1, nnz, nnzp1, 1, 1, nker, nquad, npatches, 1, npatp1, 12, nptso, nptso, 1, 1, maxitp1, 1, 3, npts, npts);
      
      densities = complex(zeros(4,npts));
      densities(1:3,:) = zjvec;
      densities(4,:) = rho;
    end

    errs = errs(1:niter);
    varargout{1} = errs;
    varargout{2} = rres;
    varargout{3} = Q;
end    
%
%
%-------------
%   POLYNYA routs
%----------------------------------
%
