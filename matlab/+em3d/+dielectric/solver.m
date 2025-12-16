function [densities, varargout] = solver(S, einc, hinc, eps, om, rep_params, opts)
%
%  em3d.dielectric.solver
%    Solve the Maxwell dielectric boundary value problem
%
%  Syntax
%   [densities] = em3d.dielectric.solver(S, einc, hinc, eps, om, rep_params)
%   [densities] = em3d.dielectric.solver(S, einc, hinc, eps, om, rep_params, opts)
%
%  This routine will support the following representations:
%  * muller   
%
%  For notes on the specific representations, boundary integral equations,
%  and order of kernels returned by this routine, checkout
%  em3d.dielectric.Contents.m
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * einc, hinc: incident electric and magnetic fields 
%    * eps: precision requested
%    * om: wave number
%    * rep_params: parameters for integral representation 
%                  for muller, it should be a 4 vector
%                  consisting of [\ep0, \mu0, \ep1, \mu1]
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

    rep = 'muller';
    if isfield(opts, 'rep')
      rep = opts.rep;
    end

    if strcmpi(rep, 'muller')
      nker = 16;
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
    mex_id_ = 'get_patch_id_uvs(i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], io int64_t[x], io double[xx])';
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
         Q = em3d.dielectric.get_quadrature_correction(S,eps,om,rep_params);
      end
    else
      opts_qcorr = [];
      opts_qcorr.type = 'complex';
      opts_qcorr.nker = nker;
      Q = init_empty_quadrature_correction(S,opts_qcorr);
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

    ep0 = rep_params(1);
    mu0 = rep_params(2);
    ep1 = rep_params(3);
    mu1 = rep_params(4);

    rhs = complex(zeros(npts,4));

    rhs(:,1) = -sum(einc.*S.drv, 1).'/(mu0 + mu1);
    rhs(:,2) = sum(einc.*S.dru, 1).'/(mu0 + mu1);
    rhs(:,3) = -sum(hinc.*S.drv, 1).'/(ep0 + ep1);
    rhs(:,4) = sum(hinc.*S.dru, 1).'/(ep0 + ep1);

    densities = complex(zeros(npts,4));
    zpars = complex([om, rep_params(:).']);


% Call the solver 
    mex_id_ = 'em_muller_trans_solver_guru(i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i double[xx], i double[x], i dcomplex[x], i int64_t[x], i dcomplex[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i dcomplex[xx], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i double[x], i double[x], io int64_t[x], io double[x], io double[x], io dcomplex[xx])';
[niter, errs, rres, densities] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, maxit, rhs, nnz, row_ptr, col_ind, iquad, nquad, wnear, novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, densities, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 5, 1, npts, 4, 1, nptsp1, nnz, nnzp1, 1, nquad, nker, npatches, 1, npatp1, 12, nptso, nptso, 1, 1, maxitp1, 1, npts, 4);

    errs = errs(1:niter);
    varargout{1} = errs;
    varargout{2} = rres;
    varargout{3} = Q;
end    
%
%
%
