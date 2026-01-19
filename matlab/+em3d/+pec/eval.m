function [E, H] = eval(S, densities, targinfo, eps, zk, rep_params, varargin)
%
%  em3d.pec.eval
%
%    This subroutine evaluates the electric and magnetic
%    field at a colelction of targets given the solution
%    to the corresponding integral equation
%
%
%  Notes for this routine:
%  The PDE takes the form
%  1v) \nabla \times E =  ik H
%  2v) \nabla \cdot  E =     0
%  3v) \nabla \times H = -ik E
%  4v) \nabla \cdot  H =     0
%  
%  where E is the electric field, H is the magnetic field, 
%  and k is the wavenumber
%
%  The PEC boundary conditions are given by
%  1b) n \times (E + E_in) = 0
%  2b) n \cdot  (E + E_in) = \rho
%  3b) n \times (H + H_in) = J
%  4b) n \cdot  (H + H_in) = 0
%
%  where (E_in, H_in) are the incoming electric and magnetic 
%  fields, \rho is the surface current, and J is the 
%  surface current.
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
%  Syntax
%   [E, H] = em3d.pec.eval(S, densities, targinfo, eps, zk, rep_params)
%   [E, H] = em3d.pec.eval(S, densities, targinfo, eps, zk, rep_params, opts)
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * densities: layer potential densities, of size (ndim, npts)
%        where ndim depends on the integral representation used
%    * targinfo: target info 
%       targinfo.r = (3,nt) target locations
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%          patch if on-surface (optional)
%    * eps: precision requested
%    * zk : wave number
%    * rep_params: parameters for integral representation 
%                  for nrccie, it should be a scalar
%    * opts: options struct
%        opts.rep - integral representation being used
%                         Supported representations
%                         'nrccie'
%        opts.nonsmoothonly - use smooth quadrature rule for
%                             evaluating layer potential (false)
%        opts.precomp_quadrature: precomputed quadrature corrections struct 
%           currently only supports quadrature corrections
%           computed in rsc format 
%
%
%
%

    if(nargin < 7) 
      opts = [];
    else
      opts = varargin{1};
    end

    nonsmoothonly = false;
    if(isfield(opts,'nonsmoothonly'))
      nonsmoothonly = opts.nonsmoothonly;
    end

    isprecompq = false;
    if isfield(opts, 'precomp_quadrature')
      Q = opts.precomp_quadrature;
      isprecompq = true;
    end

    if(isprecompq)
      if ~(strcmpi(Q.format,'rsc'))
        fprintf('Invalid precomputed quadrature format\n');
        fprintf('Ignoring quadrature corrections\n');
        opts_qcorr = [];
        opts_qcorr.type = 'complex';

        opts_qcorr.nker = nker;
        Q = init_empty_quadrature_correction(targinfo,opts_qcorr);
      end
    end

    
    rep = 'nrccie';

    if isfield(opts, 'rep')
      rep = opts.rep;
    end

    if strcmpi(rep, 'nrccie')
      nker = 4;
      ndim_s = 4;
      [nn, ~] = size(densities);
      zpars = complex(zeros(2,1));
      zpars(1) = zk;
      zpars(2) = rep_params;
      
      if nn ~= ndim_s
        error('EM3D.PEC.EVAL: number of densities not consistent with representation\n'); 
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
        opts_quad.rep = [rep '-eval'];
%
%  For now Q is going to be a struct with 'quad_format', 
%  'nkernels', 'pde', 'bc', 'kernel', 'ker_order',
%  and either, 'wnear', 'row_ind', 'col_ptr', or
%  with 'spmat' as a sparse matrix or a cell array of wnear/spmat
%  if nkernel is >1
%

        [Q] = em3d.pec.get_quadrature_correction(S, eps, zk, rep_params, targinfo, opts_quad);
      else
        opts_qcorr = [];
        opts_qcorr.type = 'complex';
        opts_qcorr.nker = nker;
        Q = init_empty_quadrature_correction(targinfo,opts_qcorr);
      end
    end
    nquad = Q.iquad(end)-1;
    nnz = length(Q.col_ind);
    nnzp1 = nnz+1; 

    [novers] = get_oversampling_parameters(S, Q, eps);
    Sover = oversample(S,novers);


% Extract oversampled arrays

    [srcover,~,~,ixyzso,~,wover] = extract_arrays(Sover);
    nptso = Sover.npts; 

% Extract quadrature arrays
    row_ptr = Q.row_ptr;
    col_ind = Q.col_ind;
    iquad = Q.iquad;
    wnear = Q.wnear;

    p = complex(zeros(6,ntarg));

    ipotflag = 3;


    ndd = 0;
    dpars = [];
    ndz = length(zpars);
    ndi = 0;
    ipars = [];
    lwork = 0;
    work = [];
    ndim_p = 6;
    idensflag = 1;
% Call the layer potential evaluator
    
    if strcmpi(rep, 'nrccie')
      mex_id_ = 'em_nrccie_eval_addsub(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i dcomplex[xx], c i int64_t[x], c i int64_t[x], c io dcomplex[xx])';
[p] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, ixyzso, srcover, wover, lwork, work, idensflag, ndim_s, densities, ipotflag, ndim_p, p, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, 1, 1, ndd, 1, ndz, 1, ndi, 1, ntargp1, nnz, nnzp1, 1, 1, nker, nquad, npatches, 1, npatp1, 12, nptso, nptso, 1, lwork, 1, 1, ndim_s, npts, 1, 1, ndim_p, ntarg);
    end
    E = p(1:3,:);
    H = p(4:6,:);
end    
%
%
%
%----------------------------------
%
