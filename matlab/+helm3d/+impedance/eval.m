function p = eval(S,densities,eps,zk,alpha,varargin)
%
%  helm3d.impedance.eval
%    Evaluates the helmholtz impedance layer potential at a collection 
%    of targets
%
%  Syntax
%   pot = helm3d.impedance.eval(S,densities,eps,zk,alpha)
%   pot = helm3d.impedance.eval(S,densities,eps,zk,alpha,targinfo)
%   pot = helm3d.impedance.eval(S,densities,eps,zk,alpha,targinfo,Q)
%   pot = helm3d.impedance.eval(S,densities,eps,zk,alpha,targinfo,Q,opts)
%
%  Integral representation
%     pot = S_{k} [\sigma] + i \alpha D_{k} S_{i|k|}[\sigma]
%
%  S_{k}, D_{k}: helmholtz single and double layer potential
%  
%  densities(2,1:npts) for rpcomb, with densities(1,:) = sigma, 
%  and densities(2,:) = S_{i|k|} sigma
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * densities: layer potential densities
%    * eps: precision requested
%    * zk: wave number
%    * alpha: alpha above
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

    if(nargin < 8) 
      opts = [];
    else
      opts = varargin{3};
    end

    isprecompq = true;
    if(nargin < 7)
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
        opts_qcorr.type = 'complex';
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

    if(nargin < 6)
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
        opts_quad.rep = 'rpcomb-eval';
%
%  For now Q is going to be a struct with 'quad_format', 
%  'nkernels', 'pde', 'bc', 'kernel', 'ker_order',
%  and either, 'wnear', 'row_ind', 'col_ptr', or
%  with 'spmat' as a sparse matrix or a cell array of wnear/spmat
%  if nkernel is >1
%

        [Q] = helm3d.impedance.get_quadrature_correction(S,eps,zk,alpha,targinfo,opts_quad);
      else
        opts_qcorr = [];
        opts_qcorr.type = 'complex';
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

    p = complex(zeros(ntarg,1));

    ndd = 0;
    dpars = [];
    ndz = 2;
    zpars = complex(zeros(2,1));
    zpars(1) = zk;
    zpars(2) = alpha;
    ndi = 0;
    ipars = [];
    nker = 2;
    lwork = 0;
    work = [];
    ndim_s = 2;
    ndim_p = 1;
    idensflag = 1;
    ipotflag = 1;
% Call the layer potential evaluator
    mex_id_ = 'helm_rpcomb_eval_addsub(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i double[x], i int[x], i double[x], i int[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i dcomplex[xx], i int[x], i int[x], i int[x], i double[xx], i double[x], i int[x], i double[x], i int[x], i int[x], i dcomplex[xx], i int[x], i int[x], io dcomplex[x])';
[p] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, ixyzso, srcover, wover, lwork, work, idensflag, ndim_s, densities, ipotflag, ndim_p, p, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, 1, 1, ndd, 1, ndz, 1, ndi, 1, ntargp1, nnz, nnzp1, 1, 1, nker, nquad, npatches, 1, npatp1, 12, nptso, nptso, 1, lwork, 1, 1, 2, npts, 1, 1, ntarg);
end    
%
%
%
%----------------------------------
%
