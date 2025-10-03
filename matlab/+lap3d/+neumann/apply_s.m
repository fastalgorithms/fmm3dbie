function p = apply_s(S,sigma,eps,opts)
%
%  lap3d.neumann.apply_s
%    Evaluates the Laplace neumann layer potential at a collection 
%    of targets
%
%  Syntax

%
%  Integral representation
%     pot = S_{0} [\sigma] 
%
%  S_{0}: Laplace single layer potential
%  
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * sigma: layer potential density
%    * opts: options struct
%    * eps: precision requested
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.precomp_quadrature: precomputed quadrature corrections 
%           struct currently only supports quadrature corrections
%           computed in rsc format 
%    

    if(nargin < 4) 
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
        opts_quad.rep = 'eval';
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


    p = zeros(npts,1);

    ndd = 2;
    dpars = [1,0];
    ndz = 0;
    zpars = [];
    ndi = 0;
    ipars = [];
    nker = 1;
    lwork = 0;
    work = [];
    ndim_s = 1;
    ndim_p = 1;
    idensflag = 1;
    ipotflag = 1;


mex_id_ = 'lpcomp_lap_comb_dir_addsub(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[x], i int[x], i double[x], i int[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[x], i int[x], i int[x], i int[x], i double[xx], i double[x], i int[x], i double[x], i int[x], i double[x], i double[x])';
fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, ixyzso, srcover, wover, lwork, work, ndim_s, sigma, p, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndd, 1, ndz, 1, ndi, 1, nptsp1, nnz, nnzp1, 1, 1, nquad, npatches, 1, npatp1, 12, nptso, nptso, 1, lwork, 1, npts, npts);

end
%
%
%
%----------------------------------
%
