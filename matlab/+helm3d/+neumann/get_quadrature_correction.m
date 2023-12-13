function Q = get_quadrature_correction(S,zpars,eps,targinfo,opts)
%
%  helm3d.neumann.get_quadrature_correction
%    This subroutine returns the near quadrature correction
%    for the chosen representation to solve the Helmholtz Neumann 
%    problem, with densities supported on the surface, and targets 
%    given by targinfo, as a cell array of sparse matrices or an array of matrices
%    in the rsc format, where each cell of the matrix is the
%    representation of the sparse matrix corresponding to one
%    of the kernels
%
%  The routine supports the following representations:
%  * s       (single layer)
%  * rpcomb  (right preconditioned combined field representation)
%
%  The default quadrature corrections are for targets on surface, so  
%  they assume that if targinfo is provided, then it comes equipped
%  with normal information.
%
%  The representation must be set via setting opts.int_rep. 
%  If unspecified, the default is rpcomb
%
%  To obtain the corresponding post-processing/off-surface 
%  quadrature corrections, append '-eval'
%  to the end of the integral representation. For example, for
%  the off-surface S corrections, use opts.int_rep='s-eval'
%
%  In the input array zpars, the first parameter must always be 
%  the wavenumber k
%
%
%  Syntax
%   Q = helm3d.neumann.get_quadrature_correction(S,zpars,eps)
%   Q = helm3d.neumann.get_quadrature_correction(S,zpars,eps,targinfo)
%   Q = helm3d.neumann.get_quadrature_correction(S,zpars,eps,targinfo,opts)
%
%  Integral representation             
%   's':       pot = S_{k} [\sigma] 
%   'rpcomb':  pot = S_{k} [\sigma] + 1i*\alpha D_{k} S_{i|k|}[\sigma]
%
%   opts.int_rep       kernels returned
%     's'               S'_{k} 
%     's-eval'          S_{k} 
%     'rpcomb'          S_{k}', S_{i|k|}, S_{i|k|}', D_{k}' - D_{i|k|}'
%     'rpcomb-eval'     S_{k}, D_{k}
%
%  S_{k}, D_{k}: helmholtz single and double layer potential
%  
%  k, \alpha = zpars(1:2)
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
%        opts.int_rep - integral representation being used
%                         Supported representations
%                         's', 'rpcomb'
%        opts.format  - Storage format for sparse matrices
%                         'rsc' - row sparse compressed format
%                         'csc' - column sparse compressed format
%                         'sparse' - sparse matrix format
%        opts.quadtype - quadrature type, currently only 'ggq' supported
%
%  TODO: fix behavior for eval
%

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;
    npp1 = npatches+1;
    n3 = 3;


    if(nargin < 4)
      targinfo = [];
      targinfo.r = S.r;
      targinfo.du = S.du;
      targinfo.dv = S.dv;
      targinfo.n = S.n;
      patch_id = S.patch_id;
      uvs_targ = S.uvs_targ;
      opts = [];
    end

    ff = 'rsc';
    if(isfield(opts,'format'))
       ff = opts.format;
    end

    if(~(strcmpi(ff,'rsc') || strcmpi(ff,'csc') || strcmpi(ff,'sparse')))
       fprintf('invalid quadrature format, reverting to rsc format\n');
       ff = 'rsc';
    end


    targs = extract_targ_array(targinfo); 
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;
    
    if(isfield(targinfo,'patch_id') || isprop(targinfo,'patch_id'))
      patch_id = targinfo.patch_id;
    else
      patch_id = zeros(ntarg,1);
    end

    if(isfield(targinfo,'uvs_targ') || isprop(targinfo,'uvs_targ'))
      uvs_targ = targinfo.uvs_targ;
    else
      uvs_targ = zeros(2,ntarg);
    end

    if(length(patch_id)~=ntarg)
      fprintf('Incorrect size of patch id in target info struct. Aborting! \n');
    end

    [n1,n2] = size(uvs_targ);
    if(n1 ~=2 && n2 ~=ntarg)
      fprintf('Incorrect size of uvs_targ array in targinfo struct. Aborting! \n');
    end


    
    iptype0 = iptype(1);
    norder0 = norders(1);
    rfac = 0.0;
    rfac0 = 0.0;
    mex_id_ = 'get_rfacs(i int[x], i int[x], io double[x], io double[x])';
[rfac, rfac0] = fmm3dbie_routs(mex_id_, norder0, iptype0, rfac, rfac0, 1, 1, 1, 1);
    

    cms = zeros(3,npatches);
    rads = zeros(npatches,1);
    mex_id_ = 'get_centroid_rads(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], io double[xx], io double[x])';
[cms, rads] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, cms, rads, 1, npatches, npp1, npatches, 1, n9, npts, n3, npatches, npatches);

    rad_near = rads*rfac;
    nnz = 0;
    mex_id_ = 'findnearmem(i double[xx], i int[x], i double[x], i int[x], i double[xx], i int[x], io int[x])';
[nnz] = fmm3dbie_routs(mex_id_, cms, npatches, rad_near, ndtarg, targs, ntarg, nnz, n3, npatches, 1, npatches, 1, ndtarg, ntarg, 1, 1);

    row_ptr = zeros(ntarg+1,1);
    col_ind = zeros(nnz,1);
    ntp1 = ntarg+1;
    nnzp1 = nnz+1;
    mex_id_ = 'findnear(i double[xx], i int[x], i double[x], i int[x], i double[xx], i int[x], io int[x], io int[x])';
[row_ptr, col_ind] = fmm3dbie_routs(mex_id_, cms, npatches, rad_near, ndtarg, targs, ntarg, row_ptr, col_ind, n3, npatches, 1, npatches, 1, ndtarg, ntarg, 1, ntp1, nnz);

    iquad = zeros(nnz+1,1);
    mex_id_ = 'get_iquad_rsc(i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], io int[x])';
[iquad] = fmm3dbie_routs(mex_id_, npatches, ixyzs, npts, nnz, row_ptr, col_ind, iquad, 1, npp1, 1, 1, ntp1, nnz, nnzp1);

    nquad = iquad(nnz+1)-1;
    iquadtype = 1;
    if(isfield(opts,'quadtype'))
      if(strcmpi(opts.quadtype,'ggq'))
         iquadtype = 1;
      else
        fprintf('Unsupported quadrature type, reverting to ggq\n');
        iquadtype = 1;
      end
    end

    int_rep = 'rpcomb';
    if isfield(opts, 'int_rep')
      int_rep = opts.int_rep;
    end

    if strcmpi(int_rep, 's')
       nker = 1;
       wnear = complex(zeros(nquad,nker));
       mex_id_ = 'getnearquad_helm_s_neu(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[x], i int[x], io dcomplex[xx])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad, nker);
    elseif strcmpi(int_rep', 'rpcomb')
       nker = 4;
       wnear = complex(zeros(nquad,nker));
       mex_id_ = 'getnearquad_helm_rpcomb_neu(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[x], i int[x], io dcomplex[xx])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 2, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad, nker);
    end
    
    Q = [];
    Q.targinfo = targinfo;
    Q.ifcomplex = 1;
    Q.wavenumber = zpars(1);
    Q.kernel_order = -1;
    Q.rfac = rfac;
    Q.nquad = nquad;
    Q.format = ff;
    Q.kernel_order = 0;

    if(strcmpi(ff,'rsc'))
        Q.iquad = iquad;
        Q.wnear = wnear;
        Q.row_ptr = row_ptr;
        Q.col_ind = col_ind;
    elseif(strcmpi(ff,'csc'))
        col_ptr = zeros(npatches+1,1);
        row_ind = zeros(nnz,1);
        iper = zeros(nnz,1);
        npatp1 = npatches+1;
        mex_id_ = 'rsc_to_csc(i int[x], i int[x], i int[x], i int[x], i int[x], io int[x], io int[x], io int[x])';
[col_ptr, row_ind, iper] = fmm3dbie_routs(mex_id_, npatches, ntarg, nnz, row_ptr, col_ind, col_ptr, row_ind, iper, 1, 1, 1, ntp1, nnz, npatp1, nnz, nnz);
        Q.iquad = iquad;
        Q.iper = iper;
        Q.wnear = wnear;
        Q.col_ptr = col_ptr;
        Q.row_ind = row_ind;
    else
        spmat = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear);
        Q.spmat = spmat;
    end
    
end
%
%

