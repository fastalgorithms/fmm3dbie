function Q = get_quadrature_correction(S, eps, zk, alpha, targinfo, opts)
%
%  helm3d.neumann.get_quadrature_correction
%    This subroutine returns the near quadrature correction
%    for the combined field representation with density supported
%    on the surface S, and targets given by targinfo 
%    as a sparse matrix/rsc format 
%
%  Syntax
%   Q = helm3d.neumann.get_quadrature_correction(S,eps,zk,alpha)
%   Q = helm3d.neumann.get_quadrature_correction(S,eps,zk,alpha,targinfo)
%   Q = helm3d.neumann.get_quadrature_correction(S,eps,zk,alpha,targinfo,opts)
%
%  Integral representation
%     pot = S_{k} [\sigma] + 1i \alpha D_{k} S_{i|k|} [\sigma]
%
%  S_{k}, D_{k}: helmholtz single and double layer potential
%  
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
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
%    * opts: options struct
%        opts.format - Storage format for sparse matrices
%           'rsc' - row sparse compressed format
%           'csc' - column sparse compressed format
%           'sparse' - sparse matrix format
%        opts.quadtype - quadrature type, currently only 'ggq' supported
%        opts.rep - whether to compute correction for on surface
%                   integral equation or post-processing
%                   'rpcomb-bc' - on surface integral equation. If option
%                   is 'rpcomb-bc' then targinfo is ignored
%                   'rpcomb-eval' - post-processing
%

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;
    npp1 = npatches+1;
    nsp1 = npts + 1;
    n3 = 3;

    if nargin < 5
      targinfo = [];
      targinfo.r = S.r;
      targinfo.du = S.du;
      targinfo.dv = S.dv;
      targinfo.n = S.n;
      targinfo.patch_id = S.patch_id;
      targinfo.uvs_targ = S.uvs_targ;
      opts = [];
    end

    qtype = 'rpcomb-bc';

    if nargin < 6
      opts = [];
    end

    if isfield(opts, 'rep')
      qtype = opts.rep;
    end

    ff = 'rsc';
    if(isfield(opts,'format'))
       ff = opts.format;
    end

    if(~(strcmpi(ff,'rsc') || strcmpi(ff,'csc') || strcmpi(ff,'sparse')))
       fprintf('invalid quadrature format, reverting to rsc format\n');
       ff = 'rsc';
    end

    tinfouse = [];
    if strcmpi(qtype, 'rpcomb-bc')
      tinfouse.r = S.r;
      tinfouse.du = S.du;
      tinfouse.dv = S.dv;
      tinfouse.n = S.n;
      tinfouse.patch_id = S.patch_id;
      tinfouse.uvs_targ = S.uvs_targ;
    else
      tinfouse = targinfo;
    end




    targs = extract_targ_array(tinfouse); 
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;
    
    if(isfield(tinfouse,'patch_id') || isprop(tinfouse,'patch_id'))
      patch_id = tinfouse.patch_id;
    else
      patch_id = zeros(ntarg,1);
    end

    if(isfield(tinfouse,'uvs_targ') || isprop(tinfouse,'uvs_targ'))
      uvs_targ = tinfouse.uvs_targ;
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
    mex_id_ = 'get_rfacs(c i int[x], c i int[x], c io double[x], c io double[x])';
[rfac, rfac0] = fmm3dbie_routs(mex_id_, norder0, iptype0, rfac, rfac0, 1, 1, 1, 1);
    

    cms = S.cms; 
    rads = S.rads;

    rad_near = rads*rfac;
    nnz = 0;
    mex_id_ = 'findnearmem(c i double[xx], c i int[x], c i double[x], c i int[x], c i double[xx], c i int[x], c io int[x])';
[nnz] = fmm3dbie_routs(mex_id_, cms, npatches, rad_near, ndtarg, targs, ntarg, nnz, n3, npatches, 1, npatches, 1, ndtarg, ntarg, 1, 1);

    row_ptr = zeros(ntarg+1,1);
    col_ind = zeros(nnz,1);
    ntp1 = ntarg+1;
    nnzp1 = nnz+1;
    mex_id_ = 'findnear(c i double[xx], c i int[x], c i double[x], c i int[x], c i double[xx], c i int[x], c io int[x], c io int[x])';
[row_ptr, col_ind] = fmm3dbie_routs(mex_id_, cms, npatches, rad_near, ndtarg, targs, ntarg, row_ptr, col_ind, n3, npatches, 1, npatches, 1, ndtarg, ntarg, 1, ntp1, nnz);

    iquad = zeros(nnz+1,1);
    mex_id_ = 'get_iquad_rsc(c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c io int[x])';
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

    zpars = complex(zeros(2,1));
    zpars(1) = zk;
    zpars(2) = alpha;
    zkuse = complex(zk);

    if strcmpi(qtype, 'rpcomb-bc')
      nker = 4;
      wnear = complex(zeros(nker,nquad));
      mex_id_ = 'getnearquad_helm_rpcomb_neu(c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[x], c i int[x], c io dcomplex[xx])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 2, 1, 1, ntp1, nnz, nnzp1, 1, 1, nker, nquad);
    elseif strcmpi(qtype, 'rpcomb-eval')
      nker = 2;
      wnear = complex(zeros(nker,nquad));
      mex_id_ = 'getnearquad_helm_rpcomb_eval(c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[xx], c i double[xx], c i int[x], c i int[x], c i double[xx], c i int[x], c i double[xx], c i double[x], c i dcomplex[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[x], c i int[x], c io dcomplex[xx])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, zkuse, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 1, 1, 1, ntp1, nnz, nnzp1, 1, 1, nker, nquad);
    else
      error('HELM3D.neumann.GET_QUADRATURE_CORRECTION:Unsupported quadrature correction');
    end
    
    Q = [];
    Q.targinfo = targinfo;
    Q.ifcomplex = 1;
    Q.wavenumber = zk;
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
        mex_id_ = 'rsc_to_csc(c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c io int[x], c io int[x], c io int[x])';
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
%
%
%-------------------------------------------------

