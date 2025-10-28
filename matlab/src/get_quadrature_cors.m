function Q = get_quadrature_cors(S, eps, ipars,dpars,zpars, fort_quad_handle, kern,zk,kernel_order, targinfo)
%
%  get_quadrature_cors
%    This subroutine returns the near quadrature corrections for the kernel
%    kern. fort_handle is a handle for the mexified quadrature routine.
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * eps: precision requested
%    * ipars: integer parameters for kernel, typically vector indices
%    * dpars: real parameters for kern
%    * zpars: complex parameters for kern
%    * fort_quad_handle: handle for fortran quadrature routine, see ....
%    * kern: kernel evaluator, if empty, return the smooth rule
%    * zk: typical problem wavelength, used to estimate oversampling
%    * kernel_order: kernel singularity strength, used to estimate
%       oversampling order, set to -1 if empty
%       kernel_order = -1 -> single layer type operator
%       kernel_order = 0 -> Double layer type operator
%       kernel_order = 1 -> derivative of double layer
%    * targinfo: target info (optional), default is to use S
%       targinfo.r = (3,nt) target locations
%       targinfo.du = u tangential derivative info
%       targinfo.dv = v tangential derivative info
%       targinfo.n = normal info
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%          patch if on-surface (optional)
%

    [srcvals,srccoefs,norders,ixyzs,iptype,~] = extract_arrays(S);
    npatches = S.npatches;

    if nargin < 11
      targinfo = [];
      targinfo.r = S.r;
      targinfo.du = S.du;
      targinfo.dv = S.dv;
      targinfo.n = S.n;
      targinfo.patch_id = S.patch_id;
      targinfo.uvs_targ = S.uvs_targ;
    end

    if isempty(kernel_order), kernel_order = -1; end

    targs = extract_targ_array(targinfo); 
    [ndtarg,ntarg] = size(targs);
    
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

     ff = 'sparse';
    % if(isfield(opts,'format'))
    %    ff = opts.format;
    % end

    if(~(strcmpi(ff,'rsc') || strcmpi(ff,'csc') || strcmpi(ff,'sparse')))
       fprintf('invalid quadrature format, reverting to rsc format\n');
       ff = 'rsc';
    end

    %
    % start processing
    %

    iptype_avg = floor(sum(iptype)/(npatches+0.0d0));
    norder_avg = floor(sum(norders)/(npatches+0.0d0));

    % get nearfield definition
    [rfac, rfac0] = get_rfacs(norder_avg,iptype_avg);
    [cms, rads] = get_centroid_rads(npatches,norders,ixyzs,iptype,npts, ... 
       srccoefs);

    rad_near = rads*rfac;

    % find near quadrature correction interactions
    nnz = findnearmem(cms, npatches, rad_near, ndtarg, targs,ntargs);
    [row_ptr,col_ind] = findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,nnz);

    iquad = get_iquad_rsc(npatches,ixyzs,npts,ntarg,nnz,row_ptr,col_ind);
    nquad = iquad(nnz+1)-1;

    % get quadratures
    iquadtype = 1;
    wnear = fort_quad_handle(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,eps,iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,targinfo,nquad);

    %% tidy up
    Q = [];
    Q.targinfo = targinfo;
    Q.ifcomplex = ~isreal(wnear);
    Q.wavenumber = zk;
    Q.kernel_order = kernel_order;
    Q.rfac = rfac;
    Q.nquad = nquad;
    Q.format = ff;

    if icor
        novers = get_oversampling_parameters(eps,S,Q);
        novers = min(novers(:),S.norders(:));
        wnear_smth = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,novers); 
        wnear = wnear - wnear_smth;
    else
        novers = S.norders;
    end
    Q.novers = novers;

    if(strcmpi(ff,'rsc'))
        Q.iquad = iquad;
        Q.wnear = wnear;
        Q.row_ptr = row_ptr;
        Q.col_ind = col_ind;
    elseif(strcmpi(ff,'csc'))
        [col_ptr, row_ind, iper] = rsc_to_csc(npatches,ntarg,row_ptr,col_ind,nnz);
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
