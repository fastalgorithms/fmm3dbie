function Q = get_quadrature_correction(S, eps, om, rep_params)
%
%  em3d.dielectric.get_quadrature_correction
%    This subroutine returns the near quadrature correction
%    for the chosen maxwell representation, with densities supported
%    on the surface, and targets given by targinfo 
%    as a cell array of sparse matrices or an array of matrices
%    in the rsc format, where each folumn of the matrix is the
%    representation of the sparse matrix corresponding to one
%    of the kernels
%
%  Notes for this routine:
%  The PDE takes the form
%  1v) \nabla \times E =  i\om \mu H
%  2v) \nabla \cdot  E =     0
%  3v) \nabla \times H = -i\om \ep E
%  4v) \nabla \cdot  H =     0
%  
%  where E is the electric field, H is the magnetic field, 
%  and \om is the wavenumber, \ep is the permittivity, 
%  and \mu is the permeability
%
%  The dielectric boundary conditions are given by
%  1b) n \times (E0 + E_in) = n \times E1
%  2b) n \times (H0 + H_in) = n \times H1
%
%  where (E_in, H_in) are the incoming electric and magnetic 
%  fields, E0, H0 are the fields in the exterior, and
%  (E1, H1) are the fields in the interior
%
%  Syntax
%   Q = em3d.pec.get_quadrature_correction(S,eps,zk,rep_params)
%
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * eps: precision requested
%    * om: wave number
%    * rep_params: parameters for integral representation 
%                  for muller, it should be a 4 vector
%                  consisting of [\ep0, \mu0, \ep1, \mu1]
%

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;
    npp1 = npatches+1;
    nsp1 = npts + 1;
    n3 = 3;


    ff = 'rsc';

    tinfouse = [];
    tinfouse.r = S.r;
    tinfouse.du = S.du;
    tinfouse.dv = S.dv;
    tinfouse.n = S.n;
    tinfouse.patch_id = S.patch_id;
    tinfouse.uvs_targ = S.uvs_targ;


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

    nker = 16;
    zpars = complex([om rep_params(:).']);
    wnear = complex(zeros(nquad,nker));
    mex_id_ = 'getnearquad_em_muller_trans(c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[xx], c i double[xx], c i int[x], c i int[x], c i double[xx], c i int[x], c i double[xx], c i double[x], c i dcomplex[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[x], c i int[x], c io dcomplex[xx])';
[wnear] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 5, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad, nker);

    Q = [];
    Q.targinfo = tinfouse;
    Q.ifcomplex = 1;
    Q.wavenumber = om;
    Q.kernel_order = -1;
    Q.rfac = rfac;
    Q.nquad = nquad;
    Q.format = ff;
    Q.kernel_order = 0;

    Q.iquad = iquad;
    Q.wnear = wnear;
    Q.row_ptr = row_ptr;
    Q.col_ind = col_ind;
    
end
%
%
%
%-------------------------------------------------

