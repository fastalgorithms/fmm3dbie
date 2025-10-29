c  This file contains the following user callable routine
c
c    getnear_quad_str_guru - string specification of kernel
c
c   Note that in getnear_quad_str_guru, the real and imaginary parts
c   of the kernel are interlaced in a real array, and it is the
c   users responsibility to send in nker to be 2* the number of kernels
c   desired
      
      subroutine getnearquad_str_guru(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targvals,
     2   ipatch_id,uvs_targ,eps,str_pde,str_ker,ndd,dpars,ndz,zpars,
     3   ndi,ipars,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,
     4   nker,wnear)
c
c------------------------
c  This subroutine generates the near field quadrature
c  for a given kernel which is assumed to be
c  a compact/principal value integral operator
c  where the near field is specified by the user 
c  in row sparse compressed format.
c
c  Input arguments:
c
c    - npatches: integer
c        number of patches
c    - norders: integer(npatches)
c        order of discretization on each patch
c    - ixyzs: integer(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to patch i
c    - iptype: integer(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c        *  iptype = 11, quadrangular patch discretized using GL nodes,
c                        full degree polynomials
c        *  iptype = 12, quadrangular patch discretized using cheb nodes,
c                        full degree polynomials
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: double precision (9,npts)
c        basis expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch. For each patch
c        * srccoefs(1:3,i) is xyz info
c        * srccoefs(4:6,i) is dxyz/du info
c        * srccoefs(7:9,i) is dxyz/dv info
c    - srcvals: double precision (12,npts)
c        xyz(u,v) and derivative info sampled at the 
c        discretization nodes on the surface
c        * srcvals(1:3,i) - xyz info
c        * srcvals(4:6,i) - dxyz/du info
c        * srcvals(7:9,i) - dxyz/dv info
c        * srcvals(10:12,i) - normals info
c    - ndtarg: integer
c        leading dimension of target array
c    - ntarg: integer
c        number of targets
c    - targvals: double precision (ndtarg,ntarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - ipatch_id: integer(ntarg)
c        patch id of target, patch_id = -1, if target off-surface
c    - uvs_targ: double precision (2,ntarg)
c        local uv coordinates on patch if target on surface
c    - eps: double precision
c        precision requested
c    - str_pde: character(*)
c        string for specifying pde
c        * str_pde = 'laplace'
c        * str_pde = 'helmholtz'
c        * str_pde = 'yukawa'
c        * str_pde = 'maxwell'
c    - str_ker: character(*) 
c        string for specifying kernel
c        * str_ker = 'single'
c        * str_ker = 'sprime', normal derivative of single layer
c        * str_ker = 'double'
c        * str_ker = 'dprime' normal derivative of double layer
c        * str_ker = 'comb' combined field representation
c        * str_ker = 'sgradx' x component of gradient of single layer
c        * str_ker = 'sgrady' y component of gradient of single layer
c        * str_ker = 'sgradz' z component of gradient of single layer
c    - ipv: integer
c        Flag for choosing type of self-quadrature
c        * ipv = 0, for compact/weakly singular operators
c        * ipv = 1, for singular operators
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence 
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c    - ndd: integer
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer
c        number of integer parameters
c    - ipars: integer(ndi)
c        integer parameters
c    - iquadtype: integer
c        quadrature type
c        * iquadtype = 1, use ggq for self + adaptive integration for rest
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(ntarg+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - rfac0: integer
c        radius parameter for near field
c    - nquad: integer
c        number of entries in wnear
c  Output parameters:
c    - nker: integer
c        number of kernels
c    - wnear: double precision(:,:)
c        near field quadrature corrections
c----------------------------------------------------               
c
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: ndi,ndd,ndz,iquadtype
      integer, intent(in) :: ipars(ndi)
      integer, intent(in) :: ndtarg
      integer, intent(in) :: npatches,norders(npatches),npts
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer, intent(in) :: ntarg,ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: targvals(ndtarg,ntarg)
      real *8, intent(in) :: dpars(ndd)
      character (len=*), intent(in)  :: str_pde, str_ker
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(out) :: wnear(*)
c
c  temporary variables
c
c
      integer ipv, ifcomplex
      procedure (), pointer :: fker

c
c  helmholtz kernels
c
      external h3d_slp, h3d_dlp, h3d_sprime, h3d_comb, h3d_qlp
      external h3d_dprime_diff, h3d_slp_diff, h3d_dlp_diff
      external h3d_sprime_diff, h3d_sgradx, h3d_sgrady, h3d_sgradz
c
c  laplace kernels
c
      external l3d_slp, l3d_dlp, l3d_sprime, l3d_comb, l3d_qlp
      external l3d_sgradx, l3d_sgrady, l3d_sgradz
c
c  yukawa kernels
c
      external y3d_slp, y3d_dlp, y3d_sprime, y3d_comb
      external y3d_sgradx, y3d_sgrady, y3d_sgradz

      nker = 0
      
      if(iquadtype.eq.1) then
        ier = 0
        nker = 1
        if (lle(trim(str_pde), "laplace")) then
          ifcomplex = 0
          if (lle(trim(str_ker), "single")) then
            fker => l3d_slp
            ipv = 0
          elseif (lle(trim(str_ker), 'sprime')) then
            fker => l3d_sprime
            ipv = 1
          elseif (lle(trim(str_ker), 'double')) then
            fker => l3d_dlp
            ipv = 1
          elseif (lle(trim(str_ker), 'comb')) then
            fker => l3d_comb
            ipv = 1
          elseif (lle(trim(str_ker), 'sgradx')) then
            fker => l3d_sgradx
            ipv = 1
          elseif (lle(trim(str_ker), 'sgrady')) then
            fker => l3d_sgrady
            ipv = 1
          elseif (lle(trim(str_ker), 'sgradz')) then
            fker => l3d_sgradz
            ipv = 1
          else
            print *, "kernel not found for laplace"
            print *, "returning without computing anything"
            ier = 2
        endif
        elseif (lle(trim(str_pde), 'yukawa')) then
          ifcomplex = 0
          if (lle(trim(str_ker), 'single')) then
            fker => y3d_slp
            ipv = 0
          elseif (lle(trim(str_ker), 'sprime')) then
            fker => y3d_sprime
            ipv = 1
          elseif (lle(trim(str_ker), 'double')) then
            fker => y3d_dlp
            ipv = 1
          elseif (lle(trim(str_ker), 'comb')) then
            fker => y3d_comb
            ipv = 1
          elseif (lle(trim(str_ker), 'sgradx')) then
            fker => y3d_sgradx
            ipv = 1
          elseif (lle(trim(str_ker), 'sgrady')) then
            fker => y3d_sgrady
            ipv = 1
          elseif (lle(trim(str_ker), 'sgradz')) then
            fker => y3d_sgradz
            ipv = 1
          else
            print *, "kernel not found for Yukawa"
            print *, "returning without computing anything"
            ier = 2
          endif
        elseif (lle(trim(str_pde), 'helmholtz')) then
          ifcomplex = 1
          if (lle(trim(str_ker), 'single')) then
            fker => h3d_slp
            ipv = 0
          elseif (lle(trim(str_ker), 'sprime')) then
            fker => h3d_sprime
            ipv = 1
          elseif (lle(trim(str_ker), 'double')) then
            fker => h3d_dlp
            ipv = 1
          elseif (lle(trim(str_ker), 'comb')) then
            fker => h3d_comb
            ipv = 1
          elseif (lle(trim(str_ker), 'sgradx')) then
            fker => h3d_sgradx
            ipv = 1
          elseif (lle(trim(str_ker), 'sgrady')) then
            fker => h3d_sgrady
            ipv = 1
          elseif (lle(trim(str_ker), 'sgradz')) then
            fker => h3d_sgradz
            ipv = 1
          else
            print *, "kernel not found for Helmholtz"
            print *, "returning without computing anything"
            ier = 2
          endif
        endif

        
        if (ier.ne.0) return

        if(nker.eq.1) then 
          if(ifcomplex.eq.1) then
            call zgetnearquad_ggq_guru(npatches,norders,
     1       ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targvals,
     2       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,
     3       ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
          else
            call dgetnearquad_ggq_guru(npatches,norders,
     1       ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targvals,
     2       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,
     3       ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
          endif
        endif
      endif

      return
      end



