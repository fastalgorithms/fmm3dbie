
      subroutine getnearquad_polynya_lps_neu_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = (Enter representation here) 
!
!  On imposing the boundary condition, we get the following operator
!
!  du/dn + ik \lambda u =  
!    z \sigma + S_{k}'[\sigma] + i\alpha S_{i|k|}'^2 [\sigma] + i \alpha 
!       (D_{k}' - D_{i|k|}') S_{i|k|}[\sigma]  + 
!       ik \lambda (S_{k} + i \alpha D_{k} S_{i|k|} + 
!       i \alpha w S_{i|k|}) = f
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a patch centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!  
! 
!  Input arguments:
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev 
!                     nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - ndtarg: integer
!        leading dimension of target information array
!    - ntarg: integer
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information 
!    - ipatch_id: integer(ntarg)
!        ipatch_id(i) indicates the patch on which target i
!        is, if it is on surface. ipatch_id(i) should be 0 
!        otherwise
!    - uvs_targ: real *8(2,ntarg)
!        if ipatch_id(i) > 0, then uvs_targ(1:2,i) are the
!        local uv coordinates of the target on the patch,
!        uvs_targ(1:2,i) is unused otherwise
!    - eps: real *8
!        precision requested
!    - zpars: complex *16 (1)
!        kernel parameters 
!        zpars(1) = k, Helmholtz wavenumber
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(ntarg+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: real *8
!        radius parameter for switching to predetermined quadarature
!        rule        
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair
!
!  Output arguments
!    - wnear: complex *16(2,nquad)
!        The desired near field quadrature
!        wnear(1,:) - stores the quadrature corrections for ker1 
!        wnear(2,:) - stores the quadrature correction for ker2
!               
!

      implicit none 
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      integer, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg, ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars(1)
      integer, intent(in) :: iquadtype, nnz
      integer, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer, intent(in) :: nquad
      complex *16, intent(out) :: wnear(2,nquad)

!  Temporary variables      
      complex *16 zk, ima, zpars_tmp(3)
      integer ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer ipars
      complex *16, allocatable :: wneartmp(:)
  
      procedure (), pointer :: fker
! replace kernels here      
      external h3d_slp, h3d_dlp

      data ima/(0.0d0,1.0d0)/
    
      ndd = 0
      ndz = 1
      ndi = 0


      allocate(wneartmp(nquad))
    
    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => h3d_slp
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(1,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO

      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0 
        enddo
!$OMP END PARALLEL DO

        ipv = 1
        fker => h3d_dlp
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(2,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO
      endif


      return
      end
!


