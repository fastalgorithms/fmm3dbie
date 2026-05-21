!
!  Near-field quadrature for capillary (surface tension) wave problems.
!
!  PDE: Laplace equation in the fluid with a capillary free-surface
!  boundary condition (surface tension modifies the dispersion relation).
!  The velocity potential phi satisfies:
!    Delta phi = 0  in the fluid
!    d phi/dn = f   on the body boundary
!
!  Representation: single layer potential using the capillary wave
!  Green's functions G_s and G_phi built from 3 cubic dispersion roots:
!    phi(x) = int G_s(x,y) sigma(y) dS(y)
!
!  zpars(1:3) = cubic dispersion roots
!  zpars(4:6) = corresponding residues
!
!  iker selects which kernel to assemble:
!    iker = 0 -> G_s         (gshelmkern)
!    iker = 1 -> G_phi       (gphihelmkern)
!    iker = 2 -> Delta G_s   (lapgshelmkern)
!    iker = 3 -> Delta G_phi (lapgphihelmkern)
!    iker = 5 -> G_{3d,phi}  (s3dgphihelmkern)
!

      subroutine getnearquad_capillary_all(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        eps, zpars, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, iker, wnear)
!
!  This subroutine generates the near field quadrature for the
!  capillary wave single layer representation:
!
!    phi(x) = int G_s(x,y) sigma(y) dS(y)
!
!  On imposing the Neumann boundary condition d phi/dn = f, we get:
!
!    sigma/2 + int d G_s/dn(x,y) sigma(y) dS(y) = f
!
!  The quadrature is computed by the following strategy:
!  targets within a sphere of radius rfac0*rs of a patch centroid
!  are handled using adaptive integration (ggq_guru);
!  all other near-field targets use oversampled quadrature.
!
!  The recommended parameter for rfac0 is 1.25d0
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
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du, and dxyz/dv
!    - srcvals: real *8 (12,npts)
!        xyz, derivatives, and normals at discretization nodes
!    - eps: real *8
!        precision requested
!    - zpars: complex *16(6)
!        zpars(1:3) = cubic dispersion roots
!        zpars(4:6) = corresponding residues
!    - iquadtype: integer
!        quadrature type (iquadtype = 1: ggq + adaptive)
!    - nnz: integer
!        number of source patch -> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) points to start of source patch list for target i
!    - col_ind: integer(nnz)
!        list of source patches for all targets
!    - iquad: integer(nnz+1)
!        iquad(i) is start of quadrature in wnear for col_ind(i)
!    - rfac0: real *8
!        radius parameter for switching to predetermined quadrature rule
!    - nquad: integer
!        number of near field quadrature entries
!    - iker: integer
!        kernel selector:
!        iker = 0 -> G_s (gshelmkern)
!        iker = 1 -> G_phi (gphihelmkern)
!        iker = 2 -> Delta G_s (lapgshelmkern)
!        iker = 3 -> Delta G_phi (lapgphihelmkern)
!        iker = 5 -> G_{3d,phi} (s3dgphihelmkern)
!
!  Output arguments:
!    - wnear: complex *16(nquad)
!        near field quadrature corrections for the selected kernel

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars(6)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      complex *16, intent(out) :: wnear(nquad)
      
      complex *16 zpars_tmp(3)
      integer *8 ipars(2)
      real *8 dpars(1)
      

      real *8, allocatable :: uvs_targ(:,:)
      integer *8, allocatable :: ipatch_id(:)


      integer *8 ipv, i, ndi, ndd, ndz
      
      integer *8 ndtarg, ntarg

      integer *8 iker

      procedure (), pointer :: fker
      external gphihelmkern, p3d_log, gshelmkern, lapgshelmkern
      external s3dgphihelmkern, lapgphihelmkern

      ndz=6
      ndd=0
      ndi=0
      ndtarg = 12
      ntarg = npts

      allocate(ipatch_id(npts),uvs_targ(2,npts))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO      

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)
      if (iquadtype.eq.1) then
        ipv=0

        if (iker.eq.0) then
          fker => gshelmkern 
        elseif (iker.eq.1) then
          fker => gphihelmkern
        elseif (iker.eq.2) then
          fker => lapgshelmkern
        elseif (iker.eq.3) then
          fker => lapgphihelmkern
        elseif (iker.eq.5) then
          fker => s3dgphihelmkern        
        endif
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wnear)

      endif

      return
      end subroutine getnearquad_capillary_all
!
!
!
      subroutine getnearquad_capillary_all_vpp(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        eps, zpars0, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, iker, maxdist, wnear)
!
!  _vpp variant of getnearquad_capillary_all: uses the variable-point-pair
!  (vpp) quadrature scheme to handle the logarithmic singularity of G_s and
!  G_phi, precomputing a kernel table via vpp_buildkern before calling
!  zgetnearquad_ggq_guru.
!
!  Input arguments: same as getnearquad_capillary_all, plus:
!    - maxdist: real *8
!        maximum source-target distance for the vpp kernel table
!
!  Output arguments:
!    - wnear: complex *16(nquad)
!        near field quadrature corrections for the selected kernel
  
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars0(6)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0, maxdist
      integer *8, intent(in) :: nquad
      complex *16, intent(out) :: wnear(nquad)
      complex *16 zpars_tmp(3), zpars
      integer *8 ipars(5), ipars0
      integer *8, parameter :: ldpars = 10000000
      real *8 :: dpars0
      real *8, allocatable :: dpars(:)
      real *8 :: tol,pi,done,a,b
      real *8, allocatable :: uvs_targ(:,:)
      integer *8, allocatable :: ipatch_id(:)
      integer *8 ipv, i, ndi, ndd, ndz, ndi0, ndd0, ndz0
      integer *8 ndtarg, ntarg, n, nf, maxsub, maxdepth, ietype
      integer *8 iker, ier
      procedure (), pointer :: fker
      external gphihelmkern, gshelmkern, lapgshelmkern
      external s3dgphihelmkern, lapgphihelmkern, vpp_kern

      allocate(dpars(ldpars))

!      parameters defining the true kernel
      ndz0=6
      ndd0=0
      ndi0=0
      
      ndtarg = 12
      ntarg = npts

      allocate(ipatch_id(npts),uvs_targ(2,npts))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO      

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)
      if (iquadtype.eq.1) then
        ipv=0

        done = 1
        pi = atan(done)*4
        nf = 2
        a = 1d-12
        b = maxdist
        n = 16
        tol = 1d-12
        maxsub = 1000
        maxdepth = 50
        ietype = 1
        ier = 0

        if (iker.eq.0) then
          call vpp_buildkern(gshelmkern,ndd0,dpars0,ndz0,zpars0,ndi0, &
          ipars0,nf,a,b,n,tol,maxsub,maxdepth,ietype, &
          ipars,ndd,ldpars,dpars,ier)
        elseif (iker.eq.1) then
          call vpp_buildkern(gphihelmkern,ndd0,dpars0,ndz0,zpars0,ndi0, &
          ipars0,nf,a,b,n,tol,maxsub,maxdepth,ietype, &
          ipars,ndd,ldpars,dpars,ier)
        elseif (iker.eq.2) then
          call vpp_buildkern(lapgshelmkern,ndd0,dpars0,ndz0,zpars0,ndi0, &
          ipars0,nf,a,b,n,tol,maxsub,maxdepth,ietype, &
          ipars,ndd,ldpars,dpars,ier)
        elseif (iker.eq.3) then
          call vpp_buildkern(lapgphihelmkern,ndd0,dpars0,ndz0,zpars0,ndi0, &
          ipars0,nf,a,b,n,tol,maxsub,maxdepth,ietype, &
          ipars,ndd,ldpars,dpars,ier)
        elseif (iker.eq.5) then
          call vpp_buildkern(s3dgphihelmkern,ndd0,dpars0,ndz0,zpars0,ndi0, &
          ipars0,nf,a,b,n,tol,maxsub,maxdepth,ietype, &
          ipars,ndd,ldpars,dpars,ier)
        endif

        fker => vpp_kern

        ndi = 5
        ndz = 0
        zpars = 0

        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wnear)

      endif

      return
      end subroutine getnearquad_capillary_all_vpp
!
!
!
