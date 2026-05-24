cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Combined prime representation for acoustic scattering
c
c  PDE:
c    (\Delta + k^2) u = 0
c
c  This file provides quadrature and evaluation for the
c  combined prime kernel:
c
c    K[\sigma] = (\alpha S_{k}' + \beta D_{k}')[\sigma]
c
c  where
c    S_{k}' = d/dn_x G_k   (normal derivative at target)
c    D_{k}' = d/dn_x d/dn_y G_k  (hypersingular)
c
c  zpars(1) = k
c  zpars(2) = alpha
c  zpars(3) = beta
c
c  User callable routines:
c    - getnearquad_helm_comb_cprime: compute near quadrature for
c        on-surface targets
c
c    - getnearquad_helm_comb_cprime_eval: compute near quadrature for
c        general (on or off surface) targets
c
c    - helm_comb_cprime_eval: evaluate the layer potential at a
c        collection of targets using the combined prime kernel
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine getnearquad_helm_comb_cprime(npatches, norders,
     1   ixyzs, iptype, npts, srccoefs, srcvals,
     2   eps, zpars, iquadtype, nnz, row_ptr, col_ind,
     3   iquad, rfac0, nquad, wnear)
c
c       this subroutine generates the near field quadrature
c       for the combined prime representation
c       K = (\alpha S_{k}' + \beta D_{k}')[\sigma]
c       where the near field is specified by the user
c       in row sparse compressed format, and targets are on-surface.
c
c       The quadrature is computed by the following strategy:
c        targets within a sphere of radius rfac0*rs of a patch centroid
c        are handled using adaptive integration, where rs is the radius
c        of the bounding sphere for the patch.
c        All other targets in the near field are handled via
c        oversampled quadrature.
c
c       The recommended parameter for rfac0 is 1.25d0
c
c       input:
c         npatches - integer *8
c            number of patches
c         norders - integer *8(npatches)
c            order of discretization on each patch
c         ixyzs - integer *8(npatches+1)
c            starting location of data on patch i
c         iptype - integer *8(npatches)
c            type of patch
c         npts - integer *8
c            total number of discretization points
c         srccoefs - real *8 (9,npts)
c            koornwinder expansion coefficients of xyz, dxyz/du, dxyz/dv
c         srcvals - real *8 (12,npts)
c            xyz(u,v), dxyz/du, dxyz/dv, normals
c         eps - real *8
c            precision requested
c         zpars - complex *16(3)
c            zpars(1) = k, zpars(2) = alpha, zpars(3) = beta
c         iquadtype - integer *8
c            quadrature type (1 = ggq)
c         nnz - integer *8
c            number of source patch -> target interactions in near field
c         row_ptr - integer *8(npts+1)
c         col_ind - integer *8(nnz)
c         iquad - integer *8(nnz+1)
c         rfac0 - real *8
c            radius parameter for near field
c         nquad - integer *8
c            number of entries in wnear
c
c        output:
c            wnear - complex *16(nquad)
c

      implicit none
      integer *8, intent(in) :: npatches,norders(npatches),npts,nquad
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: iquadtype
      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      complex *16, intent(in) :: zpars(3)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(npts+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      complex *16, intent(out) :: wnear(nquad)
      integer *8 i
      integer *8 ndtarg,ntarg

      ndtarg = 12
      ntarg = npts
      allocate(ipatch_id(npts), uvs_targ(2,npts))

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
C$OMP END PARALLEL DO

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)

      call getnearquad_helm_comb_cprime_eval(npatches, norders, ixyzs,
     1  iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals,
     2  ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr,
     3  col_ind, iquad, rfac0, nquad, wnear)

      return
      end
c
c
c


      subroutine getnearquad_helm_comb_cprime_eval(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear)
c
c       this subroutine generates the near field quadrature
c       for the combined prime representation
c       K = (\alpha S_{k}' + \beta D_{k}')[\sigma]
c       where the near field is specified by the user
c       in row sparse compressed format.
c       Targets can be on-surface or off-surface.
c
c       ipv = 2 is used for hypersingular quadrature.
c
c       input:
c         npatches - integer *8
c         norders - integer *8(npatches)
c         ixyzs - integer *8(npatches+1)
c         iptype - integer *8(npatches)
c         npts - integer *8
c         srccoefs - real *8 (9,npts)
c         srcvals - real *8 (12,npts)
c         ndtarg - integer *8
c            leading dimension of target array (must be >= 12)
c         ntarg - integer *8
c            number of targets
c         targs - real *8 (ndtarg,ntarg)
c            target information (must include normals in rows 10:12)
c         ipatch_id - integer *8(ntarg)
c            patch id of target (-1 if off-surface)
c         uvs_targ - real *8 (2,ntarg)
c         eps - real *8
c         zpars - complex *16(3)
c            zpars(1) = k, zpars(2) = alpha, zpars(3) = beta
c         iquadtype - integer *8 (1 = ggq)
c         nnz - integer *8
c         row_ptr - integer *8(ntarg+1)
c         col_ind - integer *8(nnz)
c         iquad - integer *8(nnz+1)
c         rfac0 - real *8
c         nquad - integer *8
c
c        output:
c            wnear - complex *16(nquad)
c

      implicit none
      integer *8, intent(in) :: npatches,norders(npatches),npts,nquad
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: ndtarg,ntarg
      integer *8, intent(in) :: iquadtype
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      complex *16, intent(in) :: zpars(3)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      complex *16, intent(out) :: wnear(nquad)

      integer *8 ipars
      integer *8 ndd,ndz,ndi
      real *8 dpars

      complex *16 alpha,beta
      integer *8 i
      integer *8 ipv

      procedure (), pointer :: fker
      external h3d_combprime

      alpha = zpars(2)
      beta  = zpars(3)

      ndz = 3
      ndi = 0
      ndd = 0

      if(iquadtype.eq.1) then
        ipv = 2
        fker => h3d_combprime

        call zgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1     ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,
     1     ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      endif

      return
      end
c
c
c


      subroutine helm_comb_cprime_eval(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,sigma,pot)
c
cf2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
cf2py intent(in) ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps,zpars
cf2py intent(in) sigma
cf2py intent(out) pot
c
c------------------------------
c  This subroutine evaluates the layer potential for the representation
c
c    K[\sigma] = (\alpha S_{k}' + \beta D_{k}')[\sigma]
c
c  where
c    S_{k}' = d/dn_x G_k   (normal derivative at target x)
c    D_{k}' = d/dn_x d/dn_y G_k  (hypersingular)
c
c  The far field is computed via FMM using gradient output at targets,
c  then dotted with the target normal.
c  The near field uses the precomputed hypersingular quadrature (ipv=2).
c
c  Note: For targets on the boundary, this routine computes only the
c  principal value part.
c
c  Input arguments:
c    - npatches: integer *8
c    - norders: integer *8(npatches)
c    - ixyzs: integer *8(npatches+1)
c    - iptype: integer *8(npatches)
c    - npts: integer *8
c    - srccoefs: real *8(9,npts)
c    - srcvals: real *8(12,npts)
c    - ndtarg: integer *8
c        leading dimension of target array (must be >= 12)
c    - ntarg: integer *8
c    - targs: real *8(ndtarg,ntarg)
c        target info; rows 10:12 must contain target normals
c    - ipatch_id: integer *8(ntarg)
c        -1 if off-surface
c    - uvs_targ: real *8(2,ntarg)
c    - eps: real *8
c    - zpars: complex *16(3)
c        zpars(1) = k, zpars(2) = alpha, zpars(3) = beta
c    - sigma: complex *16(npts)
c
c  Output arguments:
c    - pot: complex *16(ntarg)
c
c-----------------------------------

      implicit none
      integer *8, intent(in) :: npatches,npts
      integer *8, intent(in) :: ndtarg,ntarg
      integer *8, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      complex *16, intent(in) :: zpars(3)
      complex *16, intent(in) :: sigma(npts)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      complex *16, intent(out) :: pot(ntarg)

      integer *8 nptso,nnz,nquad

      integer *8 nover,npolso
      integer *8 norder,npols
      integer *8, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer *8, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)

      integer *8 i,j,jpatch,jquadstart,jstart

      integer *8 ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime,work

      real *8 ttot,done
      real *8 rfac,rfac0
      integer *8 iptype_avg,norder_avg
      integer *8 ikerorder, iquadtype,npts_over
      integer *8 ndd,ndz,ndi,nker,lwork,ndim

      complex *16 alpha,beta

      integer *8 ns
      complex *16, allocatable :: sigmaover(:)
      real *8, allocatable :: sources(:,:),srctmp(:,:),srctmp2(:,:)
      complex *16, allocatable :: charges(:),dipvec(:,:)
      complex *16, allocatable :: ctmp2(:),dtmp2(:,:)
      complex *16, allocatable :: pot_aux(:),grad_aux(:,:)

      integer *8 ifcharge,ifdipole,ifpgh,ifpghtarg
      integer *8 ier
      complex *16 pottmp,gradtmp(3)
      integer *8 nss,l
      real *8 thresh,ra
      integer *8 nmax
      integer *8 nd,ntarg0
      integer *8 int8_3,int8_12

      parameter (nd=1,ntarg0=1)
      int8_3  = 3
      int8_12 = 12
      done = 1

      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)

      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,
     1     srccoefs,cms,rads)

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
C$OMP END PARALLEL DO

      call findnearmem(cms,npatches,rad_near,ndtarg,targs,ntarg,nnz)

      allocate(row_ptr(ntarg+1),col_ind(nnz))

      call findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,row_ptr,
     1        col_ind)

      allocate(iquad(nnz+1))
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind,
     1         iquad)

c     ikerorder: 0 for S' only (alpha only), 1 if D' is present
      ikerorder = 0
      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 1

      allocate(novers(npatches),ixyzso(npatches+1))

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,ntarg,targs,ikerorder,zpars(1),
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)

      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad))

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO

      iquadtype = 1

      call getnearquad_helm_comb_cprime_eval(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)

c
c   compute layer potential via FMM
c   S'[sigma](x) = grad_x S[sigma](x) . n(x)
c   D'[sigma](x) = grad_x D[sigma](x) . n(x)
c
c   Far field: call hfmm3d_t_cd_g with
c     charges = alpha * sigma_over * wts  (for S')
c     dipoles = beta  * sigma_over * wts * n_src  (for D')
c   then dot grad_aux with target normal
c
      ns = npts_over
      alpha = zpars(2)
      beta  = zpars(3)

      allocate(sources(3,ns),srctmp(3,ntarg))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(ns))
      allocate(pot_aux(ntarg),grad_aux(3,ntarg))

      call oversample_fun_surf(2_8,npatches,norders,ixyzs,iptype,
     1    npts,sigma,novers,ixyzso,ns,sigmaover)

      ra = 0

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
        charges(i) = alpha*sigmaover(i)*wover(i)
        dipvec(1,i) = beta*sigmaover(i)*wover(i)*srcover(10,i)
        dipvec(2,i) = beta*sigmaover(i)*wover(i)*srcover(11,i)
        dipvec(3,i) = beta*sigmaover(i)*wover(i)*srcover(12,i)
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        srctmp(1,i) = targs(1,i)
        srctmp(2,i) = targs(2,i)
        srctmp(3,i) = targs(3,i)
        pot(i)      = 0
        pot_aux(i)  = 0
        grad_aux(1,i) = 0
        grad_aux(2,i) = 0
        grad_aux(3,i) = 0
      enddo
C$OMP END PARALLEL DO

      ifcharge = 1
      ifdipole = 1
      if(abs(alpha).lt.1.0d-16) ifcharge = 0
      if(abs(beta).lt.1.0d-16)  ifdipole = 0

      ier = 0
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        call hfmm3d_t_cd_g(eps,zpars(1),ns,sources,charges,dipvec,
     1    ntarg,srctmp,pot_aux,grad_aux,ier)
      elseif(ifcharge.eq.1.and.ifdipole.eq.0) then
        call hfmm3d_t_c_g(eps,zpars(1),ns,sources,charges,
     1    ntarg,srctmp,pot_aux,grad_aux,ier)
      elseif(ifcharge.eq.0.and.ifdipole.eq.1) then
        call hfmm3d_t_d_g(eps,zpars(1),ns,sources,dipvec,
     1    ntarg,srctmp,pot_aux,grad_aux,ier)
      endif

c     dot gradient with target normal
C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        pot(i) = grad_aux(1,i)*targs(10,i) +
     1           grad_aux(2,i)*targs(11,i) +
     2           grad_aux(3,i)*targs(12,i)
      enddo
C$OMP END PARALLEL DO

      call get_fmm_thresh(int8_12,ns,srcover,int8_12,ntarg,srctmp,
     1  thresh)

c     Add precomputed near-field quadrature
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,pottmp,npols,l)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            pot(i) = pot(i) + wnear(jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO

c     Remove near contribution from FMM (add-subtract)
      call get_near_corr_max(ntarg,row_ptr,nnz,col_ind,npatches,
     1  ixyzso,nmax)
      allocate(srctmp2(3,nmax),ctmp2(nmax),dtmp2(3,nmax))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,pottmp,gradtmp)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            if(ifcharge.eq.1) ctmp2(nss)   = charges(l)
            if(ifdipole.eq.1) then
              dtmp2(1,nss) = dipvec(1,l)
              dtmp2(2,nss) = dipvec(2,l)
              dtmp2(3,nss) = dipvec(3,l)
            endif
          enddo
        enddo

        pottmp = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          call h3ddirectcdg(nd,zpars(1),srctmp2,ctmp2,dtmp2,
     1      nss,srctmp(1,i),ntarg0,pottmp,gradtmp,thresh)
        elseif(ifcharge.eq.1.and.ifdipole.eq.0) then
          call h3ddirectcg(nd,zpars(1),srctmp2,ctmp2,
     1      nss,srctmp(1,i),ntarg0,pottmp,gradtmp,thresh)
        elseif(ifcharge.eq.0.and.ifdipole.eq.1) then
          call h3ddirectdg(nd,zpars(1),srctmp2,dtmp2,
     1      nss,srctmp(1,i),ntarg0,pottmp,gradtmp,thresh)
        endif

        pot(i) = pot(i) - (gradtmp(1)*targs(10,i) +
     1                     gradtmp(2)*targs(11,i) +
     2                     gradtmp(3)*targs(12,i))
      enddo
C$OMP END PARALLEL DO

      return
      end
c
c
c


      subroutine helm_comb_cprime_eval_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   eps,ndd,dpars,ndz,zpars,ndi,ipars,
     3   nnz,row_ptr,col_ind,iquad,nquad,nker,wnear,novers,
     4   nptso,ixyzso,srcover,whtsover,lwork,work,idensflag,
     5   ndim_s,sigma,ipotflag,ndim_p,pot)
c
c------------------------------
c  This subroutine evaluates the combined prime layer potential
c
c    K[\sigma] = (\alpha S_{k}' + \beta D_{k}')[\sigma]
c
c  where the near field quadrature is precomputed and stored
c  in row sparse compressed format, and the oversampled geometry
c  is provided by the caller.
c
c  For targets on the boundary, this routine returns only the
c  principal value part.
c
c  The FMM computes grad_x (alpha*S_k + beta*D_k)[sigma] at targets,
c  dotted with the target normal n_x, using add-subtract to remove
c  the near contribution and replace it with the precomputed wnear.
c
c  Input arguments:
c    - npatches: integer *8
c    - norders: integer *8(npatches)
c    - ixyzs: integer *8(npatches+1)
c    - iptype: integer *8(npatches)
c    - npts: integer *8
c    - srccoefs: real *8(9,npts)
c    - srcvals: real *8(12,npts)
c    - ndtarg: integer *8
c        leading dimension of target array (must be >= 12)
c    - ntarg: integer *8
c    - targs: real *8(ndtarg,ntarg)
c        rows 10:12 must contain target normals
c    - eps: real *8
c    - ndd: integer *8 (unused, pass 0)
c    - dpars: real *8(ndd) (unused)
c    - ndz: integer *8
c        must be 3
c    - zpars: complex *16(ndz)
c        zpars(1) = k, zpars(2) = alpha, zpars(3) = beta
c    - ndi: integer *8 (unused, pass 0)
c    - ipars: integer *8(ndi) (unused)
c    - nnz: integer *8
c    - row_ptr: integer *8(ntarg+1)
c    - col_ind: integer *8(nnz)
c    - iquad: integer *8(nnz+1)
c    - nquad: integer *8
c    - nker: integer *8 (must be 1)
c    - wnear: complex *16(nker,nquad)
c        precomputed near quadrature for K = alpha*S_k' + beta*D_k'
c    - novers: integer *8(npatches)
c    - nptso: integer *8
c    - ixyzso: integer *8(npatches+1)
c    - srcover: real *8(12,nptso)
c    - whtsover: real *8(nptso)
c    - lwork: integer *8 (unused)
c    - work: real *8(lwork) (unused)
c    - idensflag: integer *8 (unused)
c    - ndim_s: integer *8 (must be 1)
c    - sigma: complex *16(npts)
c    - ipotflag: integer *8 (unused)
c    - ndim_p: integer *8 (must be 1)
c
c  Output arguments:
c    - pot: complex *16(ntarg)
c

      implicit none
      integer *8, intent(in) :: npatches,npts
      integer *8, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer *8, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer *8, intent(in) :: ndd,ndz,ndi
      real *8, intent(in) :: dpars(max(ndd,1))
      complex *16, intent(in) :: zpars(ndz)
      integer *8, intent(in) :: ipars(max(ndi,1))
      integer *8, intent(in) :: nnz,nquad
      integer *8, intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      integer *8, intent(in) :: nker
      complex *16, intent(in) :: wnear(nker,nquad)
      integer *8, intent(in) :: nptso
      integer *8, intent(in) :: ixyzso(npatches+1),novers(npatches)
      real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
      integer *8, intent(in) :: lwork
      real *8, intent(in) :: work(max(lwork,1))
      integer *8, intent(in) :: ndim_s,ndim_p,idensflag,ipotflag
      complex *16, intent(in) :: sigma(npts)
      complex *16, intent(out) :: pot(ntarg)

      integer *8 ns
      complex *16 alpha,beta
      integer *8 ifcharge,ifdipole

      real *8, allocatable :: sources(:,:),srctmp(:,:),srctmp2(:,:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      complex *16, allocatable :: ctmp2(:),dtmp2(:,:)
      complex *16, allocatable :: pot_aux(:),grad_aux(:,:)

      integer *8 i,j,jpatch,jquadstart,jstart
      integer *8 nss,l,npols,nmax
      complex *16 pottmp,gradtmp(3)
      real *8 thresh

      integer *8 ier
      integer *8 nd,ntarg0
      integer *8 int8_12

      parameter (nd=1,ntarg0=1)
      int8_12 = 12

      ns = nptso
      alpha = zpars(2)
      beta  = zpars(3)

      ifcharge = 1
      ifdipole = 1
      if(abs(alpha).lt.1.0d-16) ifcharge = 0
      if(abs(beta).lt.1.0d-16)  ifdipole = 0

      allocate(sources(3,ns),srctmp(3,ntarg))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(ns))
      allocate(pot_aux(ntarg),grad_aux(3,ntarg))

c     oversample density
      call oversample_fun_surf(2_8,npatches,norders,ixyzs,iptype,
     1    npts,sigma,novers,ixyzso,ns,sigmaover)

c     set FMM source/target arrays
C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
        charges(i) = alpha*sigmaover(i)*whtsover(i)
        dipvec(1,i) = beta*sigmaover(i)*whtsover(i)*srcover(10,i)
        dipvec(2,i) = beta*sigmaover(i)*whtsover(i)*srcover(11,i)
        dipvec(3,i) = beta*sigmaover(i)*whtsover(i)*srcover(12,i)
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        srctmp(1,i) = targs(1,i)
        srctmp(2,i) = targs(2,i)
        srctmp(3,i) = targs(3,i)
        pot(i)        = 0
        pot_aux(i)    = 0
        grad_aux(1,i) = 0
        grad_aux(2,i) = 0
        grad_aux(3,i) = 0
      enddo
C$OMP END PARALLEL DO

c     call FMM: compute gradient of (alpha*S_k + beta*D_k)[sigma]
      ier = 0
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        call hfmm3d_t_cd_g(eps,zpars(1),ns,sources,charges,dipvec,
     1    ntarg,srctmp,pot_aux,grad_aux,ier)
      elseif(ifcharge.eq.1.and.ifdipole.eq.0) then
        call hfmm3d_t_c_g(eps,zpars(1),ns,sources,charges,
     1    ntarg,srctmp,pot_aux,grad_aux,ier)
      elseif(ifcharge.eq.0.and.ifdipole.eq.1) then
        call hfmm3d_t_d_g(eps,zpars(1),ns,sources,dipvec,
     1    ntarg,srctmp,pot_aux,grad_aux,ier)
      endif

c     dot gradient with target normal
C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        pot(i) = grad_aux(1,i)*targs(10,i) +
     1           grad_aux(2,i)*targs(11,i) +
     2           grad_aux(3,i)*targs(12,i)
      enddo
C$OMP END PARALLEL DO

      call get_fmm_thresh(int8_12,ns,srcover,int8_12,ntarg,srctmp,
     1  thresh)

c     add precomputed near quadrature
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,pottmp,npols,l)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            pot(i) = pot(i) + wnear(1,jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO

c     subtract near FMM contribution (add-subtract correction)
      call get_near_corr_max(ntarg,row_ptr,nnz,col_ind,npatches,
     1  ixyzso,nmax)
      allocate(srctmp2(3,nmax),ctmp2(nmax),dtmp2(3,nmax))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,pottmp,gradtmp)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            if(ifcharge.eq.1) ctmp2(nss)   = charges(l)
            if(ifdipole.eq.1) then
              dtmp2(1,nss) = dipvec(1,l)
              dtmp2(2,nss) = dipvec(2,l)
              dtmp2(3,nss) = dipvec(3,l)
            endif
          enddo
        enddo

        pottmp = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          call h3ddirectcdg(nd,zpars(1),srctmp2,ctmp2,dtmp2,
     1      nss,srctmp(1,i),ntarg0,pottmp,gradtmp,thresh)
        elseif(ifcharge.eq.1.and.ifdipole.eq.0) then
          call h3ddirectcg(nd,zpars(1),srctmp2,ctmp2,
     1      nss,srctmp(1,i),ntarg0,pottmp,gradtmp,thresh)
        elseif(ifcharge.eq.0.and.ifdipole.eq.1) then
          call h3ddirectdg(nd,zpars(1),srctmp2,dtmp2,
     1      nss,srctmp(1,i),ntarg0,pottmp,gradtmp,thresh)
        endif

        pot(i) = pot(i) - (gradtmp(1)*targs(10,i) +
     1                     gradtmp(2)*targs(11,i) +
     2                     gradtmp(3)*targs(12,i))
      enddo
C$OMP END PARALLEL DO

      return
      end
