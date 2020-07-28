      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:),targs(:,:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear1(:),wnear2(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:),rhs(:),sigma2(:),sigma3(:)
      real *8, allocatable :: errs(:)
      real *8 eps_gmres
      real *8 uvs_tmp(2)
      complex * 16 zpars(3),zpars_tmp(3),zpars_fp(3)
      complex *16 zk
      complex *16, allocatable :: pottmp(:)
      integer numit,niter


      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima,z1,z2

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


c
c       select geometry type
c       igeomtype = 1 => sphere
c 
      igeomtype = 1
      ipars(1) = 1
      npatches = 12*(4**ipars(1))

      zk = 1.11d0+ima*0.0d0
      zpars(1) = zk 
      zpars(2) = 1.0d0

      xyz_out(1) = 3.17d0
      xyz_out(2) = -0.03d0
      xyz_out(3) = 3.15d0

      xyz_in(1) = 0.17d0
      xyz_in(2) = 0.23d0
      xyz_in(3) = -0.11d0

      norder = 4 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 0


      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

      isout0 = .false.
      isout1 = .false.
      call test_exterior_pt(npatches,norder,npts,srcvals,srccoefs,
     1  wts,xyz_in,isout0)
      
      call test_exterior_pt(npatches,norder,npts,srcvals,srccoefs,
     1   wts,xyz_out,isout1)

       print *, isout0,isout1


      allocate(sigma(npts),rhs(npts),sigma2(npts),sigma3(npts))

      do i=1,npts
        call h3d_sprime(xyz_in,12,srcvals(1,i),0,dpars,1,zpars,0,
     1     ipars,rhs(i))
        sigma(i) = 0 
      enddo


      numit = 200
      ifinout = 0
      niter = 0
      allocate(errs(numit+1))

      eps = 0.51d-4

      eps_gmres = 0.5d-6

      ipars(1) = 2
      ipars(2) = 2

      ipt = 3
      jpt = 20

      ndd = 0
      ndz = 2
      ndi = 2

      call fker_h_neumanns(srcvals(1,ipt),12,srcvals(1,jpt),ndd,
     1   dpars,ndz,zpars,ndi,ipars,z1)
      
      zpars_tmp(1) = zpars(1)
      zpars_tmp(2) = zpars(1)*ima
      z2 = 0
      call h3d_dprime_diff(srcvals(1,ipt),12,srcvals(1,jpt),ndd,
     1  dpars,ndz,zpars_tmp,ndi,ipars,z2)

      z2 = z2/4/pi

      call prin2('z1=*',z1,2)
      call prin2('z2=*',z2,2)

      erra = abs(z1-z2)/abs(z1)
      call prin2('error =*',erra,1)

      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      
      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,cms,rads)
      
      allocate(ipatch_id(npts),uvs_targ(2,npts))
      allocate(targs(3,npts))
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
        targs(1:3,i) = srcvals(1:3,i)
      enddo

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1 ipatch_id,uvs_targ)

c
c    find near field
c
      iptype = 1
      call get_rfacs(norder,iptype,rfac,rfac0)
      do i=1,npatches 
        rad_near(i) = rads(i)*rfac
      enddo
      

      call findnearmem(cms,npatches,rad_near,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)

      nquad = iquad(nnz+1)-1

      goto 1111      

      allocate(wnear1(4*nquad),wnear2(4*nquad))

      rfac0 = 1.25d0
      iquadtype = 1
      call getnearquad_h_neumann(npatches,norders,
     1 ixyzs,iptype,npts,srccoefs,srcvals,12,npts,srcvals,
     2 ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3 iquad,rfac0,nquad,wnear1,wnear1(nquad+1),wnear1(2*nquad+1),
     4 wnear1(3*nquad+1))

      call getnearquad_helm_rpcomb_neu(npatches,norders,ixyzs,
     1 iptype,npts,srccoefs,srcvals,eps,zpars,iquadtype,nnz,
     2 row_ptr,col_ind,iquad,rfac0,nquad,wnear2)
      
      wnear2 = wnear2/4/pi
      
      erra1 = 0
      erra2 = 0
      erra3 = 0
      erra4 = 0
      ra = 0
      do i=1,nquad
        erra1 = erra1 + abs(wnear1(i)-wnear2(i))
        erra2 = erra2 + abs(wnear1(nquad+i)-wnear2(nquad+i))
        erra3 = erra3 + abs(wnear1(2*nquad+i)-wnear2(2*nquad+i))
        erra4 = erra4 + abs(wnear1(3*nquad+i)-wnear2(3*nquad+i))
        ra = ra + abs(wnear1(i))+abs(wnear1(nquad+i))+
     1     abs(wnear1(2*nquad+i))+abs(wnear1(3*nquad+i))
      enddo

      erra1 = erra1/ra
      erra2 = erra2/ra
      erra3 = erra3/ra
      erra4 = erra4/ra
      call prin2('erra1=*',erra1,1)
      call prin2('erra2=*',erra2,1)
      call prin2('erra3=*',erra3,1)
      call prin2('erra4=*',erra4,1)

 1111 continue      

      niter = 0

      call helm_rpcomb_neu_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,zpars,numit,rhs,eps_gmres,niter,errs,rres,
     2  sigma)
      sigma = sigma*4*pi

      ipatch_tmp = -1
      uvs_tmp(1) = 0
      uvs_tmp(2) = 0

      allocate(pottmp(npts))

      call lpcomp_helm_rpcomb_dir(npatches,norders,ixyzs,iptype,
     1 npts,srccoefs,srcvals,3,1,xyz_out,ipatch_tmp,uvs_tmp,
     2 eps,zpars,sigma,pot,pottmp)

      pottmp = pottmp/4/pi


      zpars_fp(1) = zk
      zpars_fp(2) = zpars(2)*4*pi
      ifinout = 1

      call h_neumann_solver(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1 srcvals,eps,zpars_fp,numit,ifinout,rhs,eps_gmres,niter,errs,rres,
     2 sigma2,sigma3)

      erra = 0
      ra = 0
      erra2 = 0
      ra2 = 0
      do i=1,npts
        erra = erra + abs(sigma2(i)-sigma(i))
        ra = ra + abs(sigma2(i))
        erra2 = erra2 + abs(sigma3(i)-pottmp(i))
        ra2 = ra2 + abs(sigma3(i))
      enddo
      erra = erra/ra
      erra2 = erra2/ra2

      call prin2('error in solve=*',erra,1)
      call prin2('error in second density=*',erra2,1)
      call prin2('sigma3=*',sigma3,12)
      call prin2('pottmp=*',pottmp,12)




      stop
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri, 
     1      triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         iptr3,iptr4, norder,
     2         'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end


      subroutine test_exterior_pt(npatches,norder,npts,srcvals,
     1   srccoefs,wts,xyzout,isout)
c
c
c  this subroutine tests whether the pt xyzin, is
c  in the exterior of a surface, and also estimates the error
c  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
c  centered at the interior point. Whether a point 
c  is in the interior or not is tested using Gauss' 
c  identity for the flux due to a point charge
c
c
c  input:
c    npatches - integer
c       number of patches
c    norder - integer
c       order of discretization
c    npts - integer
c       total number of discretization points on the surface
c    srccoefs - real *8 (9,npts)
c       koornwinder expansion coefficients of geometry info
c    xyzout -  real *8 (3)
c       point to be tested
c
c  output: 
c    isout - boolean
c      whether the target is in the interior or not
c

      implicit none
      integer npatches,norder,npts,npols
      real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer ipatch,j,i
      real *8 ra,ds
      logical isout

      done = 1
      pi = atan(done)*4

      npols = (norder+1)*(norder+2)/2


      zk = 0

      ra = 0



      do ipatch=1,npatches
        do j=1,npols
          i = (ipatch-1)*npols + j
          call h3d_sprime(xyzout,12,srcvals(1,i),0,dpars,1,zk,0,ipars,
     1       val)

          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-1) isout = .false.
      if(abs(ra).le.1.0d-1) isout = .true.

      return
      end

   




