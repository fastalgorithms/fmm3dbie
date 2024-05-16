      implicit real *8 (a-h,o-z)

      character *100 fname
      integer ipars(2),igeomtype,npatches,norder,npols
      integer ntargs,idx,ndim,maxiter,nerrlabels
      real *8 radius,pi,thet,phi,rad,r2min,r2temp,t1,t2,tolnear

      integer, allocatable :: lpatchidxvec(:),itvec(:)
      integer, allocatable :: newtonflags(:),flags(:),flagstrue(:)

      
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)
      real *8, allocatable :: targs(:,:)
      real *8, allocatable :: surftargstrue(:,:),surftargs(:,:)      
      real *8, allocatable :: sxyz(:,:),dists(:),uvslocvec(:,:)

      call prini(6,13)            
      
      done = 1
      pi = atan(done)*4

      ntargs = 100000
      ndim = 3
c     hard coded radius for now
      radius = 1.0d0 
      maxiter = 100
      tolnear = 1.0d-12
c
c     select geometry type, only one option atm
c     igeomtype = 1 => sphere
      
      igeomtype = 1
      
      ipars(1) = 2
      npatches = 12*(4**ipars(1))

      norder = 8 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))

      ifplot = 1
      fname = 'surface-geometry.vtk'

      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(targs(ndim,ntargs),flagstrue(ntargs))
      allocate(surftargstrue(ndim,ntargs),surftargs(ndim,ntargs))
      allocate(lpatchidxvec(ntargs),itvec(ntargs),dists(ntargs))
c     initialize targets
      do i=1,ntargs
         thet = hkrand(0)*pi
         phi = hkrand(0)*2.0d0*pi
         rad = 0.8 + 0.4d0*hkrand(0)
         targs(1,i) = rad*radius*dsin(thet)*dcos(phi)
         targs(2,i) = rad*radius*dsin(thet)*dsin(phi)
         targs(3,i) = rad*radius*dcos(thet)
         surftargstrue(1,i) = radius*dsin(thet)*dcos(phi)
         surftargstrue(2,i) = radius*dsin(thet)*dsin(phi)
         surftargstrue(3,i) = radius*dcos(thet)
c     outside
         if((targs(1,i)**2+targs(2,i)**2+targs(3,i)**2).gt.radius) then
            flagstrue(i) = 1
         endif
c     inside
         if((targs(1,i)**2+targs(2,i)**2+targs(3,i)**2).lt.radius) then
            flagstrue(i) = 2
         endif
      enddo


c     find closest surface discretization point
      
      do i=1,ntargs
         r2min = 1.0e6
         do j=1,npts
            rtemp = (targs(1,i)-srcvals(1,j))**2
     1           + (targs(2,i)-srcvals(2,j))**2
     2           + (targs(3,i)-srcvals(3,j))**2
            
            if(rtemp.lt.r2min) then
                        lpatchidxvec(i) = mod(j+(npols-1),npols) + 1
c     Find the starting index for the patch
                        itvec(i) = j-lpatchidxvec(i)+1
                        r2min = rtemp
            endif
         enddo
         
      enddo


      allocate(sxyz(3,ntargs),uvslocvec(2,ntargs),newtonflags(ntargs))
      allocate(flags(ntargs))

      print *," Find closest surface point"      
      call cpu_time(t1)
      call findnearsrfcpnt(ntargs,targs,
     1     nptssurf,norder,npols,maxiter,srcvals,
     2     srccoefs,tolnear,lpatchidxvec,itvec,sxyz,
     3     uvslocvec,dists,newtonflags)

      call inoutsurf3d(npatches,norder,npols,nptssurf,
     1     srccoefs,ntargs,targs,sxyz,itvec,
     2     uvslocvec,flags)

      call cpu_time(t2)
      
      call prin2('speed in points per second=*',(ntargs+0.0d0)/
     1     (t2-t1),1)

      nerrlabels = 0
      do i=1,ntargs
         if(flags(i).ne.flagstrue(i)) then
            nerrlabels = nerrlabels + 1
            if (flagstrue(i).eq.1) then
               print *," Correct label: ", outside
            endif
            if (flagstrue(i).eq.2) then
               print *," Correct label: ", inside
            endif
            if (flags(i).eq.1) then
               print *," Estimated label: ", outside
            endif
            if (flags(i).eq.2) then
               print *," Estimated label: ", inside
            endif
            
         endif
      enddo

      print *," Number of incorrectly labeled points: ", nerrlabels
      
      print *,"DONE"
      stop
      end



      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1     srcvals,srccoefs,ifplot,fname)
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
     1        triaskel, isides)

         xtri_geometry => xtri_sphere_eval
         ptr1 => triaskel(1,1,1)
         ptr2 => p2(1)
         ptr3 => p3(1)
         ptr4 => p4(1)


         if(ifplot.eq.1) then
            call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1           iptr3,iptr4, norder,
     2           'Triangulated surface of the sphere')
         endif


         call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,
     1        npols,uvs,umatr,srcvals,srccoefs)
      endif

      
     
      return  
      end
