
      implicit real *8 (a-h,o-z) 

      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      character *100 fname
      integer ipars(2)

      integer, allocatable :: norders(:), ixyzs(:), iptype(:)
      integer, allocatable :: ixys2d(:), iptype2d(:)
      real *8, allocatable :: srccoefs2d(:,:), srcvals2d(:,:)

      real *8, allocatable :: ts(:), ws(:), umat(:,:), vmat(:,:)

      real *8 abc(3), c0(3)
      integer nabc(3), nuv(2)
      real *8 radii(3)

      integer ipatch_id
      real *8 uvs_targ(2)
      real *8 pars(2)
      real *8, allocatable :: tchse(:)
      integer, allocatable :: nrts(:,:)

      external funcurve_oocyte_riemann

      call prini(6,13)


      done = 1
      pi = atan(done)*4

      dlam = 0.5d0

      abc(1) = 2.1d0
      abc(2) = 1.0d0
      abc(3) = 4.0d0

      nabc(1) = 5
      nabc(2) = 3
      nabc(3) = 10

      c0(1) = 0
      c0(2) = 5
      c0(3) = -3

      npatches = 0
      npts = 0

      norder = 6

      ifellip = 0
      ifsphere = 0
      ifstartorus = 0
      ifstell = 0
      ifoocyte = 0
      ifoocyte_chunks = 1

      if (ifellip.eq.1) then

        print *, "========================================="
        print *, "Testing ellipsoids"    
        print *, ""
        print *, "" 
        iptype0 = 1

        call get_ellipsoid_npat_mem(abc, nabc, c0, norder, iptype0, 
     1    npatches, npts)
        print *, "npatches tri=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

        call get_ellipsoid_npat(abc, nabc, c0, norder, iptype0, 
     1    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'ellip_tri.vtk','a')

        call surf_quadratic_msh_vtk_plot(npatches, norders, ixyzs, 
     1     iptype, npts, srccoefs, srcvals, 'ellip_tri_msh.vtk','a')

        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)

        npatches = 0
        npts = 0
        iptype0 = 11

        call get_ellipsoid_npat_mem(abc, nabc, c0, norder, iptype0, 
     1    npatches, npts) 
        print *, "npatches quad=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

        call get_ellipsoid_npat(abc, nabc, c0, norder, iptype0, 
     1    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'ellip_quad.vtk','a')

        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)
      endif

      if (ifsphere.eq.1) then

        print *, "========================================="
        print *, "Testing sphere"    
        print *, ""
        print *, "" 

        a = abc(1)
        na = nabc(1)

        iptype0 = 1
        call get_sphere_npat_mem(a, na, c0, norder, iptype0, 
     1    npatches, npts)
        print *, "npatches tri=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

        call get_sphere_npat(a, na, c0, norder, iptype0, 
     1    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'sphere_tri.vtk','a')

        call surf_quadratic_msh_vtk_plot(npatches, norders, ixyzs,
     1     iptype, npts, srccoefs, srcvals, 'sphere_tri_msh.vtk','a')


        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)

        npatches = 0
        npts = 0
        iptype0 = 12

        call get_sphere_npat_mem(a, na, c0, norder, iptype0, 
     1    npatches, npts) 
        print *, "npatches quad=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

        call get_sphere_npat(a, na, c0, norder, iptype0, 
     1    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'sphere_quad.vtk','a')

        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)
      endif

      if (ifstartorus.eq.1) then

        print *, "========================================="
        print *, "Testing star tourus"    
        print *, ""
        print *, ""

        radii(1) = 2
        radii(2) = 0.75d0
        radii(3) = 0.25d0

        abc(1) = 1
        abc(2) = 1
        abc(3) = 1
        
        nosc = 5

        nuv(1) = 10
        nuv(2) = 15

        iptype0 = 1

        call get_startorus_npat_mem(radii, nosc, abc, nuv, norder, 
     1    iptype0, npatches, npts)
        print *, "npatches tri=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

        call get_startorus_npat(radii, nosc, abc, nuv, norder, iptype0, 
     1    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'star_torus_tri.vtk','a')

        call surf_quadratic_msh_vtk_plot(npatches, norders, ixyzs,
     1     iptype, npts, srccoefs, srcvals, 'star_torus_tri_msh.vtk',
     2     'a')


        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)

        npatches = 0
        npts = 0
        iptype0 = 12

        call get_startorus_npat_mem(radii, nosc, abc, nuv, norder, 
     1    iptype0, npatches, npts)
        print *, "npatches quad=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

        call get_startorus_npat(radii, nosc, abc, nuv, norder, iptype0, 
     1    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'star_torus_quad.vtk','a')

        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)
      endif

      if (ifstell.eq.1) then

        print *, "========================================="
        print *, "Testing stellarator"    
        print *, ""
        print *, ""


        nuv(1) = 10
        nuv(2) = 30
        npatches = 0

        iptype0 = 1

        call get_stellarator_npat_mem(nuv, norder, iptype0, 
     1     npatches, npts)
        print *, "npatches tri=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

        call get_stellarator_npat(nuv, norder, iptype0, npatches, npts, 
     1     norders, ixyzs, iptype, srccoefs, srcvals)

        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'stellarator_tri.vtk','a')

        call surf_quadratic_msh_vtk_plot(npatches, norders, ixyzs,
     1     iptype, npts, srccoefs, srcvals, 'stellarator_tri_msh.vtk',
     2     'a')

        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)

        npatches = 0
        npts = 0
        iptype0 = 11

        call get_stellarator_npat_mem(nuv, norder, iptype0, 
     1     npatches, npts)
        print *, "npatches quad=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

        call get_stellarator_npat(nuv, norder, iptype0, npatches, npts, 
     1     norders, ixyzs, iptype, srccoefs, srcvals)

        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'stellarator_quad.vtk','a')

        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)
      endif


      if (ifoocyte.eq.1) then

        print *, "========================================="
        print *, "Testing oocyte"    
        print *, ""
        print *, ""
        
        a = 0.25d0
        b = 0.1d0
        
        pars(1) = a
        pars(2) = b

        nmid = 5
        iref = 0
        nch2d = 2*(iref+1) + nmid
        k = 16

        allocate(tchse(nch2d+1), nrts(2,nch2d))
        call get_oocyte3d_tchse(iref,nmid,nch2d,tchse)

        rmid = 0.5d0
        iort = 1

        nmid = 3

        nrts(1,1:nch2d) = 3
        nrts(2,1:nch2d) = 4

c        nrts(1:2,1) = 0
c        nrts(1:2,nch2d) = 0

        npatches = 0
        iptype0 = 1

        np = 2
        iort = 1

        call get_axissym_fcurve_npat_mem(nch2d, tchse,  
     1     funcurve_oocyte_riemann, np, pars, nrts, rmid, nmid, iort, 
     2     norder, iptype0, npatches, npts)
        print *, "npatches tri=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

        call get_axissym_fcurve_npat(nch2d, tchse,  
     1     funcurve_oocyte_riemann, np, pars, nrts, rmid, nmid, iort, 
     2     norder, iptype0, npatches, npts, norders, ixyzs, iptype, 
     3     srccoefs, srcvals)
        

        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'oocyte_tri.vtk','a')

        call surf_quadratic_msh_vtk_plot(npatches, norders, ixyzs,
     1     iptype, npts, srccoefs, srcvals, 'oocyte_tri_msh.vtk',
     2     'a')
        call vtk_scatter_plot_scalar(npts, srcvals(1:3,1:npts),
     1      srcvals(1,1:npts), 'oocyte_tri_scatter.vtk','a')
        call surf_vtk_plot_vec(npatches, norders, ixyzs, iptype, 
     1     npts, srccoefs, srcvals, srcvals(10:12,1:npts),
     2     'oocyte_tri_normals.vtk','a')

        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)

        npatches = 0
        npts = 0
        iptype0 = 11

        call get_axissym_fcurve_npat_mem(nch2d, tchse,  
     1     funcurve_oocyte_riemann, np, pars, nrts, rmid, nmid, iort, 
     2     norder, iptype0, npatches, npts)
        print *, "npatches quad=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

        call get_axissym_fcurve_npat(nch2d, tchse,  
     1     funcurve_oocyte_riemann, np, pars, nrts, rmid, nmid, iort, 
     2     norder, iptype0, npatches, npts, norders, ixyzs, iptype, 
     3     srccoefs, srcvals)
        
        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'oocyte_quad.vtk','a')

        call surf_quadratic_msh_vtk_plot(npatches, norders, ixyzs,
     1     iptype, npts, srccoefs, srcvals, 'oocyte_quad_msh.vtk',
     2     'a')
        call surf_vtk_plot_vec(npatches, norders, ixyzs, iptype, 
     1     npts, srccoefs, srcvals, srcvals(10:12,1:npts),
     2     'oocyte_quad_normals.vtk','a')
        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)
      endif

      if (ifoocyte_chunks.eq.1) then

        print *, "========================================="
        print *, "Testing oocyte chunks"     
        print *, ""
        print *, ""
        
        a = 0.25d0
        b = 0.1d0

        np = 2
        
        pars(1) = a
        pars(2) = b

        nmid = 20
        iref = 0
        nch2d = 2*(iref+1) + nmid
        k = 16

        allocate(tchse(nch2d+1), nrts(2,nch2d))
        call get_oocyte3d_tchse(iref,nmid,nch2d,tchse)

        npts2d = nch2d*k

        allocate(srcvals2d(8,npts2d), srccoefs2d(6,npts2d))
        allocate(ixys2d(nch2d+1), iptype2d(nch2d))

        k = 16
        itype = 2
        allocate(ts(k), ws(k), umat(k,k), vmat(k,k))
        call legeexps(itype, k, ts, umat, vmat, ws)

        do i=1,nch2d
          ixys2d(i) = (i-1)*k + 1
          iptype2d(i) = 1
        enddo

        ixys2d(nch2d+1) = npts2d+1
        do ich = 1,nch2d
          tchstart = tchse(ich)
          tchend = tchse(ich+1)
          
          hh = (tchend - tchstart)/2.0d0
          
          do j=1,k
            ipt = (ich-1)*k + j
            tt = tchstart + (ts(j)+1.0d0)/2.0d0*(tchend-tchstart)

            call funcurve_oocyte_riemann(tt,np,pars,x,y,dx,dy,dx2,dy2)
            srcvals2d(1,ipt) = x
            srcvals2d(2,ipt) = y
            srcvals2d(3,ipt) = dx*hh
            srcvals2d(4,ipt) = dy*hh
            srcvals2d(5,ipt) = dx2*hh**2
            srcvals2d(6,ipt) = dy2*hh**2
            ds = sqrt(dx**2 + dy**2)
            srcvals2d(7,ipt) = dy/ds
            srcvals2d(8,ipt) = -dx/ds
          enddo
          do j=1,k
            jpt = (ich-1)*k + j
            srccoefs2d(1:6,jpt) = 0
            do l=1,k
              lpt = (ich-1)*k + l
              srccoefs2d(1:6,jpt) = srccoefs2d(1:6,jpt) + 
     1            umat(j,l)*srcvals2d(1:6,lpt)
            enddo
          enddo
        enddo


        rmid = 0.5d0
        iort = 1

        nmid = 3

        nrts(1,1:nch2d) = 3
        nrts(2,1:nch2d) = 4

c        nrts(1:2,1) = 0
c        nrts(1:2,nch2d) = 0

        npatches = 0
        iptype0 = 1

        iort = 1

        call get_axissym_npat_mem(nch2d, npts2d, ixys2d, iptype2d,  
     1     srcvals2d, srccoefs2d, nrts, rmid, nmid, iort, 
     2     norder, iptype0, npatches, npts)
        print *, "npatches tri=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))


        call get_axissym_npat(nch2d, npts2d, ixys2d, iptype2d,   
     1     srcvals2d, srccoefs2d, nrts, rmid, nmid, iort, 
     2     norder, iptype0, npatches, npts, norders, ixyzs, iptype, 
     3     srccoefs, srcvals)
        

        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'oocyte_c_tri.vtk','a')

        call surf_quadratic_msh_vtk_plot(npatches, norders, ixyzs,
     1     iptype, npts, srccoefs, srcvals, 'oocyte_c_tri_msh.vtk',
     2     'a')
        call vtk_scatter_plot_scalar(npts, srcvals(1:3,1:npts),
     1      srcvals(1,1:npts), 'oocyte_c_tri_scatter.vtk','a')
        call surf_vtk_plot_vec(npatches, norders, ixyzs, iptype, 
     1     npts, srccoefs, srcvals, srcvals(10:12,1:npts),
     2     'oocyte_c_tri_normals.vtk','a')

        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)

        npatches = 0
        npts = 0
        iptype0 = 11

        call get_axissym_npat_mem(nch2d, npts2d, ixys2d, iptype2d,  
     1     srcvals2d, srccoefs2d, nrts, rmid, nmid, iort, 
     2     norder, iptype0, npatches, npts)
        print *, "npatches quad=", npatches

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))


        call get_axissym_npat(nch2d, npts2d, ixys2d, iptype2d,   
     1     srcvals2d, srccoefs2d, nrts, rmid, nmid, iort, 
     2     norder, iptype0, npatches, npts, norders, ixyzs, iptype, 
     3     srccoefs, srcvals)
        
        call plot_surface_info_all(dlam, npatches, norders, ixyzs, 
     1    iptype, npts, srccoefs, srcvals, 'oocyte_c_quad.vtk','a')

        call surf_quadratic_msh_vtk_plot(npatches, norders, ixyzs,
     1     iptype, npts, srccoefs, srcvals, 'oocyte_c_quad_msh.vtk',
     2     'a')
        call vtk_scatter_plot_scalar(npts, srcvals(1:3,1:npts),
     1      srcvals(1,1:npts), 'oocyte_c_quad_scatter.vtk','a')
        call surf_vtk_plot_vec(npatches, norders, ixyzs, iptype, 
     1     npts, srccoefs, srcvals, srcvals(10:12,1:npts),
     2     'oocyte_c_quad_normals.vtk','a')
        deallocate(srcvals, srccoefs, norders, ixyzs, iptype)
      endif

      stop
      end


   





      subroutine funcurve_oocyte_riemann(t,np,pars,x,y,dxdt, 
     1   dydt,d2xdt2,d2ydt2)
      implicit real *8 (a-h,o-z)
      real *8 t,pars(np),pi


      done = 1.0d0
      pi = atan(done)*4

      a = pars(1)
      b = pars(2)
      w = sqrt(1.0d0 + a**2 + b**2)

      tt = pi*t
      x = w/2*sin(tt) - a/2*sin(tt) - b/2/sqrt(2.0d0)*sin(2*tt)
      dxdt = pi*(w/2*cos(tt) - a/2*cos(tt) - b/sqrt(2.0d0)*cos(2*tt))
      d2xdt2 = -pi*pi*(w/2*sin(tt) - a/2*sin(tt) -  
     1    b*sqrt(2.0d0)*sin(2*tt))
      
      y = w*cos(tt) + a*cos(tt) + b/sqrt(2.0d0)*cos(2*tt)
      dydt = -pi*(w*sin(tt) + a*sin(tt) + b*sqrt(2.0d0)*sin(2*tt))
      d2ydt2 = -pi*pi*(w*cos(tt) + a*cos(tt) + 
     1   b*sqrt(2.0d0)*2*cos(2*tt))

      
      return
      end subroutine funcurve_oocyte_riemann
      


      subroutine get_oocyte3d_tchse(iref,nmid,nch2d,tchse)
      implicit real *8 (a-h,o-z)
      integer iref,nmid,nch2d
      real *8 tchse(nch2d+1)


      tchse(1) = 0.0d0
      hpan = 1.0d0/(nmid+2)
      rpan = hpan*(0.5d0)**(iref)
      tchse(2) = tchse(1) + rpan
      do i=2,iref+1
        tchse(i+1) = tchse(i) + rpan
        rpan = rpan*2
      enddo

      do i=iref+2,iref+1+nmid
        tchse(i+1) = tchse(i) + hpan
      enddo

      rpan = hpan/2
      do i=iref+2+nmid,nch2d-1
        tchse(i+1) = tchse(i) + rpan
        rpan = rpan/2
      enddo
      tchse(nch2d+1) = 1.0d0


      return
      end subroutine get_oocyte3d_tchse


