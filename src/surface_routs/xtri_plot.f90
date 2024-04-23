
!
!
! xtri plotting library (c) 2018 Mike O'Neil
! Contact: oneil@cims.nyu.edu
!
! This file contains plotting routines to accompany the xtri
! surface triangulation routines
!
!




subroutine xtri_vtk_surf(fname, ntri, xeval, &
    par1, par2, par3, par4, norder0,title)
  implicit real *8 (a-h,o-z)
  real *8 :: xtriainfo(3,1)
  character (len=*) :: title,fname
  external xeval

  real *8 :: uvs(2,1000), xyz(10), dxyzduv(3,20)
  real *8, allocatable :: sigma(:,:)

  !
  ! Dump out a .vtk ASCII file for reading in paraview or other
  !
  ! Merely plot the triangulated surface using the triangle evalation
  ! routine xeval, which should have the calling sequence:
  !
  !      xeval(itri, u, v, xyz, dxyzduv, par1, par2, par3, par4)
  !
  ! This calling sequence is consistent with legacy "patchmatc"
  ! codes. When using the xtri data structure, usually
  ! xeval=xtri_eval, and the calling sequence and parameters will be:
  !
  !      xtri_eval(itri, u, v, xyz, dxyzduv, korder, xtriainfo, &
  !          par3, par4)
  !
  ! Each triangle is oversampled kover times into flat triangles.
  ! Setting kover=1 results in 4 times as many triangles, kover=2
  ! results in 16 times as many triangles, etc.
  !
  ! input:
  !   fname - file name
  !   ntri - number of triangles
  !   xeval - subroutine evaluating points on the surface
  !   par1, par2, par3, par4 - parameters for xeval
  !   norder0 - order of discretization 
  !   title - title of the plot, example: 'plot title'
  !
  !

  vmax = -1000
  vmin = 1000

  itype=1
  norder = 1
  npols = (norder+1)*(norder+2)/2
  call get_vioreanu_nodes(norder,npols,uvs)


  allocate(sigma(npols,ntri))

  !
  ! evaluate z height
  !
  do itri = 1,ntri
    do j = 1,npols
      call xeval(itri, uvs(1,j), uvs(2,j), xyz, dxyzduv, &
          par1, par2, par3, par4)
      sigma(j,itri) = xyz(3)
      if (xyz(3) .gt. vmax) vmax = xyz(3)
      if (xyz(3) .lt. vmin) vmin = xyz(3)
    end do
  end do

  !
  ! now call the regular plotting routine...
  !
  call xtri_vtk_plot(fname, ntri, xeval, par1, par2, &
      par3, par4, norder, norder0,sigma, title)

  return
end subroutine xtri_vtk_surf






subroutine xtri_vtk_plot(fname, ntri, xeval, par1, par2, &
    par3, par4, norder, norder0, sigma, title)
  implicit real *8 (a-h,o-z)
  real *8 :: sigma(*)
  character (len=*) :: title,fname
  external xeval

  real *8, allocatable :: uvs(:,:),umatr(:,:)

  real *8, allocatable :: sigmaout(:,:), xtriout(:,:,:)

  !
  ! Plot the real-valued function sigma on eaah of the triangles using
  ! the triangle evalation routine xeval, which should have the
  ! calling sequence:
  !
  !      xeval(itri, u, v, xyz, dxyzduv, par1, par2, par3, par4)
  !
  ! This calling sequence is consistent with legacy "patchmatc"
  ! codes. When using the xtri data structure, usually
  ! xeval=xtri_eval, and the calling sequence and parameters will be:
  !
  !      xtri_eval(itri, u, v, xyz, dxyzduv, korder, xtriainfo, &
  !          par3, par4)
  !
  ! Each triangle (and function sigma) is oversampled kover times into
  ! flat triangles.  Setting kover=1 results in 4 times as many
  ! triangles, kover=2 results in 16 times as many triangles, etc.
  ! contruct a script to plot the function sigma on the triangles
  !
  !
 npols0 = (norder0+1)*(norder0+2)/2

!
!   set kover based on order
!
  kover = 1
  if(npols0.gt.4**0) kover = 1
  if(npols0.gt.4**1) kover = 2
  if(npols0.gt.4**2) kover = 3
  if(npols0.gt.4**3) kover = 4
  if(npols0.gt.4**4) kover = 5

   
  len = ntri*(4**kover)
  allocate(sigmaout(3,len))
  allocate(xtriout(3,3,len))

  npols = (norder+1)*(norder+2)/2
  allocate(uvs(2,npols),umatr(npols,npols))

  itype=1
  call get_vioreanu_nodes(norder,npols,uvs)
  call koorn_vals2coefs(norder,npols,uvs,umatr)

  call xtri_flatten(ntri, xeval, par1, par2, &
    par3, par4, norder, sigma, npols, umatr, &
    kover, ntriout, xtriout, sigmaout)
  

  m = 1
  call xtri_vtk_flat_scalars(fname, ntriout, xtriout, m, sigmaout, &
      title)

  return
end subroutine xtri_vtk_plot






subroutine xtri_vtk_flat_scalars(fname, ntri, xtri1s, m, sigma, title)
  implicit real *8 (a-h,o-z)
  real *8 :: xtri1s(3,3,ntri), sigma(m,3,ntri)
  character(len=*) :: title,fname
  character(len=25) :: fmt2

  !
  ! This routine plots a sequence of FLAT TRIANGLES with surface
  ! color vals.
  !
  ! Input:
  !   iw - plot number, controls the filenames
  !   ntri - number of flat triangles
  !   xtri1s - full triangle information
  !   sigma - function to plot, tabulated at corners of triangles
  !
  ! Output:
  !   a vtk file plotIW.vtk that can be plotted in paraview
  !
  !


  !
  ! write the vtk plotting script
  !
  iunit1 = 877
  open(unit = iunit1, file=trim(fname), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i8,a)') "POINTS ", ntri*3, " float"

  fmt2 = "(E11.5,2X,E11.5,2X,E11.5)"
  do i = 1,ntri
    write(iunit1,fmt2) xtri1s(1,1,i), xtri1s(2,1,i), xtri1s(3,1,i)
    write(iunit1,fmt2) xtri1s(1,2,i), xtri1s(2,2,i), xtri1s(3,2,i)
    write(iunit1,fmt2) xtri1s(1,3,i), xtri1s(2,3,i), xtri1s(3,3,i)
  end do


  write(iunit1,'(a,i8,i8)') "CELLS ", ntri, ntri*4

  do i = 1,ntri
    i1 = 3*(i-1) + 1
    write(iunit1,'(a,i8,i8,i8)') "3 ", i1-1, i1, i1+1
  end do

  write(iunit1,'(a,i8)') "CELL_TYPES ", ntri
  do i = 1,ntri
    write(iunit1,'(a)') "5"
  end do

  write(iunit1,'(a)') ""
  write(iunit1,'(a,i8)') "POINT_DATA ", ntri*3
  write(iunit1,'(a,i4)') "SCALARS scalar_function float ", m
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntri
    do j = 1,3
      do k = 1,m
        write(iunit1,'(E11.5)') sigma(k,j,i)
      end do
    end do
  end do

  close(iunit1)

  return
end subroutine xtri_vtk_flat_scalars





subroutine xtri_flatten(ntri, xeval, par1, par2, &
    par3, par4, norder, sigma, npols, umatr, &
    kover, ntriout, xtriout, sigmaout)
  implicit real *8 (a-h,o-z)
  real *8 :: sigma(*), xtriout(3,3,1)
  real *8 :: sigmaout(3,*), umatr(npols,npols)
  external xeval

  real *8 :: uv1(2), uv2(2), uv3(2), xyz(3), dxyzduv(3,20)
  real *8 :: coefs(10000), pols(10000)
  real *8, allocatable :: uvs(:,:,:), centers(:,:)

  !
  ! Oversample each curvilinear triangle to FLAT ones, korder times,
  ! and evaluate sigma at the corners of each triangle - this routine
  ! is mainly for plotting purposes.
  !
  ! Input:
  !   ntri -
  !   xeval - subroutine evaluating the triangles, must have the
  !       calling sequence:
  !
  !          xeval(itri, u, v, xyz, dxyzduv, par1, par2, par3, par4)
  !   norder - order of sigma on each triangle
  !   npols - number of point at which sigma is sampled on each
  !       triangle
  !   umatr - matrix converting values of sigma to expansion coefs
  !   kover - number of times to oversample each triangle
  !
  ! Output:
  !   xtriout - output flat triangles
  !   sigmaout - values of sigma at the corners of each triangle in
  !       xtriout
  !
  !

  neach = 4**kover
  allocate(uvs(2,3,neach))
  allocate(centers(2,neach))

  k1 = 1

  kout = 3
  uvs(1,1,1) = 0
  uvs(2,1,1) = 0
  uvs(1,2,1) = 1
  uvs(2,2,1) = 0
  uvs(1,3,1) = 0
  uvs(2,3,1) = 1

  !
  ! if necessary, recursively subdivide the triangle - first construct
  ! all the uv points
  !
  if (kover .gt. 0) then

    ntot = 1
    do i = 1,kover

      ltot = ntot

      do j = 1,ltot
        uv1(1) = uvs(1,1,j)
        uv1(2) = uvs(2,1,j)
        uv2(1) = uvs(1,2,j)
        uv2(2) = uvs(2,2,j)
        uv3(1) = uvs(1,3,j)
        uv3(2) = uvs(2,3,j)

        uvs(1,1,j) = uv1(1)
        uvs(2,1,j) = uv1(2)
        uvs(1,2,j) = (uv1(1) + uv2(1))/2
        uvs(2,2,j) = (uv1(2) + uv2(2))/2
        uvs(1,3,j) = (uv1(1) + uv3(1))/2
        uvs(2,3,j) = (uv1(2) + uv3(2))/2

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv1(1) + uv2(1))/2
        uvs(2,1,ntot) = (uv1(2) + uv2(2))/2
        uvs(1,2,ntot) = uv2(1)
        uvs(2,2,ntot) = uv2(2)
        uvs(1,3,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,3,ntot) = (uv2(2) + uv3(2))/2

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv1(1) + uv3(1))/2
        uvs(2,1,ntot) = (uv1(2) + uv3(2))/2
        uvs(1,2,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,2,ntot) = (uv2(2) + uv3(2))/2
        uvs(1,3,ntot) = uv3(1)
        uvs(2,3,ntot) = uv3(2)

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,1,ntot) = (uv2(2) + uv3(2))/2
        uvs(1,2,ntot) = (uv1(1) + uv3(1))/2
        uvs(2,2,ntot) = (uv1(2) + uv3(2))/2
        uvs(1,3,ntot) = (uv1(1) + uv2(1))/2
        uvs(2,3,ntot) = (uv1(2) + uv2(2))/2

      end do
    end do

  end if


  !
  ! calculate the centers
  !
  do i = 1,neach
    x = (uvs(1,1,i) + uvs(1,2,i) + uvs(1,3,i))/3
    y = (uvs(2,1,i) + uvs(2,2,i) + uvs(2,3,i))/3
    centers(1,i) = x
    centers(2,i) = y
  end do


  !
  ! evaluate the triangles and the interpolated sigma
  !
  ntriout = 0
  do itri = 1,ntri
    do j = 1,neach
      ntriout = ntriout + 1
      do l = 1,3
        u = uvs(1,l,j)
        v = uvs(2,l,j)
        call xeval(itri, u, v, xyz, dxyzduv, par1, par2, &
            par3, par4)
        xtriout(1,l,ntriout) = xyz(1)
        xtriout(2,l,ntriout) = xyz(2)
        xtriout(3,l,ntriout) = xyz(3)
      end do
    end do
  end do


  !call prinf('npols = *', npols, 1)

  ijk = 0
  !call prinf('ntri = *', ntri, 1)

  do itri = 1,ntri
    !call prinf('itri = *', itri, 1)
    ind = (itri-1)*npols + 1
    call dmatvec(npols, npols, umatr, sigma(ind), coefs)
    !call prin2('coefs = *', coefs, npols)

    do j = 1,neach
      ijk = ijk + 1

      do k = 1,3
        !u = centers(1,j)
        !v = centers(2,j)
        u = uvs(1,k,j)
        v = uvs(2,k,j)
        call koorn_pols(uvs(1,k,j),norder,npols,pols)
        !call prin2('pols = *', pols, npols)

        sigmaout(k,ijk) = 0
        do l = 1,npols
          sigmaout(k,ijk) = sigmaout(k,ijk) + coefs(l)*pols(l)
        end do
      end do

    end do
  end do

  return
end subroutine xtri_flatten


