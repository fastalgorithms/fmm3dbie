!
!  FIX this file for quad patches
!
!
!      surf_quadratic_msh_vtk_plot - generate a vtk file to
!        plot a quadratic mesh corresponding to discretization
!      
!      surf_flat_msh_vtk_plot - generate a vtk file to plot
!        a flat mesh corresponding to discretization
!
!      surf_vtk_plot - generate a vtk file to plot the surface, 
!         scalar plotted is z coordinate
!
!      surf_vtk_plot_scalar - generate a vtk file to plot the surface,
!         along with prescribed scalar
!
!      vtk_write_plane - write a structured grid of data defined on
!       a plane
!
!      vtk_write_plane_vec - write a structured grid of data defined on
!       a plane
!
!      vtk_scatter_plot_scalar - write a scatter plot of points with function
!      values
!
!
subroutine surf_quadratic_msh_vtk_plot(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,fname,title)
!
! This subroutine writes a vtk file to plot quadratic mesh 
! corresponding to discretization
!
! Currently only supports triangulation
!
!  Input arguments:
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location in srccoefs,srcvals array where data for
!        patch i begins
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch with RV nodes
!    - npts: integere
!        number of points in discretization
!    - srccoefs: real *8 (9,npts)
!        xyz,dxyz/du,dxyz/dv koornwinder expansion coeffs for all patches
!    - srcvals: real *8 (12,npts)
!        xyz,dxyz/du,dxyz/dv,normals patch info 
!    - fname: character (len=*)
!        file name where vtk output should be written
!
  implicit none
  integer, intent(in) :: npatches,norders(npatches),ixyzs(npatches+1)
  integer, intent(in) :: iptype(npatches),npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  character (len=*), intent(in) :: fname,title

  real *8 uvs(2,6)
  real *8, allocatable :: xyzs(:,:),pols(:)

  integer i,j,k,l,n0,npuv,ipatch,ipt,i1,m,norder,npols,iunit1

  uvs(1,1) = 0.0d0
  uvs(2,1) = 0.0d0
  
  uvs(1,2) = 1.0d0
  uvs(2,2) = 0.0d0

  uvs(1,3) = 0.0d0
  uvs(2,3) = 1.0d0

  uvs(1,4) = 0.5d0
  uvs(2,4) = 0.0d0
  
  uvs(1,5) = 0.5d0
  uvs(2,5) = 0.5d0

  uvs(1,6) = 0.0d0
  uvs(2,6) = 0.5d0

  npuv = 6

  n0 = npatches*npuv
  allocate(xyzs(3,n0))

  do ipatch=1,npatches
    npols = ixyzs(ipatch+1)-ixyzs(ipatch)
    norder = norders(ipatch)
    allocate(pols(npols))
    do j=1,npuv
      call koorn_pols(uvs(1,j),norder,npols,pols)
      ipt = (ipatch-1)*npuv + j
      do m=1,3
        xyzs(m,ipt) = 0
      enddo

      do l=1,npols
        do m=1,3
          xyzs(m,ipt) = xyzs(m,ipt) + &
             pols(l)*srccoefs(m,ixyzs(ipatch)+l-1)
        enddo
      enddo
    enddo
    deallocate(pols)
  enddo
  
  iunit1 = 877
  open(unit = iunit1, file=trim(fname), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i9,a)') "POINTS ", n0, " float"

  do i = 1,n0
    write(iunit1,"(E11.5,2(2x,e11.5))") xyzs(1,i), xyzs(2,i), xyzs(3,i)
  end do

  write(iunit1,'(a,i9,i9)') "CELLS ", npatches, npatches*7

  do ipatch=1,npatches
    i1 = 6*(ipatch-1) 
    write(iunit1,'(a,i9,i9,i9,i9,i9,i9)') "6 ", i1, i1+1,i1+2, &
      i1+3,i1+4,i1+5
  enddo

  write(iunit1,'(a,i9)') "CELL_TYPES ", npatches
  do ipatch = 1,npatches
    if(iptype(ipatch).eq.1) then
      write(iunit1,'(a)') "22"
    endif
  end do

  write(iunit1,'(a)') ""
  write(iunit1,'(a,i9)') "CELL_DATA ", npatches
  write(iunit1,'(a)') "SCALARS material int"
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,npatches
    write(iunit1,'(a)') "1"
  end do

  close(iunit1)

end subroutine surf_quadratic_msh_vtk_plot
!
!
!
!
!
subroutine surf_vtk_plot(npatches,norders,ixyzs,iptype,npts,srccoefs,&
   srcvals,fname,title)
  !
  ! This subroutine writes a vtk file to plot the surface as a
  ! collection of flat triangles, independent of the underlying order
  ! of the discretization
  !
  !  
  !f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs
  !f2py intent(in) srcvals,fname,title
  !
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),npts
  integer iptype(npatches)
  real *8 srccoefs(9,npts),srcvals(12,npts)
  real *8, allocatable :: sigma(:)
  character (len=*) fname,title

  integer i
   
  allocate(sigma(npts))
  !$OMP PARALLEL DO DEFAULT(SHARED)   
  do i=1,npts
     sigma(i) = srcvals(3,i)
  enddo
  !$OMP END PARALLEL DO

  call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,npts,  &
       srccoefs,srcvals,sigma,fname,title)

end subroutine surf_vtk_plot
!
!
!
!
!
!
subroutine surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sigma,fname,title)
!
!   This subroutine writes a vtk to plot the surface along
!   with a scalar. Currently only supports triangular patches
!
!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs
!f2py intent(in) srcvals,sigma,fname,title
!
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),npts
  integer iptype(npatches)
  real *8 srccoefs(9,npts),srcvals(12,npts),sigma(npts)
  real *8, allocatable :: sigma_coefs(:),pols(:),rtmp(:)
  character (len=*) fname,title

  real *8, allocatable :: xyzs(:,:),uvs(:,:,:),splot(:)
  real *8, allocatable :: uvs_quad(:,:,:)
  integer, allocatable :: kovers(:),nps(:),ipstart(:)

  integer i,j,k,l,ipatch,npout,kover,npols
  integer itrip,itric1,nb,nlmax,nuv,istart,iend,nd
  integer ilstart,itri,iunit1,m,ncell,ncsize,norder,nuvl,i1
  integer imul,nordermax,npolmax

  real *8 ra,erra

!
!  get the coefs of the density
!
   nd = 1
   allocate(sigma_coefs(npts))
   call surf_vals_to_coefs(nd,npatches,norders,ixyzs,iptype,npts, &
     sigma,sigma_coefs)
 
!
!   estimate kovers, nps 
!

  npout = 0
  allocate(kovers(npatches),nps(npatches),ipstart(npatches+1))
  do i=1,npatches
    npols = ixyzs(i+1)-ixyzs(i)
    kover = 0
    if(npols.gt.4**0) kover = 1
    if(npols.gt.4**1) kover = 2
    if(npols.gt.4**2) kover = 3
    if(npols.gt.4**3) kover = 4
    kovers(i) = kover
    nps(i) = 4**kover
    if(iptype(i).eq.1) nps(i) = nps(i)*3
    if(iptype(i).eq.11) nps(i) = nps(i)*4
    if(iptype(i).eq.12) nps(i) = nps(i)*4
    npout = npout + nps(i) 
  enddo


  ipstart(1) = 1
  nps(1) = nps(1) + 1
  call cumsum(npatches,nps,ipstart(2))
  nps(1) = nps(1) - 1

  allocate(xyzs(3,npout),splot(npout))

!
!   get uvs of all patches of type = 1
!
  
  nlmax = 5
  nuv = (4**(nlmax+1)-1)/3
  allocate(uvs(2,3,nuv))
  allocate(uvs_quad(2,4,nuv))


  uvs(1,1,1) = 0
  uvs(2,1,1) = 0
  uvs(1,2,1) = 1
  uvs(2,2,1) = 0
  uvs(1,3,1) = 0
  uvs(2,3,1) = 1

  uvs_quad(1,1,1) = -1
  uvs_quad(2,1,1) = -1
  uvs_quad(1,2,1) = 1
  uvs_quad(2,2,1) = -1
  uvs_quad(1,3,1) = -1
  uvs_quad(2,3,1) = 1

  uvs_quad(1,4,1) = 1
  uvs_quad(2,4,1) = 1

  do i=0,nlmax-1
    istart = (4**(i)-1)/3+1
    nb = 4**i
    iend = istart + nb-1
    do itrip = istart,iend
      itric1 = (itrip-istart)*4 + iend
      call gettrichildren(uvs(1,1,itrip),uvs(1,1,itric1+1), &
       uvs(1,1,itric1+2),uvs(1,1,itric1+3),uvs(1,1,itric1+4))   
      call getquadchildren(uvs_quad(1,1,itrip),uvs_quad(1,1,itric1+1), &
       uvs_quad(1,1,itric1+2),uvs_quad(1,1,itric1+3),uvs_quad(1,1,itric1+4))  
      
      do j=1,4
        uvs_quad(1,4,itric1+j) = uvs_quad(1,2,itric1+j) + & 
          uvs_quad(1,3,itric1+j) - uvs_quad(1,1,itric1+j)
        uvs_quad(2,4,itric1+j) = uvs_quad(2,2,itric1+j) + &
          uvs_quad(2,3,itric1+j) - uvs_quad(2,1,itric1+j)
      enddo
    enddo
  enddo


  nordermax = maxval(norders)
  npolmax = (nordermax+1)**2
  allocate(pols(npolmax))
  do ipatch=1,npatches
    istart = ipstart(ipatch)
    npols = ixyzs(ipatch+1)-ixyzs(ipatch)
    norder = norders(ipatch)

    nuvl = ipstart(ipatch+1)-ipstart(ipatch)
    ilstart = (4**(kovers(ipatch))-1)/3+1
    nb = 4**(kovers(ipatch))

    if(iptype(ipatch).eq.1) imul = 3
    if(iptype(ipatch).eq.11) imul = 4
    if(iptype(ipatch).eq.12) imul = 4
      
    do i=1,nb
      itri = i+ilstart-1
      do j=1,imul
        if(iptype(ipatch).eq.1) then
          call get_basis_pols(uvs(1,j,itri),norder,npols, &
            iptype(ipatch),pols)
        else if (iptype(ipatch).eq.11.or.iptype(ipatch).eq.12) then
          call get_basis_pols(uvs_quad(1,j,itri),norder,npols, &
            iptype(ipatch),pols)
        endif
         
          
        do m=1,3
          xyzs(m,istart+imul*(i-1)+j-1) = 0
        enddo
        splot(istart+imul*(i-1)+j-1) = 0

        do l=1,npols
          do m=1,3
            xyzs(m,istart+imul*(i-1)+j-1) = & 
              xyzs(m,istart+imul*(i-1)+j-1) + &
              pols(l)*srccoefs(m,ixyzs(ipatch)+l-1)
          enddo
          splot(istart+imul*(i-1)+j-1) = &
           splot(istart+imul*(i-1)+j-1)+ &
           pols(l)*sigma_coefs(ixyzs(ipatch)+l-1)
        enddo
      enddo
    enddo


  enddo
  
  iunit1 = 877
  open(unit = iunit1, file=trim(fname), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i9,a)') "POINTS ", npout, " float"

  do i = 1,npout
    write(iunit1,"(E11.5,2(2x,e11.5))") xyzs(1,i), xyzs(2,i), xyzs(3,i)
  end do

  ncell = 0
  ncsize = 0
  do i=1,npatches
    if(iptype(i).eq.1) ncell = ncell + 4**kovers(i)
    if(iptype(i).eq.11) ncell = ncell + 5**kovers(i)
    if(iptype(i).eq.12) ncell = ncell + 5**kovers(i)
    ncsize = ncsize + ncell 
  enddo

  write(iunit1,'(a,i9,i9)') "CELLS ", ncell, ncsize

  do ipatch=1,npatches
    nb = 4**kovers(ipatch)
    istart = ipstart(ipatch)
    if(iptype(ipatch).eq.1) imul = 3
    if(iptype(ipatch).eq.11) imul = 4
    if(iptype(ipatch).eq.12) imul = 4

    do i = 1,nb
        i1 = istart + imul*(i-1) 
        if(iptype(ipatch).eq.1) &
          write(iunit1,'(a,i9,i9,i9)') "3 ", i1-1, i1, i1+1
        if(iptype(ipatch).eq.11.or.iptype(ipatch).eq.12) &
          write(iunit1,'(a,i9,i9,i9,i9)') "3 ", i1-1, i1, i1+1,i1+2
    enddo


   
  end do

  write(iunit1,'(a,i9)') "CELL_TYPES ", ncell
  do ipatch = 1,npatches
    nb = 4**kovers(ipatch)
    do i=1,nb
      if(iptype(ipatch).eq.1) write(iunit1,'(a)') "5"
      if(iptype(ipatch).eq.11.or.iptype(ipatch).eq.12) write(iunit1,'(a)') "9"
    enddo
  end do

  write(iunit1,'(a)') ""
  write(iunit1,'(a,i9)') "POINT_DATA ", npout
  write(iunit1,'(a,i4)') "SCALARS scalar_function float ", 1
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,npout
    write(iunit1,'(E11.5)') splot(i)
  end do

  close(iunit1)



end subroutine surf_vtk_plot_scalar
!
!
!
!
!
!
!
!
subroutine surf_vtk_plot_scalar_many(nd,npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sigma,nsc,scalar_names,fname,title)
!
!   This subroutine writes a vtk to plot the surface along
!   with many scalar. Currently only supports triangular patches
!
!
  implicit none
  integer nd,nsc
  integer npatches,norders(npatches),ixyzs(npatches+1),npts
  integer iptype(npatches)
  real *8 srccoefs(9,npts),srcvals(12,npts),sigma(nd,npts)
  real *8, allocatable :: sigma_coefs(:,:),pols(:),rtmp(:)
  character (len=*) fname,title
  character (len=nsc), dimension (nd) :: scalar_names

  real *8, allocatable :: xyzs(:,:),uvs(:,:,:),splot(:,:)
  integer, allocatable :: kovers(:),nps(:),ipstart(:)

  integer i,j,k,l,ipatch,npout,kover,npols,idim
  integer itrip,itric1,nb,nlmax,nuv,istart,iend
  integer ilstart,itri,iunit1,m,ncell,ncsize,norder,nuvl,i1

  real *8 ra,erra

!
!  get the coefs of the density
!
   allocate(sigma_coefs(nd,npts))
   call surf_vals_to_coefs(nd,npatches,norders,ixyzs,iptype,npts, &
     sigma,sigma_coefs)
 
!
!   estimate kovers, nps 
!

  npout = 0
  allocate(kovers(npatches),nps(npatches),ipstart(npatches+1))
  do i=1,npatches
    npols = ixyzs(i+1)-ixyzs(i)
    kover = 0
    if(npols.gt.4**0) kover = 1
    if(npols.gt.4**1) kover = 2
    if(npols.gt.4**2) kover = 3
    if(npols.gt.4**3) kover = 4
    kovers(i) = kover
    nps(i) = 4**kover
    if(iptype(i).eq.1) nps(i) = nps(i)*3
    npout = npout + nps(i) 
  enddo


  ipstart(1) = 1
  nps(1) = nps(1) + 1
  call cumsum(npatches,nps,ipstart(2))
  nps(1) = nps(1) - 1

  allocate(xyzs(3,npout),splot(nd,npout))

!
!   get uvs of all patches of type = 1
!
  
  nlmax = 5
  nuv = (4**(nlmax+1)-1)/3
  allocate(uvs(2,3,nuv))

  uvs(1,1,1) = 0
  uvs(2,1,1) = 0
  uvs(1,2,1) = 1
  uvs(2,2,1) = 0
  uvs(1,3,1) = 0
  uvs(2,3,1) = 1

  do i=0,nlmax-1
    istart = (4**(i)-1)/3+1
    nb = 4**i
    iend = istart + nb-1
    do itrip = istart,iend
      itric1 = (itrip-istart)*4 + iend
      call gettrichildren(uvs(1,1,itrip),uvs(1,1,itric1+1), &
       uvs(1,1,itric1+2),uvs(1,1,itric1+3),uvs(1,1,itric1+4))   
    enddo
  enddo



  do ipatch=1,npatches
    istart = ipstart(ipatch)
    npols = ixyzs(ipatch+1)-ixyzs(ipatch)
    norder = norders(ipatch)
    allocate(pols(npols))
    if(iptype(ipatch).eq.1) then

      nuvl = ipstart(ipatch+1)-ipstart(ipatch)
      ilstart = (4**(kovers(ipatch))-1)/3+1
      nb = 4**(kovers(ipatch))
      
      do i=1,nb
        itri = i+ilstart-1
        do j=1,3
          call koorn_pols(uvs(1,j,itri),norder,npols,pols)
          
          do m=1,3
            xyzs(m,istart+3*(i-1)+j-1) = 0
          enddo
          do idim=1,nd
            splot(idim,istart+3*(i-1)+j-1) = 0
          enddo
         
          
          do l=1,npols
            do m=1,3
              xyzs(m,istart+3*(i-1)+j-1) = & 
                xyzs(m,istart+3*(i-1)+j-1) + &
                pols(l)*srccoefs(m,ixyzs(ipatch)+l-1)
            enddo
            do idim=1,nd
              splot(idim,istart+3*(i-1)+j-1) = &
               splot(idim,istart+3*(i-1)+j-1)+ &
               pols(l)*sigma_coefs(idim,ixyzs(ipatch)+l-1)
            enddo
          enddo
        enddo
      enddo
    endif
    deallocate(pols)
  enddo
  
  iunit1 = 877
  open(unit = iunit1, file=trim(fname), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i9,a)') "POINTS ", npout, " float"

  do i = 1,npout
    write(iunit1,"(E11.5,2(2x,e11.5))") xyzs(1,i), xyzs(2,i), xyzs(3,i)
  end do

  ncell = 0
  ncsize = 0
  do i=1,npatches
    ncell = ncell + 4**kovers(i)
    if(iptype(i).eq.1) ncsize = ncsize + 4*(4**kovers(i))
  enddo

  write(iunit1,'(a,i9,i9)') "CELLS ", ncell, ncsize

  do ipatch=1,npatches
    nb = 4**kovers(ipatch)
    if(iptype(ipatch).eq.1) then
      istart = ipstart(ipatch) 
      do i = 1,nb
        i1 = istart + 3*(i-1) 
        write(iunit1,'(a,i9,i9,i9)') "3 ", i1-1, i1, i1+1
      enddo
    endif
  end do

  write(iunit1,'(a,i9)') "CELL_TYPES ", ncell
  do ipatch = 1,npatches
    nb = 4**kovers(ipatch)
    if(iptype(ipatch).eq.1) then
      do i=1,nb
        write(iunit1,'(a)') "5"
      enddo
    endif
  end do

  write(iunit1,'(a)') ""
  write(iunit1,'(a,i9)') "POINT_DATA ", npout
  
  do idim=1,nd

    write(iunit1,'(a,a,a,i4)') "SCALARS ",trim(scalar_names(idim))," float ", 1
    write(iunit1,'(a)') "LOOKUP_TABLE default"
    do i = 1,npout
      write(iunit1,'(E11.5)') splot(idim,i)
    end do
  enddo

  close(iunit1)



end subroutine surf_vtk_plot_scalar_many
!
!
!
!
!
!
!
subroutine vtk_write_plane(ndims,ntarg,xyz,dxyz,f,title,fname)
!
!   writes data for a structured grid on a plane in vtk format
!
  implicit none
  integer ndims(3),ntarg,i,iunit1
  real *8 xyz(3),dxyz(3),f(ntarg)
  character (len=*) fname,title

  print *, "Here - in vtk_write_plane"
  

  iunit1 = 877
  open(unit = iunit1, file=trim(fname), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET STRUCTURED_POINTS"
  write(iunit1,'(a,i4,i4,i4)') "DIMENSIONS ", ndims(1),ndims(2),ndims(3)
  write(iunit1,'(a,e11.5,1x,e11.5,1x,e11.5)') "ORIGIN ", &
    xyz(1),xyz(2),xyz(3)
  write(iunit1,'(a,e11.5,1x,e11.5,1x,e11.5)') "SPACING ", &
    dxyz(1),dxyz(2),dxyz(3)
  write(iunit1,'(a)') ""
  write(iunit1,'(a,i9)') "POINT_DATA ",ntarg 
  write(iunit1,'(a,i4)') "SCALARS scalar_function float ", 1
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntarg
    write(iunit1,'(E11.5)') f(i)
  end do
  close(iunit1)



end subroutine vtk_write_plane
!
!
!
!
!
!
!
subroutine vtk_write_plane_vec(ndims,ntarg,xyz,dxyz,f,title,fname)
!
!   writes data for a structured grid on a plane in vtk format
!
  implicit none
  integer ndims(3),ntarg,i,iunit1
  real *8 xyz(3),dxyz(3),f(3,ntarg)
  character (len=*) fname,title

  print *, "Here - in vtk_write_plane"
  

  iunit1 = 877
  open(unit = iunit1, file=trim(fname), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET STRUCTURED_POINTS"
  write(iunit1,'(a,i6,i6,i6)') "DIMENSIONS ", ndims(1),ndims(2),ndims(3)
  write(iunit1,'(a,e11.5,1x,e11.5,1x,e11.5)') "ORIGIN ", &
    xyz(1),xyz(2),xyz(3)
  write(iunit1,'(a,e11.5,1x,e11.5,1x,e11.5)') "SPACING ", &
    dxyz(1),dxyz(2),dxyz(3)
  write(iunit1,'(a)') ""
  write(iunit1,'(a,i9)') "POINT_DATA ",ntarg 
  write(iunit1,'(a,i4)') "SCALARS vec_comps float ", 3
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntarg
    write(iunit1,'(3(E11.5,2x))') f(1,i),f(2,i),f(3,i)
  end do
  close(iunit1)



end subroutine vtk_write_plane_vec
!
!
!
!
!



subroutine surf_vtk_plot_vec(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sigma,fname,title)
!
!   This subroutine writes a vtk to plot the surface along
!   with a vector field. Currently only supports triangular patches
!
!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs
!f2py intent(in) srcvals,sigma,fname,title
!
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),npts
  integer iptype(npatches)
  real *8 srccoefs(9,npts),srcvals(12,npts),sigma(3,npts)
  real *8, allocatable :: sigma_coefs(:,:),pols(:),rtmp(:)
  character (len=*) fname,title

  real *8, allocatable :: xyzs(:,:),uvs(:,:,:),splot(:,:)
  integer, allocatable :: kovers(:),nps(:),ipstart(:)

  integer i,j,k,l,ipatch,npout,kover,npols
  integer itrip,itric1,nb,nlmax,nuv,istart,iend,nd
  integer ilstart,itri,iunit1,m,ncell,ncsize,norder,nuvl,i1
  integer idim

  real *8 ra,erra

!
!  get the coefs of the density
!
   nd = 3
   allocate(sigma_coefs(nd,npts))
   call surf_vals_to_coefs(nd,npatches,norders,ixyzs,iptype,npts, &
     sigma,sigma_coefs)


 
!
!   estimate kovers, nps 
!

  npout = 0
  allocate(kovers(npatches),nps(npatches),ipstart(npatches+1))
  do i=1,npatches
    npols = ixyzs(i+1)-ixyzs(i)
    kover = 0
    if(npols.gt.4**0) kover = 1
    if(npols.gt.4**1) kover = 2
    if(npols.gt.4**2) kover = 3
    if(npols.gt.4**3) kover = 4
    kovers(i) = kover
    nps(i) = 4**kover
    if(iptype(i).eq.1) nps(i) = nps(i)*3
    if(iptype(i).eq.11) nps(i) = nps(i)*4
    if(iptype(i).eq.12) nps(i) = nps(i)*4
    npout = npout + nps(i) 
  enddo


  ipstart(1) = 1
  nps(1) = nps(1) + 1
  call cumsum(npatches,nps,ipstart(2))
  nps(1) = nps(1) - 1

  allocate(xyzs(3,npout),splot(nd,npout))

!
!   get uvs of all patches of type = 1
!
  
  nlmax = 5
  nuv = (4**(nlmax+1)-1)/3
  allocate(uvs(2,3,nuv))

  uvs(1,1,1) = 0
  uvs(2,1,1) = 0
  uvs(1,2,1) = 1
  uvs(2,2,1) = 0
  uvs(1,3,1) = 0
  uvs(2,3,1) = 1

  do i=0,nlmax-1
    istart = (4**(i)-1)/3+1
    nb = 4**i
    iend = istart + nb-1
    do itrip = istart,iend
      itric1 = (itrip-istart)*4 + iend
      call gettrichildren(uvs(1,1,itrip),uvs(1,1,itric1+1), &
       uvs(1,1,itric1+2),uvs(1,1,itric1+3),uvs(1,1,itric1+4))   
    enddo
  enddo


  do ipatch=1,npatches
    istart = ipstart(ipatch)
    npols = ixyzs(ipatch+1)-ixyzs(ipatch)
    norder = norders(ipatch)
    allocate(pols(npols))
    if(iptype(ipatch).eq.1) then

      nuvl = ipstart(ipatch+1)-ipstart(ipatch)
      ilstart = (4**(kovers(ipatch))-1)/3+1
      nb = 4**(kovers(ipatch))
      do i=1,nb
        itri = i+ilstart-1
        do j=1,3
          call koorn_pols(uvs(1,j,itri),norder,npols,pols)
          
          do m=1,3
            xyzs(m,istart+3*(i-1)+j-1) = 0
          enddo
          do idim=1,nd
            splot(idim,istart+3*(i-1)+j-1) = 0
          enddo

          do l=1,npols
            do m=1,3
              xyzs(m,istart+3*(i-1)+j-1) = & 
                xyzs(m,istart+3*(i-1)+j-1) + &
                pols(l)*srccoefs(m,ixyzs(ipatch)+l-1)
            enddo
            do idim=1,nd
              splot(idim,istart+3*(i-1)+j-1) = &
               splot(idim,istart+3*(i-1)+j-1)+ &
               pols(l)*sigma_coefs(idim,ixyzs(ipatch)+l-1)
            enddo
          enddo
        enddo
      enddo
    endif
    deallocate(pols)
  enddo
  
  iunit1 = 877
  open(unit = iunit1, file=trim(fname))

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i9,a)') "POINTS ", npout, " float"

  do i = 1,npout
    write(iunit1,"(E11.5,2(2x,e11.5))") xyzs(1,i), xyzs(2,i), xyzs(3,i)
  end do

  ncell = 0
  ncsize = 0
  do i=1,npatches
    ncell = ncell + 4**kovers(i)
    if(iptype(i).eq.1) ncsize = ncsize + 4*(4**kovers(i))
    if(iptype(i).eq.11) ncsize = ncsize + 5*(4**kovers(i))
    if(iptype(i).eq.12) ncsize = ncsize + 5*(4**kovers(i))
  enddo

  write(iunit1,'(a,i9,i9)') "CELLS ", ncell, ncsize

  do ipatch=1,npatches
    nb = 4**kovers(ipatch)
    if(iptype(ipatch).eq.1) then
      istart = ipstart(ipatch) 
      do i = 1,nb
        i1 = istart + 3*(i-1) 
        write(iunit1,'(a,i9,i9,i9)') "3 ", i1-1, i1, i1+1
      enddo
    endif
  end do

  write(iunit1,'(a,i9)') "CELL_TYPES ", ncell
  do ipatch = 1,npatches
    nb = 4**kovers(ipatch)
    if(iptype(ipatch).eq.1) then
      do i=1,nb
        write(iunit1,'(a)') "5"
      enddo
    endif
  end do

  write(iunit1,'(a)') ""
  write(iunit1,'(a,i9)') "POINT_DATA ", npout
  write(iunit1,'(a,i4)') "SCALARS vec_comps float ", 3
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,npout
    write(iunit1,'(E11.5,2x,E11.5,2x,e11.5)') &
      splot(1,i),splot(2,i),splot(3,i)
  end do

  close(iunit1)



end subroutine surf_vtk_plot_vec
!
!


subroutine surf_vtk_plot_zvec(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sigma,fname,title)
!
!   This subroutine writes a vtk to plot the surface along
!   with a complex vector field. 
!   Currently only supports triangular patches
!
!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs
!f2py intent(in) srcvals,sigma,fname,title
!
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),npts
  integer iptype(npatches)
  real *8 srccoefs(9,npts),srcvals(12,npts)
  complex *16 sigma(3,npts)
  real *8, allocatable :: pols(:),rtmp(:)
  complex *16, allocatable :: sigma_coefs(:,:)
  character (len=*) fname,title

  real *8, allocatable :: xyzs(:,:),uvs(:,:,:)
  complex *16, allocatable :: splot(:,:)
  integer, allocatable :: kovers(:),nps(:),ipstart(:)

  integer i,j,k,l,ipatch,npout,kover,npols
  integer itrip,itric1,nb,nlmax,nuv,istart,iend,nd,nd2
  integer ilstart,itri,iunit1,m,ncell,ncsize,norder,nuvl,i1
  integer idim

  real *8 ra,erra

!
!  get the coefs of the density
!
   nd = 3
   nd2 = 6
   allocate(sigma_coefs(nd,npts))
   call surf_vals_to_coefs(nd2,npatches,norders,ixyzs,iptype,npts, &
     sigma,sigma_coefs)


 
!
!   estimate kovers, nps 
!

  npout = 0
  allocate(kovers(npatches),nps(npatches),ipstart(npatches+1))
  do i=1,npatches
    npols = ixyzs(i+1)-ixyzs(i)
    kover = 0
    if(npols.gt.4**0) kover = 1
    if(npols.gt.4**1) kover = 2
    if(npols.gt.4**2) kover = 3
    if(npols.gt.4**3) kover = 4
    kovers(i) = kover
    nps(i) = 4**kover
    if(iptype(i).eq.1) nps(i) = nps(i)*3
    npout = npout + nps(i) 
  enddo


  ipstart(1) = 1
  nps(1) = nps(1) + 1
  call cumsum(npatches,nps,ipstart(2))
  nps(1) = nps(1) - 1

  allocate(xyzs(3,npout),splot(nd,npout))

!
!   get uvs of all patches of type = 1
!
  
  nlmax = 5
  nuv = (4**(nlmax+1)-1)/3
  allocate(uvs(2,3,nuv))

  uvs(1,1,1) = 0
  uvs(2,1,1) = 0
  uvs(1,2,1) = 1
  uvs(2,2,1) = 0
  uvs(1,3,1) = 0
  uvs(2,3,1) = 1

  do i=0,nlmax-1
    istart = (4**(i)-1)/3+1
    nb = 4**i
    iend = istart + nb-1
    do itrip = istart,iend
      itric1 = (itrip-istart)*4 + iend
      call gettrichildren(uvs(1,1,itrip),uvs(1,1,itric1+1), &
       uvs(1,1,itric1+2),uvs(1,1,itric1+3),uvs(1,1,itric1+4))   
    enddo
  enddo


  do ipatch=1,npatches
    istart = ipstart(ipatch)
    npols = ixyzs(ipatch+1)-ixyzs(ipatch)
    norder = norders(ipatch)
    allocate(pols(npols))
    if(iptype(ipatch).eq.1) then

      nuvl = ipstart(ipatch+1)-ipstart(ipatch)
      ilstart = (4**(kovers(ipatch))-1)/3+1
      nb = 4**(kovers(ipatch))
      do i=1,nb
        itri = i+ilstart-1
        do j=1,3
          call koorn_pols(uvs(1,j,itri),norder,npols,pols)
          
          do m=1,3
            xyzs(m,istart+3*(i-1)+j-1) = 0
          enddo
          do idim=1,nd
            splot(idim,istart+3*(i-1)+j-1) = 0
          enddo

          do l=1,npols
            do m=1,3
              xyzs(m,istart+3*(i-1)+j-1) = & 
                xyzs(m,istart+3*(i-1)+j-1) + &
                pols(l)*srccoefs(m,ixyzs(ipatch)+l-1)
            enddo
            do idim=1,nd
              splot(idim,istart+3*(i-1)+j-1) = &
               splot(idim,istart+3*(i-1)+j-1)+ &
               pols(l)*sigma_coefs(idim,ixyzs(ipatch)+l-1)
            enddo
          enddo
        enddo
      enddo
    endif
    deallocate(pols)
  enddo
  
  iunit1 = 877
  open(unit = iunit1, file=trim(fname))

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i9,a)') "POINTS ", npout, " float"

  do i = 1,npout
    write(iunit1,"(E11.5,2(2x,e11.5))") xyzs(1,i), xyzs(2,i), xyzs(3,i)
  end do

  ncell = 0
  ncsize = 0
  do i=1,npatches
    ncell = ncell + 4**kovers(i)
    if(iptype(i).eq.1) ncsize = ncsize + 4*(4**kovers(i))
  enddo

  write(iunit1,'(a,i9,i9)') "CELLS ", ncell, ncsize

  do ipatch=1,npatches
    nb = 4**kovers(ipatch)
    if(iptype(ipatch).eq.1) then
      istart = ipstart(ipatch) 
      do i = 1,nb
        i1 = istart + 3*(i-1) 
        write(iunit1,'(a,i9,i9,i9)') "3 ", i1-1, i1, i1+1
      enddo
    endif
  end do

  write(iunit1,'(a,i9)') "CELL_TYPES ", ncell
  do ipatch = 1,npatches
    nb = 4**kovers(ipatch)
    if(iptype(ipatch).eq.1) then
      do i=1,nb
        write(iunit1,'(a)') "5"
      enddo
    endif
  end do

  write(iunit1,'(a)') ""
  write(iunit1,'(a,i9)') "POINT_DATA ", npout
  write(iunit1,'(a,i4)') "SCALARS vec_comps float ", 6
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,npout
    write(iunit1,'(5(E11.5,2x),E11.5)') &
      real(splot(1,i)),imag(splot(1,i)),real(splot(2,i)), &
      imag(splot(2,i)),real(splot(3,i)),imag(splot(3,i))
  end do

  close(iunit1)



end subroutine surf_vtk_plot_zvec
!
!


subroutine vtk_scatter_plot_scalar(n,xyzs,f,fname,title)
!
!  this subroutine creates a scatter plot of points with
!  given function values
!
!  Input arguments:
!    - n: integer
!        number of points
!    - xyzs: real *8 (3,n)
!         xyz coordinates of the points
!    - f: real *8(n)
!         function values
!    - fname: character*
!         file name for writing out vtk file
!    - title: character*
!         title of plot
!
  implicit none
  integer n
  character (len=*) fname,title
  real *8 xyzs(3,n),f(n)

  integer i
  integer iunit1

  
  iunit1 = 877
  open(unit = iunit1, file=trim(fname), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET POLYDATA"
  write(iunit1,'(a,i9,a)') "POINTS ", n, " float"

  do i = 1,n
    write(iunit1,"(E11.5,2(2x,e11.5))") xyzs(1,i), xyzs(2,i), xyzs(3,i)
  end do

  write(iunit1,'(a)') ""
  write(iunit1,'(a,i9)') "POINT_DATA ", n
  
  write(iunit1,'(a,i4)') "SCALARS scat float ", 1
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,n
    write(iunit1,'(E11.5)') f(i)
  enddo

  close(iunit1)

end subroutine vtk_scatter_plot_scalar
