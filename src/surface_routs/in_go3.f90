!
!
!  This file contains the input routines for handling Felipe's
!  go3 format for input - a format meant to handle 
!  triangulation of surfaces with RV nodes, where each patch
!  is discretized using the same order RV nodes
!
!
!  The file has 2 user callable routines 
!    open_gov3_geometry_mem - estimates the memory requirements
!    open_gov3_geometry - initializes the various arrays 
!


subroutine open_gov3_geometry(filename,npatches,norders,ixyzs, &
  iptype,npoints,srcvals,srccoefs,wts)
implicit none

    !List of calling arguments
    character (len=*), intent(in) :: filename
	integer ( kind = 4 ), intent(in) :: npoints
	integer ( kind = 4 ), intent(in) :: npatches
    integer ( kind = 4), intent(out) :: norders(npatches)
    integer ( kind = 4), intent(out) :: ixyzs(npatches+1)
    integer ( kind = 4), intent(out) :: iptype(npatches)
	real ( kind = 8 ), intent(out) :: srcvals(12,npoints)
	real ( kind = 8 ), intent(out) :: srccoefs(9,npoints)
	real ( kind = 8 ), intent(out) :: wts(npoints)
	
    !List of local variables
    integer ( kind = 4 ) umio,count1,count2,flag,aux,npols,icount
    integer ( kind = 4) norder
    integer :: ierror
	real ( kind = 8 ) aux_real,aux_vect(3)
	real ( kind = 8 ), allocatable :: h_points(:),h_coefs(:)
	real ( kind = 8 ), allocatable :: uv(:,:),umatr(:,:),vmatr(:,:),w(:)
    integer i

    open(UNIT=18, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)

        read(18,*) norder
        read(18,*) count1

        do count1=1,npoints
            read(18,*) srcvals(1,count1) !points(1,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(2,count1) !points(2,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(3,count1) !points(3,count1)
        enddo

        do count1=1,npoints
            read(18,*) srcvals(4,count1) !du(1,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(5,count1) !du(2,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(6,count1) !du(3,count1)
        enddo


        do count1=1,npoints
            read(18,*) srcvals(7,count1) !dv(1,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(8,count1) !dv(2,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(9,count1) !dv(3,count1)
        enddo

        do count1=1,npoints
            read(18,*) srcvals(10,count1) !normals(1,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(11,count1) !normals(2,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(12,count1) !normals(3,count1)
        enddo

        close (18)

		npols = (norder+1)*(norder+2)/2
        do i = 1,npatches
          norders(i) = norder
          iptype(i) = 1
          ixyzs(i) = (i-1)*npols + 1
        enddo
        ixyzs(npatches+1) = npoints+1

		allocate(h_points(npols))
		allocate(h_coefs(npols))
		allocate(uv(2,npols),umatr(npols,npols),vmatr(npols,npols),w(npols))
		call vioreanu_simplex_quad(norder,npols,uv,umatr,vmatr,w)
        
        call get_qwts(npatches,norders,ixyzs,iptype,npoints,srcvals,wts)


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(count1,h_points,h_coefs)
		do count1=1,npatches
			h_points(:)=srcvals(1,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(1,(count1-1)*npols+1:count1*npols)=h_coefs(:)					
			
			h_points(:)=srcvals(2,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(2,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

			h_points(:)=srcvals(3,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(3,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

			h_points(:)=srcvals(4,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(4,(count1-1)*npols+1:count1*npols)=h_coefs(:)					
			
			h_points(:)=srcvals(5,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(5,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

			h_points(:)=srcvals(6,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(6,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

			h_points(:)=srcvals(7,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(7,(count1-1)*npols+1:count1*npols)=h_coefs(:)					
			
			h_points(:)=srcvals(8,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(8,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

			h_points(:)=srcvals(9,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(9,(count1-1)*npols+1:count1*npols)=h_coefs(:)		
		enddo
!$OMP END PARALLEL DO        

return
end subroutine open_gov3_geometry



subroutine open_gov3_geometry_mem(filename,ntri,npoints)
implicit none

    !List of calling arguments
    character (len=*), intent(in) :: filename
	integer ( kind = 4 ), intent(out) :: npoints,ntri
	
    !List of local variables
    integer ( kind = 4 ) umio,count1,count2,flag,aux
    integer :: ierror,norder
	
        open(UNIT=18, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)

        read(18,*) norder
        read(18,*) ntri
		npoints=(norder+1)*(norder+2)/2*ntri

        close (18)

return
end subroutine open_gov3_geometry_mem
