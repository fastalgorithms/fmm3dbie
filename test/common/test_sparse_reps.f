      implicit real *8 (a-h,o-z)

      call test_conv_to_csc()
      call test_csc_to_rsc()

      stop
      end



      subroutine test_conv_to_csc()
      implicit real *8 (a-h,o-z)

      return
      end

      

      
      subroutine test_csc_to_rsc()
      
      implicit real *8 (a-h,o-z)
      integer(8), allocatable :: row_ptr(:),col_ind(:),col_ptr(:),
     1   row_ind(:),iper(:)

      call prini(6,13)
      
      ns = 301
      nt = 101

      allocate(row_ptr(nt+1),col_ptr(ns+1))
      

      nnz_max = ns*nt
      allocate(col_ind(nnz_max),iper(nnz_max),row_ind(nnz_max))

      
      row_ptr(1) = 1
      do i=1,nt
        nrel = hkrand(0)*0.5*ns
        do j=1,nrel
          col_ind(row_ptr(i)+j-1) = ceiling(hkrand(0)*ns)
        enddo
        row_ptr(i+1) = row_ptr(i) + nrel
      enddo
      nnz = row_ptr(nt+1)-1

      call rsc_to_csc(ns,nt,nnz,row_ptr,col_ind,col_ptr,row_ind,iper)

      do i=1,ns
        do j=col_ptr(i),col_ptr(i+1)-1
           itarg = row_ind(j)
           ii1 = iper(j)
           isrc = col_ind(ii1)
           if(ii1.lt.row_ptr(itarg).or.ii1.ge.row_ptr(itarg+1)) then
              call prinf('error in permutation array*',i,0)
              print *, ii1,row_ptr(itarg),row_ptr(itarg+1)
              stop
           endif

           if(isrc.ne.i) then
             call prinf('error in converting rsc to csc*',i,0)
             print *, "No match found for interaction between ",isrc,
     1          itarg
             stop

           endif
        enddo
      enddo

      print *, "csc to rsc Test completed successfully"
       

      return
      end
