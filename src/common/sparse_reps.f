c
c              rsc_to_csc - convert a row sparse compressed data
c               structure to a column sparse compressed data 
c               structure or vice vrsa
c
c
c              conv_to_csc - convert a given set of (i,j) indices
c                 to csc format
c
c-----------------------------------------------------
c
c
c
c
c

      subroutine rsc_to_csc(ncol,nrow,nnz,row_ptr,col_ind,
     1    col_ptr,row_ind,iper)
c
c
c
c       convert a row sparse compressed representation
c        to a column sparse compressed representation
c        or vice-vrsa
c
c     input
c       ncol - integer
c          number of columns
c       nrow - integer
c          number of rows
c       nnz - integer
c          number of non-zero entries in row sparse
c          compressed format
c       row_ptr - integer (nrow+1)
c          row_ptr(i) indicates location in col_ind
c          where relevant entries corresponding
c          to row i begin
c       col_ind - integer(nnz)
c          list of column indices. 
c          (i,col_ind(row_ptr(i):row_ptr(i+1)-1)) are 
c            all the non-zero entries
c       
c      output
c        col_ptr - integer (nsrc+1)
c          col_ptr(i) indicates location in row_ind
c           where relevant entries for column i 
c           begin
c
c       row_ind - integer(nnz)
c          list of row indices. 
c          (row_ind(row_ptr(i):row_ptr(i+1)-1),i) are 
c            all the non-zero entries
c
c       iper - integer(nnz)
c          (irow,jcol) corresponding to row_ind(iper(i))
c            is the same as (irow,jcol) corresponding to
c            col_ind(i)
c



      implicit real *8 (a-h,o-z)      
       integer ncol,nnz,row_ptr(nrow+1),col_ind(nnz)
       integer col_ptr(ncol+1),row_ind(nnz),iper(nnz)
       integer, allocatable :: nslr(:)

      allocate(nslr(ncol))
      lflg = nnz
      do i=1,ncol
         nslr(i) = 0
      enddo

      do i=1,lflg
        nslr(col_ind(i)) = nslr(col_ind(i)) + 1
      enddo


      col_ptr(1) = 1
      do i=2,ncol+1
         col_ptr(i) = col_ptr(i-1)+nslr(i-1)
      enddo


      do itarg=1,nrow
         do ictr=row_ptr(itarg),row_ptr(itarg+1)-1
           jsrc = col_ind(ictr)

           iper(col_ptr(jsrc)) = ictr

           row_ind(col_ptr(jsrc)) = itarg
           col_ptr(jsrc) = col_ptr(jsrc) + 1 
         enddo
      enddo

      col_ptr(1) = 1
      do i=2,ncol+1
         col_ptr(i) = col_ptr(i-1)+nslr(i-1)
      enddo

      return
      end
c
c
c
c
c-----------------------------------
      subroutine conv_to_csc(nent,m,iind,jind,col_ptr,row_ind)
      implicit none
      integer nent,m,iind(nent),jind(nent),col_ptr(m+1),row_ind(nent)
      integer, allocatable :: jsort(:),jper(:),icnt(:)
      integer i,icur


      allocate(jper(nent),jsort(nent),icnt(m))

      call sorti(nent,jind,jper)

      do i=1,nent
        jsort(i) = jind(jper(i))
        row_ind(i) = iind(jper(i))
      enddo

      icur = 1

      do i=1,nent
        if(jsort(i).gt.icur) then
          icnt(icur) = i-1
          icur = icur+1
        endif
      enddo

      icnt(1) = icnt(1) + 1
      col_ptr(1) = 1
      call cumsum(m,icnt,col_ptr(2))


      return
      end

      
