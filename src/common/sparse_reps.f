c
c
c This file contains the following user callable routines:
c - conv_to_csc: convert a given set of (i,j) indices
c              to csc format
c - rsc_to_csc: convert a row sparse compressed data
c                 structure to a column sparse compressed data 
c                 structure or vice vrsa
c-----------------------------------------------------
c
c
c
c
c
c-------------------------------------------------
c> @brief
c> Convert a row sparse compressed representation to
c> a column sparse compressed representation or vice-vrsa
c>
c> @param[in] ncol: number of columns
c> @param[in] nrow: number of rows
c> @param[in] nnz: number of non-zero entries in row sparse format
c> @param[in] row_ptr: row_ptr(i) indiciates location in col_ind
c>               where relevant entries corresponding to row i begin
c> @param[in] col_ind: list of column indices. 
c>              (i,col_ind(row_ptr(i):row_ptr(i+1)-1)) are the non-zero
c>              entries
c> @param[out] col_ptr: col_ptr(i) indicates location in row_ind where
c>              relevant entries for column i begin
c> @param[out] row_ind: list of row_indices
c>          (row_ind(row_ptr(i):row_ptr(i+1)-1),i) are 
c>            all the non-zero entries
c> @param[out] iper: (irow,jcol) corresponding to row_ind(iper(i))
c>          is the same as (irow,jcol) corresponding to col_ind(i)
c> 
c> @todo openmp version, python interface
c>-------------------------------------------------
      subroutine rsc_to_csc(ncol,nrow,nnz,row_ptr,col_ind,
     1    col_ptr,row_ind,iper)
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: ncol,nnz,nrow
      integer, intent(in) :: row_ptr(nrow+1),col_ind(nnz)
      integer, intent(out) :: col_ptr(ncol+1),row_ind(nnz)
      integer, intent(out) :: iper(nnz)

      integer, allocatable :: nslr(:)
      integer i,itarg,ictr

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
c
c----------------------------------------------------------
c> @brief
c> Convert a list of (i,j) matrix entries to column sparse
c> compressed format
c>
c> @author manas rachh
c> 
c> @param[in] nent: number of entries (i,j)
c> @param[in] m: number of columns in the matrix
c> @param[in] iind: list of row indices of the (i,j) tuples
c> @param[in] jind: list of column indices of the (i,j) tuples
c> @param[out] col_ptr: col_ptr(i) indicates location in row_ind where
c>              relevant entries for column i begin
c> @param[out] iper: (irow,jcol) corresponding to row_ind(iper(i))
c>          is the same as (irow,jcol) corresponding to col_ind(i)
c>
c> @todo 
c> openmp version
c------------------------------------------------------------
      subroutine conv_to_csc(nent,m,iind,jind,col_ptr,row_ind)
      implicit none
c
c
c
cf2py  intent(in) nent,m,iind,jind
cf2py  intent(out) col_ptr,row_ind
      integer, intent(in) :: nent,m,iind(nent),jind(nent)
      integer, intent(out) :: col_ptr(m+1),row_ind(nent)
      integer, allocatable :: jsort(:),jper(:),icnt(:)
      integer i,icur,icur0,icur1


      allocate(jper(nent),jsort(nent))

      call sorti(nent,jind,jper)

      do i=1,nent
        jsort(i) = jind(jper(i))
        row_ind(i) = iind(jper(i))
      enddo

      icur = 1
      col_ptr(1) = 1

      do i=1,nent
        if(jsort(i).gt.icur) then
          icur0 = icur
          icur = jsort(i)

          do icur1 = icur0+1,icur
            col_ptr(icur1) = i
          enddo
        endif
      enddo

      do i=icur+1,m+1
        col_ptr(i) = nent+1
      enddo



      return
      end

      
