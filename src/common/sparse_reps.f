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
c> @todo python interface
c-------------------------------------------------
      subroutine rsc_to_csc(ncol,nrow,nnz,row_ptr,col_ind,
     1    col_ptr,row_ind,iper)
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: ncol,nnz,nrow
      integer, intent(in) :: row_ptr(nrow+1),col_ind(nnz)
      integer, intent(out) :: col_ptr(ncol+1),row_ind(nnz)
      integer, intent(out) :: iper(nnz)

      integer, allocatable :: row_ind_exp(:)

      integer, allocatable :: nslr(:)
      integer i,itarg,ictr

      allocate(row_ind_exp(nnz))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,nrow
        do j=row_ptr(i),row_ptr(i+1)-1
          row_ind_exp(j) = i
        enddo
      enddo
C$OMP END PARALLEL DO     

      if(nrow.lt.1000.or.nnz.lt.10000) then
        call sorti(nnz,col_ind,iper)
        call conv_to_csc(nnz,ncol,row_ind_exp,col_ind,col_ptr,row_ind)
      else
        call sorti(nnz,col_ind,iper)
        call conv_to_csc(nnz,ncol,row_ind_exp,col_ind,col_ptr,
     1     row_ind)
      endif

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
c> @todo openmp version
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
      integer i,icur,icur0,icur1,j


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
c
c
c
c
c
c----------------------------------------------------------
c> @brief
c> Convert a list of (i,j) matrix entries to column sparse
c> compressed format.
c> This is openmp version
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
c------------------------------------------------------------
      subroutine conv_to_csc_para(nent,m,iind,jind,col_ptr,row_ind)
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
      integer, external :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
      integer, external :: OMP_GET_MAX_THREADS
      integer nthreads,ithread,first,mend,istart,iend
      integer, allocatable :: split(:),ind(:)

      allocate(jper(nent),jsort(nent))
      allocate(ind(nent))

      call sorti_para(nent,jind,jper)

      nthreads=1
C$    nthreads=OMP_GET_MAX_THREADS()

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i)
      do i=1,nent
        ind(i) = i
        jsort(i) = jind(jper(i))
        row_ind(i) = iind(jper(i))
      enddo
C$OMP END PARALLEL DO

      allocate(split(nthreads+1))
      split(nthreads+1) = nent+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ithread,i,first)
      do ithread = 1,nthreads
        i = ((ithread-1)*nent)/nthreads+1
        call lower_bound(jsort,ind,nent,jsort(i),first)
        split(ithread) = first
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ithread,istart,iend,icur,i,icur0,icur1,mend)
      do ithread=1,nthreads
        istart = split(ithread)
        iend = split(ithread+1)
        icur = jsort(istart)
        if(icur.eq.0) then
          icur=1
        endif
        col_ptr(icur) = istart
        do i=istart,iend-1
          if(jsort(i).gt.icur) then
            icur0 = icur
            icur = jsort(i)
            do icur1 = icur0+1,icur
              col_ptr(icur1) = i
            enddo
          endif
        enddo
        if(ithread.eq.nthreads) then
          mend=m+1
        else
          mend=jsort(iend)
        endif
        do i=icur+1,mend
          col_ptr(i) = iend
        enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
