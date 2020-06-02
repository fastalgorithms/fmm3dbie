c      PROGRAM driver
c      IMPLICIT NONE
c      INTEGER xx(50)
c      INTEGER n, index(50), i, index2(50)
c
c      n = 12
c
c      index(01) = 07
c      index(02) = 03
c      index(03) = 01
c      index(04) = 12
c      index(05) = 11
c      index(06) = 02
c      index(07) = 05
c      index(08) = 10
c      index(09) = 04
c      index(10) = 08
c      index(11) = 06
c      index(12) = 09
c
c      DO 1000, i = 1,n
c        xx(index(i)) = i
c 1000 CONTINUE
c
c      DO 2000, i = 1,n
c        write(*,*) 'xx(',i,') = ', xx(i)
c 2000 CONTINUE
c
c      CALL SORTI(n,xx,index2)
c
c      DO 3000, i = 1,n
c        write(*,*) index(i), index2(i)
c 3000 CONTINUE
c
c      STOP
c      END PROGRAM driver
c
ccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SORTI(N,DATA,INDEX)

C===================================================================
C
C     SORTRX -- SORT, integer input, indeX output
C
C
C     Input:  N     INTEGER
C             DATA  INTEGER
C
C     Output: INDEX INTEGER (DIMENSION N)
C
C This routine performs an in-memory sort of the first N elements of
C array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================
C
C SORTRX uses a hybrid QuickSort algorithm, based on several
C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
C "pivot key" [my term] for dividing each subsequence is chosen to be
C the median of the first, last, and middle values of the subsequence;
C and the QuickSort is cut off when a subsequence has 9 or fewer
C elements, and a straight insertion sort of the entire array is done
C at the end.  The result is comparable to a pure insertion sort for
C very short arrays, and very fast for very large arrays (of order 12
C micro-sec/element on the 3081K for arrays of 10K elements).  It is
C also not subject to the poor performance of the pure QuickSort on
C partially ordered data.
C
C Created:  15 Jul 1986  Len Moss
C
C===================================================================
 
      INTEGER   N,INDEX(N)
      INTEGER   DATA(N)
 
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      integer      DATAP
 
C     QuickSort Cutoff
C
C     Quit QuickSort-ing when a subsequence contains M or fewer
C     elements and finish off at end with straight insertion sort.
C     According to Knuth, V.3, the optimum value of M is around 9.
 
      INTEGER   M
      PARAMETER (M=9)
 
C===================================================================
C
C     Make initial guess for INDEX
 
      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE
 
C     If array is short, skip QuickSort and go directly to
C     the straight insertion sort.
 
      IF (N.LE.M) GOTO 900
 
C===================================================================
C
C     QuickSort
C
C     The "Qn:"s correspond roughly to steps in Algorithm Q,
C     Knuth, V.3, PP.116-117, modified to select the median
C     of the first, last, and middle elements as the "pivot
C     key" (in Knuth's notation, "K").  Also modified to leave
C     data in place and produce an INDEX array.  To simplify
C     comments, let DATA[I]=DATA(INDEX(I)).
 
C Q1: Initialize
      ISTK=0
      L=1
      R=N
 
  200 CONTINUE
 
C Q2: Sort the subsequence DATA[L]..DATA[R].
C
C     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
C     r > R, and L <= m <= R.  (First time through, there is no
C     DATA for l < L or r > R.)
 
      I=L
      J=R
 
C Q2.5: Select pivot key
C
C     Let the pivot, P, be the midpoint of this subsequence,
C     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
C     so the corresponding DATA values are in increasing order.
C     The pivot key, DATAP, is then DATA[P].
 
      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)
 
      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
C     Now we swap values between the right and left sides and/or
C     move DATAP until all smaller values are on the left and all
C     larger values are on the right.  Neither the left or right
C     side will be internally ordered yet; however, DATAP will be
C     in its final position.
 
  300 CONTINUE
 
C Q3: Search for datum on left >= DATAP
C
C     At this point, DATA[L] <= DATAP.  We can therefore start scanning
C     up from L, looking for a value >= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300
 
  400 CONTINUE
 
C Q4: Search for datum on right <= DATAP
C
C     At this point, DATA[R] >= DATAP.  We can therefore start scanning
C     down from R, looking for a value <= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400
 
C Q5: Have the two scans collided?
 
      IF (I.LT.J) THEN
 
C Q6: No, interchange DATA[I] <--> DATA[J] and continue
 
         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE
 
C Q7: Yes, select next subsequence to sort
C
C     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
C     for all L <= l < I and J < r <= R.  If both subsequences are
C     more than M elements long, push the longer one on the stack and
C     go back to QuickSort the shorter; if only one is more than M
C     elements long, go back and QuickSort it; otherwise, pop a
C     subsequence off the stack and QuickSort it.
 
         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
C Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
 
  900 CONTINUE
 
C===================================================================
C
C Q9: Straight Insertion sort
 
      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE
 
      END
c
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SORTR(N,DATA,INDEX)

C===================================================================
C
C     SORTRX -- SORT, integer input, indeX output
C
C
C     Input:  N     INTEGER
C             DATA  REAL
C
C     Output: INDEX INTEGER (DIMENSION N)
C
C This routine performs an in-memory sort of the first N elements of
C array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================
C
C SORTRX uses a hybrid QuickSort algorithm, based on several
C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
C "pivot key" [my term] for dividing each subsequence is chosen to be
C the median of the first, last, and middle values of the subsequence;
C and the QuickSort is cut off when a subsequence has 9 or fewer
C elements, and a straight insertion sort of the entire array is done
C at the end.  The result is comparable to a pure insertion sort for
C very short arrays, and very fast for very large arrays (of order 12
C micro-sec/element on the 3081K for arrays of 10K elements).  It is
C also not subject to the poor performance of the pure QuickSort on
C partially ordered data.
C
C Created:  15 Jul 1986  Len Moss
C
C===================================================================
 
      INTEGER   N,INDEX(N)
      REAL *8   DATA(N)
 
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      REAL *8      DATAP
 
C     QuickSort Cutoff
C
C     Quit QuickSort-ing when a subsequence contains M or fewer
C     elements and finish off at end with straight insertion sort.
C     According to Knuth, V.3, the optimum value of M is around 9.
 
      INTEGER   M
      PARAMETER (M=9)
 
C===================================================================
C
C     Make initial guess for INDEX
 
      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE
 
C     If array is short, skip QuickSort and go directly to
C     the straight insertion sort.
 
      IF (N.LE.M) GOTO 900
 
C===================================================================
C
C     QuickSort
C
C     The "Qn:"s correspond roughly to steps in Algorithm Q,
C     Knuth, V.3, PP.116-117, modified to select the median
C     of the first, last, and middle elements as the "pivot
C     key" (in Knuth's notation, "K").  Also modified to leave
C     data in place and produce an INDEX array.  To simplify
C     comments, let DATA[I]=DATA(INDEX(I)).
 
C Q1: Initialize
      ISTK=0
      L=1
      R=N
 
  200 CONTINUE
 
C Q2: Sort the subsequence DATA[L]..DATA[R].
C
C     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
C     r > R, and L <= m <= R.  (First time through, there is no
C     DATA for l < L or r > R.)
 
      I=L
      J=R
 
C Q2.5: Select pivot key
C
C     Let the pivot, P, be the midpoint of this subsequence,
C     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
C     so the corresponding DATA values are in increasing order.
C     The pivot key, DATAP, is then DATA[P].
 
      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)
 
      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
C     Now we swap values between the right and left sides and/or
C     move DATAP until all smaller values are on the left and all
C     larger values are on the right.  Neither the left or right
C     side will be internally ordered yet; however, DATAP will be
C     in its final position.
 
  300 CONTINUE
 
C Q3: Search for datum on left >= DATAP
C
C     At this point, DATA[L] <= DATAP.  We can therefore start scanning
C     up from L, looking for a value >= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300
 
  400 CONTINUE
 
C Q4: Search for datum on right <= DATAP
C
C     At this point, DATA[R] >= DATAP.  We can therefore start scanning
C     down from R, looking for a value <= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400
 
C Q5: Have the two scans collided?
 
      IF (I.LT.J) THEN
 
C Q6: No, interchange DATA[I] <--> DATA[J] and continue
 
         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE
 
C Q7: Yes, select next subsequence to sort
C
C     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
C     for all L <= l < I and J < r <= R.  If both subsequences are
C     more than M elements long, push the longer one on the stack and
C     go back to QuickSort the shorter; if only one is more than M
C     elements long, go back and QuickSort it; otherwise, pop a
C     subsequence off the stack and QuickSort it.
 
         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
C Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
 
  900 CONTINUE
 
C===================================================================
C
C Q9: Straight Insertion sort
 
      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE
 
      END
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sorti_para(n,data,index)
C===================================================================
C
C     SORTRX -- SORT, integer input, indeX output
C
C
C     Input:  N     INTEGER
C             DATA  INTEGER
C
C     Output: INDEX INTEGER (DIMENSION N)
C
C This openmp threaded routine performs an sort of the first N elements
C of array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================

      implicit none
      integer n,data(n),index(n)
      integer nthreads,nelems,ielems,istart,iend,i,j,swap
      integer(8) ithread
      integer istart2,iend2,nelems2
      integer, allocatable :: split(:)
      integer, allocatable :: index2(:)
      integer, external :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
      integer, external :: OMP_GET_MAX_THREADS

      nthreads=1
C$    nthreads=OMP_GET_MAX_THREADS()
c      write(*,*) "max num threads: ", nthreads

      if(N<2*nthreads .or. nthreads==1) then
c        write(*,*) "call sorti"
        call sorti(n,data,index)
        return
      endif

c      write(*,*) "call omp sort"
ccc   split array evenly
      allocate(split(nthreads+1))
      allocate(index2(n))
      split(nthreads+1)=n+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ithread)
      do ithread = 1,nthreads
        split(ithread) = ((ithread-1)*n)/nthreads+1
      enddo
C$OMP END PARALLEL DO
c      write(*,*) nthreads+1,split(1),split(2),split(3)

ccc   sort each part
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ithread,nelems,istart,iend)
      do ithread = 1,nthreads
        istart = split(ithread)
        iend = split(ithread+1)
        nelems = iend-istart
c        write(*,*) "nelems: ", nelems
        call sorti(nelems,data(istart),index(istart))
        do ielems = istart,iend-1
          index(ielems) = index(ielems) + istart - 1
        enddo
      enddo
C$OMP END PARALLEL DO

ccc   merge two parts at a time
      swap=1
      j=1
      do while (j<nthreads)
        i=1
        do while (i<nthreads+1)
          if(i+j<nthreads+1) then
            istart=split(i)
            iend=split(i+j)
            istart2=split(i+j)
            if(i+2*j>nthreads+1) then
              iend2=split(nthreads+1)
            else
              iend2=split(i+2*j)
            endif
            nelems=iend-istart
            nelems2=iend2-istart2
            if(swap==1) then
              call merge_arr_para(data,index(istart),index(istart2),
     1             nelems,nelems2,index2(istart),nthreads)
            else
              call merge_arr_para(data,index2(istart),index2(istart2),
     1             nelems,nelems2,index(istart),nthreads)
            endif
          else
            if(swap==0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(istart)
              do istart=split(i),split(nthreads+1)-1
                index(istart) = index2(istart)
              enddo
C$OMP END PARALLEL DO
            else
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(istart)
              do istart=split(i),split(nthreads+1)-1
                index2(istart) = index(istart)
              enddo
C$OMP END PARALLEL DO
            endif
          endif
          i=i+2*j
        end do
        j=j*2
        if(swap==1) then
          swap=0
        else
          swap=1
        endif
      end do

      if(swap==0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(istart)
        do istart=1,n
          index(istart) = index2(istart)
        enddo
C$OMP END PARALLEL DO
      endif

      deallocate(split)
      deallocate(index2)
      END
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine merge_arr_para(data,ind1,ind2,n1,n2,indout,nthreads)
C===================================================================
C
C This openmp threaded routine performs merge of two arrays.
C ind1 contains the index of array1 of data, ind2 contains the index
C of array2 of data.
C The merged index set is in indout of increaing order, it contains
C the index respect to data, the original data is not physically
C rearranged.
C
C===================================================================
      implicit none
      integer n1,n2,nthreads,i,j,k,ns,indx,first,val,split1,split2
      integer(8) tmpk
      integer split_size1,split_size2,n1tmp,n2tmp
      integer istart1,istart2,istart3
      integer data(*)
      integer ind1(n1),ind2(n2),indout(n1+n2)
      integer, allocatable :: split(:), split_size(:)
      integer, allocatable :: split_ind1(:), split_ind2(:),indtmp(:)
      if(n1==0 .and. n2==0) return
      if(n1==0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,k,tmpk)
        do k=0,nthreads-1
          tmpk = k
          i = (tmpk*n2)/nthreads+1
          j = ((tmpk+1)*n2)/nthreads
          indout(i:j) = ind2(i:j)
        enddo
C$OMP END PARALLEL DO
        return
      endif

      if(n2==0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,k,tmpk)
        do k=0,nthreads-1
          tmpk = k
          i = (tmpk*n1)/nthreads+1
          j = ((tmpk+1)*n1)/nthreads
          indout(i:j) = ind1(i:j)
        enddo
C$OMP END PARALLEL DO
        return
      endif

      ns=10
      allocate(split(nthreads*ns*2))
      allocate(split_size(nthreads*ns*2))
      allocate(indtmp(nthreads*ns*2))
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,k,indx,first,val,tmpk)
      do i=0,nthreads-1
        do j=0,ns-1
          k = i*ns+j
          tmpk = k
          indtmp(k+1)=k+1
          indx = (tmpk*n1)/(nthreads*ns)
          split(k+1) = indx+1
          val = data(ind1(split(k+1)))
          call lower_bound(data,ind2,n2,val,first)
          split_size(k+1) = indx+first-1

          indx = (tmpk*n2)/(nthreads*ns)
          k = k+nthreads*ns
          split(k+1) = indx+1
          val = data(ind2(split(k+1)))
          call lower_bound(data,ind1,n1,val,first)
          split_size(k+1) = indx+first-1
        enddo
      enddo
C$OMP END PARALLEL DO

      allocate(split_ind1(nthreads+1),split_ind2(nthreads+1))
      split_ind1(1)=1
      split_ind2(1)=1
      split_ind1(nthreads+1)=n1+1
      split_ind2(nthreads+1)=n2+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,k,split1,split2,split_size1,split_size2,val,first)
C$OMP$PRIVATE(tmpk)
      do i=1,nthreads-1
        tmpk = i
        k = (tmpk*(n1+n2))/nthreads

        call lower_bound(split_size(1),indtmp,nthreads*ns,
     1       k,first)
        j=first
        if(j>nthreads*ns) j=nthreads*ns
        split1=split(j)
        split_size1=split_size(j)

        call lower_bound(split_size(1+nthreads*ns),indtmp,nthreads*ns,
     1       k,first)
        j=first+nthreads*ns
        if(j>2*nthreads*ns) j=2*nthreads*ns
        if(abs(k-split_size(j)) < abs(k-split_size1)) then
          val = data(ind2(split(j)))
        else
          val = data(ind1(split1))
        endif
        call lower_bound(data,ind1,n1,val,first)
        split_ind1(i+1) = first
        call lower_bound(data,ind2,n2,val,first)
        split_ind2(i+1) = first
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,n1tmp,n2tmp,istart1,istart2,istart3)
      do i=1,nthreads
        n1tmp = split_ind1(i+1) - split_ind1(i)
        n2tmp = split_ind2(i+1) - split_ind2(i)
        istart1 = split_ind1(i)
        istart2 = split_ind2(i)
        istart3 = split_ind1(i)+split_ind2(i)-1
c        write(*,*) "i:",i,"n1:",n1tmp,"n2:",n2tmp
        call merge_arr(data,ind1(istart1),ind2(istart2),n1tmp,n2tmp,
     1       indout(istart3),nthreads)
      enddo
C$OMP END PARALLEL DO

      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine merge_arr(data,ind1,ind2,n1,n2,indout,nthreads)
      implicit none
      integer n1,n2,nthreads,i,j,k
      integer data(*)
      integer ind1(n1),ind2(n2),indout(n1+n2)

      i=1
      j=1
      k=1

      do while(i<=n1.and.j<=n2)
        if(data(ind1(i))<data(ind2(j))) then
          indout(k)=ind1(i)
          i=i+1
        else
          indout(k)=ind2(j)
          j=j+1
        endif
        k=k+1
      enddo

      do while(i<=n1)
        indout(k) = ind1(i)
        i=i+1
        k=k+1
      enddo

      do while(j<=n2)
        indout(k) = ind2(j)
        j=j+1
        k=k+1
      enddo

      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lower_bound(data,ind,n,val,first)
        integer data(*),ind(*),n
        integer val
        integer it,first
        integer cnt,step

        first=1
        cnt=n
        do while(cnt>0)
          it = first
          step = cnt/2
          it = it + step
          if(data(ind(it))<val) then
            first = it+1
            cnt = cnt-step-1
          else
            cnt = step
          endif
        enddo

      return
      end
