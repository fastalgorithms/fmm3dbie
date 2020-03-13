c
c
c    This file contains routines for estimating unique
c    set of indices from a given collection of arrays
c
c
c 

      subroutine get_iuni1_omp(n,a,nuni,iuni,iuniind)
c
c        given an array a, find unique set of indices
c        in the array and a mapping from where the indices
c        are contained in the shorter array
c
      integer n, a(n),nuni,iuni(n),iuniind(n)
      integer, allocatable :: iind(:),asort(:),iunisort(:)

      allocate(iind(n),asort(n),iunisort(n))
      call sorti(n,a,iind)

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,n
        asort(i) = a(iind(i))
      enddo
C$OMP END PARALLEL DO      

      nuni = 1
      iuni(1) = asort(1)
      iunisort(1) = 1
      do i=2,n
        if(asort(i).gt.iuni(nuni)) then
          nuni = nuni + 1
          iuni(nuni) = asort(i)
        endif
        iunisort(i) = nuni
      enddo

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,n
        iuniind(iind(i)) = iunisort(i)
      enddo
C$OMP END PARALLEL DO

      return
      end


      subroutine get_iuni3(n,a,b,c,nuni,iuni,iuniind)
c
c        given an array a, find unique set of indices
c        in the array and a mapping from where the indices
c        are contained in the shorter array
c
      integer n, a(n),b(n),c(n),nuni,iuni(3,n),iuniind(n)
      integer, allocatable :: iind(:),asort(:),iunisort(:)
      integer, allocatable :: bsort(:),csort(:)

      allocate(iind(n),asort(n),iunisort(n),bsort(n),csort(n))
      call sorti(n,c,iind)

      do i=1,n
        asort(i) = a(iind(i))
        bsort(i) = b(iind(i))
        csort(i) = c(iind(i))
      enddo

      nuni = 1
      iuni(1,1) = asort(1)
      iuni(2,1) = bsort(1)
      iuni(3,1) = csort(1)
      iunisort(1) = 1
 
      nuni0 = 1
      do i=2,n
        if(csort(i).gt.iuni(3,nuni)) then
          nuni = nuni + 1
          iuni(1,nuni) = asort(i)
          iuni(2,nuni) = bsort(i)
          iuni(3,nuni) = csort(i)
          iunisort(i) = nuni
          nuni0 = nuni
        else
          do j=nuni,nuni0,-1
            if(asort(i).eq.iuni(1,j).and.bsort(i).eq.iuni(2,j)) then 
              iunisort(i) = j
              goto 1111
            endif
          enddo

          nuni = nuni + 1
          iuni(1,nuni) = asort(i)
          iuni(2,nuni) = bsort(i)
          iuni(3,nuni) = csort(i)
          iunisort(i) = nuni
        endif
 1111 continue

      enddo

      do i=1,n
        iuniind(iind(i)) = iunisort(i)
      enddo

      return
      end

