c----------------------------
c  This file has the following user callable routines:
c
c    - setsub:
c        given two integer *8 sets a,b, compute a\ (a \cap b)^{c}, 
c        and b\ (a \cap b)^{c}
c
c    - setdecomp: 
c        given two integer *8 sets a,b, compute a \cap b and
c        a \cap b^{c} and mappings of the set a to various
c        elements of either sets
c
c    - get_iuni1:
c        compute list of unique integer *8s from a given
c        array of integer *8s
c 
c    - get_iuni1_omp: 
c        open mp version of get_iuni1 (note the
c        sorting routine needs to be updated to be openmped)
c
c    - get_iuni3: 
c        compute list of unique integer *8 triplets
c        from a given collection of 3 integer *8 arrays
c

      subroutine setsub(a,n,b,m,amb,namb,bma,nbma)
c
c-------------------------------
c  Compute a \cap b^{c} and b \cap a^{c} for two integer *8
c  sets with unique elements
c
c  Input arguments:
c
c    - a: integer *8(n)
c        integer *8 array 1
c    - n: integer *8
c        number of input elements for array 1
c    - b: integer *8(m)
c        integer *8 array 2
c    - m: integer *8 
c        number of input elements for array 2
c
c  Output arguments:
c
c    - amb: integer *8(n) 
c        elements which are in a but not in b.
c        On input of size (n), but only first namb
c        entries are relevant.
c    - namb: integer *8
c        number of elements in amb
c    - bma: integer *8(m)
c        elements which are in b but not in a.
c        On input of size (n), but only first nbma
c        entries are relevant.
c    - nbma: integer *8
c        number of elements in bma
c---------------------------

      
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: n,m
      integer *8, intent(in) :: a(n),b(m)
      integer *8, intent(out) :: amb(n),bma(m)
      integer *8, allocatable :: asort(:),bsort(:),w(:),abint(:)
      integer *8, allocatable :: iabinta(:),iabintb(:)

      namb = 0
      nbma = 0
      
      allocate(asort(n),bsort(m))
      nmax = max(n,m)

      
      allocate(w(2*(nmax+100)))

      do i=1,n
        asort(i) = a(i)
      enddo

      do i=1,m
        bsort(i) = b(i)
      enddo
      
      call sorti(n,a,w)
      do i=1,n
        asort(i) = a(w(i))
      enddo
      call sorti(m,b,w)
      do i=1,m
        bsort(i) = b(w(i))
      enddo


c
c      now compute the set intersection
c
      iaind = 1
      ibind = 1
      
      nmin = min(n,m)
      allocate(abint(nmin),iabinta(nmin),iabintb(nmin))
      iabint = 0

      if(n.le.m) then
        do i=1,n
          do j=ibind,m
            if(bsort(j).eq.asort(i)) then
              iabint = iabint + 1
              abint(iabint) = asort(i)
              iabinta(iabint) = i
              iabintb(iabint) = j 
              ibind = j+1
              goto 1000
            endif

            if(bsort(j).gt.asort(i)) then
              ibind = max(j-1,1)
              goto 1000
            endif
          enddo
 1000 continue          
        enddo
      endif

      

      
      if(m.lt.n) then
        do i=1,m
          do j=iaind,n
            if(bsort(i).eq.asort(j)) then
              iabint = iabint + 1
              abint(iabint) = bsort(i)
              iabinta(iabint) = j
              iabintb(iabint) = i 
              iaind = j+1
              goto 1010
            endif

            if(asort(j).gt.bsort(i)) then
              iaind = max(j-1,1)
              goto 1010
            endif
          enddo
 1010 continue          
        enddo
      endif



c
c   now find set subtractions
c
      iaind = 1
      namb = 0
      do i=1,n
        if(iaind.le.iabint) then
          if(i.eq.iabinta(iaind)) then
            iaind = iaind + 1
          else
            namb = namb + 1
            amb(namb) = asort(i)
          endif
        else
          namb = namb + 1
          amb(namb) = asort(i)
        endif
      enddo

      ibind = 1
      nbma = 0
      do i=1,m
        if(ibind.le.iabint) then
          if(i.eq.iabintb(ibind)) then
            ibind = ibind + 1
          else
            nbma = nbma + 1
            bma(nbma) = bsort(i) 
          endif
        else
          nbma = nbma + 1
          bma(nbma) = bsort(i)
        endif
      enddo


      return
      end
c
c
c
c
c
      subroutine setdecomp(n,a,m,b,naintb,aintb,iaintba,iaintbb,naintbc,
     1   aintbc,iaintbc)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 n,a(n),m,b(m),aintb(n),iaintba(n),aintbc(n),iaintbc(n)
      integer *8 iaintbb(n)
      integer *8, allocatable :: asort(:),bsort(:),iasort(:),ibsort(:)
      integer *8 nmin

      
      allocate(asort(n),bsort(m),iasort(n),ibsort(m))

      
      call sorti(n,a,iasort)
      do i=1,n
        asort(i) = a(iasort(i))
      enddo
      call sorti(m,b,ibsort)
      do i=1,m
        bsort(i) = b(ibsort(i))
      enddo




c
c      now compute the set intersection
c
      iaind = 1
      ibind = 1
      
      nmin = min(n,m)
      iabint = 0

      if(n.le.m) then
        do i=1,n
          do j=ibind,m
            if(bsort(j).eq.asort(i)) then
              iabint = iabint + 1
              aintb(iabint) = asort(i)
              iaintba(iabint) = i
              iaintbb(iabint) = j
              ibind = j+1
              goto 1000
            endif

            if(bsort(j).gt.asort(i)) then
              ibind = max(j-1,1)
              goto 1000
            endif
          enddo
 1000 continue          
        enddo
      endif

      

      
      if(m.lt.n) then
        do i=1,m
          do j=iaind,n
            if(bsort(i).eq.asort(j)) then
              iabint = iabint + 1
              aintb(iabint) = bsort(i)
              iaintba(iabint) = j
              iaintbb(iabint) = i
              iaind = j+1
              goto 1010
            endif

            if(asort(j).gt.bsort(i)) then
              iaind = max(j-1,1)
              goto 1010
            endif
          enddo
 1010 continue          
        enddo
      endif


      naintb = iabint

      if(naintb.eq.0) then
        naintbc = n
        do i=1,n
          aintbc(i) = asort(i)
          iaintbc(i) = i
        enddo
      else
        iabint = 1
        iabintc = 0
        do i=1,n
          if(i.ne.iaintba(iabint)) then
            iabintc = iabintc + 1
            aintbc(iabintc) = asort(i)
            iaintbc(iabintc) = i
          else
            iabint = iabint+1
            if(iabint.gt.naintb) iabint=naintb
          endif
        enddo
        naintbc = iabintc
      endif



c
c   unsort iaintb, and iaintbc
c
      do i=1,naintb
        iaintba(i) = iasort(iaintba(i))
        iaintbb(i) = ibsort(iaintbb(i))
      enddo

      do i=1,naintbc
        iaintbc(i) = iasort(iaintbc(i))
      enddo


      return
      end
c
c
c
c
c
c
c
c
      subroutine get_iuni1(n,a,nuni,iuni,iuniind)
c
c        given an array a, find unique set of indices
c        in the array and a mapping from where the indices
c        are contained in the shorter array
c
      implicit integer *8 (i-n)
      integer *8 n, a(n),nuni,iuni(n),iuniind(n)
      integer *8, allocatable :: iind(:),asort(:),iunisort(:)

      if(n.eq.0) then
        nuni = 0
        return
      endif

      allocate(iind(n),asort(n),iunisort(n))
      call sorti(n,a,iind)

      do i=1,n
        asort(i) = a(iind(i))
      enddo

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

      do i=1,n
        iuniind(iind(i)) = iunisort(i)
      enddo

      return
      end





      subroutine get_iuni1_omp(n,a,nuni,iuni,iuniind)
c
c        given an array a, find unique set of indices
c        in the array and a mapping from where the indices
c        are contained in the shorter array
c
      implicit integer *8 (i-n)
      integer *8 n, a(n),nuni,iuni(n),iuniind(n)
      integer *8, allocatable :: iind(:),asort(:),iunisort(:)

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
      implicit integer *8 (i-n)
      integer *8 n, a(n),b(n),c(n),nuni,iuni(3,n),iuniind(n)
      integer *8, allocatable :: iind(:),asort(:),iunisort(:)
      integer *8, allocatable :: bsort(:),csort(:)

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

