c
c
c   this subroutine computes the set subtraction
c      between two sets of integers   
c

      subroutine setsub(a,n,b,m,amb,namb,bma,nbma)
c
c       find set subtraction for two integer arrays 
c        with unique elements
c
c      input:
c       a - integer(n)
c         integer array 1
c       n - integer
c         number of input elements for array 1
c       b - integer(m)
c         integer array 2
c       m - number of input elements for array 2
c
c      output
c       amb - integer(n) 
c         a \ (a \cap b)^{c} elements which are in a but not in b
c       namb - integer
c         number of elements in amb
c       bma - integer(m)
c         b \ (a \cap b)^{c} elements which are in b but not in a
c       nbma - integer
c         number of elements in bma
c

      
      implicit real *8 (a-h,o-z)
      integer a(n),b(m),amb(n),bma(m)
      integer, allocatable :: asort(:),bsort(:),w(:),abint(:)
      integer, allocatable :: iabinta(:),iabintb(:)

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



cc      call sortanyn(asort,n,w)
cc      call sortanyn(bsort,m,w)


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
