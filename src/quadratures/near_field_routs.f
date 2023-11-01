c              this file contains the near field routines 
c              to be used for quadrature generation
c
c              findnear - find near field of a collection
c                         of sources with extent and targets
c 
c              findnearmem - memory management routine for find
c                            near
c
c              O(n^2) analogs of findnear and findnearmem
c
c              findnearslow - find near field of a collection
c                         of sources with extent and targets
c 
c              findnearslowmem - memory management routine for find
c                            near
c
c              get_iquad_rsc - given a row sparse compressed
c                  format for quadrature from targets
c                  to patches, and order of discretization
c                  of patches, determine a pointer
c                  array for where quadrature for an interaction
c                  starts
c
c              get_rfacs - get default 
c                factors for defining the near-field
c                and region inside near-field for using
c                adaptive integration as a function of order
c
c              get_quadparams_adap - get various parameters
c                when computing nearly-singular integrals
c                using adaptive integration
c
c              get_near_corr_max - estimate max size
c                of array needed for near correction
c
c-------------------------------------------------------------------
      subroutine findnearslow(xyzs,ns,rads,ndt,targets,nt,row_ptr,
     1  col_ind)
c
cc      identify targets which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i)
c       This is a dense O(ns \times nt) algorithm
c
c       Calling sequence variables
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer(8)
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       ndt     in: integer(8)
c               leading dimension of targets array
c
c       targets in: real *8(ndt,nt)
c               target info.
c               The first three components must be xyz coordinates
c
c       nt      in: integer(8)
c               number of targets
c
c
c       OUTPUT
c       row_ptr    out: integer(8)(nt+1)
c                rowptr(i) is the starting point in the iflg
c                array for the list of sources
c                relevant for target i
c
c       col_ind  out: integer(8)(nnz) (nnz is computed using mem routine)
c                col_ind(row_ptr(i):row_ptr(i+1)) is the list
c                of sources which satisfy
c                d(s_{j},t_{i}) <= rs_{j}
c               
c-------------------------------

      implicit none
      integer(8) ns,ndt,nt,i,j
      real *8 xyzs(3,ns),targets(ndt,nt),rads(ns),rr
      integer(8) row_ptr(*),col_ind(*)
      integer(8), allocatable :: nlst(:)

      allocate(nlst(nt))

      do i=1,nt
        nlst(i) = 0
        do j=1,ns
          rr = (xyzs(1,j)-targets(1,i))**2 + 
     1          (xyzs(2,j)-targets(2,i))**2 +
     2          (xyzs(3,j)-targets(3,i))**2
          if(rr.le.rads(j)**2) then
            nlst(i) = nlst(i) + 1
          endif
        enddo
      enddo

      row_ptr(1) = 1
      do i=1,nt
        row_ptr(i+1) = row_ptr(i)+nlst(i)
      enddo

      do i=1,nt
        nlst(i) = 0
      enddo

      do i=1,nt
        do j=1,ns
          rr = (xyzs(1,j)-targets(1,i))**2 + 
     1         (xyzs(2,j)-targets(2,i))**2 +
     2        (xyzs(3,j)-targets(3,i))**2

      
          if(rr.le.rads(j)**2) then
            col_ind(row_ptr(i)+nlst(i)) = j
            nlst(i) = nlst(i) + 1
          endif
        enddo
      enddo



      return
      end
c
c
c
c
c
c------------------------------------------------------------------      

      subroutine findnearslowmem(xyzs,ns,rads,ndt,targets,nt,nnz)
c
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i)
c       This is a dense O(ns \times nt) algorithm.
c
c       Calling sequence variables
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer(8)
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c      
c       ndt     in: integer(8)
c               leading dimension of target array
c
c       targets in: real *8(ndt,nt)
c               target info
c               the first three components must be
c               the xyz coordinates
c
c       nt      in: integer(8)
c               number of targets
c
c
c       OUTPUT
c       nnz     out: integer(8)
c               number of elements in the flag array
c-------------------------------

      implicit real *8 (a-h,o-z)
      implicit integer(8) (i-n)
      real *8 xyzs(3,*),targets(ndt,*),rads(*)
      integer(8) nnz

      nnz = 0

      do i=1,nt
        do j=1,ns
          rr = (xyzs(1,j)-targets(1,i))**2 + 
     1         (xyzs(2,j)-targets(2,i))**2 +
     2         (xyzs(3,j)-targets(3,i))**2

          if(rr.le.rads(j)**2) then
            nnz = nnz+1
          endif
        enddo
      enddo

      return
      end
c
c
c
c
c
c------------------------------------------------------------------      
      subroutine findnearmem(xyzs,ns,rads,ndt,targets,nt,nnz)
c
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i).
c       We use hung lists to identify a nearly minimal set of
c       sources to loop over relevant for each target.
c
c       This is a memory management routine for findnearmem
c       
c
c       Calling sequence variables
c
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer(8)
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       ndt     in: integer(8)
c               leading dimension of target array
c
c       targets in: real *8(ndt,nt)
c               target info
c               the first three components must be
c               the xyz coordinates
c
c       nt      in: integer(8)
c               number of targets
c
c       radt    in: real *8(ns)
c               radii associated with the targets
c
c       targets in: real *8(3,nt)
c               location of the targets
c
c       nt      in: integer(8)
c               number of targets
c
c       OUTPUT
c       nnz     out: integer(8)
c               number of elements in the flag array
c
c               
c-------------------------------

       implicit real *8 (a-h,o-z)
       implicit integer(8) (i-n)
       real *8 xyzs(3,*),rads(*),targets(ndt,*)
       real *8, allocatable :: rstmp(:)

       integer(8), allocatable :: nlst(:)
c
cc       tree variables
c 
       real *8, allocatable :: centers(:,:),boxsize(:),targs(:,:)
       integer(8), allocatable :: itree(:)
       integer(8) ipointer(8), ltree
       integer(8), allocatable :: ilevel(:)
       integer(8), allocatable :: nlist1(:),list1(:,:)
       integer(8), allocatable :: nlist2(:),list2(:,:)
       integer(8), allocatable :: nlist3(:),list3(:,:)
       integer(8), allocatable :: nlist4(:),list4(:,:)
       integer(8) mnlist1,mnlist2,mnlist3,mnlist4
       integer(8), allocatable :: isrcper(:),isrcse(:,:)
       integer(8), allocatable :: itargper(:),itargse(:,:)

       allocate(rstmp(ns))

C$OMP PARALLEL DO DEFAULT(SHARED)
       do i=1,ns
         rstmp(i) = 2*rads(i)
       enddo
C$OMP END PARALLEL DO       
       allocate(targs(3,nt))

C$OMP PARALLEL DO DEFAULT(SHARED)       
       do i=1,nt
         targs(1,i) = targets(1,i)
         targs(2,i) = targets(2,i)
         targs(3,i) = targets(3,i)
       enddo
C$OMP END PARALLEL DO       

       idivflag = 0
       ndiv = 2
       isep = 1
       nlmax = 51
       nlmin = 0
       iper = 0
       ifunif = 0
       nbmax = 0

       nlevels = 0
       nboxes = 0
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       mnbors = 27

       call pts_tree_mem(xyzs,ns,targs,nt,idivflag,ndiv,nlmin,nlmax,
     1  iper,ifunif,nlevels,nboxes,ltree) 

       allocate(centers(3,nboxes),itree(ltree),boxsize(0:nlevels))

       call pts_tree_build(xyzs,ns,targs,nt,idivflag,ndiv,nlmin,
     1   nlmax,iper,ifunif,nlevels,nboxes,ltree,itree,ipointer,centers,
     2   boxsize)
       
       allocate(isrcse(2,nboxes),itargse(2,nboxes))
       allocate(isrcper(ns),itargper(nt))

       call pts_tree_sort(ns,xyzs,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,isrcper,isrcse)
      
       call pts_tree_sort(nt,targs,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,itargper,itargse)

c
c   initialize various tree lists
c
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
      mnbors = 27

      isep = 1
      
      call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     2  itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))
      allocate(list2(mnlist2,nboxes),nlist2(nboxes))
      allocate(list3(mnlist3,nboxes),nlist3(nboxes))
      allocate(list4(mnlist4,nboxes),nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     3  itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2,
     4  mnlist2,list2,nlist3,mnlist3,list3,
     4  nlist4,mnlist4,list4)

      nnz = 0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,nchild,ilev,iss)
C$OMP$PRIVATE(isrc,rtest,ilevup,jbox,i,ncoll,l,lbox,xdis,ydis,zdis)
C$OMP$PRIVATE(llev,distest,itt,itarg,rr,jlev) REDUCTION(+:nnz)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox = 1,nboxes
        nchild = itree(ipointer(4)+ibox-1)
        ilev = itree(ipointer(2)+ibox-1)
        if(nchild.eq.0) then
          do iss=isrcse(1,ibox),isrcse(2,ibox)
             isrc = isrcper(iss)
             rtest = rads(isrc)**2
             ilevup = ceiling(log(rstmp(isrc)/boxsize(ilev))/log(2.0d0))
             jbox = ibox
             ilevup = min(ilevup,ilev)

             do i=1,ilevup
               jbox = itree(ipointer(3)+jbox-1)
             enddo
c
c   loop over colleagues
c
             ncoll = itree(ipointer(6)+jbox-1)
             do l=1,ncoll
               lbox = itree(ipointer(7)+mnbors*(jbox-1)+l-1)
c
c  check if lbox does not intersect with region of interest
c 
               xdis = abs(centers(1,lbox)-xyzs(1,isrc))
               ydis = abs(centers(2,lbox)-xyzs(2,isrc))
               zdis = abs(centers(3,lbox)-xyzs(3,isrc))

               llev = itree(ipointer(2)+lbox-1)
               distest = rads(isrc)+boxsize(llev)/2
               if((xdis.le.distest).and.(ydis.le.distest).and.
     1             (zdis.le.distest)) then
                 
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     nnz = nnz + 1
                   endif
                 enddo
               endif
             enddo

c
c   loop over list1 boxes which are larger
c
             jlev = itree(ipointer(2)+jbox-1)
             do l=1,nlist1(jbox)
               lbox = list1(l,jbox)
               llev = itree(ipointer(2)+lbox-1)
               if(llev.lt.jlev) then
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     nnz = nnz + 1
                   endif
                 enddo
               endif
             enddo
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      


       return
       end
c
c
c
c
c
c-----------------------------------------------------
      subroutine findnear(xyzs,ns,rads,ndt,targets,nt,row_ptr,
     1       col_ind) 
c     
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i).
c       We use hung lists to identify a nearly minimal set of
c       sources to loop over relevant for each target.
c       
c
c       Calling sequence variables
c
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer(8)
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       ndt     in: integer(8)
c               leading dimension of target array
c
c       targets in: real *8(ndt,nt)
c               target info
c               the first three components must be
c               the xyz coordinates
c
c       nt      in: integer(8)
c               number of targets
c
c       row_ptr    out: integer(8)(nt+1)
c                row_ptr(i) is the starting point in the col_ind
c                array for the list of sources
c                relevant for target i
c
c       col_ind     out: integer(8)(nnz) (nnz is computed using
c                                      findnearmem)
c                col_ind(row_ptr(i):row_ptr(i+1)-1) is the list
c                of sources which satisfy
c                d(s_{j},t_{i}) <= rs_{j}
c               
c               
c-------------------------------
       implicit real *8 (a-h,o-z)
       implicit integer(8) (i-n)
       real *8 xyzs(3,*),rads(*),targets(ndt,*)
       real *8, allocatable :: rstmp(:)

       integer(8) row_ptr(*),col_ind(*)
       integer(8), allocatable :: nlst(:)
c
cc       tree variables
c 
       real *8, allocatable :: centers(:,:),boxsize(:),targs(:,:)
       integer(8), allocatable :: itree(:)
       integer(8) ipointer(8), ltree
       integer(8), allocatable :: ilevel(:)
       integer(8), allocatable :: nlist1(:),list1(:,:)
       integer(8), allocatable :: nlist2(:),list2(:,:)
       integer(8), allocatable :: nlist3(:),list3(:,:)
       integer(8), allocatable :: nlist4(:),list4(:,:)
       integer(8) mnlist1,mnlist2,mnlist3,mnlist4
       integer(8), allocatable :: isrcper(:),isrcse(:,:)
       integer(8), allocatable :: itargper(:),itargse(:,:)
       integer(8), allocatable :: nlistsrc(:)
       integer(8), allocatable :: col_ind2(:),row_ind2(:),col_ptr(:)

       allocate(rstmp(ns))

C$OMP PARALLEL DO DEFAULT(SHARED)
       do i=1,ns
         rstmp(i) = 2*rads(i)
       enddo
C$OMP END PARALLEL DO       
       allocate(targs(3,nt))

C$OMP PARALLEL DO DEFAULT(SHARED)       
       do i=1,nt
         targs(1,i) = targets(1,i)
         targs(2,i) = targets(2,i)
         targs(3,i) = targets(3,i)
       enddo
C$OMP END PARALLEL DO       

       idivflag = 0
       ndiv = 2
       isep = 1
       nlmax = 51
       nbmax = 0

       nlevels = 0
       nboxes = 0
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       mnbors = 27
       ifunif = 0
       nlmin = 0
       iper = 0

       call pts_tree_mem(xyzs,ns,targs,nt,idivflag,ndiv,nlmin,
     1  nlmax,iper,ifunif,nlevels,nboxes,ltree) 

       allocate(centers(3,nboxes),itree(ltree),boxsize(0:nlevels))

       call pts_tree_build(xyzs,ns,targs,nt,idivflag,ndiv,nlmin,
     1   nlmax,iper,ifunif,nlevels,nboxes,ltree,itree,ipointer,centers,
     2   boxsize)
       
       allocate(isrcse(2,nboxes),itargse(2,nboxes))
       allocate(isrcper(ns),itargper(nt))

       call pts_tree_sort(ns,xyzs,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,isrcper,isrcse)
      
       call pts_tree_sort(nt,targs,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,itargper,itargse)

c
c   initialize various tree lists
c
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
      mnbors = 27

      isep = 1
      
      call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     2  itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))
      allocate(list2(mnlist2,nboxes),nlist2(nboxes))
      allocate(list3(mnlist3,nboxes),nlist3(nboxes))
      allocate(list4(mnlist4,nboxes),nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     3  itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2,
     4  mnlist2,list2,nlist3,mnlist3,list3,
     4  nlist4,mnlist4,list4)

c
c   estimate list of tethered sources
c
c
      allocate(nlistsrc(ns),nlst(nt))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,ns
        nlistsrc(i) = 0
      enddo
C$OMP END PARALLEL DO      

      nnz = 0
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,nchild,ilev,iss)
C$OMP$PRIVATE(isrc,rtest,ilevup,jbox,i,ncoll,l,lbox,xdis,ydis,zdis)
C$OMP$PRIVATE(llev,distest,itt,itarg,rr,jlev) REDUCTION(+:nnz)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox = 1,nboxes
        nchild = itree(ipointer(4)+ibox-1)
        ilev = itree(ipointer(2)+ibox-1)
        if(nchild.eq.0) then
          do iss=isrcse(1,ibox),isrcse(2,ibox)
             isrc = isrcper(iss)
             rtest = rads(isrc)**2
             ilevup = ceiling(log(rstmp(isrc)/boxsize(ilev))/log(2.0d0))
             jbox = ibox
             ilevup = min(ilevup,ilev)

             do i=1,ilevup
               jbox = itree(ipointer(3)+jbox-1)
             enddo
c
c   loop over colleagues
c
             ncoll = itree(ipointer(6)+jbox-1)
             do l=1,ncoll
               lbox = itree(ipointer(7)+mnbors*(jbox-1)+l-1)
c
c  check if lbox does not intersect with region of interest
c 
               xdis = abs(centers(1,lbox)-xyzs(1,isrc))
               ydis = abs(centers(2,lbox)-xyzs(2,isrc))
               zdis = abs(centers(3,lbox)-xyzs(3,isrc))

               llev = itree(ipointer(2)+lbox-1)
               distest = rads(isrc)+boxsize(llev)/2
               if((xdis.le.distest).and.(ydis.le.distest).and.
     1             (zdis.le.distest)) then
                 
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     nlistsrc(isrc) = nlistsrc(isrc)+1
                     nnz = nnz + 1
                   endif
                 enddo
               endif
             enddo

c
c   loop over list1 boxes which are larger
c
             jlev = itree(ipointer(2)+jbox-1)
             do l=1,nlist1(jbox)
               lbox = list1(l,jbox)
               llev = itree(ipointer(2)+lbox-1)
               if(llev.lt.jlev) then
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     nlistsrc(isrc) = nlistsrc(isrc)+1
                     nnz = nnz + 1
                   endif
                 enddo
               endif
             enddo
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      allocate(col_ind2(nnz),row_ind2(nnz))
      allocate(col_ptr(ns+1))

      call cumsum(ns,nlistsrc,col_ptr(2))
      col_ptr(1) = 1
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=2,ns+1
        col_ptr(i) = col_ptr(i) + 1
      enddo
C$OMP END PARALLEL DO     

      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,nchild,ilev,iss)
C$OMP$PRIVATE(isrc,rtest,ilevup,jbox,i,ncoll,l,lbox,xdis,ydis,zdis)
C$OMP$PRIVATE(llev,distest,itt,itarg,rr,jlev,nsloc) 
C$OMP$SCHEDULE(DYNAMIC)
      do ibox = 1,nboxes
        nchild = itree(ipointer(4)+ibox-1)
        ilev = itree(ipointer(2)+ibox-1)
        if(nchild.eq.0) then
          do iss=isrcse(1,ibox),isrcse(2,ibox)
             isrc = isrcper(iss)
             nsloc = 0
             rtest = rads(isrc)**2
             ilevup = ceiling(log(rstmp(isrc)/boxsize(ilev))/log(2.0d0))
             jbox = ibox
             ilevup = min(ilevup,ilev)

             do i=1,ilevup
               jbox = itree(ipointer(3)+jbox-1)
             enddo
c
c   loop over colleagues
c
             ncoll = itree(ipointer(6)+jbox-1)
             do l=1,ncoll
               lbox = itree(ipointer(7)+mnbors*(jbox-1)+l-1)
c
c  check if lbox does not intersect with region of interest
c 
               xdis = abs(centers(1,lbox)-xyzs(1,isrc))
               ydis = abs(centers(2,lbox)-xyzs(2,isrc))
               zdis = abs(centers(3,lbox)-xyzs(3,isrc))

               llev = itree(ipointer(2)+lbox-1)
               distest = rads(isrc)+boxsize(llev)/2
               if((xdis.le.distest).and.(ydis.le.distest).and.
     1             (zdis.le.distest)) then
                 
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     col_ind2(col_ptr(isrc)+nsloc) = isrc
                     row_ind2(col_ptr(isrc)+nsloc) = itarg
                     nsloc = nsloc + 1
                   endif
                 enddo
               endif
             enddo

c
c   loop over list1 boxes which are larger
c
             jlev = itree(ipointer(2)+jbox-1)
             do l=1,nlist1(jbox)
               lbox = list1(l,jbox)
               llev = itree(ipointer(2)+lbox-1)
               if(llev.lt.jlev) then
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     col_ind2(col_ptr(isrc)+nsloc) = isrc
                     row_ind2(col_ptr(isrc)+nsloc) = itarg
                     nsloc = nsloc + 1
                   endif
                 enddo
               endif
             enddo
          enddo
        endif
      enddo
C$OMP END PARALLEL DO     
      

      if(ns.lt.1000.or.nnz.lt.10000) then
        call conv_to_csc(nnz,nt,col_ind2,row_ind2,row_ptr,col_ind)
      else
        call conv_to_csc(nnz,nt,col_ind2,row_ind2,row_ptr,col_ind)
      endif


      return
      end
c
c
c
c
c
c
c-----------------------------------------------------------
      subroutine get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,
     1    col_ind,iquad)
c
c       given an row sparse compressed format from targets
c       to patches, and number of discretization nodes on patches
c       determine location in quadrature array where particular
c       interactions starts
c
c       input:
c         npatches - integer(8)
c           number of patches
c         ixyzs - integer(8)(npatches+1)
c           location in source array where data for patch i starts
c            (array for determining number of points per patch)
c         npts - integer(8)
c           number of targets
c         nnz - integer(8)
c           number of near-field interactions
c         row_ptr - integer(8)(npts+1)
c           pointer in col_ind array for interactions corresponding
c           to target i
c         col_ind - integer(8)(nnz)
c            list of source patches corresponding to near field
c            interactions
c       
c
c       output
c         iquad - integer(8)(nnz)
c            location in quadrature array where quadrature
c            corresponding to interaction of col_ind(i) starts
c            
c
      implicit none
      integer(8) npatches,ixyzs(npatches+1),npts,nnz
      integer(8) row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      integer(8) i,ipatch,npols
      integer(8), allocatable :: iqtmp(:)
      
      allocate(iqtmp(nnz))

      do i=1,nnz
        ipatch = col_ind(i)
        npols = ixyzs(ipatch+1)-ixyzs(ipatch)
        iqtmp(i) = npols
      enddo

      iquad(1) = 1
      call cumsum(nnz,iqtmp,iquad(2))
      do i=2,nnz+1
        iquad(i) = iquad(i) + 1
      enddo

      

      return
      end
c
c
c
c
c
      subroutine get_rfacs(norder,iptype,rfac,rfac0)
c
c
c       this subroutine gets the factors for defining
c       the two regions of the near field as function
c       of patch type and order
c
c       If h is the radius of the bounding sphere,
c       centered at the centroid of the patch,
c       the near-field for storing quadrature
c       should be defined by h*rfac,
c       the subsection for which adaptive integration
c       should be used is defined by h*rfac0 
c
c        Note that these are recommended parameters 
c        based on empirical testing,
c        and should be appropriately used when
c        calling findnear (rfac), and in the 
c        getnearquad routines (rfac0)
c
c      input:
c        norder - order of discretization
c        iptype - type of patch
c      output:
c        rfac - factor for defining near field for precomputing 
c               quadrature
c        rfac0 - factor for defining near field for using
c              adaptive integration 
c

      implicit none
      integer(8) norder,iptype
      real *8 rfac,rfac0

      rfac = 1.25d0
      rfac0 = 1.25d0



      if(norder.le.2) rfac = 2.75d0
      if(norder.le.6.and.norder.gt.2) rfac = 2.0d0
      if(norder.gt.6) rfac = 1.25d0
      rfac0 = 1.25d0

      

      return
      end

c
c
c
c
c
c
c

      subroutine get_quadparams_adap(eps,iptype,nqorder,eps_adap,nlev,
     1   nqorder_f)
c
c
c
c         This subroutine returns the quadrature parameters
c         for computing integrals using adaptive integration
c
c        input:
c          eps - real *8
c             tolerance
c          iptype - integer(8)
c             patch type
c        outputs:
c          nqorder - integer(8)
c             order of XG nodes to use on each triangle in the
c             adaptive integration strategy
c          eps_adap - real *8
c             stopping criterion for adaptive integration
c          nlev - integer(8)
c             number of uniform levels for using oversampled 
c             quadrature in the near-field
c          nqorder_f - order of XG nodes to use in each of
c            triangles when using oversampled quadrature in
c            the near field
c
      implicit none
      real *8 eps,eps_adap,eps0
      integer(8) norder,nqorder,i,iprec
      integer(8) nlev,nqorder_f,iptype
      

      iprec = 0
      if(eps.lt.0.5d-2) iprec = 1
      if(eps.lt.0.5d-3) iprec = 2
      if(eps.lt.0.5d-6) iprec = 3
      if(eps.lt.0.5d-9) iprec = 4


      nqorder = 10
      eps_adap = eps
      nqorder_f = 8 
      nlev = 1

      if(norder.ge.4) then
        if(iprec.eq.0) then
          nqorder = 12
          eps_adap = 1.0d-2
          nqorder_f = 8
        endif

        if(iprec.eq.1) then
          nqorder= 12
          eps_adap = 0.2d-2
          nqorder_f = 10
        endif

        if(iprec.eq.2) then
          nqorder = 16
          eps_adap = 0.1d-4
          nqorder_f = 16
        endif

        if(iprec.eq.3) then
          nqorder = 22
          eps_adap = 7.0d-8
          nqorder_f = 22
        endif

        if(iprec.eq.4) then
          nqorder = 30
          eps_adap = 3.0d-10
          nqorder_f = 30
        endif
      else
        if(iprec.eq.0) then
          nqorder = 8
          eps_adap = 1.0d-2
          nqorder_f = 4
        endif

        if(iprec.eq.1) then
          nqorder= 8
          eps_adap = 0.2d-2
          nqorder_f = 6
        endif

        if(iprec.eq.2) then
          nqorder = 12
          eps_adap = 0.1d-4
          nqorder_f = 12
        endif

        if(iprec.eq.3) then
          nqorder = 16
          eps_adap = 3.0d-8
          nqorder_f = 20
        endif

        if(iprec.eq.4) then
          nqorder = 24
          eps_adap = 2.0d-10
          nqorder_f = 26
        endif
      endif



      return
      end
c
c
c
c
c
      subroutine get_near_corr_max(ntarg,row_ptr,nnz,col_ind,npatches,
     1   ixyzs,nmax)
c------------------------------------
c
c  This subroutine estimates the max size of sources
c  (in terms of discretization points) contained
c  in the near field of a given target.
c
c  Note that one could easily call this routine by the oversampled
c  sources by sending in ixyzs corresponding to oversampled
c  discretization
c
c
c  Input arguments:
c  - ntarg: integer(8)
c      number of targets
c  - row_ptr: integer(8)(ntarg+1)
c      rowptr(i) is the starting point in the col_ind
c      array for the list of sources relevant for target i
c  - nnz: integer(8)
c      number of entries in the col_ind array
c  - col_ind: integer(8)(nnz) 
c      col_ind(row_ptr(i):row_ptr(i+1)) is the list
c      of sources in the near field of target i
c  - npatches: integer(8)
c      number of patches in the surface discretization
c  - ixyzs: integer(8)(npatches+1)
c      ixyzs(i+1)-ixyzs(i) is the number of discretization
c      nodes on patch i
c  
c  Output arguments:
c  
c  - nmax: integer(8)
c      max number of sources in the near field of any target
c               
c-------------------------------
      implicit none
      integer(8), intent(in) :: npatches,nnz,ntarg
      integer(8), intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer(8), intent(in) :: ixyzs(npatches+1)
      integer(8), intent(out) :: nmax

      integer(8) i,j,jpatch,ntmp

      nmax = 0

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,ntmp)
C$OMP$REDUCTION(max:nmax) 
      do i=1,ntarg
        ntmp = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          ntmp = ntmp + ixyzs(jpatch+1)-ixyzs(jpatch)
        enddo
        if(ntmp.gt.nmax) nmax = ntmp
      enddo
C$OMP END PARALLEL DO

      return
      end


      
