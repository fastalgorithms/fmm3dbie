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
c              rsc_to_csc - convert a row sparse compressed data
c               structure to a column sparse compressed data 
c               structure or vice vrsa
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
c
c-------------------------------------------------------------------
      subroutine findnearslow(xyzs,ns,rads,targets,nt,row_ptr,col_ind)
c
cc      identify targets which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i)
c       This is a dense O(ns \times nt) algorithm
c
c       Calling sequence variables
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(3,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c
c       OUTPUT
c       row_ptr    out: integer(nt+1)
c                rowptr(i) is the starting point in the iflg
c                array for the list of sources
c                relevant for target i
c
c       col_ind  out: integer(nnz) (nnz is computed using mem routine)
c                col_ind(row_ptr(i):row_ptr(i+1)) is the list
c                of sources which satisfy
c                d(s_{j},t_{i}) <= rs_{j}
c               
c-------------------------------

      implicit real *8 (a-h,o-z)
      real *8 xyzs(3,ns),targets(3,nt),rads(ns)
      integer row_ptr(*),col_ind(*)
      integer, allocatable :: nlst(:)

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

      subroutine findnearslowmem(xyzs,ns,rads,targets,nt,nnz)
c
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i)
c       This is a dense O(ns \times nt) algorithm.
c
c       Calling sequence variables
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(3,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c
c       OUTPUT
c       nnz     out: integer
c               number of elements in the flag array
c-------------------------------

      implicit real *8 (a-h,o-z)
      real *8 xyzs(3,*),targets(3,*),rads(*)
      integer nnz

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
      subroutine findnearmem(xyzs,ns,rads,targets,nt,nnz)
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
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(3,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c       radt    in: real *8(ns)
c               radii associated with the targets
c
c       targets in: real *8(3,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c       OUTPUT
c       nnz     out: integer
c               number of elements in the flag array
c
c               
c-------------------------------

       implicit real *8 (a-h,o-z)
       real *8 xyzs(3,ns),rads(ns),targets(3,nt)
       real *8 xtmp(3)
       real *8, allocatable :: rstmp(:)

       integer nnz

cc       tree variables
c 
       real *8, allocatable :: centers(:,:),boxsize(:)
       integer, allocatable :: ilevel(:)
       integer, allocatable :: itree(:)
       integer *8 ipointer(32),ltree

       logical res



       allocate(rstmp(ns))



       do i=1,ns
         rstmp(i) = 2*rads(i)
       enddo


       mnbors = 27
       mnlist2 = 7*mnbors

       rttmp = 0


       idivflag = 0
       ndiv = 2
       isep = 1
       nlmax = 200
       nbmax = 0

       nlevels = 0
       nboxes = 0
       mhung = 0
       mnbors = 0
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       ntmp = 0


       call mklraptreemem(ier,xyzs,ns,rstmp,targets,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,isep,nlmax,nbmax,nlevels,nboxes,mnbors,mnlist1,
     2    mnlist2,mnlist3,mnlist4,mhung,ltree)



       allocate(centers(3,nboxes),itree(ltree),boxsize(0:nlevels))
  
       call mklraptree(xyzs,ns,rstmp,targets,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,isep,mhung,mnbors,mnlist1,mnlist2,mnlist3,
     2    mnlist4,nlevels,nboxes,centers,boxsize,
     2    itree,ltree,ipointer)


       allocate(ilevel(nboxes))

       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           ilevel(ibox) = ilev
         enddo
       enddo

       nnz = 0


       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(ipointer(3)+ibox-1)
           if(nchild.eq.0) then 
             itstart = itree(ipointer(12)+ibox-1)
             itend = itree(ipointer(13)+ibox-1)
             do itt = itstart,itend
               itarg = itree(ipointer(6)+itt-1)
               nhunglistsrc = itree(ipointer(30)+ibox-1)
               do ii=1,nhunglistsrc
                 iss = itree(ipointer(31)+mhung*(ibox-1)+ii-1)
                 rr = (targets(1,itarg)-xyzs(1,iss))**2+ 
     1                 (targets(2,itarg)-xyzs(2,iss))**2+
     2                 (targets(3,itarg)-xyzs(3,iss))**2
                 if(rr.le.rads(iss)**2) nnz = nnz + 1              
               enddo

               nlist1 = itree(ipointer(20)+ibox-1)
               do j=1,nlist1
                 jbox = itree(ipointer(21)+mnlist1*(ibox-1)+j-1)
                 isstart = itree(ipointer(10)+jbox-1)
                 isend = itree(ipointer(11)+jbox-1)
                 do ii = isstart,isend
                   iss = itree(ipointer(5)+ii-1)
                   rr = (targets(1,itarg)-xyzs(1,iss))**2+ 
     1                    (targets(2,itarg)-xyzs(2,iss))**2+
     2                    (targets(3,itarg)-xyzs(3,iss))**2
                   if(rr.le.rads(iss)**2) nnz = nnz + 1
                 enddo
               enddo
             enddo
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
c-----------------------------------------------------
      subroutine findnear(xyzs,ns,rads,targets,nt,row_ptr,
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
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(3,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c       row_ptr    out: integer(nt+1)
c                row_ptr(i) is the starting point in the col_ind
c                array for the list of sources
c                relevant for target i
c
c       col_ind     out: integer(nnz) (nnz is computed using
c                                      findnearmem)
c                col_ind(row_ptr(i):row_ptr(i+1)-1) is the list
c                of sources which satisfy
c                d(s_{j},t_{i}) <= rs_{j}
c               
c               
c-------------------------------
       implicit real *8 (a-h,o-z)
       real *8 xyzs(3,*),rads(*),targets(3,*)
       real *8, allocatable :: rstmp(:)

       integer row_ptr(*),col_ind(*)
       integer, allocatable :: nlst(:)
c
cc       tree variables
c 
       real *8, allocatable :: centers(:,:),boxsize(:)
       integer, allocatable :: itree(:)
       integer *8 ipointer(32), ltree
       integer, allocatable :: ilevel(:)

       allocate(rstmp(ns))

       do i=1,ns
         rstmp(i) = 2*rads(i)
       enddo

       rttmp = 0.0d0


       mnbors = 27
       mnlist2 = 7*mnbors

       idivflag = 0
       ndiv = 2
       isep = 1
       nlmax = 200
       nbmax = 0

       nlevels = 0
       nboxes = 0
       mhung = 0
       ntmp = 0
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       mnbors = 0

       call mklraptreemem(ier,xyzs,ns,rstmp,targets,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,isep,nlmax,nbmax,nlevels,nboxes,mnbors,mnlist1,
     2    mnlist2,mnlist3,mnlist4,mhung,ltree)

       allocate(centers(3,nboxes),itree(ltree),boxsize(0:nlevels))

       
       call mklraptree(xyzs,ns,rstmp,targets,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,isep,mhung,mnbors,mnlist1,mnlist2,mnlist3,
     2    mnlist4,nlevels,nboxes,centers,boxsize,
     2    itree,ltree,ipointer)

       allocate(ilevel(nboxes))

       do ilev=0,nlevels
          do ibox=itree(2*ilev+1),itree(2*ilev+2)
              ilevel(ibox) = ilev
          enddo
       enddo

       allocate(nlst(nt))


       do i=1,nt
          nlst(i) = 0
       enddo

       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(ipointer(3)+ibox-1)
           if(nchild.eq.0) then
             itstart = itree(ipointer(12)+ibox-1)
             itend = itree(ipointer(13)+ibox-1)

             do it = itstart,itend
               itarg = itree(ipointer(6)+it-1)

               nhunglistsrc = itree(ipointer(30)+ibox-1)
               do ii=1,nhunglistsrc
                 is = itree(ipointer(31)+mhung*(ibox-1)+ii-1)
                 rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                 (targets(2,itarg)-xyzs(2,is))**2+
     2                 (targets(3,itarg)-xyzs(3,is))**2
                 if(rr.le.rads(is)**2) 
     1               nlst(itarg) = nlst(itarg) + 1              
               enddo

               nlist1 = itree(ipointer(20)+ibox-1)
               do j=1,nlist1
                 jbox = itree(ipointer(21)+mnlist1*(ibox-1)+j-1)
                 isstart = itree(ipointer(10)+jbox-1)
                 isend = itree(ipointer(11)+jbox-1)
                 do ii = isstart,isend
                   is = itree(ipointer(5)+ii-1)
                   rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                     (targets(2,itarg)-xyzs(2,is))**2+
     2                     (targets(3,itarg)-xyzs(3,is))**2
                   if(rr.le.rads(is)**2) 
     1                    nlst(itarg) = nlst(itarg) + 1
                 enddo
               enddo
             enddo
           endif
         enddo
       enddo

       row_ptr(1) = 1
       do i=1,nt
          row_ptr(i+1) = row_ptr(i) + nlst(i)
          nlst(i) = 0
       enddo



       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(ipointer(3)+ibox-1)
           if(nchild.eq.0) then
             itstart = itree(ipointer(12)+ibox-1)
             itend = itree(ipointer(13)+ibox-1)

             do it = itstart,itend
               itarg = itree(ipointer(6)+it-1)
               nhunglistsrc = itree(ipointer(30)+ibox-1)
               do ii=1,nhunglistsrc
                 is = itree(ipointer(31)+mhung*(ibox-1)+ii-1)
                 rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                 (targets(2,itarg)-xyzs(2,is))**2+
     2                 (targets(3,itarg)-xyzs(3,is))**2
                 if(rr.le.rads(is)**2) then
                   col_ind(row_ptr(itarg)+nlst(itarg)) =is
                   nlst(itarg) = nlst(itarg)+1
                 endif
               enddo

               nlist1 = itree(ipointer(20)+ibox-1)
               do j=1,nlist1
                 jbox = itree(ipointer(21)+mnlist1*(ibox-1)+j-1)
                 isstart = itree(ipointer(10)+jbox-1)
                 isend = itree(ipointer(11)+jbox-1)
                 do ii = isstart,isend
                   is = itree(ipointer(5)+ii-1)
                   rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                  (targets(2,itarg)-xyzs(2,is))**2+
     2                  (targets(3,itarg)-xyzs(3,is))**2
                   if(rr.le.rads(is)**2) then 
                     col_ind(row_ptr(itarg)+nlst(itarg)) =is
                     nlst(itarg) = nlst(itarg)+1
                   endif
                 enddo
               enddo
             enddo
           endif
         enddo
       enddo


       return
       end
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
c         npatches - integer
c           number of patches
c         ixyzs - integer(npatches+1)
c           location in source array where data for patch i starts
c            (array for determining number of points per patch)
c         npts - integer
c           number of targets
c         nnz - integer
c           number of near-field interactions
c         row_ptr - integer(npts+1)
c           pointer in col_ind array for interactions corresponding
c           to target i
c         col_ind - integer(nnz)
c            list of source patches corresponding to near field
c            interactions
c       
c
c       output
c         iquad - integer(nnz)
c            location in quadrature array where quadrature
c            corresponding to interaction of col_ind(i) starts
c            
c
      implicit none
      integer npatches,ixyzs(npatches+1),npts,nnz
      integer row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      integer i,ipatch,npols
      integer, allocatable :: iqtmp(:)
      
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
      integer norder,iptype
      real *8 rfac,rfac0

      rfac = 1.25d0
      rfac0 = 1.25d0


      if(iptype.eq.1) then

        if(norder.le.3) rfac = 2.75d0
        if(norder.le.6.and.norder.gt.3) rfac = 2.0d0
        if(norder.gt.6) rfac = 1.25d0
        rfac0 = 1.25d0
      endif

      

      return
      end

c
c
c
c
c
c
c

      subroutine get_quadparams_adap(eps,nqorder,eps_adap,nlev,
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
c        outputs:
c          nqorder - integer
c             order of XG nodes to use on each triangle in the
c             adaptive integration strategy
c          eps_adap - real *8
c             stopping criterion for adaptive integration
c          nlev - integer
c             number of uniform levels for using oversampled 
c             quadrature in the near-field
c          nqorder_f - order of XG nodes to use in each of
c            triangles when using oversampled quadrature in
c            the near field
c
      implicit none
      real *8 eps,eps_adap,eps0
      integer norder,nqorder,i,iprec
      integer nlev,nqorder_f
      

      iprec = 0
      if(eps.lt.0.5d-2) iprec = 1
      if(eps.lt.0.5d-3) iprec = 2
      if(eps.lt.0.5d-6) iprec = 3
      if(eps.lt.0.5d-9) iprec = 4


      nqorder = 10
      eps_adap = eps
      
      nlev = 1

      if(iprec.eq.0) then
        nqorder = 6
        eps_adap = 0.9d-2
        nqorder_f = 4
      endif

      if(iprec.eq.1) then
        nqorder= 7
        eps_adap = 0.5d-2
        nqorder_f = 4
      endif

      if(iprec.eq.2) then
        nqorder = 12
        eps_adap = 0.5d-4
        nqorder_f = 7
      endif

      if(iprec.eq.3) then
        nqorder = 25
        eps_adap = 3.0d-7
        nqorder_f = 12
      endif

      if(iprec.eq.4) then
        nqorder = 30
        eps_adap = 3.0d-10
        nqorder_f = 15
      endif
      
      return
      end
c
c



      
