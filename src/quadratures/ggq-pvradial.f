

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the
c       quadrature code.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The following subroutines are user-callable:
c
c   radial_init - read a precomputed table of quadraturs from a
c       text file on the disk.
c
c   raddiag - return a quadrature for integrating a function of the
c        form (1) over a triangle containing the origin.
c
c   radquad - return a quadrature for integrating a function of the
c        form (1) over a triangle
c
c
c       JOINT UTILITY ROUTINE FOR DIVIDING THE PROBLEM??
c
c
c       FETCH PVQUAD IN BEGINNING OF RADDIAG???
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine pv_radial_init_mem(norder,lkeep)
        implicit real *8 (a-h,o-z)
        character *56 filename
        dimension rs(2,1000),ts(2,1000)

 0050 format ("../../src/quadratures/ggq-self-quads/pvradquads",
     1  I2.2,".txt")
 0100 format (I3.3)
 0200 format (D44.36)
c
        ier   = 0
        lkeep = 0
c
        max   = 1000
c
        write(filename,0050) norder
c
        iw = 101
        open(iw,FILE=filename,STATUS='old')

c
c       Grab the header data.
c
        read(iw,0100) nr
        read(iw,0100) nt
c
c       Read the parameter intervals.
c
        read (iw,0200) (rs(1,j),rs(2,j),j=1,nr)
        read (iw,0200) (ts(1,j),ts(2,j),j=1,nt)
c
c       Fill out the header of the rad structure.
c
        nlege = (norder)/2+1
c
        irs = 100
        lrs = 2*nr
c
        its = irs+lrs
        lts = 2*nt
c
        ixslege = its+lts
        lxslege = nlege
c
        iwhtslege = ixslege+lxslege
        lwhtslege = nlege
c
        ins = iwhtslege+lwhtslege
        lns = nr*nt

        ixs = ins+lns
        lxs = max*nt*nr
c
        iwhts = ixs+lxs
        lwhts = max*nr*nt
c
        lkeep = iwhts+lwhts
        close(iw)
c

        return
        end
c
c
c


        subroutine pv_radial_init(ier,norder,rad,lrad,lkeep)
        implicit double precision (a-h,o-z)
        dimension rad(1)
        dimension rs(2,1000),ts(2,1000)
        character*2 str
        character*56 filename
c
c       Read a table of precomputed quadrature rules from a file
c       on the disk into the user-supplied array.  The table is
c       stored in the file ``pvradquads??.txt'' where ?? is the order
c       of the quadrature rules.
c
c                            Input Parameters:
c
c   norder - the order of the quadrature rules; at the present time,
c       norder can be 4, 8, 12, 16 or 20
c   lrad - length of the user-supplied array rad which will hold the
c       quadrature table
c
c                            Output Parameters:
c
c   ier - an error return code;
c       ier = 0    indicates successful execution
c       ier = 4    mean that the user-supplied array rad is of
c                  insufficient length
c       ier = 128  means that the input file could not be opened
c       ier = 1024 means that a formatting error was encountered while
c                  attempting to read the table of radial quadratures
c                  from the disk
c
c   rad - on return, this user-supplied array will contain a quadrature
c       table and a structure header describing the quadratures
c
c
 0050 format ("../../src/quadratures/ggq-self-quads/pvradquads",
     1  I2.2,".txt")

 0100 format (I3.3)
 0200 format (D44.36)
c
        ier   = 0
        lkeep = 0
c
        max   = 1000
c
        write(filename,0050) norder
c
        iw = 101
        open(iw,FILE=filename,STATUS='old')

c
c       Grab the header data.
c
        read(iw,0100) nr
        read(iw,0100) nt
c
c       Read the parameter intervals.
c
        read (iw,0200) (rs(1,j),rs(2,j),j=1,nr)
        read (iw,0200) (ts(1,j),ts(2,j),j=1,nt)
c
c       Fill out the header of the rad structure.
c
        nlege = (norder)/2+1
c
        irs = 100
        lrs = 2*nr
c
        its = irs+lrs
        lts = 2*nt
c
        ixslege = its+lts
        lxslege = nlege
c
        iwhtslege = ixslege+lxslege
        lwhtslege = nlege
c
        ins = iwhtslege+lwhtslege
        lns = nr*nt

        ixs = ins+lns
        lxs = max*nt*nr
c
        iwhts = ixs+lxs
        lwhts = max*nr*nt
c
        lkeep = iwhts+lwhts
c
        if (lkeep .gt. lrad) then
        ier = 4
        return
        endif
c
        rad(1)  = nr
        rad(2)  = nt
        rad(3)  = irs
        rad(4)  = its
        rad(5)  = ins
        rad(6)  = ixs
        rad(7)  = iwhts
        rad(8)  = max
c
        rad(10) = nlege
        rad(11) = ixslege
        rad(12) = iwhtslege
c
        rad(21) = norder
c
c       Construct the Legendre quadrature.
c
cccc        call legequad(nlege,rad(ixslege),rad(iwhtslege))
         itype111 = 1
         call legerts(itype111, nlege, rad(ixslege), rad(iwhtslege))
c
c       Copy the as and bs into the array.
c
        call pv_radmove(nr*2,rs,rad(irs))
        call pv_radmove(nt*2,ts,rad(its))
c
c       Call the auxillary routine to read the quadrature rules.
c
        call pv_radinit0(ier,iw,max,nr,nt,rad(ins),rad(ixs),rad(iwhts),
     -    rad(irs),rad(its),norder)
        if (ier .ne. 0) return
c
        close(iw)
        return
c
c       We arrive here on IO errors.
c
 1000 continue
        ier = 128
        return
        end
c
c
c
        subroutine pv_radinit0(ier,iw,max,nr,nt,ns,xs,whts,rs,ts,norder)
        implicit double precision (a-h,o-z)
        dimension ns(nt,nr),rs(2,1),ts(2,1)
        dimension xs(max,nt,nr),whts(max,nt,nr)
        character*2 str
c
        ier = 0
 0100 format (I3.3)
 0200 format (D44.36)
c
c       Read each of the quadrature rules from the text file.
c
        do 1100 ir=1,nr
        do 1200 it=1,nt
c
c       Read the quadrature.
c
        read (iw,0100) nn
        read (iw,0200) (xs(i,it,ir),i=1,nn)
        read (iw,0200) (whts(i,it,ir),i=1,nn)
c

c$$$ 3000 format("$",I2,"$ &","$",E13.6,"$ & $",E13.6,"$ & $",E13.6,"$ & $",
c$$$     -    E13.6," $ & $",I3,"$ \\")
c$$$
c$$$        call corrand3(1,dd)
c$$$        if (dd .gt. .93d0) then
c$$$        write (*,3000) norder,rs(1,ir),rs(2,ir),ts(1,it),ts(2,it),nn
c$$$        endif
c
        ns(it,ir)=nn
c
 1200 continue
 1100 continue
c

        end



        subroutine pv_linedist(x1,y1,x2,y2,r)
        implicit double precision (a-h,o-z)
c
c       Find the point on the line segment connecting (x1,y1)
c       and (x2,y2) which is closed to the origin.
c
        dx = (x2-x1)
        dy = (y2-y1)
c
        t0 = 0
        t1 = 1
        t2 = -(x1*dx+y1*dy)/(dx**2+dy**2)
c
        if (0 .lt. t2 .AND. t2 .lt. 1) then
        r2 = (x1+dx*t2)**2+(y1+dy*t2)**2
        else
        r2 = 10
        endif

        r0 = (x1+dx*t0)**2+(y1+dy*t0)**2
        r1 = (x1+dx*t1)**2+(y1+dy*t1)**2
c
        r = min(r0,r1)
        r = min(r,r2)
        r = sqrt(r)
c
        end


        subroutine pv_raddiag(ier,rad,verts,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension rad(*),xs(*),ys(*),whts(*),verts(2,3)
        real *8, allocatable :: zs(:,:)
c
c
c       Return a quadrature for the evaluation of radially singular
c       functions of the form (1) over an arbitrary triangle
c       containing the origin.
c
c       Note that this routine will *fail* if the origin is too close
c       to the boundary; the precise tolerances depend on the
c       precomputed quadrature table.
c
c                            Input Parameters:
c
c   rad - the radial quadrature structure generated by radial_init
c   verts - a (2,3) array each column of which gives the coordinates
c       of one vertex of the user-supplied triangle
c
c                           Output Parameters:
c
c   ier - an error return code;
c       ier = 0     indicates succesful execution
c
c       ier = 64    means that one of the necessary quadrature rules is
c                   missing from the table in the rad structure; typically,
c                   this  indicates that the triangle is very close to
c                   degenerate
c
c       ier = 128   means that one of more of the vertices is *very*
c                   close to the origin
c
c       ier = 256   means the origin is *very* close to the boundary
c                   of the triangle
c
c       ier = 512   means that the vertices are nearly colinear
c
c
c   nquad - the size of the quadrature rule
c   (xs,ys) - the coordinates of the quadrature
c
c
        ier   = 0
        nquad = 0

        allocate(zs(2,10000))
c
        norder = rad(21)

c
c
c       Extract the coordinates of the vertices.
c
        x1 = verts(1,1)
        y1 = verts(2,1)
c
        x2 = verts(1,2)
        y2 = verts(2,2)
c
        x3 = verts(1,3)
        y3 = verts(2,3)
c
        call pv_linedist(x1,y1,x2,y2,dd1)
        call pv_linedist(x2,y2,x3,y3,dd2)
        call pv_linedist(x1,y1,x3,y3,dd3)
c
        rcut = min(dd1,dd2)
        rcut = min(dd3,rcut)
        rcut = rcut/2
c
c        call prin2("rcut = *",rcut,1)
c
c       Build quadratures for each subtriangle.
c
        nquad = 0
c
        call pv_raddiag0(ier,rad,x1,y1,x2,y2,nquad0,xs(nquad+1),
     -    ys(nquad+1),whts(nquad+1),rcut,norder)
        if (ier .ne. 0) return
        nquad = nquad+nquad0
c
        call pv_raddiag0(ier,rad,x1,y1,x3,y3,nquad0,xs(nquad+1),
     1    ys(nquad+1),whts(nquad+1),rcut,norder)
        if (ier .ne. 0) return
        nquad = nquad+nquad0
c
        call pv_raddiag0(ier,rad,x2,y2,x3,y3,nquad0,xs(nquad+1),
     1    ys(nquad+1),whts(nquad+1),rcut,norder)
        if (ier .ne. 0) return
        nquad = nquad+nquad0
c
        nlege = (norder/2)+1
        n     = 3*(norder+1)+2
c
        call pv_raddiag_disk0(n,nlege,rcut,
     -    nquad0,xs(nquad+1),ys(nquad+1),whts(nquad+1))
c
        nquad = nquad+nquad0
        end


        subroutine pv_raddiag_disk0(n,nlege,rcut,
     -    nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xslege(1000),whtslege(1000)
        dimension xs(1),ys(1),whts(1)
c
        data pi      / 3.14159265358979323846264338327950288d0 /
c
c       Construct a quadrature for a circle of a given radius.
c

        nquad = 0
c
        dd = n
        dd = 2*pi/dd
c
cccc        call legequad(nlege,xslege,whtslege)
         itype111 = 1
         call legerts(itype111, nlege, xslege, whtslege)
c
        r0 = 0
        r1 = rcut
c
        do 1000 i=1,nlege
        r    = (r1-r0)/2*xslege(i) + (r1+r0)/2
        rwht = (r1-r0)/2*whtslege(i)
c
        do 1100 j=1,n
        t    = dd*(j-1)
        twht = dd
c
        x   = r*cos(t)
        y   = r*sin(t)
        wht = twht*rwht*r
c
        nquad = nquad+1
        xs(nquad)   = x
        ys(nquad)   = y
        whts(nquad) = wht
 1100 continue
 1000 continue
c
        end

c
c
c
        subroutine pv_raddiag0(ier,rad,x1,y1,x2,y2,nquad,xs,ys,whts,
     1    rcut,norder)
        implicit double precision (a-h,o-z)
        dimension rad(*),xs(*),ys(*),whts(*),verts(2,3)
c
        dimension amatr(2,2),ainv(2,2)
        real *8, allocatable :: ts(:),twhts(:)
c
c       Return a quadrature for evaluating radial singular functions
c       of the form (1) over a triangle one of whose vertices is
c       the origin.
c
c       Triangles of this type can be transformed via rotation, scaling,
c       and, if necessary, reflection into the canonical form shown
c       below:
c
c                              z3 = a e^{ib}
c                              0
c                             / \
c                            /   \         0<a<1
c                           /     \        0<b<Pi
c                          /       \
c                         /         \
c                        /           \
c                       /             \
c                      0---------------0
c                   (0,0)            (1,0)
c
c
c       Note: the user is responsible for ensure that the ouput arrays
c       xs,ys, and whts are sufficiently large to hold the resulting
c       quadrature formula.
c
c                            Input Parameters:
c
c   rad - the radial structure as returned by radial_init
c   (x1,y1) - the coordinates of one of the vertices which is not the
c       origin
c   (x2,y2) - the coordinates of one of the other vertex which is not
c       the origin
c
c                           Output Parameters:
c
c   ier - an error return code;
c       ier = 0     indicates successful execution
c
c       ier = 64    means that one of the necessary quadrature rules is
c                   missing from the table in the rad structure; typically,
c                   this  indicates that the triangle is very close to
c                   degenerate
c
c       ier = 128   means that one of the user-supplied vertices is
c                   the origin
c
c   nquad - the number of nodes in the resulting quadrature
c   (xs,ys) - the coordinates of the quadrature nodes
c   whts - the quadrature weights
c
        ier   = 0
        nquad = 0

        allocate(ts(10000),twhts(10000))
c
        call pv_radtransform(x1,y1,x2,y2,alpha,phi,amatr,ainv,det)
        rcut0 = rcut*sqrt(det)
c
c       Fetch the theta quadrature.
c
c
        call pv_radfetch(ier,rad,alpha,phi,nts,its,itwhts)
c
        if (ier .ne. 0) then
        print *,"radfetch failed"
        print *,"alpah = ",alpha
        print *,"phi = ",phi
        ier = 64
        stop
        return
        endif
c
c       Call an auxillary routine to build the product quadrature.
c
        call pv_raddiag1(nts,rad(its),rad(itwhts),
     -    nquad,xs,ys,whts,alpha,phi,rcut0,norder)
c
        do 1000 i=1,nquad
        x = xs(i)
        y = ys(i)
        wht = whts(i)
c
        u = ainv(1,1)*x+ainv(1,2)*y
        v = ainv(2,1)*x+ainv(2,2)*y
        wht=wht/det
c
        xs(i) = u
        ys(i) = v
        whts(i) = wht
 1000 continue
        end


        subroutine pv_raddiag1(nts,ts,twhts,
     1     nquad,xs,ys,whts,alpha,phi,rcut,norder)
        implicit double precision (a-h,o-z)
        dimension ts(1),twhts(1),amatr(2,2),ainv(2,2)
        dimension xs(1),ys(1),whts(1)
c
        double precision, allocatable :: rs(:),rwhts(:)
c
c        dimension rs(10 000),rwhts(10 000)
c
        allocate(rs(10 000),rwhts(10 000))
c
        nquad = 0
        do 1000 j=1,nts
c
        t    = ts(j)*phi
        twht = twhts(j)*phi
c
        r0   = rcut
        r1   = alpha*sin(phi)/(alpha*sin(phi-t)+sin(t))
c
        if (r0 .gt. r1) then
           print *,"r1 < r0: ",r0,r1
           stop
        endif
c
c$$$        nrs = norder/2+1
c$$$        call legequad(nrs,rs,rwhts)
c$$$        do 0100 i=1,nrs
c$$$        rs(i)    = (r1-r0)/2*rs(i) + (r1+r0)/2
c$$$        rwhts(i) = (r1-r0)/2*rwhts(i)
c$$$ 0100 continue
c
        ier = 1024

        if (norder .eq. 4)  call pvquad4(ier,r0,r1,nrs,rs,rwhts)
        if (norder .eq. 8)  call pvquad8(ier,r0,r1,nrs,rs,rwhts)
        if (norder .eq. 12) call pvquad12(ier,r0,r1,nrs,rs,rwhts)
        if (norder .eq. 16) call pvquad16(ier,r0,r1,nrs,rs,rwhts)
c$$$c
c$$$        if (ier .ne. 0) return
c
        do 1100 i=1,nrs
c
        r    = rs(i)
        rwht = rwhts(i)
c
        x    = r*cos(t)
        y    = r*sin(t)
        wht  = rwht*twht*r
c
        nquad = nquad+1
        xs(nquad)   = x
        ys(nquad)   = y
        whts(nquad) = wht
c
 1100 continue
 1000 continue

c
        end



        subroutine pv_radfetch(ier,rad,r,t,n,ixs0,iwhts0)
        implicit double precision (a-h,o-z)
        dimension rad(1)
c
c       Return pointers to the nodes and weights of one of the
c       precomputed quadrature formulae residing in the rad structure.
c
c                          Input Parameters:
c
c   rad - the structure returned by radial_init
c   r - the value of the radial parameter
c   t - the value of the angular parameter
c
c                         Output Parameters:
c
c   ier - an error return code;
c       ier = 0   indicates succcessful execution
c       ier = 64  means that the requested quadrature rule is mising
c                 from the table stored in the rad structure
c
c   n - the number of quadrature nodes
c   ixs0 - a pointer into the rad structure to the quadrature nodes
c   iwhts0 - a pointer into the rad structure to the quadrature weights
c
c
        ier = 0
        n   = 0
c
c       Fetch data from the structure's header.
c
        nr     =  rad(1)
        nt     =  rad(2)
        irs    =  rad(3)
        its    =  rad(4)
        ins    =  rad(5)
        ixs    =  rad(6)
        iwhts  =  rad(7)
        max    =  rad(8)
c
c       Call an auxillary routine to find the quadrature.
c
        call pv_radfetch0(ier,max,nr,nt,rad(irs),rad(its),rad(ins),
     1     rad(ixs),rad(iwhts),n,ii,r,t)
c
        ixs0   = ixs+ii-1
        iwhts0 = iwhts+ii-1
c
        end
c
c
c
        subroutine pv_radfetch0(ier,max,nr,nt,rs,ts,ns,xs,xwhts,
     1    n,ii,r,t)
        implicit double precision (a-h,o-z)
        dimension rs(2,nr),ts(2,nt)
        dimension ns(nt,nr),xs(max,nt,nr),whts(max,nt,nr)
c
        data eps / 1.0d-15 /
        ier = 0
c
c       Figure out which quadrature rule to return.
c
        do 1000 ir=1,nr
        if (rs(1,ir) .le. r .AND. r .le. rs(2,ir)+eps) goto 1100
 1000 continue
        ier = 64
        return
 1100 continue
c
        do 2000 it=1,nt
        if (ts(1,it) .le. t .AND. t .le. ts(2,it)+eps) goto 2100
 2000 continue
        ier = 64
        return
 2100 continue
c
c       Get the number of nodes and set the pointer.
c
        n   = ns(it,ir)
        ii  = 1+max*(it-1)+max*nt*(ir-1)
c
        end
c
c
c
        subroutine pv_rad2x2inv(amatr,ainv,det)
        implicit double precision (a-h,o-z)
        dimension amatr(2,2),ainv(2,2)
c
        det = amatr(1,1)*amatr(2,2)-amatr(1,2)*amatr(2,1)
        ainv(1,1) =  amatr(2,2)/det
        ainv(2,2) =  amatr(1,1)/det
        ainv(1,2) = -amatr(1,2)/det
        ainv(2,1) = -amatr(2,1)/det
        end
c
c
c
        subroutine pv_radinsort2(k,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
c
        do 1000 i=2,k
        val=a(i)
        val2=b(i)
        j=i-1
        do 1100 while (j .ge. 1 .AND. a(j) .gt. val)
        a(j+1)=a(j)
        b(j+1)=b(j)
        j=j-1
 1100 continue
        a(j+1)=val
        b(j+1)=val2
 1000 continue
        end
c
c
c
        subroutine pv_radmove(k,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
        do 1000 j=1,k
        b(j)=a(j)
 1000 continue
        end





        subroutine pv_radadap_2x2inv(amatr,ainv,det)
        implicit double precision (a-h,o-z)
        dimension amatr(2,2),ainv(2,2)
c
        det = amatr(1,1)*amatr(2,2)-amatr(1,2)*amatr(2,1)
        ainv(1,1) =  amatr(2,2)/det
        ainv(2,2) =  amatr(1,1)/det
        ainv(1,2) = -amatr(1,2)/det
        ainv(2,1) = -amatr(2,1)/det
        end




        subroutine pv_radtransform(x1,y1,x2,y2,alpha,phi,amatr,ainv,det)
        implicit double precision (a-h,o-z)
        dimension amatr(2,2),ainv(2,2)
c
c       Construct a transformation taking a user-specified triangle T
c       with vertices
c
c         (0,0), (x1,y1) and (x2,y2)
c
c       to a triangle T0 with vertices
c
c         (0,0), (1,0), and (alpha cos(phi), alpha sin(phi))
c
c            0 < alpha < 1,  0 < phi   < pi.
c
c       Also, return the inverse mapping and the determinant of
c       the mapping.
c
c
c                           Input Parameters:
c
c   (x1,y1) - coordinates of one of the nonzero vertices
c   (x2,y2) - coordinates of the other nonzero vertex
c
c                          Output Parameters:
c
c   alpha - the radius of the
c   phi - the ngle
c
c   amatr - the matrix mapping the triangle T to T0
c   ainv - the matrix specifying the inverse of the mapping amatr
c   det - the absolute value of the determinant of the matrix amatr
c
c
        r1 = (x1**2+y1**2)
        r2 = (x2**2+y2**2)
c
        if (r1 .ge. r2) then
        u1 = x1
        v1 = y1
        u2 = x2
        v2 = y2
        else
        u1 = x2
        v1 = y2
        u2 = x1
        v2 = y1
        endif
c
        r1 = sqrt(u1**2+v1**2)
c
        theta      = atan2(v1,u1)
        amatr(1,1) = cos(theta)/r1
        amatr(1,2) = sin(theta)/r1
        amatr(2,1) =-sin(theta)/r1
        amatr(2,2) = cos(theta)/r1
c
        ainv(1,1) = cos(theta)*r1
        ainv(1,2) =-sin(theta)*r1
        ainv(2,1) = sin(theta)*r1
        ainv(2,2) = cos(theta)*r1
c
        xx = amatr(1,1)*u2 + amatr(1,2)*v2
        yy = amatr(2,1)*u2 + amatr(2,2)*v2
        phi = atan2(yy,xx)
        alpha = sqrt(xx**2+yy**2)
c
        if (phi .lt. 0) then
        amatr(2,1)=-amatr(2,1)
        amatr(2,2)=-amatr(2,2)
        ainv(1,2) =-ainv(1,2)
        ainv(2,2) =-ainv(2,2)
        phi=-phi
        endif
c
        det = (1/r1)**2
c
        end


        subroutine pvquad4(ier,r0,r1,nquad,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(*),whts(*)
        dimension deltas(2,23)
        dimension xs01(19),whts01(19)
        dimension xs02(16),whts02(16)
        dimension xs03(15),whts03(15)
        dimension xs04(15),whts04(15)
        dimension xs05(15),whts05(15)
        dimension xs06(15),whts06(15)
        dimension xs07(15),whts07(15)
        dimension xs08(15),whts08(15)
        dimension xs09(15),whts09(15)
        dimension xs10(19),whts10(19)
        dimension xs11(15),whts11(15)
        dimension xs12(15),whts12(15)
        dimension xs13(15),whts13(15)
        dimension xs14(15),whts14(15)
        dimension xs15(10),whts15(10)
        dimension xs16(09),whts16(09)
        dimension xs17(08),whts17(08)
        dimension xs18(08),whts18(08)
        dimension xs19(07),whts19(07)
        dimension xs20(07),whts20(07)
        dimension xs21(07),whts21(07)
        dimension xs22(06),whts22(06)
        dimension xs23(05),whts23(05)
        data deltas /
     -   0.1000000000000000D-14, 0.1000000000000000D-13,
     -   0.1000000000000000D-13, 0.1000000000000000D-12,
     -   0.1000000000000000D-12, 0.1000000000000000D-11,
     -   0.1000000000000000D-11, 0.1000000000000000D-10,
     -   0.1000000000000000D-10, 0.1000000000000000D-09,
     -   0.1000000000000000D-09, 0.1000000000000000D-08,
     -   0.1000000000000000D-08, 0.1000000000000000D-07,
     -   0.1000000000000000D-07, 0.1000000000000000D-06,
     -   0.1000000000000000D-06, 0.1000000000000000D-05,
     -   0.1000000000000000D-05, 0.1000000000000000D-04,
     -   0.1000000000000000D-04, 0.1000000000000000D-03,
     -   0.1000000000000000D-03, 0.1000000000000000D-02,
     -   0.1000000000000000D-02, 0.1000000000000000D-01,
     -   0.1000000000000000D-01, 0.1000000000000000D+00,
     -   0.1000000000000000D+00, 0.2000000000000000D+00,
     -   0.2000000000000000D+00, 0.3000000000000000D+00,
     -   0.3000000000000000D+00, 0.4000000000000000D+00,
     -   0.4000000000000000D+00, 0.5000000000000000D+00,
     -   0.5000000000000000D+00, 0.6000000000000000D+00,
     -   0.6000000000000000D+00, 0.7000000000000000D+00,
     -   0.7000000000000000D+00, 0.8000000000000000D+00,
     -   0.8000000000000000D+00, 0.9000000000000000D+00,
     -   0.9000000000000000D+00, 0.1000000000000000D+01/
        data xs01 /
     -   0.769091598778523001625263270569D-18,
     -   0.456860782877733950317484996940D-17,
     -   0.196391961647774506744817468515D-15,
     -   0.215153127444011396380924664005D-15,
     -   0.275370316896263262547979365670D-15,
     -   0.753712542822363617235328451697D-15,
     -   0.145716938195559275967250507066D-14,
     -   0.263815025162449596039157254647D-14,
     -   0.471213844622937558367983116378D-14,
     -   0.864425933455017367681439473417D-14,
     -   0.171726443599427013709763756601D-13,
     -   0.411980115903614578737294011452D-13,
     -   0.195991871784523771920834551409D-12,
     -   0.696331896854274817269179009639D-12,
     -   0.447365942877085283493379900656D-11,
     -   0.449621197250403511924568301512D-01,
     -   0.186984085265610791046557512130D+00,
     -   0.537098852213924477602377127069D+00,
     -   0.895389206134401080987022284092D+00/
        data whts01 /
     -  -0.188180051033170313286454741289D-15,
     -   0.246359468804765093494966864820D-15,
     -   0.248157818234164041258911800414D-14,
     -  -0.318111334619158184155516479088D-14,
     -   0.114440332416224482784398272974D-14,
     -   0.544315582321979274060745918861D-15,
     -   0.896181146816516414540264827701D-15,
     -   0.153031075747963805547942905994D-14,
     -   0.276335867696583068015911739686D-14,
     -   0.549222703511957809098418448140D-14,
     -   0.129399001501734157351984757456D-13,
     -   0.430729867376371307825369216371D-13,
     -   0.608136108748541609467256240016D-12,
     -  -0.465804873752699255471812892080D-11,
     -   0.131927190107297695105294751273D-09,
     -   0.902514741098702574964568730731D-01,
     -   0.242837039252224084225055309382D+00,
     -   0.409193002550336920565104139597D+00,
     -   0.257718483959623717907199835800D+00/
        data xs02 /
     -   0.343246704906304435380993559773D-17,
     -   0.205983518106915275457657287228D-16,
     -   0.128247319329778780726608377097D-14,
     -   0.369828195237443001207049358988D-14,
     -   0.787723568999288586474639438400D-14,
     -   0.148508684156047855204816099734D-13,
     -   0.265864350238834449934104527226D-13,
     -   0.471054077739949203707904738642D-13,
     -   0.856536796170082725299428991539D-13,
     -   0.167697845623125964626593576603D-12,
     -   0.387305016664552544653921562654D-12,
     -   0.136677851520535670985412641410D-11,
     -   0.112919890734595799623888367717D-09,
     -   0.112701667655567016051838390475D+00,
     -   0.500000001282718919925774994785D+00,
     -   0.887298334909870805791458569834D+00/
        data whts02 /
     -  -0.369682731458339430278763463424D-14,
     -   0.420503490456039957240373770674D-14,
     -   0.174577749142814645995402998356D-14,
     -   0.317497154646816133406881100998D-14,
     -   0.534804453015867294917097244081D-14,
     -   0.891094404648614971360682257354D-14,
     -   0.151868626195478924801305708046D-13,
     -   0.272529814695031049623180667078D-13,
     -   0.535223048258529652058246491441D-13,
     -   0.122869330775142766600893705964D-12,
     -   0.378274876394164690955254232179D-12,
     -   0.232501251136540388579918446008D-11,
     -   0.256249596700863598417227326261D-08,
     -   0.277777777065156219567635568075D+00,
     -   0.444444443304249851622756069975D+00,
     -   0.277777777065156154988318244337D+00/
        data xs03 /
     -   0.267714616075980688386156698705D-14,
     -   0.147170888971074497538348758219D-13,
     -   0.390620246383088041566947311342D-13,
     -   0.811864094241679196368208431474D-13,
     -   0.151441588174456921371931934325D-12,
     -   0.269595485610787650878630133301D-12,
     -   0.476065442824089193496140496796D-12,
     -   0.863777648444423684776791363386D-12,
     -   0.168864971591400630648624911170D-11,
     -   0.389555586003543230448880888406D-11,
     -   0.137250713119947693695524652184D-10,
     -   0.101449587862137838988506729643D-08,
     -   0.112701683758073379399893190676D+00,
     -   0.500000010356615994340436805882D+00,
     -   0.887298336955157448070465442397D+00/
        data whts03 /
     -   0.693852965502075463770020233599D-14,
     -   0.175623589153803799853279256842D-13,
     -   0.320094549361144371575122642370D-13,
     -   0.538968277850157544570582507269D-13,
     -   0.897459044723680837780150619147D-13,
     -   0.152860198977624676912281419280D-12,
     -   0.274169811134323067600239855045D-12,
     -   0.538225217317308485828157845994D-12,
     -   0.123512755183062975198246323148D-11,
     -   0.380037621694917258075023889116D-11,
     -   0.233072574220496836588380270145D-10,
     -   0.206837195614156023547721662638D-07,
     -   0.277777772024106350342123137075D+00,
     -   0.444444435238563732617497386182D+00,
     -   0.277777772024102186130754480338D+00/
        data xs04 /
     -   0.267702764679122295913302214897D-13,
     -   0.147163953636828916487178574245D-12,
     -   0.390599700418954512816101567444D-12,
     -   0.811814260072202773831205595241D-12,
     -   0.151430259645768774543873986886D-11,
     -   0.269569749748881579017474608035D-11,
     -   0.476004127667788003907679007848D-11,
     -   0.863615613675219041168164064684D-11,
     -   0.168813045316436313151237334001D-10,
     -   0.389313689383971416777693208648D-10,
     -   0.136975969623498276075068289713D-09,
     -   0.893813683095098116760420515074D-08,
     -   0.112701809106940540290376813172D+00,
     -   0.500000080991785859373038019778D+00,
     -   0.887298352876561133878510576313D+00/
        data whts04 /
     -   0.693821799939438151584340823931D-13,
     -   0.175614711743092128134673366959D-12,
     -   0.320074988905591992433539408396D-12,
     -   0.538926421247919827168316646497D-12,
     -   0.897367001927433386920133083384D-12,
     -   0.152838675296177374411560094306D-11,
     -   0.274114079965194114987160039542D-11,
     -   0.538055289429762941273970141586D-11,
     -   0.123444834162362090735600090563D-10,
     -   0.379592522710198245660786534694D-10,
     -   0.232158023857226965006163210649D-09,
     -   0.161689201683634157617855728050D-06,
     -   0.277777732782590008775659179643D+00,
     -   0.444444372451756279768717750337D+00,
     -   0.277777732782338822526253128062D+00/
        data xs05 /
     -   0.267687083224776253123216698411D-12,
     -   0.147154777054591197685551403709D-11,
     -   0.390572515124292144320309207807D-11,
     -   0.811748323658678075281742042717D-11,
     -   0.151415271166962575704246927128D-10,
     -   0.269535700977341729086245891436D-10,
     -   0.475923013596118101344459074904D-10,
     -   0.863401287393257365075348997662D-10,
     -   0.168744381473755063956104356890D-09,
     -   0.388994066914227440710660158876D-09,
     -   0.136614176919756853185915885724D-08,
     -   0.773142466972843854572455018977D-07,
     -   0.112702751385698995323084092276D+00,
     -   0.500000611975921986646565535323D+00,
     -   0.887298472562217255030355007409D+00/
        data whts05 /
     -   0.693780562561950472865037927234D-12,
     -   0.175602965572639955662156807279D-11,
     -   0.320049108057575639423276995143D-11,
     -   0.538871041822979771025720710099D-11,
     -   0.897245228410722694460624107605D-11,
     -   0.152810201351590618153774046768D-10,
     -   0.274040361255149950564194065525D-10,
     -   0.537830563981043923135142155008D-10,
     -   0.123355049604305190806538485557D-09,
     -   0.379004826576261026037690021775D-09,
     -   0.230956614209194638704487990633D-08,
     -   0.122100903685339936670035277143D-05,
     -   0.277777437805106799495332578151D+00,
     -   0.444443900466429014074087322509D+00,
     -   0.277777437791021738098721214834D+00/
        data xs06 /
     -   0.267665359055813048408543297565D-11,
     -   0.147142064482834413888841482183D-10,
     -   0.390534855298708932075999054360D-10,
     -   0.811656984547187953811896525389D-10,
     -   0.151394509132670848316896282279D-09,
     -   0.269488539791208209465121705613D-09,
     -   0.475810674088281957797124586029D-09,
     -   0.863104511873312517716778448581D-09,
     -   0.168649341170556808384724019904D-08,
     -   0.388552129842391828261643109574D-08,
     -   0.136116246957376354000818299387D-07,
     -   0.652495449942902296713527207535D-06,
     -   0.112709506466158820599895375399D+00,
     -   0.500004418632195746679555360630D+00,
     -   0.887299330598432565355009933150D+00/
        data whts06 /
     -   0.693723434842957967842023803976D-11,
     -   0.175586693461324885931485824090D-10,
     -   0.320013256177756762207637281034D-10,
     -   0.538794330162985152295927501998D-10,
     -   0.897076558902379644329366402733D-10,
     -   0.152770765819831494511606859146D-09,
     -   0.273938279689884158287105098775D-09,
     -   0.537519466939938957877367018239D-09,
     -   0.123230831542646554666232265786D-08,
     -   0.378193025673305674656864557638D-08,
     -   0.229308580073338366757192647280D-07,
     -   0.880742242733034159285454644324D-05,
     -   0.277775323691817485706642845167D+00,
     -   0.444440516800983828519950394980D+00,
     -   0.277775322975361945269926101525D+00/
        data xs07 /
     -   0.267633284901739214358832917336D-10,
     -   0.147123295559674782747190666507D-09,
     -   0.390479255556978174512811748768D-09,
     -   0.811522140046040338492379607531D-09,
     -   0.151363859923570004093079798547D-08,
     -   0.269418926433045056763694022438D-08,
     -   0.475644878529380726432850368102D-08,
     -   0.862666639524434676952917129105D-08,
     -   0.168509196568593584578689298403D-07,
     -   0.387901452450548716544915336988D-07,
     -   0.135388037997030545675815702298D-06,
     -   0.531934811495854813761121868043D-05,
     -   0.112754840867812460847776540625D+00,
     -   0.500029969912373779532757661434D+00,
     -   0.887305090084386416681186263917D+00/
        data whts07 /
     -   0.693639090191540604400160244747D-10,
     -   0.175562669585258906223716202220D-09,
     -   0.319960327479804549949854492099D-09,
     -   0.538681087079068696467009928938D-09,
     -   0.896827589914244409694098653362D-09,
     -   0.152712564386429962427446905955D-08,
     -   0.273787656612699541659772266653D-08,
     -   0.537060631035274538915058373352D-08,
     -   0.123047782427913969111528812533D-07,
     -   0.376999491938184755065420330699D-07,
     -   0.226909876516521721537558452070D-06,
     -   0.596187549587133853556317292382D-04,
     -   0.277761159330479405694885544271D+00,
     -   0.444417805837453074838698970672D+00,
     -   0.277761127526500749527894845276D+00/
        data xs08 /
     -   0.267581352710097989068048041765D-09,
     -   0.147092906829099630388150607950D-08,
     -   0.390389237375595056258918079742D-08,
     -   0.811303834175995889197603060399D-08,
     -   0.151314244990995280030994788083D-07,
     -   0.269306252282478766193627926567D-07,
     -   0.475376589922121640705972940512D-07,
     -   0.861958384430951183078587915306D-07,
     -   0.168282719177022491630226678556D-06,
     -   0.386852474699013584596044836469D-06,
     -   0.134226424816070428899330210441D-05,
     -   0.411794960163005068439653978671D-04,
     -   0.113031223331215681452226851728D+00,
     -   0.500185893720693765923750699535D+00,
     -   0.887340241001086653416026193726D+00/
        data whts08 /
     -   0.693502525904390741820613587977D-09,
     -   0.175523773392987233659479095262D-08,
     -   0.319874637978106853396168286594D-08,
     -   0.538497768405407930178508582957D-08,
     -   0.896424614081353834233440021619D-08,
     -   0.152618381349827239608540369547D-07,
     -   0.273544002660937600627212979970D-07,
     -   0.536318888345251846797809973077D-07,
     -   0.122752276905827813115932647736D-06,
     -   0.375079705347023881798886977604D-06,
     -   0.223112304115308945590924170094D-05,
     -   0.367760465431869272000262634697D-03,
     -   0.277675648204166470683961972110D+00,
     -   0.444279253596244204431258531747D+00,
     -   0.277674492534296349587010449595D+00/
        data xs09 /
     -   0.267485050282545992878543822320D-08,
     -   0.147036555411058226558832806687D-07,
     -   0.390222318548139356899409131237D-07,
     -   0.810899062684960095576985491996D-07,
     -   0.151222262630587466969542082488D-06,
     -   0.269097406363804251006956078312D-06,
     -   0.474879492976993733979401730569D-06,
     -   0.860647051203975264379493155610D-06,
     -   0.167864061863555903721113847880D-05,
     -   0.384921708275276258265397753205D-05,
     -   0.132128216489844121599917642837D-04,
     -   0.293720843260262106101076409018D-03,
     -   0.114491376009441023767273372548D+00,
     -   0.501013526912565837221051359208D+00,
     -   0.887526932942236725261984685571D+00/
        data whts09 /
     -   0.693249283791104340564629074060D-08,
     -   0.175451647727887048397878151678D-07,
     -   0.319715754532937629108776801636D-07,
     -   0.538157905287941574515880007463D-07,
     -   0.895677670481720094139708674060D-07,
     -   0.152443866434136834749704809649D-06,
     -   0.273092800905177989047738709050D-06,
     -   0.534946892790425370944040308189D-06,
     -   0.122207020833025833430996957667D-05,
     -   0.371560436391727854549185427444D-05,
     -   0.216344436397630290750636339461D-04,
     -   0.196766733260875579411670173412D-02,
     -   0.277245345775017266828782454200D+00,
     -   0.443544847868157303893687357142D+00,
     -   0.277214406589653892217585858112D+00/
        data xs10 /
     -   0.181134949563862657104686417057D-07,
     -   0.104138547438182378938601074629D-06,
     -   0.287628740187050745582254354688D-06,
     -   0.612102415891435411300463156912D-06,
     -   0.115166422550280977107117705631D-05,
     -   0.204031032426298884727071698856D-05,
     -   0.353353888715807058327317533718D-05,
     -   0.616006679342554701682447535686D-05,
     -   0.111547635522663053274414474119D-04,
     -   0.219549653959447922663746246713D-04,
     -   0.510825890435454402104691836004D-04,
     -   0.172394101493932665765103136845D-03,
     -   0.175204170066048322668049347385D-02,
     -   0.376651148320314562492078635828D-01,
     -   0.174075227838425183642339186136D+00,
     -   0.441501707312664710405344582662D+00,
     -   0.667287752333306255540098632098D+00,
     -   0.868505540191903804570193961381D+00,
     -   0.988425613425481612026438688181D+00/
        data whts10 /
     -   0.476043166920279250351514195640D-07,
     -   0.129092084295051893632660624220D-06,
     -   0.244876874185989346379544023748D-06,
     -   0.415973568637047552429510881102D-06,
     -   0.684567337002892181417915175277D-06,
     -   0.113312218625303504058401151810D-05,
     -   0.193550866318024097190033943544D-05,
     -   0.350432938094027036489413747638D-05,
     -   0.698169230108708808107162332212D-05,
     -   0.162764707119575830318095137542D-04,
     -   0.499170107186733991002461023120D-04,
     -   0.267768853761357694767844821455D-03,
     -   0.594407277556010160411679404403D-02,
     -   0.758421018555833234396994572469D-01,
     -   0.212923791794636521850566038255D+00,
     -   0.280101005864622772653079679677D+00,
     -   0.188848455279169221494660993401D+00,
     -   0.190737309670327960084782762447D+00,
     -   0.452542236581958365528370295972D-01/
        data xs11 /
     -   0.266691615662893581399460818606D-06,
     -   0.146572238711896595463585422167D-05,
     -   0.388846833295177552362814954597D-05,
     -   0.807563507268464204584900467162D-05,
     -   0.150464430173607175988413030408D-04,
     -   0.267377934139432818631915071380D-04,
     -   0.470793715287157139809104765726D-04,
     -   0.849910645763920985859833714913D-04,
     -   0.164467993711297964767131830528D-03,
     -   0.369659536055964690602028699061D-03,
     -   0.117220836449107359083592531654D-02,
     -   0.957705990396998517856203368807D-02,
     -   0.142230693364726021748844505642D+00,
     -   0.517573788041043016239196983367D+00,
     -   0.891293195772994289107958056566D+00/
        data whts11 /
     -   0.691162780770025355064209500872D-06,
     -   0.174857310199561657146564204456D-05,
     -   0.318406381350668731028589739137D-05,
     -   0.535357444506430243371991086450D-05,
     -   0.889526642599173440760151611937D-05,
     -   0.151008888408226594121967367963D-04,
     -   0.269394455245499456282983202660D-04,
     -   0.523775451650026886244472870660D-04,
     -   0.117832588967373640430540622356D-03,
     -   0.344418634693514810580091774966D-03,
     -   0.171744684804814672691874795789D-02,
     -   0.277506495350903038909296859545D-01,
     -   0.272891337333495530501518716311D+00,
     -   0.429118420789174868179702427745D+00,
     -   0.267945603750432558590176710114D+00/
        data xs12 /
     -   0.264939566067450881957475953553D-05,
     -   0.145546562764522131802033090774D-04,
     -   0.385806884965920119039955277460D-04,
     -   0.800189195534323375055758855426D-04,
     -   0.148789482986468062723589557250D-03,
     -   0.263583988664278811231915638887D-03,
     -   0.461820042118383186610316331803D-03,
     -   0.826586133108247567003391111190D-03,
     -   0.157279660045842385735165084538D-02,
     -   0.339533254102254245621241817254D-02,
     -   0.943115136888298066312233324723D-02,
     -   0.407844144582554012201291597482D-01,
     -   0.201546386770886399876805965543D+00,
     -   0.555761689003584610169217532854D+00,
     -   0.900145825113832590760923022504D+00/
        data whts12 /
     -   0.686554971870684497940951067366D-05,
     -   0.173543936965497042017347083923D-04,
     -   0.315511115864756621808857024067D-04,
     -   0.529165283482862968394322949389D-04,
     -   0.875945422175423123703325063782D-04,
     -   0.147853264845200907937256655991D-03,
     -   0.261333965460530001310227281110D-03,
     -   0.499883573919252486325759746068D-03,
     -   0.108859440120590779791150995136D-02,
     -   0.294342675988074248250820081195D-02,
     -   0.114546704950561276893746859765D-01,
     -   0.686292490837690675492072692907D-01,
     -   0.271354177897905858328894103264D+00,
     -   0.397189286397534225548970328630D+00,
     -   0.246235242034855526386988863669D+00/
        data xs13 /
     -   0.261942590818831193488577752657D-04,
     -   0.143749911642469848397201547749D-03,
     -   0.380258639823720392903263318721D-03,
     -   0.785961986011668318770087210798D-03,
     -   0.145338963690044216373782603370D-02,
     -   0.255206611408983165756244930710D-02,
     -   0.440692354932290769543264522218D-02,
     -   0.769142626797452237470324191434D-02,
     -   0.139568231273766240039673508505D-01,
     -   0.272559365244730987621785015252D-01,
     -   0.594572098903391301817532529228D-01,
     -   0.144418830091860976637541147258D+00,
     -   0.341473416734274999654985037409D+00,
     -   0.648508707403705185634820868945D+00,
     -   0.922175765148802770103043225638D+00/
        data whts13 /
     -   0.678629665764754405168664483912D-04,
     -   0.171184443472138690207341994684D-03,
     -   0.309955713908422298893393595939D-03,
     -   0.516344949415243708444159586916D-03,
     -   0.845595752194089838903703442874D-03,
     -   0.140310889753205708021487723176D-02,
     -   0.241141087586592464025900375180D-02,
     -   0.439402521265070422772028586962D-02,
     -   0.873949720877644039593396687028D-02,
     -   0.196020989713085857424529611714D-01,
     -   0.500867731691565516833948006917D-01,
     -   0.131307518941576625506784158836D+00,
     -   0.264658184964160253223384603051D+00,
     -   0.322984308312455820376965236607D+00,
     -   0.192502129620950667145924640851D+00/
        data xs14 /
     -   0.253621214435242793416900931891D-03,
     -   0.138793238458954842773515963196D-02,
     -   0.365155045246457660523890976849D-02,
     -   0.748067921177647173109372839542D-02,
     -   0.136425709311977851818191197341D-01,
     -   0.234407841807348371433139249780D-01,
     -   0.390909747524159142707099670427D-01,
     -   0.643915651326541632524463498579D-01,
     -   0.105814509593165914430627615908D+00,
     -   0.173771499205894996743893840459D+00,
     -   0.282167423942024813958207026784D+00,
     -   0.441175365943549919972369943195D+00,
     -   0.640142559085109965151094874259D+00,
     -   0.835250140043408216052613047415D+00,
     -   0.966913768310051391338734324492D+00/
        data whts14 /
     -   0.656654111901597741738973584084D-03,
     -   0.164722388079351519544881329252D-02,
     -   0.295095807718051237531199850621D-02,
     -   0.483192621420179345395536209569D-02,
     -   0.770315421217678504386494045750D-02,
     -   0.122517672976052567093470880181D-01,
     -   0.196633825258910746194588246893D-01,
     -   0.319919059367178783558041547855D-01,
     -   0.525837927621045786546961361081D-01,
     -   0.857517167214613480262306150457D-01,
     -   0.133053698812117608137216082545D+00,
     -   0.183480717565311451744031725394D+00,
     -   0.206838225644249638922955934058D+00,
     -   0.172692495643047711600429606462D+00,
     -   0.839023805952392494195097449586D-01/
        data xs15 /
     -   0.325594294393931817733148353156D-02,
     -   0.179519024329837740295501401941D-01,
     -   0.478159033934114649046110420938D-01,
     -   0.993266651390469566325445780409D-01,
     -   0.182484243488806283089606890699D+00,
     -   0.308696114036277436386473016783D+00,
     -   0.482201608651301636538974975894D+00,
     -   0.684780380344215155810960787979D+00,
     -   0.868063885241499842867451415439D+00,
     -   0.976713891668238289438111164153D+00/
        data whts15 /
     -   0.844499518865578257455161958144D-02,
     -   0.214889038011645927448584068641D-01,
     -   0.393083763886414358553856650101D-01,
     -   0.654038185852817993087398200173D-01,
     -   0.102948435520487144619715322924D+00,
     -   0.150471477042641293636116985172D+00,
     -   0.193701799006545728314211950325D+00,
     -   0.203026312158710663591709667355D+00,
     -   0.153262104246320004030955896414D+00,
     -   0.619437780615515553237546663368D-01/
        data xs16 /
     -   0.642560938143521123729295363351D-02,
     -   0.350145564070644271582972013533D-01,
     -   0.910665239263104368573957453437D-01,
     -   0.181643474154590019892479399823D+00,
     -   0.313317336622328020572161985033D+00,
     -   0.484922539863381961526179187397D+00,
     -   0.678087597365131433693670045414D+00,
     -   0.854852795172797850630229380205D+00,
     -   0.970909939112262146843862178336D+00/
        data whts16 /
     -   0.166228908412124461672085186743D-01,
     -   0.413119512556262336067555826297D-01,
     -   0.720145427423706281382036288639D-01,
     -   0.110357381433320728366522294303D+00,
     -   0.152978975728159419596841074779D+00,
     -   0.187239457279254291964225824304D+00,
     -   0.192521510296623189270810875396D+00,
     -   0.153221039691634136597345747408D+00,
     -   0.737322507317989262920864536410D-01/
        data xs17 /
     -   0.999034260038183139688187790030D-02,
     -   0.538925638754004449300766693198D-01,
     -   0.137440038415819341946635857289D+00,
     -   0.265702900377124750846261014804D+00,
     -   0.438251200168473176693182462131D+00,
     -   0.639629237046671117738518871701D+00,
     -   0.832719909723786491753428710756D+00,
     -   0.965702818018701069612259100736D+00/
        data whts17 /
     -   0.257859427414674474788898494678D-01,
     -   0.628009400753683953313646597693D-01,
     -   0.105239783829018422775218359605D+00,
     -   0.151387715644374964608263356508D+00,
     -   0.191275213125812115248827231543D+00,
     -   0.205191027358669648381891394382D+00,
     -   0.171923986012769819363475471071D+00,
     -   0.863953912125191868120696776536D-01/
        data xs18 /
     -   0.136398141755597247404010907409D-01,
     -   0.727140980583600249054552337821D-01,
     -   0.181368573651882474988251332482D+00,
     -   0.339097917218204285360472074758D+00,
     -   0.534982848682218743691368070328D+00,
     -   0.740185596710822264911967993933D+00,
     -   0.908647279912208656559831145140D+00,
     -   0.992662746832692360559806675797D+00/
        data whts18 /
     -   0.351102760016077152790072707410D-01,
     -   0.835092179151601446778205479833D-01,
     -   0.133848436242637533699813007758D+00,
     -   0.179937640911372268983790869087D+00,
     -   0.207034695783883123591810429863D+00,
     -   0.195428348692009946512587181218D+00,
     -   0.133030364959700348375235246250D+00,
     -   0.321010194936289188799354470992D-01/
        data xs19 /
     -   0.181954843109814057323834128755D-01,
     -   0.957101001537606611777746990136D-01,
     -   0.232965154505592046625387539413D+00,
     -   0.420428552922705506381443449048D+00,
     -   0.634056301498665736091396836046D+00,
     -   0.832568508069361488888575628802D+00,
     -   0.965987496448592861403060792688D+00/
        data whts19 /
     -   0.466916692596755828061611693218D-01,
     -   0.108113594710982450028786384328D+00,
     -   0.164938838073281936443174081612D+00,
     -   0.206046737957294408247293132053D+00,
     -   0.214189128806311030879607202988D+00,
     -   0.174144178090884000672741013947D+00,
     -   0.858758531015705909222370157507D-01/
        data xs20 /
     -   0.147412661728270259924435421028D-01,
     -   0.841252834972853142857797760213D-01,
     -   0.218246247242424600120216821236D+00,
     -   0.408412127443121779438736655835D+00,
     -   0.627086471743663136995636579672D+00,
     -   0.829842465189431966216450572125D+00,
     -   0.965505654754910392358817426332D+00/
        data whts20 /
     -   0.388772185910518438328053091338D-01,
     -   0.101613220406601642905633916932D+00,
     -   0.165164291000713715759966954724D+00,
     -   0.210570972858467976313921792568D+00,
     -   0.219207803105492161104888512105D+00,
     -   0.177426611084379457282052751122D+00,
     -   0.871398829532932028007307634156D-01/
        data xs21 /
     -   0.238062692914894025795140936060D-01,
     -   0.122375038138908984017206839786D+00,
     -   0.286464701553976975768136812801D+00,
     -   0.491379535239204904473285697228D+00,
     -   0.700643669879715600413171329065D+00,
     -   0.872841681082328333036733608966D+00,
     -   0.975664665760870360740982474428D+00/
        data whts21 /
     -   0.607538504049195367560893858498D-01,
     -   0.134372130970792257618942149129D+00,
     -   0.189694136311683440505043247750D+00,
     -   0.213947482145661476567034009454D+00,
     -   0.197423258840452394092974115690D+00,
     -   0.141339410184732360578147731807D+00,
     -   0.624697311417585338817693603207D-01/
        data xs22 /
     -   0.311300241088926504908259035138D-01,
     -   0.157984882705996762242277258712D+00,
     -   0.361452091545795283812566772344D+00,
     -   0.599938809293803471741671864246D+00,
     -   0.818975238665459717754268811795D+00,
     -   0.963527925274846154256742474813D+00/
        data whts22 /
     -   0.792036305404581774496491398559D-01,
     -   0.170647780077384856753305024805D+00,
     -   0.229271004267373577917507742965D+00,
     -   0.238414547339674950277479579037D+00,
     -   0.190179992634112337788269453532D+00,
     -   0.922830451409960998137890598054D-01/
        data xs23 /
     -   0.454798283439560270719406096257D-01,
     -   0.225083682330710095994910499302D+00,
     -   0.491992167820785033027040511732D+00,
     -   0.763544600599379416415132390406D+00,
     -   0.951656116556948377114040485337D+00/
        data whts23 /
     -   0.115028857431326825677712230660D+00,
     -   0.235185869253024002272452572515D+00,
     -   0.284438202006767096216540590981D+00,
     -   0.243441432252069583477375411375D+00,
     -   0.121905639056812492355919194468D+00/

        ier = 0

        delta = r0/r1

        if( deltas(1,01) .le. delta .AND. delta .le. deltas(2,01)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs01(l)
        whts(l) = whts01(l)
        end do
        goto 1000
        end if


        if( deltas(1,02) .le. delta .AND. delta .le. deltas(2,02)) then
        nquad = 16
        do l=1,16
        xs(l)   = xs02(l)
        whts(l) = whts02(l)
        end do
        goto 1000
        end if


        if( deltas(1,03) .le. delta .AND. delta .le. deltas(2,03)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs03(l)
        whts(l) = whts03(l)
        end do
        goto 1000
        end if


        if( deltas(1,04) .le. delta .AND. delta .le. deltas(2,04)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs04(l)
        whts(l) = whts04(l)
        end do
        goto 1000
        end if


        if( deltas(1,05) .le. delta .AND. delta .le. deltas(2,05)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs05(l)
        whts(l) = whts05(l)
        end do
        goto 1000
        end if


        if( deltas(1,06) .le. delta .AND. delta .le. deltas(2,06)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs06(l)
        whts(l) = whts06(l)
        end do
        goto 1000
        end if


        if( deltas(1,07) .le. delta .AND. delta .le. deltas(2,07)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs07(l)
        whts(l) = whts07(l)
        end do
        goto 1000
        end if


        if( deltas(1,08) .le. delta .AND. delta .le. deltas(2,08)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs08(l)
        whts(l) = whts08(l)
        end do
        goto 1000
        end if


        if( deltas(1,09) .le. delta .AND. delta .le. deltas(2,09)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs09(l)
        whts(l) = whts09(l)
        end do
        goto 1000
        end if


        if( deltas(1,10) .le. delta .AND. delta .le. deltas(2,10)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs10(l)
        whts(l) = whts10(l)
        end do
        goto 1000
        end if


        if( deltas(1,11) .le. delta .AND. delta .le. deltas(2,11)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs11(l)
        whts(l) = whts11(l)
        end do
        goto 1000
        end if


        if( deltas(1,12) .le. delta .AND. delta .le. deltas(2,12)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs12(l)
        whts(l) = whts12(l)
        end do
        goto 1000
        end if


        if( deltas(1,13) .le. delta .AND. delta .le. deltas(2,13)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs13(l)
        whts(l) = whts13(l)
        end do
        goto 1000
        end if


        if( deltas(1,14) .le. delta .AND. delta .le. deltas(2,14)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs14(l)
        whts(l) = whts14(l)
        end do
        goto 1000
        end if


        if( deltas(1,15) .le. delta .AND. delta .le. deltas(2,15)) then
        nquad = 10
        do l=1,10
        xs(l)   = xs15(l)
        whts(l) = whts15(l)
        end do
        goto 1000
        end if


        if( deltas(1,16) .le. delta .AND. delta .le. deltas(2,16)) then
        nquad = 09
        do l=1,09
        xs(l)   = xs16(l)
        whts(l) = whts16(l)
        end do
        goto 1000
        end if


        if( deltas(1,17) .le. delta .AND. delta .le. deltas(2,17)) then
        nquad = 08
        do l=1,08
        xs(l)   = xs17(l)
        whts(l) = whts17(l)
        end do
        goto 1000
        end if


        if( deltas(1,18) .le. delta .AND. delta .le. deltas(2,18)) then
        nquad = 08
        do l=1,08
        xs(l)   = xs18(l)
        whts(l) = whts18(l)
        end do
        goto 1000
        end if


        if( deltas(1,19) .le. delta .AND. delta .le. deltas(2,19)) then
        nquad = 07
        do l=1,07
        xs(l)   = xs19(l)
        whts(l) = whts19(l)
        end do
        goto 1000
        end if


        if( deltas(1,20) .le. delta .AND. delta .le. deltas(2,20)) then
        nquad = 07
        do l=1,07
        xs(l)   = xs20(l)
        whts(l) = whts20(l)
        end do
        goto 1000
        end if


        if( deltas(1,21) .le. delta .AND. delta .le. deltas(2,21)) then
        nquad = 07
        do l=1,07
        xs(l)   = xs21(l)
        whts(l) = whts21(l)
        end do
        goto 1000
        end if


        if( deltas(1,22) .le. delta .AND. delta .le. deltas(2,22)) then
        nquad = 06
        do l=1,06
        xs(l)   = xs22(l)
        whts(l) = whts22(l)
        end do
        goto 1000
        end if


        if( deltas(1,23) .le. delta .AND. delta .le. deltas(2,23)) then
        nquad = 05
        do l=1,05
        xs(l)   = xs23(l)
        whts(l) = whts23(l)
        end do
        goto 1000
        end if

        ier = 4
        return

 1000 continue

        do i=1,nquad
        xs(i)   = (r1-r0)*xs(i)+r0
        whts(i) = (r1-r0)*whts(i)
        end do

        end


        subroutine pvquad8(ier,r0,r1,nquad,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(*),whts(*)
        dimension deltas(2,23)
        dimension xs01(18),whts01(18)
        dimension xs02(17),whts02(17)
        dimension xs03(17),whts03(17)
        dimension xs04(21),whts04(21)
        dimension xs05(17),whts05(17)
        dimension xs06(17),whts06(17)
        dimension xs07(17),whts07(17)
        dimension xs08(17),whts08(17)
        dimension xs09(21),whts09(21)
        dimension xs10(17),whts10(17)
        dimension xs11(17),whts11(17)
        dimension xs12(17),whts12(17)
        dimension xs13(17),whts13(17)
        dimension xs14(17),whts14(17)
        dimension xs15(11),whts15(11)
        dimension xs16(11),whts16(11)
        dimension xs17(10),whts17(10)
        dimension xs18(09),whts18(09)
        dimension xs19(08),whts19(08)
        dimension xs20(08),whts20(08)
        dimension xs21(07),whts21(07)
        dimension xs22(07),whts22(07)
        dimension xs23(06),whts23(06)
        data deltas /
     -   0.1000000000000000D-14, 0.1000000000000000D-13,
     -   0.1000000000000000D-13, 0.1000000000000000D-12,
     -   0.1000000000000000D-12, 0.1000000000000000D-11,
     -   0.1000000000000000D-11, 0.1000000000000000D-10,
     -   0.1000000000000000D-10, 0.1000000000000000D-09,
     -   0.1000000000000000D-09, 0.1000000000000000D-08,
     -   0.1000000000000000D-08, 0.1000000000000000D-07,
     -   0.1000000000000000D-07, 0.1000000000000000D-06,
     -   0.1000000000000000D-06, 0.1000000000000000D-05,
     -   0.1000000000000000D-05, 0.1000000000000000D-04,
     -   0.1000000000000000D-04, 0.1000000000000000D-03,
     -   0.1000000000000000D-03, 0.1000000000000000D-02,
     -   0.1000000000000000D-02, 0.1000000000000000D-01,
     -   0.1000000000000000D-01, 0.1000000000000000D+00,
     -   0.1000000000000000D+00, 0.2000000000000000D+00,
     -   0.2000000000000000D+00, 0.3000000000000000D+00,
     -   0.3000000000000000D+00, 0.4000000000000000D+00,
     -   0.4000000000000000D+00, 0.5000000000000000D+00,
     -   0.5000000000000000D+00, 0.6000000000000000D+00,
     -   0.6000000000000000D+00, 0.7000000000000000D+00,
     -   0.7000000000000000D+00, 0.8000000000000000D+00,
     -   0.8000000000000000D+00, 0.9000000000000000D+00,
     -   0.9000000000000000D+00, 0.1000000000000000D+01/
        data xs01 /
     -   0.269008882881993457996121059870D-16,
     -   0.147928768565714405560368654727D-15,
     -   0.392868264803355134371063932916D-15,
     -   0.817327812802575599555049314567D-15,
     -   0.152687652699605418717773019266D-14,
     -   0.272441079790837758588824745671D-14,
     -   0.482905540764252738006811385271D-14,
     -   0.882156258200355766793464719865D-14,
     -   0.174971666583041653029380363507D-13,
     -   0.421195713123985847121392386264D-13,
     -   0.257343559763077173403625053363D-12,
     -   0.337627945649472431555722136200D-12,
     -   0.562094321039745355315131904203D-11,
     -   0.469100771717386934509337124918D-01,
     -   0.230765345061015983509664162194D+00,
     -   0.500000000074007020028341344883D+00,
     -   0.769234655086998056530494305294D+00,
     -   0.953089922976275346419944009967D+00/
        data whts01 /
     -   0.697256979164316566974280613159D-16,
     -   0.176594393184151845705226838195D-15,
     -   0.322238334574689926322020410057D-15,
     -   0.543571292959663738893963395278D-15,
     -   0.907633588634313536390905828087D-15,
     -   0.155258991484075094329873175503D-14,
     -   0.280467215217175051492993779003D-14,
     -   0.557924602451745003457604079103D-14,
     -   0.131838772424337271056930816741D-13,
     -   0.444507260180501873684739922258D-13,
     -   0.285810175596394973524873792585D-11,
     -  -0.404267625500580512802236607693D-11,
     -   0.149129023094310614168248132053D-09,
     -   0.118463442510560291571732963122D+00,
     -   0.239314335214261352457481515277D+00,
     -   0.284444444402342673056298731879D+00,
     -   0.239314335214261352416110197399D+00,
     -   0.118463442510560291028448550431D+00/
        data xs02 /
     -   0.267720514740514928422271441999D-15,
     -   0.147174340834726517570935341531D-14,
     -   0.390630472638833430699791684245D-14,
     -   0.811888898337757999627956180643D-14,
     -   0.151447226870831689790349774487D-13,
     -   0.269608295860626417879021363233D-13,
     -   0.476095964527080461052042123437D-13,
     -   0.863858314037868611705398996719D-13,
     -   0.168890826810388664758203370523D-12,
     -   0.389676090876918123289838794752D-12,
     -   0.137387880096994753808579056219D-11,
     -   0.108801109400553241686022550369D-09,
     -   0.469100792929085389916768981544D-01,
     -   0.230765346773002803282183622127D+00,
     -   0.500000001186792834540328736184D+00,
     -   0.769234655600582861618564748965D+00,
     -   0.953089923080677083088412308617D+00/
        data whts02 /
     -   0.693868477209268772760072564185D-15,
     -   0.175628007653720659532366879397D-14,
     -   0.320104285218913060907936677962D-14,
     -   0.538989111670027488536593527004D-14,
     -   0.897504859807922568271387095718D-14,
     -   0.152870913078744580361306672332D-13,
     -   0.274197555215974833777825454825D-13,
     -   0.538309822272243850827805519383D-13,
     -   0.123546581693428513252375171778D-12,
     -   0.380259457125427652538030001850D-12,
     -   0.233529876150459116588310881737D-11,
     -   0.237062986189677996833585087921D-08,
     -   0.118463442246911551417210848001D+00,
     -   0.239314334681650167741718593872D+00,
     -   0.284444443769291188917302198668D+00,
     -   0.239314334681650157276575700615D+00,
     -   0.118463442246911413989911831744D+00/
        data xs03 /
     -   0.267710345650141219258658507497D-14,
     -   0.147168389949945140584582221657D-13,
     -   0.390612842991544586272874574159D-13,
     -   0.811846137232436774569164581736D-13,
     -   0.151437506069687908722345426263D-12,
     -   0.269586211864584059194036802291D-12,
     -   0.476043347836815676721952059385D-12,
     -   0.863719256774013390905387974456D-12,
     -   0.168846257683896461001652554506D-11,
     -   0.389468389604861973550059723312D-11,
     -   0.137151584037104946358507639026D-10,
     -   0.967324466230111068530353906779D-09,
     -   0.469100950272361123986844469991D-01,
     -   0.230765359472110948880596477078D+00,
     -   0.500000009441171615345694532606D+00,
     -   0.769234659410232020410802563822D+00,
     -   0.953089923855104179007783684977D+00/
        data whts03 /
     -   0.693841735578118595902568644222D-14,
     -   0.175620390333616949866875164520D-13,
     -   0.320087501045399861388707303076D-13,
     -   0.538953195309859011531684957154D-13,
     -   0.897425877722380361124770305405D-13,
     -   0.152852442907482338746340112338D-12,
     -   0.274149727724522885274701025339D-12,
     -   0.538163978164783985380418100785D-12,
     -   0.123488274651130837953862304832D-11,
     -   0.379877145473420404956230985828D-11,
     -   0.232742380354111137419236427001D-10,
     -   0.188528607432097141267162398862D-07,
     -   0.118463440291235738924354359117D+00,
     -   0.239314330730868458399348257054D+00,
     -   0.284444439073466905422859220879D+00,
     -   0.239314330730867803938189289773D+00,
     -   0.118463440291227144606284424276D+00/
        data xs04 /
     -   0.108892554806303343424706363485D-13,
     -   0.846854462634370802236199142934D-13,
     -   0.262340123339131949802549934891D-12,
     -   0.584944395485251952595479450881D-12,
     -   0.112709677614831016264641167722D-11,
     -   0.202800003425007203420253376632D-11,
     -   0.355857354686236103787835099299D-11,
     -   0.629424810798385325790578033401D-11,
     -   0.116346646242459893754089454066D-10,
     -   0.237543478899831914388478798187D-10,
     -   0.599581423485786004365134313828D-10,
     -   0.255035102181448408674468329280D-09,
     -   0.187877643006030035058142072734D-07,
     -   0.249937947657116892497199745659D-01,
     -   0.923985895042976709822162200697D-01,
     -   0.213572083989322669994366526225D+00,
     -   0.410478934058975963877514811672D+00,
     -   0.615711890822325414595401869050D+00,
     -   0.753533653221228097178732935773D+00,
     -   0.861513119308073322280498847223D+00,
     -   0.970068171505101004589685731126D+00/
        data whts04 /
     -   0.330466788879232903372898184715D-13,
     -   0.120464249382439912189122425326D-12,
     -   0.241265149204667449706575034060D-12,
     -   0.415926991560730877648847132865D-12,
     -   0.690588234361806261426056827540D-12,
     -   0.115391170855831805041197790396D-11,
     -   0.199615036784598623200398198456D-11,
     -   0.368461222753375583101639913538D-11,
     -   0.758602406285532194399428997575D-11,
     -   0.188287830961005276336705079625D-10,
     -   0.661260910784464221887255595167D-10,
     -   0.499025990994045261358092264019D-09,
     -   0.298196921664934212363328690752D-06,
     -   0.577713696578539873270292055169D-01,
     -   0.783039713630616924994555982620D-01,
     -   0.168711262711354912462804283072D+00,
     -   0.212685914041025968600818129992D+00,
     -   0.185111663334196871578590977777D+00,
     -   0.957432566377321609727854812309D-01,
     -   0.126859757194892626700213750203D+00,
     -   0.748125062630572600853070495882D-01/
        data xs05 /
     -   0.267679483682674089900577261203D-12,
     -   0.147150329930052280428934277100D-11,
     -   0.390559340834252039565834575594D-11,
     -   0.811716370745262096743832917193D-11,
     -   0.151408007921431682826080893621D-10,
     -   0.269519202022082554062924182514D-10,
     -   0.475883710938497016157519341603D-10,
     -   0.863297451125997054084394756613D-10,
     -   0.168711123618420431568192896732D-09,
     -   0.388839356423254161599988271333D-09,
     -   0.136439558723646160547501976455D-08,
     -   0.725981475610605215551977370023D-07,
     -   0.469111104943165724061364863618D-01,
     -   0.230766179057746269986120414039D+00,
     -   0.500000542169737583829456830028D+00,
     -   0.769234905280896858411216249491D+00,
     -   0.953089973835803023477467310877D+00/
        data whts05 /
     -   0.693760578144042387987986629107D-12,
     -   0.175597273222869417673881879608D-11,
     -   0.320036566113143498840588546869D-11,
     -   0.538844205564464757626862388748D-11,
     -   0.897186220789066563927495786038D-11,
     -   0.152796404637643854296743897658D-10,
     -   0.274004645332270767682777972110D-10,
     -   0.537721706548835214396029932703D-10,
     -   0.123311574000707664299796707256D-09,
     -   0.378720530809370606263586310426D-09,
     -   0.230377946487963509825760208984D-08,
     -   0.108138765753245419559950653495D-05,
     -   0.118463314100807735240847625807D+00,
     -   0.239314075753747909769713820919D+00,
     -   0.284444136010394629153912969346D+00,
     -   0.239314075751664756653312904595D+00,
     -   0.118463314073453188151389242600D+00/
        data xs06 /
     -   0.267654409058370761896603389307D-11,
     -   0.147135656810061599320073851246D-10,
     -   0.390515873478714215197082765108D-10,
     -   0.811610947775316701563070284596D-10,
     -   0.151384045056136727313375968463D-09,
     -   0.269464771944812546617827803595D-09,
     -   0.475754063724301860983701587097D-09,
     -   0.862954985620159890394035788280D-09,
     -   0.168601473256967167527970025171D-08,
     -   0.388329750330406003141105461377D-08,
     -   0.135866711412178027288161278769D-07,
     -   0.605354984853710784101803998809D-06,
     -   0.469173765074621258114990396025D-01,
     -   0.230771236645381108042044822702D+00,
     -   0.500003829613334469398032610051D+00,
     -   0.769236422540918350278211572519D+00,
     -   0.953090282265347387828295622758D+00/
        data whts06 /
     -   0.693694639859069916727933992020D-11,
     -   0.175578491712687272098750791856D-10,
     -   0.319995186027594600585633155272D-10,
     -   0.538755667400035457006270391923D-10,
     -   0.896991554375423752228533159590D-10,
     -   0.152750893233500042844124825198D-09,
     -   0.273886845536138085502613231440D-09,
     -   0.537362759774247923946260814300D-09,
     -   0.123168292976868710438407784514D-08,
     -   0.377784886738701719300480169329D-08,
     -   0.228485068803130140243950628426D-07,
     -   0.762877152842303098226561478097D-05,
     -   0.118462536513937486800171780896D+00,
     -   0.239312502386181343264383371491D+00,
     -   0.284442265834030162062091032687D+00,
     -   0.239312502285110744198370880112D+00,
     -   0.118462535187103628281231487763D+00/
        data xs07 /
     -   0.267616188514859623539113258698D-10,
     -   0.147113291363097681646897294394D-09,
     -   0.390449620623303954981650371818D-09,
     -   0.811450270264693438013630664820D-09,
     -   0.151347525435002417978946030127D-08,
     -   0.269381829462281348446234359145D-08,
     -   0.475556539198565231570448230849D-08,
     -   0.862433393351298196083982206932D-08,
     -   0.168434584509839442319457430901D-07,
     -   0.387555526387630181198347920207D-07,
     -   0.135003282458664231686720069847D-06,
     -   0.484887662118194534505870062246D-05,
     -   0.469580603129363623590686483411D-01,
     -   0.230804085432009521165269497595D+00,
     -   0.500025182385124266950747294823D+00,
     -   0.769246277662207845089311389865D+00,
     -   0.953092285630359187928964088225D+00/
        data whts07 /
     -   0.693594132370276611001488553524D-10,
     -   0.175549864602096571565733751941D-09,
     -   0.319932117224675720634977645684D-09,
     -   0.538620734106603140129898268384D-09,
     -   0.896694914165310331748652734054D-09,
     -   0.152681553102798989473462697929D-08,
     -   0.273707418407557910388826330596D-08,
     -   0.536816301477579263072370855472D-08,
     -   0.122950387974222934457882041629D-07,
     -   0.376365800717626860230229188304D-07,
     -   0.225648082649880931728252227742D-06,
     -   0.500181614678192737496775986093D-04,
     -   0.118457531022838262020066935314D+00,
     -   0.239302286357859212458028096629D+00,
     -   0.284430119049181128545698023675D+00,
     -   0.239302282162032013623207686249D+00,
     -   0.118457476034710271798263328944D+00/
        data xs08 /
     -   0.267551320181360573871760072590D-09,
     -   0.147075333437509439708609526218D-08,
     -   0.390337183781774760925149193648D-08,
     -   0.811177607408183076504398152107D-08,
     -   0.151285560418378700893522425237D-07,
     -   0.269241121489492693276720129886D-07,
     -   0.475221548218046947268274403347D-07,
     -   0.861549278812478740430309733455D-07,
     -   0.168152022232344959685867593329D-06,
     -   0.386248587804180161760266627536D-06,
     -   0.133564643655039302913598764511D-05,
     -   0.365223396872440428110605390553D-04,
     -   0.471942054874697039701979843400D-01,
     -   0.230995104264829659177186646264D+00,
     -   0.500149382403010359522139255717D+00,
     -   0.769303605105133429763826203336D+00,
     -   0.953103939565004994864711114334D+00/
        data whts08 /
     -   0.693423550907583438268019850934D-09,
     -   0.175501280948767911924240489338D-08,
     -   0.319825090838963551952432088494D-08,
     -   0.538391783627612489099016259222D-08,
     -   0.896191675718042901366658652713D-08,
     -   0.152563952841773536370089811082D-07,
     -   0.273403249770141885080521406111D-07,
     -   0.535890694688732955737818632764D-07,
     -   0.122581924545929735321550839133D-06,
     -   0.373976995026112335653501689030D-06,
     -   0.220964866898609436369435051007D-05,
     -   0.293995973521532063874263060671D-03,
     -   0.118429846886426091131700606721D+00,
     -   0.239242972743952822868692716147D+00,
     -   0.284359481780392260761625566133D+00,
     -   0.239242834028598603502634014890D+00,
     -   0.118428046201208539228748463111D+00/
        data xs09 /
     -   0.444240013329249378955037019352D-09,
     -   0.734571752501744225333549689330D-08,
     -   0.248158653654769942510699805983D-07,
     -   0.565866591505336905433029263147D-07,
     -   0.109890380537203038095564790717D-06,
     -   0.198186197805847555507086031121D-06,
     -   0.347420631877589576968842773793D-06,
     -   0.611869089300332363202335539359D-06,
     -   0.112030229664412656797389156176D-05,
     -   0.224030346471991242977418400573D-05,
     -   0.536930983469508720440276276799D-05,
     -   0.194651772904640741453942692377D-04,
     -   0.313960965049545904347449330561D-03,
     -   0.236200028826364674813405557817D-01,
     -   0.872521880734281373023357499547D-01,
     -   0.191045443569470792955080279906D+00,
     -   0.323896937508014185935953686307D+00,
     -   0.457101287167819907822818487538D+00,
     -   0.656085504143527073009699599192D+00,
     -   0.845832434493977500692548093877D+00,
     -   0.969187571783460687245085597802D+00/
        data whts09 /
     -   0.232028259685385401483029956181D-08,
     -   0.117771455277646659605067379226D-07,
     -   0.237650096420395119165997030550D-07,
     -   0.409407038156549502105802098589D-07,
     -   0.678170615531768751077448105191D-07,
     -   0.112871795001879465648096755241D-06,
     -   0.194013733714786543133993886410D-06,
     -   0.354282374645335380723863251968D-06,
     -   0.715287211139674728529621576618D-06,
     -   0.170684167493600195845756403025D-05,
     -   0.547490227696358553812653042411D-05,
     -   0.325864468064317831070983051795D-04,
     -   0.152161748957021268070598637725D-02,
     -   0.518312115783755049430002985249D-01,
     -   0.737776569641989882291495805627D-01,
     -   0.135250258139386625422378413901D+00,
     -   0.115605612420830459045139758436D+00,
     -   0.173392896393327976965622571564D+00,
     -   0.206925268931887008585643929188D+00,
     -   0.163529980473499135727053262086D+00,
     -   0.781242063428481198647272711227D-01/
        data xs10 /
     -   0.267112623805969672011378139019D-07,
     -   0.146818649824800144478180492053D-06,
     -   0.389576975276035401907011100778D-06,
     -   0.809334598525632545149764857483D-06,
     -   0.150866915855028332586186542757D-05,
     -   0.268291220567758436812156686992D-05,
     -   0.472963168005986209210648576717D-05,
     -   0.855604426653004831502450563249D-05,
     -   0.166262669170562806491610043069D-04,
     -   0.377639430383533557549305329813D-04,
     -   0.124665464207197460313299663782D-03,
     -   0.146179362279278275613553354053D-02,
     -   0.529535858968769972522837325184D-01,
     -   0.235793903254287286991770665634D+00,
     -   0.503283746529811807401662901410D+00,
     -   0.770752353403482361791075684249D+00,
     -   0.953398598832570286065249914157D+00/
        data whts10 /
     -   0.692269950729827628760737401785D-07,
     -   0.175172773426303950420953764788D-06,
     -   0.319101627610018308458329625951D-06,
     -   0.536844902255552991250916843530D-06,
     -   0.892794210849725419353878830007D-06,
     -   0.151771009713014530224713290326D-05,
     -   0.271356793285303314737400523613D-05,
     -   0.529688467126628097149213702082D-05,
     -   0.120134083645821893251326580587D-04,
     -   0.358461121911529990861793629021D-04,
     -   0.193256521566528910759084929511D-03,
     -   0.569628650164307967446907296446D-02,
     -   0.118247266161753280036450977708D+00,
     -   0.237794247070324611044721238007D+00,
     -   0.282583887638379202736152095005D+00,
     -   0.237741659240177325305367025832D+00,
     -   0.117684016042389773060815720105D+00/
        data xs11 /
     -   0.268785968415476906922960971006D-06,
     -   0.147792442103305612042514249837D-05,
     -   0.392415978695941680319544162074D-05,
     -   0.815962036590220259971960975267D-05,
     -   0.152256586091831308335087832886D-04,
     -   0.270990245675946738338298754883D-04,
     -   0.477757407130543194941273415689D-04,
     -   0.862659571706217990372319569871D-04,
     -   0.166533136500782045412960879950D-03,
     -   0.370707590391638758918726124618D-03,
     -   0.113024905256841082180292743774D-02,
     -   0.720273033408033155133669873452D-02,
     -   0.678317517029567280066717927886D-01,
     -   0.248999872809458370502385119551D+00,
     -   0.512022118773685038649976734062D+00,
     -   0.774808416017802808558468950602D+00,
     -   0.954224835728424373876352735378D+00/
        data whts11 /
     -   0.696665869295077374750691862385D-06,
     -   0.176410103306468203820875700635D-05,
     -   0.321734675059869741021493161742D-05,
     -   0.542045090612278897025974451630D-05,
     -   0.902533215264285800548070040899D-05,
     -   0.153458145022748576427777034133D-04,
     -   0.273808822312716736918386587215D-04,
     -   0.531009378747982782400670778650D-04,
     -   0.118524001315616275658522310630D-03,
     -   0.339023779236614189567257751135D-03,
     -   0.156548735796668360853356737942D-02,
     -   0.170050448696478471994676408238D-01,
     -   0.119897318732465815714694867975D+00,
     -   0.234118947844539277763802401265D+00,
     -   0.277692424625200735311837589500D+00,
     -   0.233548882705098334950345175508D+00,
     -   0.115598394553209006072719379222D+00/
        data xs12 /
     -   0.261184891259299724366726176294D-05,
     -   0.143394884212498383460527861847D-04,
     -   0.379693448777316535102211880165D-04,
     -   0.786339588268385112720260393146D-04,
     -   0.145940232829491861756544689913D-03,
     -   0.257919607999279918259717981172D-03,
     -   0.450369016631686406714998020419D-03,
     -   0.801535297109051628688580951137D-03,
     -   0.150807517259729116301087264407D-02,
     -   0.317105439681351795516369031958D-02,
     -   0.817866962971485158490907685750D-02,
     -   0.284641006714808733857237347899D-01,
     -   0.107172074500017556117826663306D+00,
     -   0.286387677100619827062643084535D+00,
     -   0.537421003985172752730884957727D+00,
     -   0.786716776779057627771358725729D+00,
     -   0.956660130672374106717204445781D+00/
        data whts12 /
     -   0.676726575612908349197366443224D-05,
     -   0.170855372886744092633141045460D-04,
     -   0.310021384873521122805936907033D-04,
     -   0.518687423554430873182122594972D-04,
     -   0.856119301406001070929827336828D-04,
     -   0.143964039866574049189787668092D-03,
     -   0.252959795126706792472716577405D-03,
     -   0.478622845470148431693068275158D-03,
     -   0.101931115724844427180618486957D-02,
     -   0.261908411344557158158324416132D-02,
     -   0.891790188591313655434580005807D-02,
     -   0.390381074210517814566619645193D-01,
     -   0.127410345662673367880786940744D+00,
     -   0.225293401358934204842609097723D+00,
     -   0.263865167844802487987556307969D+00,
     -   0.221312605731700167985098454196D+00,
     -   0.109456192529739209366749356786D+00/
        data xs13 /
     -   0.256912829304367329983713849855D-04,
     -   0.140833141996505706735961174509D-03,
     -   0.371764549497266697040801752068D-03,
     -   0.765903316258618937831180866544D-03,
     -   0.140953627366146198696274014270D-02,
     -   0.245786866053915379745312434818D-02,
     -   0.420022710190486534759834161512D-02,
     -   0.721127239475156086942475432646D-02,
     -   0.127288817780543317404685321759D-01,
     -   0.236552674995237200415849191407D-01,
     -   0.471485124672026728446139031685D-01,
     -   0.995498809721434347226705036678D-01,
     -   0.206667912033276594516233647612D+00,
     -   0.382799265543497951740141703095D+00,
     -   0.605156500301921173503504911537D+00,
     -   0.819069408491987729207650126629D+00,
     -   0.963331978058545971732320611692D+00/
        data whts13 /
     -   0.665429914035210511050460927822D-04,
     -   0.167488597379490949364108032215D-03,
     -   0.302052611047056252783684759954D-03,
     -   0.500145455667973631100317387818D-03,
     -   0.812003632682052330112920748784D-03,
     -   0.133064095685552782672450463520D-02,
     -   0.224435152533685532800798776388D-02,
     -   0.396956849565390144047116092406D-02,
     -   0.751023580749499658665439417386D-02,
     -   0.154447302062529731325853136248D-01,
     -   0.341951327437860707780726576744D-01,
     -   0.752930359638343950108076782025D-01,
     -   0.141896177141122811196657320962D+00,
     -   0.206413520513363653069762548153D+00,
     -   0.228716357367771908758548638689D+00,
     -   0.188478265240210937603666646294D+00,
     -   0.926597507501358750535750718821D-01/
        data xs14 /
     -   0.103934313543029173500112603979D-03,
     -   0.822161903972756462069049475045D-03,
     -   0.254248227785170267342540625864D-02,
     -   0.561563414351037158516745416329D-02,
     -   0.106350546260809461802160606854D-01,
     -   0.185946563910497466658496556433D-01,
     -   0.311241097217215787079438213842D-01,
     -   0.508679073579738394600906305050D-01,
     -   0.820273924423307749564822972991D-01,
     -   0.130892528175046472233293767596D+00,
     -   0.205627887030659280575961655120D+00,
     -   0.313793294680312625017008234625D+00,
     -   0.456618306214180613619629211739D+00,
     -   0.622601331866969184957102745020D+00,
     -   0.786235235187822391136637475601D+00,
     -   0.915220929302643439071831016782D+00,
     -   0.985423213473296220013084886684D+00/
        data whts14 /
     -   0.318903632287647659932294999357D-03,
     -   0.117276498608154218576163066573D-02,
     -   0.232294430192858694466568022370D-02,
     -   0.391995803054014155253602181967D-02,
     -   0.628074969295975623040344124688D-02,
     -   0.990384505730289094663375689704D-02,
     -   0.155865199925509190968160209100D-01,
     -   0.245908893030787724516307028084D-01,
     -   0.387759673545274311110619212708D-01,
     -   0.603523784839766650856754017517D-01,
     -   0.904795186746353474762103765890D-01,
     -   0.126180086046712228911521737316D+00,
     -   0.157583408935672638678546891688D+00,
     -   0.170005677532640759551447775226D+00,
     -   0.151631398756528957822640410996D+00,
     -   0.101936724295618029407187806642D+00,
     -   0.389582649229576848873281289501D-01/
        data xs15 /
     -   0.310086428166403161577596155752D-02,
     -   0.169735575452142246035237580982D-01,
     -   0.445788736327989744375428897195D-01,
     -   0.905724268307542218027482016467D-01,
     -   0.161406416556897149876201817385D+00,
     -   0.263599932281510507315155981520D+00,
     -   0.399682496225075818162618674848D+00,
     -   0.562716279932795483797155415596D+00,
     -   0.733078576654331631408368561989D+00,
     -   0.881217339281124184900660609089D+00,
     -   0.976284397524418574377323986181D+00/
        data whts15 /
     -   0.802968504492989597817656344221D-02,
     -   0.201407203732064087526833796685D-01,
     -   0.358484040779948749740563184967D-01,
     -   0.572399314298669770043903408659D-01,
     -   0.855835549384775233848712389724D-01,
     -   0.119326455838415154043895247794D+00,
     -   0.151746952223865970330782015583D+00,
     -   0.171003751325648473134402321477D+00,
     -   0.164641092351317323451431552926D+00,
     -   0.126304126713996735287630675559D+00,
     -   0.601353256822806636576803452160D-01/
        data xs16 /
     -   0.309728553049385071405276902098D-02,
     -   0.197368846265211593135931641610D-01,
     -   0.563809885687018181180482149789D-01,
     -   0.117654840339587067827946806189D+00,
     -   0.207760754930451640550228476215D+00,
     -   0.328547020571001023011078258644D+00,
     -   0.475872126843748097300446523964D+00,
     -   0.636701895560116176387506783531D+00,
     -   0.789778906001120842387753063892D+00,
     -   0.911212897316431216300422627860D+00,
     -   0.982957987262557236983257932955D+00/
        data whts16 /
     -   0.850639735202691056702864735312D-02,
     -   0.258223974653835121659565101869D-01,
     -   0.482070831497402819692210397901D-01,
     -   0.750783544417883869214377516343D-01,
     -   0.105507596686953625757030595868D+00,
     -   0.135430891180437832829325903275D+00,
     -   0.157039770705858034312092667919D+00,
     -   0.160979263731021897336106401890D+00,
     -   0.141036241997644338334703729438D+00,
     -   0.987176527251396675364884351265D-01,
     -   0.436743505640055122706083175184D-01/
        data xs17 /
     -   0.710512928838859298260084109896D-02,
     -   0.380388118576764403725060446616D-01,
     -   0.958706114169819357909711855694D-01,
     -   0.183174736628131391202202400229D+00,
     -   0.300918769287834734555851965136D+00,
     -   0.445566532100675252216278337138D+00,
     -   0.606385932807841005423475137227D+00,
     -   0.764756701015583232276811415369D+00,
     -   0.896876111633191784627009772158D+00,
     -   0.979583006774319366465495936746D+00/
        data whts17 /
     -   0.183052700103320205866335208494D-01,
     -   0.439402965971366770338213445808D-01,
     -   0.721865103234339068563458960232D-01,
     -   0.102633792396224885852552182349D+00,
     -   0.132314583326036508015330223678D+00,
     -   0.155207778179680105292688681332D+00,
     -   0.163281881953194492694634217509D+00,
     -   0.149390174986985274161300853174D+00,
     -   0.110865534300546170692083812195D+00,
     -   0.518741779264299588146092683099D-01/
        data xs18 /
     -   0.529988547651139170322156430201D-02,
     -   0.391738244772072946344358427755D-01,
     -   0.113987756121909514617828365739D+00,
     -   0.230321685900557339834206096229D+00,
     -   0.382897823066252266230631470409D+00,
     -   0.558429530687391740935317222225D+00,
     -   0.734635945213669075905682183915D+00,
     -   0.883255873595053877247486358995D+00,
     -   0.976842477429295305011769987827D+00/
        data whts18 /
     -   0.157837344040553516991666412112D-01,
     -   0.537430217981115748506598563176D-01,
     -   0.959258886798078592341077294630D-01,
     -   0.135916637340041966247476136011D+00,
     -   0.167043536008614284237675166751D+00,
     -   0.180255523880706170951514738648D+00,
     -   0.167327467753560819884688422006D+00,
     -   0.125193054860443848270475617882D+00,
     -   0.588111352746581246242356917101D-01/
        data xs19 /
     -   0.157338679484540302037218328950D-01,
     -   0.822814491242574270868269784717D-01,
     -   0.198698066625796989994567382076D+00,
     -   0.356660811776159072362028809836D+00,
     -   0.539957436664990204309908127560D+00,
     -   0.723931903688598697862061411806D+00,
     -   0.878716801128286442973987476290D+00,
     -   0.975963815059625654035763022418D+00/
        data whts19 /
     -   0.403152514296860007477334829720D-01,
     -   0.923312712904070257333893028146D-01,
     -   0.139164029934368578544516888254D+00,
     -   0.174073522065522595533290943905D+00,
     -   0.188355762680273711774536955175D+00,
     -   0.174497152235396926605172260027D+00,
     -   0.130206781586590382337782930367D+00,
     -   0.610562287777547787235772364861D-01/
        data xs20 /
     -   0.118540259026739267223236460826D-01,
     -   0.687469774979433086802551700481D-01,
     -   0.179687050459121243324329045290D+00,
     -   0.338101765842085916743327785977D+00,
     -   0.525678270239590432015996896538D+00,
     -   0.715221637544536660318433382811D+00,
     -   0.874913538929842101877802259246D+00,
     -   0.975216926579226450845454504874D+00/
        data whts20 /
     -   0.314723653784141680842935378349D-01,
     -   0.838342085988845462401970360951D-01,
     -   0.136860626511482823863195037215D+00,
     -   0.176882792829588791473555509619D+00,
     -   0.193668888495240726722947235552D+00,
     -   0.179994849448634407232423132918D+00,
     -   0.134327821195529002652610643067D+00,
     -   0.629584475422255337307778676981D-01/
        data xs21 /
     -   0.232438944075952125523273119950D-01,
     -   0.119275001264134411399567693602D+00,
     -   0.278698397896311035541100207455D+00,
     -   0.478146539665980493356721714901D+00,
     -   0.684794836412594935849002808601D+00,
     -   0.861053270277067777639951473448D+00,
     -   0.972421410075918021913685488734D+00/
        data whts21 /
     -   0.592904489028554293675022902651D-01,
     -   0.130720165875034402730327692548D+00,
     -   0.184212590241802372558745401718D+00,
     -   0.209134591870889233364119267715D+00,
     -   0.197760361851548633057002847230D+00,
     -   0.148850357677441381989709411949D+00,
     -   0.700314835804285469325930885741D-01/
        data xs22 /
     -   0.225672625163741684952350874110D-01,
     -   0.116514777567497722455076023699D+00,
     -   0.274348661123247102703434513268D+00,
     -   0.473825548536778690787127890958D+00,
     -   0.681774433079400911300178211283D+00,
     -   0.859642397222122258410963385082D+00,
     -   0.972136413479697081957413518064D+00/
        data whts22 /
     -   0.576598373465209320935814988505D-01,
     -   0.128564102672032505800271979605D+00,
     -   0.183384437969199035522482230734D+00,
     -   0.209943713285192135385956311886D+00,
     -   0.199377992872862324426936887261D+00,
     -   0.150316847043498159981050003495D+00,
     -   0.707530688106949067897210881682D-01/
        data xs23 /
     -   0.332501963354096231259294379164D-01,
     -   0.167182453283988170097494377782D+00,
     -   0.377003993473571890714971965000D+00,
     -   0.615647160513126168970727212938D+00,
     -   0.828431500347486192309133558316D+00,
     -   0.965732681155525064821812487060D+00/
        data whts23 /
     -   0.844024671040197343264225616344D-01,
     -   0.178516051152488241316264108879D+00,
     -   0.233107767157480340894875475367D+00,
     -   0.234847515202174891298475306849D+00,
     -   0.182233546378856231092137355452D+00,
     -   0.868926530049805610718251918183D-01/

        ier = 0

        delta = r0/r1

        if( deltas(1,01) .le. delta .AND. delta .le. deltas(2,01)) then
        nquad = 18
        do l=1,18
        xs(l)   = xs01(l)
        whts(l) = whts01(l)
        end do
        goto 1000
        end if


        if( deltas(1,02) .le. delta .AND. delta .le. deltas(2,02)) then
        nquad = 17
        do l=1,17
        xs(l)   = xs02(l)
        whts(l) = whts02(l)
        end do
        goto 1000
        end if


        if( deltas(1,03) .le. delta .AND. delta .le. deltas(2,03)) then
        nquad = 17
        do l=1,17
        xs(l)   = xs03(l)
        whts(l) = whts03(l)
        end do
        goto 1000
        end if


        if( deltas(1,04) .le. delta .AND. delta .le. deltas(2,04)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs04(l)
        whts(l) = whts04(l)
        end do
        goto 1000
        end if


        if( deltas(1,05) .le. delta .AND. delta .le. deltas(2,05)) then
        nquad = 17
        do l=1,17
        xs(l)   = xs05(l)
        whts(l) = whts05(l)
        end do
        goto 1000
        end if


        if( deltas(1,06) .le. delta .AND. delta .le. deltas(2,06)) then
        nquad = 17
        do l=1,17
        xs(l)   = xs06(l)
        whts(l) = whts06(l)
        end do
        goto 1000
        end if


        if( deltas(1,07) .le. delta .AND. delta .le. deltas(2,07)) then
        nquad = 17
        do l=1,17
        xs(l)   = xs07(l)
        whts(l) = whts07(l)
        end do
        goto 1000
        end if


        if( deltas(1,08) .le. delta .AND. delta .le. deltas(2,08)) then
        nquad = 17
        do l=1,17
        xs(l)   = xs08(l)
        whts(l) = whts08(l)
        end do
        goto 1000
        end if


        if( deltas(1,09) .le. delta .AND. delta .le. deltas(2,09)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs09(l)
        whts(l) = whts09(l)
        end do
        goto 1000
        end if


        if( deltas(1,10) .le. delta .AND. delta .le. deltas(2,10)) then
        nquad = 17
        do l=1,17
        xs(l)   = xs10(l)
        whts(l) = whts10(l)
        end do
        goto 1000
        end if


        if( deltas(1,11) .le. delta .AND. delta .le. deltas(2,11)) then
        nquad = 17
        do l=1,17
        xs(l)   = xs11(l)
        whts(l) = whts11(l)
        end do
        goto 1000
        end if


        if( deltas(1,12) .le. delta .AND. delta .le. deltas(2,12)) then
        nquad = 17
        do l=1,17
        xs(l)   = xs12(l)
        whts(l) = whts12(l)
        end do
        goto 1000
        end if


        if( deltas(1,13) .le. delta .AND. delta .le. deltas(2,13)) then
        nquad = 17
        do l=1,17
        xs(l)   = xs13(l)
        whts(l) = whts13(l)
        end do
        goto 1000
        end if


        if( deltas(1,14) .le. delta .AND. delta .le. deltas(2,14)) then
        nquad = 17
        do l=1,17
        xs(l)   = xs14(l)
        whts(l) = whts14(l)
        end do
        goto 1000
        end if


        if( deltas(1,15) .le. delta .AND. delta .le. deltas(2,15)) then
        nquad = 11
        do l=1,11
        xs(l)   = xs15(l)
        whts(l) = whts15(l)
        end do
        goto 1000
        end if


        if( deltas(1,16) .le. delta .AND. delta .le. deltas(2,16)) then
        nquad = 11
        do l=1,11
        xs(l)   = xs16(l)
        whts(l) = whts16(l)
        end do
        goto 1000
        end if


        if( deltas(1,17) .le. delta .AND. delta .le. deltas(2,17)) then
        nquad = 10
        do l=1,10
        xs(l)   = xs17(l)
        whts(l) = whts17(l)
        end do
        goto 1000
        end if


        if( deltas(1,18) .le. delta .AND. delta .le. deltas(2,18)) then
        nquad = 09
        do l=1,09
        xs(l)   = xs18(l)
        whts(l) = whts18(l)
        end do
        goto 1000
        end if


        if( deltas(1,19) .le. delta .AND. delta .le. deltas(2,19)) then
        nquad = 08
        do l=1,08
        xs(l)   = xs19(l)
        whts(l) = whts19(l)
        end do
        goto 1000
        end if


        if( deltas(1,20) .le. delta .AND. delta .le. deltas(2,20)) then
        nquad = 08
        do l=1,08
        xs(l)   = xs20(l)
        whts(l) = whts20(l)
        end do
        goto 1000
        end if


        if( deltas(1,21) .le. delta .AND. delta .le. deltas(2,21)) then
        nquad = 07
        do l=1,07
        xs(l)   = xs21(l)
        whts(l) = whts21(l)
        end do
        goto 1000
        end if


        if( deltas(1,22) .le. delta .AND. delta .le. deltas(2,22)) then
        nquad = 07
        do l=1,07
        xs(l)   = xs22(l)
        whts(l) = whts22(l)
        end do
        goto 1000
        end if


        if( deltas(1,23) .le. delta .AND. delta .le. deltas(2,23)) then
        nquad = 06
        do l=1,06
        xs(l)   = xs23(l)
        whts(l) = whts23(l)
        end do
        goto 1000
        end if

        ier = 4
        return

 1000 continue

        do i=1,nquad
        xs(i)   = (r1-r0)*xs(i)+r0
        whts(i) = (r1-r0)*whts(i)
        end do

        end



        subroutine pvquad12(ier,r0,r1,nquad,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(*),whts(*)
        dimension deltas(2,23)
        dimension xs01(19),whts01(19)
        dimension xs02(20),whts02(20)
        dimension xs03(19),whts03(19)
        dimension xs04(19),whts04(19)
        dimension xs05(19),whts05(19)
        dimension xs06(19),whts06(19)
        dimension xs07(19),whts07(19)
        dimension xs08(19),whts08(19)
        dimension xs09(19),whts09(19)
        dimension xs10(19),whts10(19)
        dimension xs11(21),whts11(21)
        dimension xs12(20),whts12(20)
        dimension xs13(19),whts13(19)
        dimension xs14(18),whts14(18)
        dimension xs15(13),whts15(13)
        dimension xs16(12),whts16(12)
        dimension xs17(11),whts17(11)
        dimension xs18(10),whts18(10)
        dimension xs19(10),whts19(10)
        dimension xs20(09),whts20(09)
        dimension xs21(08),whts21(08)
        dimension xs22(08),whts22(08)
        dimension xs23(07),whts23(07)
        data deltas /
     -   0.1000000000000000D-14, 0.1000000000000000D-13,
     -   0.1000000000000000D-13, 0.1000000000000000D-12,
     -   0.1000000000000000D-12, 0.1000000000000000D-11,
     -   0.1000000000000000D-11, 0.1000000000000000D-10,
     -   0.1000000000000000D-10, 0.1000000000000000D-09,
     -   0.1000000000000000D-09, 0.1000000000000000D-08,
     -   0.1000000000000000D-08, 0.1000000000000000D-07,
     -   0.1000000000000000D-07, 0.1000000000000000D-06,
     -   0.1000000000000000D-06, 0.1000000000000000D-05,
     -   0.1000000000000000D-05, 0.1000000000000000D-04,
     -   0.1000000000000000D-04, 0.1000000000000000D-03,
     -   0.1000000000000000D-03, 0.1000000000000000D-02,
     -   0.1000000000000000D-02, 0.1000000000000000D-01,
     -   0.1000000000000000D-01, 0.1000000000000000D+00,
     -   0.1000000000000000D+00, 0.2000000000000000D+00,
     -   0.2000000000000000D+00, 0.3000000000000000D+00,
     -   0.3000000000000000D+00, 0.4000000000000000D+00,
     -   0.4000000000000000D+00, 0.5000000000000000D+00,
     -   0.5000000000000000D+00, 0.6000000000000000D+00,
     -   0.6000000000000000D+00, 0.7000000000000000D+00,
     -   0.7000000000000000D+00, 0.8000000000000000D+00,
     -   0.8000000000000000D+00, 0.9000000000000000D+00,
     -   0.9000000000000000D+00, 0.1000000000000000D+01/
        data xs01 /
     -   0.267726599066937389030150986446D-16,
     -   0.147177901358101995293343923148D-15,
     -   0.390641020866043540655022206004D-15,
     -   0.811914483625624981406149767891D-15,
     -   0.151453043236234124801855464262D-14,
     -   0.269621510027113045666131289352D-14,
     -   0.476127449676234779797395203689D-14,
     -   0.863941530965045066906906619305D-14,
     -   0.168917503133443218306677750315D-13,
     -   0.389800463872957218936140831210D-13,
     -   0.137529658142791916142840303085D-12,
     -   0.117625474431579917673560804679D-10,
     -   0.254460440979516546090984334665D-01,
     -   0.129234407440950406543949756015D+00,
     -   0.297077424505563398104123420867D+00,
     -   0.500000000138181635011970445091D+00,
     -   0.702922575770799871892257455845D+00,
     -   0.870765592835412863309259000322D+00,
     -   0.974553956178411614145610885808D+00/
        data whts01 /
     -   0.693884477162844341050634885667D-16,
     -   0.175632565257954830666926746957D-15,
     -   0.320114327687958959577930428970D-15,
     -   0.539010601931464316954066337953D-15,
     -   0.897552119461410047587511145708D-15,
     -   0.152881965350827362705268178029D-14,
     -   0.274226176451239126969587926831D-14,
     -   0.538397110049438380588540707014D-14,
     -   0.123581487424877796327955899689D-13,
     -   0.380488484041673309218465104214D-13,
     -   0.234003044597809320660058044934D-12,
     -   0.276067199388101402592668380725D-09,
     -   0.647424830665424057897895032257D-01,
     -   0.139852695705988185915893021109D+00,
     -   0.190915025199797571861927941054D+00,
     -   0.208979591778980410523886485933D+00,
     -   0.190915025199797571799843606406D+00,
     -   0.139852695705988185629753876765D+00,
     -   0.647424830665424022984791283624D-01/
        data xs02 /
     -   0.367065801119752276681360681741D-17,
     -   0.220323557527590301284937463656D-16,
     -   0.128336323385496070545973347823D-14,
     -   0.369920984421067278333273102118D-14,
     -   0.787819111418487747129390534093D-14,
     -   0.148517745089188057598198616874D-13,
     -   0.265870365967280121371178393553D-13,
     -   0.471049102562717765168396493589D-13,
     -   0.856493720629611767777815836393D-13,
     -   0.167678282232558747922197988923D-12,
     -   0.387197901415236270857385731239D-12,
     -   0.136545888218428242863247869289D-11,
     -   0.105004564018832904339984599725D-09,
     -   0.254460459984852517901917032940D-01,
     -   0.129234409139080420898358467504D+00,
     -   0.297077425876373107776747190397D+00,
     -   0.500000001113260370000691997972D+00,
     -   0.702922576350147630453712078727D+00,
     -   0.870765593087440308132285083886D+00,
     -   0.974553956228035406654143057713D+00/
        data whts02 /
     -  -0.343953708508919829611657674992D-14,
     -   0.394862125534672898278987318812D-14,
     -   0.174580845933785728065244267308D-14,
     -   0.317501099769524758425469622668D-14,
     -   0.534805013548736441825099446282D-14,
     -   0.891081329575800125178114448349D-14,
     -   0.151863066124606102579148964548D-13,
     -   0.272510767639839738989801732214D-13,
     -   0.535154965087214039246884533397D-13,
     -   0.122839325058681583744097975397D-12,
     -   0.378065773739647112247837702591D-12,
     -   0.232054350884661782696339595505D-11,
     -   0.222358340277399614630726902275D-08,
     -   0.647424829402845894921625179384D-01,
     -   0.139852695433253424502767702454D+00,
     -   0.190915024827483213172064994126D+00,
     -   0.208979591371437299354442479165D+00,
     -   0.190915024827483209182709789966D+00,
     -   0.139852695433253406116315339867D+00,
     -   0.647424829402843651509523816640D-01/
        data xs03 /
     -   0.267707150930178585956757368407D-14,
     -   0.147166520427584924885045710383D-13,
     -   0.390607304521418548009679285278D-13,
     -   0.811832703689880910074798602926D-13,
     -   0.151434452295653622990377309311D-12,
     -   0.269579274377348277602527582393D-12,
     -   0.476026819411999420575934053479D-12,
     -   0.863675577832376906342532496454D-12,
     -   0.168832260147365626624285086293D-11,
     -   0.389403182385385743560256225339D-11,
     -   0.137077521356001532491052846283D-10,
     -   0.934879201199613151354754341469D-09,
     -   0.254460610512291060910707385319D-01,
     -   0.129234422588736053545938457284D+00,
     -   0.297077436733562488372025162334D+00,
     -   0.500000008836151829216868538414D+00,
     -   0.702922580938741059903063265264D+00,
     -   0.870765595083566922463623610438D+00,
     -   0.974553956621069479169321586674D+00/
        data whts03 /
     -   0.693833334437985805982289283984D-14,
     -   0.175617997295618849719867550388D-13,
     -   0.320082228233369280379787702723D-13,
     -   0.538941912237895458638766004106D-13,
     -   0.897401066207480285041314693980D-13,
     -   0.152846640874808121827870896140D-12,
     -   0.274134704540361440276100612924D-12,
     -   0.538118171476328588479022362476D-12,
     -   0.123469965455365724785747050299D-11,
     -   0.379757161339054343887198729917D-11,
     -   0.232495842373371753773727040451D-10,
     -   0.176428411980749893870087723880D-07,
     -   0.647424819402999723095817356053D-01,
     -   0.139852693273120157708038582849D+00,
     -   0.190915021878651422045115560261D+00,
     -   0.208979588143583946895817794009D+00,
     -   0.190915021878651173891162502951D+00,
     -   0.139852693273119013996868366236D+00,
     -   0.647424819402860174025113806202D-01/
        data xs04 /
     -   0.267692958568932750313587054729D-13,
     -   0.147158215220186953264742175855D-12,
     -   0.390582700515927560590652129911D-12,
     -   0.811773027566971496458937469216D-12,
     -   0.151420886724910304906773818813D-11,
     -   0.269548457409236823856988890578D-11,
     -   0.475953402267853031667538937306D-11,
     -   0.863481578809340097203216478024D-11,
     -   0.168770101811701542786710600126D-10,
     -   0.389113759342528741658621834426D-10,
     -   0.136749496938551808922272905894D-09,
     -   0.814203511353455817001905719482D-08,
     -   0.254461756292899323890825498168D-01,
     -   0.129234524964642684007234767923D+00,
     -   0.297077519376190687946568560483D+00,
     -   0.500000067621172150534743975674D+00,
     -   0.702922615866147266395681393079D+00,
     -   0.870765610277662299609229676825D+00,
     -   0.974553959612762082238702839976D+00/
        data whts04 /
     -   0.693796012891710056658810997703D-13,
     -   0.175607366466162446876156092817D-12,
     -   0.320058804654912963426089858434D-12,
     -   0.538891790191386370847203472559D-12,
     -   0.897290851084508924651961989880D-12,
     -   0.152820868887880403850615689556D-11,
     -   0.274067978212985166052562425984D-11,
     -   0.537914745644744295356478305160D-11,
     -   0.123388677440172952121151498870D-10,
     -   0.379224850994299377065912811436D-10,
     -   0.231405585001911865767974389448D-09,
     -   0.134948142995312258986054453528D-06,
     -   0.647424743293131605015284753473D-01,
     -   0.139852676830697032289735883418D+00,
     -   0.190914999432778228953194266079D+00,
     -   0.208979563573848145900583195586D+00,
     -   0.190914999432763931707449314255D+00,
     -   0.139852676830631138114837867775D+00,
     -   0.647424743285091650339106724353D-01/
        data xs05 /
     -   0.267673633014764320478072403409D-12,
     -   0.147146906229909907845744169451D-11,
     -   0.390549198435789250231935058128D-11,
     -   0.811691771642792594279734054029D-11,
     -   0.151402416366579939397022160638D-10,
     -   0.269506500733704510412880431992D-10,
     -   0.475853456000728291114930111126D-10,
     -   0.863217524185929931430937015184D-10,
     -   0.168685527391668653531217115338D-09,
     -   0.388720331797287891933293162326D-09,
     -   0.136305442691831686937756205711D-08,
     -   0.693545497096691905795641486755D-07,
     -   0.254470117881691397043140423446D-01,
     -   0.129235272085449467272292407180D+00,
     -   0.297078122488211677219778144069D+00,
     -   0.500000496624622782610245988791D+00,
     -   0.702922870760698983567347254953D+00,
     -   0.870765721161721401545287738222D+00,
     -   0.974553981445653619913774032751D+00/
        data whts05 /
     -   0.693745192730803432165093765647D-12,
     -   0.175592890881115578704395202621D-11,
     -   0.320026910621694981464312099782D-11,
     -   0.538823545863033260031523621404D-11,
     -   0.897140795301225532818805445832D-11,
     -   0.152785783981452289978966036567D-10,
     -   0.273977152967267242450271579793D-10,
     -   0.537637922211875224672870254153D-10,
     -   0.123278119358370747294866465590D-09,
     -   0.378501888769029781377469632504D-09,
     -   0.229934027413106846570294035673D-08,
     -   0.990284973682132188581476677853D-06,
     -   0.647424188214083770524758981173D-01,
     -   0.139852556839490422312917415794D+00,
     -   0.190914835627111046381656680462D+00,
     -   0.208979384268090228279281487105D+00,
     -   0.190914835626356609204107503128D+00,
     -   0.139852556836013348151503912690D+00,
     -   0.647424187789863316919385541801D-01/
        data xs06 /
     -   0.267645783065090811702762899433D-11,
     -   0.147130609115294175570037660653D-10,
     -   0.390500920550964351935037513785D-10,
     -   0.811574682881001280630429754877D-10,
     -   0.151375802303634261956208392247D-09,
     -   0.269446050227807940960682016472D-09,
     -   0.475709474824251013585492470742D-09,
     -   0.862837224020524065865592369449D-09,
     -   0.168563782163962409798205356328D-08,
     -   0.388154746928156136925164163034D-08,
     -   0.135670818119044125678264859628D-07,
     -   0.572941316375443560794851470434D-06,
     -   0.254527662189785851974052934601D-01,
     -   0.129240414213055096683076624685D+00,
     -   0.297082273513039106235598708345D+00,
     -   0.500003449327834961885562831307D+00,
     -   0.702924625126963330849154654599D+00,
     -   0.870766484345559309979866658847D+00,
     -   0.974554131715366059135584201382D+00/
        data whts06 /
     -   0.693671956288872503341423084688D-11,
     -   0.175572030759020542659984584431D-10,
     -   0.319980951457212573093752511350D-10,
     -   0.538725211982816441884018513995D-10,
     -   0.896924597043718761075689131952D-10,
     -   0.152735240608106014518201987195D-09,
     -   0.273846336972390645137839570823D-09,
     -   0.537239359241760874829137623440D-09,
     -   0.123119062105429596361851454213D-08,
     -   0.377463865984477823399884985785D-08,
     -   0.227839714397435038265131186346D-07,
     -   0.686751821643480333410709416795D-05,
     -   0.647420384306829999268656556144D-01,
     -   0.139851731109807522171553435926D+00,
     -   0.190913708230922468551571781495D+00,
     -   0.208978150166790827613947476311D+00,
     -   0.190913708195628616383723501132D+00,
     -   0.139851730947153649909173157894D+00,
     -   0.647420364471188244878297690123D-01/
        data xs07 /
     -   0.267602262379119164069325405119D-10,
     -   0.147105142354863051581515806137D-09,
     -   0.390425481612002615299705153577D-09,
     -   0.811391730514156341208275515926D-09,
     -   0.151334221084682294740544854754D-08,
     -   0.269351615924027651182069994497D-08,
     -   0.475484598134095644956369400448D-08,
     -   0.862243475697419243717032278207D-08,
     -   0.168373853385858097682754026562D-07,
     -   0.387274210189799637085883772775D-07,
     -   0.134691614447637028617121241937D-06,
     -   0.452582815290436361367762887835D-05,
     -   0.254891749804842904614283423920D-01,
     -   0.129272966413944947918242186113D+00,
     -   0.297108553509762339042365682213D+00,
     -   0.500022143228101715572105803126D+00,
     -   0.702935732329266719390258998855D+00,
     -   0.870771316219737037599623843068D+00,
     -   0.974555083105878389029404828758D+00/
        data whts07 /
     -   0.693557511328266843004962761229D-10,
     -   0.175539434287683043699941032406D-09,
     -   0.319909139170443200734424974577D-09,
     -   0.538571576886778742117991658943D-09,
     -   0.896586856622132088830199875441D-09,
     -   0.152656298199078914851979213220D-08,
     -   0.273642083310510991205705500365D-08,
     -   0.536617402284263751331511383638D-08,
     -   0.122871143763596289890503632615D-07,
     -   0.375850899879358054523412686816D-07,
     -   0.224628833917868045610875511566D-06,
     -   0.439145532369791960474063650179D-04,
     -   0.647396937100024277512999650233D-01,
     -   0.139846508494525552840706604686D+00,
     -   0.190906571696474187259887056781D+00,
     -   0.208970337198262216617884838370D+00,
     -   0.190906570306225699875473250733D+00,
     -   0.139846502089806428555999961906D+00,
     -   0.647396158213076297008205302935D-01/
        data xs08 /
     -   0.267525633923575996054690062731D-09,
     -   0.147060303482539757331613042543D-08,
     -   0.390292665338068081325552287793D-08,
     -   0.811069657913516946845625467782D-08,
     -   0.151261030991185973674712238187D-07,
     -   0.269185431206971462934992662979D-07,
     -   0.475089002046989769798764760218D-07,
     -   0.861199638532793573334769942135D-07,
     -   0.168040393038427346003548798928D-06,
     -   0.385733655894452452220054474172D-06,
     -   0.133004351037130094192701440027D-05,
     -   0.333462189373843529696579591325D-04,
     -   0.256925315164771975774084664268D-01,
     -   0.129455296100193219443464471763D+00,
     -   0.297255809765649381843966666296D+00,
     -   0.500126904804429316640532774862D+00,
     -   0.702997980912974415636375282713D+00,
     -   0.870798396377383792573149726823D+00,
     -   0.974560415217051443448027990065D+00/
        data whts08 /
     -   0.693356005316175507514827464432D-09,
     -   0.175482044222727946140610213708D-08,
     -   0.319782717679076155708636675470D-08,
     -   0.538301151224417417013140974949D-08,
     -   0.895992500707029691308826806573D-08,
     -   0.152517420914520089653971318127D-07,
     -   0.273282948536014628479137148798D-07,
     -   0.535524886722038743710591945575D-07,
     -   0.122436529758285914141091474293D-06,
     -   0.373038173428707379031140293915D-06,
     -   0.219155986049615696043943198156D-05,
     -   0.248416294865046322460586609953D-03,
     -   0.647283981090189543353975034330D-01,
     -   0.139817393525901202872301253156D+00,
     -   0.190866611928886815216712867439D+00,
     -   0.208926560844560707610153651540D+00,
     -   0.190866569294638643980786245913D+00,
     -   0.139817197498771137234205361339D+00,
     -   0.647260493473280483716951253090D-01/
        data xs09 /
     -   0.267365856239854883761462160991D-08,
     -   0.146966815202600741202345756764D-07,
     -   0.390015774509810595467807100715D-07,
     -   0.810398324866446883640796044727D-07,
     -   0.151108512286897903576440072653D-06,
     -   0.268839266823302489780818231168D-06,
     -   0.474265541395849513087798109265D-06,
     -   0.859029571310836900887877492482D-06,
     -   0.167348995944714409875028405135D-05,
     -   0.382561907226745908669125134562D-05,
     -   0.129633617293237658910123444090D-04,
     -   0.219737518238765011730592783420D-03,
     -   0.266407490242381113626853550700D-01,
     -   0.130315460187093700368733153237D+00,
     -   0.297951692629104597638699950242D+00,
     -   0.500622237573650205885846166681D+00,
     -   0.703292370949324539936702525650D+00,
     -   0.870926479440608672733091184689D+00,
     -   0.974585636120302762743726112173D+00/
        data whts09 /
     -   0.692935851041698342699940128941D-08,
     -   0.175362394176970026422656021702D-07,
     -   0.319519195260196398383980291653D-07,
     -   0.537737611531187412057622648426D-07,
     -   0.894754433348040422654990573046D-07,
     -   0.152228319094498746537218316252D-06,
     -   0.272536130657862888062959933519D-06,
     -   0.533257413162229069810083710322D-06,
     -   0.121538148604017965495070526337D-05,
     -   0.367285449943278955194642771563D-05,
     -   0.208455139697347486605653906841D-04,
     -   0.116490760215548036624415343951D-02,
     -   0.647093496064023260818253606836D-01,
     -   0.139682824065427837381260650251D+00,
     -   0.190678364826961117902699202097D+00,
     -   0.208719740549951243834318474344D+00,
     -   0.190677438682533985815325527413D+00,
     -   0.139678602669550082953221311329D+00,
     -   0.646618805584778613001240687325D-01/
        data xs10 /
     -   0.266960131887750154631237351427D-07,
     -   0.146729447432855230350875698989D-06,
     -   0.389312904900050610868330074182D-06,
     -   0.808694827210574767408223348979D-06,
     -   0.150721729173241271753499884473D-05,
     -   0.267962251907918497231212166525D-05,
     -   0.472182725172581498284218461122D-05,
     -   0.853557514757898963564628446837D-05,
     -   0.165616822271600560382309856832D-04,
     -   0.374749827891201280585843001005D-04,
     -   0.121890015499813578903667135533D-03,
     -   0.122724033608732924406292748338D-02,
     -   0.301732226812247969527371542851D-01,
     -   0.133632290942656498446085631826D+00,
     -   0.300650551890976336556528769365D+00,
     -   0.502546877194308323582419703896D+00,
     -   0.704437148224193689341800651322D+00,
     -   0.871424737733136834547082186840D+00,
     -   0.974683765081998901703748388479D+00/
        data whts10 /
     -   0.691868978470101307825800241767D-07,
     -   0.175058641044824070390787343812D-06,
     -   0.318850465769376272784476277055D-06,
     -   0.536308436689800659846036207511D-06,
     -   0.891617621565849238861867525712D-06,
     -   0.151496955868470354927236362336D-05,
     -   0.270651703298588612624966031325D-05,
     -   0.527562950202947777628317325671D-05,
     -   0.119304150316709277418588245778D-04,
     -   0.353342049888346483839094679892D-04,
     -   0.185063871635735755536152703874D-03,
     -   0.422090504254774105080485890166D-02,
     -   0.649667627342452791786456582095D-01,
     -   0.139199612109401368350460692747D+00,
     -   0.189956171547462853864225897888D+00,
     -   0.207918339637122832869501637871D+00,
     -   0.189942589176503067585813058476D+00,
     -   0.139139579850039678672887044886D+00,
     -   0.644122232728643201681747590797D-01/
        data xs11 /
     -   0.718906309928156409591955446743D-07,
     -   0.783702764665118154351383007846D-06,
     -   0.256120115668393176526319019423D-05,
     -   0.580057150660177466077551306557D-05,
     -   0.112550056578042140541974862333D-04,
     -   0.203315167697954411619659599914D-04,
     -   0.357553104755147558852578169855D-04,
     -   0.632599025882066791641852713679D-04,
     -   0.116565638318790984073519941235D-03,
     -   0.235367328585861074989938589221D-03,
     -   0.572811755894548541474404899485D-03,
     -   0.204840001390649922918794826000D-02,
     -   0.123143753889078335057547900012D-01,
     -   0.580528102216436896885294763592D-01,
     -   0.157302911412087372941233144266D+00,
     -   0.304191338447469518522022430202D+00,
     -   0.477780206500873093805009260751D+00,
     -   0.650733520859111084613660168550D+00,
     -   0.799665219404709567324251722131D+00,
     -   0.912943493006575784052286060726D+00,
     -   0.982635671246617160250983380896D+00/
        data whts11 /
     -   0.271957682625686078164928104271D-06,
     -   0.119847698840073875811828317514D-05,
     -   0.241944096015547292728370331247D-05,
     -   0.418043381463055986105579556634D-05,
     -   0.695309556865518129070421757944D-05,
     -   0.116299523312146935415463168992D-04,
     -   0.201056610724561012733019433855D-04,
     -   0.369658096310773492771887087601D-04,
     -   0.753152435917550208259283658124D-04,
     -   0.182203532815535429924234116353D-03,
     -   0.593831418450635806412214025669D-03,
     -   0.319954673483118309114138403033D-02,
     -   0.226748605853228439617270628742D-01,
     -   0.718002630193320919365173182256D-01,
     -   0.125577973695277783611651062596D+00,
     -   0.164568181833669700897972528983D+00,
     -   0.177881347775055110265400614065D+00,
     -   0.163937255781683519697623986859D+00,
     -   0.132172449677111067292808764092D+00,
     -   0.932598758666724795111644412628D-01,
     -   0.439931700081370776938230966073D-01/
        data xs12 /
     -   0.249920437579098650340064741953D-05,
     -   0.136883331118721748322993681276D-04,
     -   0.360837041179598771032029546396D-04,
     -   0.742166601135759720080057576038D-04,
     -   0.136385226393428399196531349051D-03,
     -   0.237676108979358268839945121627D-03,
     -   0.406705518934784663978943161624D-03,
     -   0.702011937974034249751539073657D-03,
     -   0.125683588274824674401808556318D-02,
     -   0.242074840478174976027215134708D-02,
     -   0.529073729889896938103613325883D-02,
     -   0.139575462991500879958356572544D-01,
     -   0.427598278665852944499129184507D-01,
     -   0.118049785857632132465406082680D+00,
     -   0.251502318490161593673248105616D+00,
     -   0.428847974067433368593988889214D+00,
     -   0.621832007160948045783901865292D+00,
     -   0.797885772004312286106348734808D+00,
     -   0.927143537800680561656760043814D+00,
     -   0.990169899546463731883298337693D+00/
        data whts12 /
     -   0.647186635230637834893613834844D-05,
     -   0.162633312195446426164120500167D-04,
     -   0.292605450308442919872208213551D-04,
     -   0.483422932325822549909951976209D-04,
     -   0.784084797294924814305189825824D-04,
     -   0.128704142229911856998563770969D-03,
     -   0.218483792341921619297674462648D-03,
     -   0.392451616843942765908815561877D-03,
     -   0.768811269989733423490962215015D-03,
     -   0.171412125443233720444961283037D-02,
     -   0.458920447156106933422024279170D-02,
     -   0.149490751252164749943629940620D-01,
     -   0.477395422709091864842586341781D-01,
     -   0.104915942699024580648598146368D+00,
     -   0.159345706909894411965330717533D+00,
     -   0.190487912434754251403568116430D+00,
     -   0.189938867858517100416137674384D+00,
     -   0.157083773996068107222320960287D+00,
     -   0.979127348646974264279894920010D-01,
     -   0.296359207779547741836933099356D-01/
        data xs13 /
     -   0.530119953270537367898332928928D-05,
     -   0.742941490987646011868327040960D-04,
     -   0.247969936212736048551726773989D-03,
     -   0.562824838349294292216657696542D-03,
     -   0.108815256518242053636902893781D-02,
     -   0.194974526567398577874716408329D-02,
     -   0.337978869109056675324463308492D-02,
     -   0.582817119473753301429838242611D-02,
     -   0.102245877163008490335622318259D-01,
     -   0.186139193127938941758077517504D-01,
     -   0.356011257796542704477494401537D-01,
     -   0.706679433043649355277111052013D-01,
     -   0.138175782082083802294429121079D+00,
     -   0.249249541568732963940129514410D+00,
     -   0.401245471330991065923441352047D+00,
     -   0.576671952212809745855471022321D+00,
     -   0.749282310607980573605468315983D+00,
     -   0.891118927462456585701326323010D+00,
     -   0.978578772307224524702740472124D+00/
        data whts13 /
     -   0.242463419308074381195816367028D-04,
     -   0.117232992631347706018049509481D-03,
     -   0.235986622514182112223227616421D-03,
     -   0.404910950194892361657206060865D-03,
     -   0.665950136376165597206474551947D-03,
     -   0.109426116466308326608848556797D-02,
     -   0.183691169104766474820391272632D-02,
     -   0.320531030712247403204889774289D-02,
     -   0.590554705421489747376358301571D-02,
     -   0.115964213582109590043843889003D-01,
     -   0.239352087097290783037325775586D-01,
     -   0.487409053500640787485055277464D-01,
     -   0.883501392944492946918662336795D-01,
     -   0.133357413136621846809279588407D+00,
     -   0.167600977354426914525394075153D+00,
     -   0.178757571268255520370078714458D+00,
     -   0.161702020526203936855502333535D+00,
     -   0.117957055636838470987736179467D+00,
     -   0.545119301045043849681909626674D-01/
        data xs14 /
     -   0.602498070991024727628516951762D-04,
     -   0.748514704645285337550265892620D-03,
     -   0.246257796140971956987789096100D-02,
     -   0.552840946416676396578521483625D-02,
     -   0.105253245050336581125534711273D-01,
     -   0.184090333549263431000655322060D-01,
     -   0.306992447925343286305182335304D-01,
     -   0.497465500337204226451983028751D-01,
     -   0.790338753240152784531617888627D-01,
     -   0.123307466614043426024571717790D+00,
     -   0.188062944147743251865462269722D+00,
     -   0.277831919412302842987427063739D+00,
     -   0.393401669602091064258237452372D+00,
     -   0.529345899019802953817751072570D+00,
     -   0.673533182546465559113941576276D+00,
     -   0.809053099443727658252891105188D+00,
     -   0.917706674796012984751076831742D+00,
     -   0.983872611775696051887828043224D+00/
        data whts14 /
     -   0.251410957906269378511357029416D-03,
     -   0.116278843080785568207127443365D-02,
     -   0.231787231259202373237185998290D-02,
     -   0.390841542069543596043356821304D-02,
     -   0.624161270125055723140517137860D-02,
     -   0.977457772418329180526449156718D-02,
     -   0.151921630238183010403333939330D-01,
     -   0.234795891751147968056768153345D-01,
     -   0.358929931072263017126507709079D-01,
     -   0.535986073802542547650515723036D-01,
     -   0.767177052249140997053885628764D-01,
     -   0.102997422284546857055495301930D+00,
     -   0.127246873585457803725006927204D+00,
     -   0.142585035193556162667553761866D+00,
     -   0.142909874234791786292417603047D+00,
     -   0.125030283664949986818457418626D+00,
     -   0.896157823329024753804929707052D-01,
     -   0.410769932450317402414171786622D-01/
        data xs15 /
     -   0.286509910043969780133332027961D-02,
     -   0.155757644604089730830601760344D-01,
     -   0.403891349645962559953226073862D-01,
     -   0.805240155839306267767415700392D-01,
     -   0.140063816228980266315763371800D+00,
     -   0.222730793812179119508246373848D+00,
     -   0.329848414650067677009020741699D+00,
     -   0.458154049674790030311855590890D+00,
     -   0.598659745749299212750383142944D+00,
     -   0.737393057152962158777512465476D+00,
     -   0.857925017029945142851030636503D+00,
     -   0.945133042516753596247196370635D+00,
     -   0.990599885386466980778716927531D+00/
        data whts15 /
     -   0.740752525124305933066103646374D-02,
     -   0.183304665211804323599205543914D-01,
     -   0.318395684072201681290075517972D-01,
     -   0.491308109596077281389077944783D-01,
     -   0.706047508219881292212118880090D-01,
     -   0.949895627947452835787877267072D-01,
     -   0.118723154784420241945073826626D+00,
     -   0.136377326755961256686681177215D+00,
     -   0.142261359927746142856120226787D+00,
     -   0.132404316675686625680875742723D+00,
     -   0.106056866648198879393193397808D+00,
     -   0.668032256774483491444687741944D-01,
     -   0.250710647745537035350903027996D-01/
        data xs16 /
     -   0.376674291496084864966347960828D-02,
     -   0.209220865216864918481219562669D-01,
     -   0.550299080797356846798127092639D-01,
     -   0.109459402652929480245300265847D+00,
     -   0.186987349915873608594970208594D+00,
     -   0.288638896301829307851533049741D+00,
     -   0.412062198521445749676506881579D+00,
     -   0.550270101837446261136044734433D+00,
     -   0.691593538337182358914478825608D+00,
     -   0.821134064097494062241455854090D+00,
     -   0.923345305327266366465297838017D+00,
     -   0.985023040954599407728629472919D+00/
        data whts16 /
     -   0.980641855009886232007836637406D-02,
     -   0.250559550961576197438194277908D-01,
     -   0.437302091565530949895953739815D-01,
     -   0.656111232696054182817151681087D-01,
     -   0.896450580958962500828875248180D-01,
     -   0.113289290766857901284782587673D+00,
     -   0.132394284163100430319732295955D+00,
     -   0.142051434431836968370154162538D+00,
     -   0.138067448718041336129932268366D+00,
     -   0.118371029049313896716751173516D+00,
     -   0.838040194842340786513396466569D-01,
     -   0.381737292183041431092120042213D-01/
        data xs17 /
     -   0.778368824980186611613114741233D-02,
     -   0.412983718034507181998277830446D-01,
     -   0.102450191654194835329712285285D+00,
     -   0.191523718118357994544062299181D+00,
     -   0.306586000992522059866922112881D+00,
     -   0.441905460737884676419719832481D+00,
     -   0.587261917021399122822154058130D+00,
     -   0.728705509589107126343533466334D+00,
     -   0.850788254856145219337440349169D+00,
     -   0.939982484603216789756530576002D+00,
     -   0.988947902650984384738685673510D+00/
        data whts17 /
     -   0.200111012543461255255434581265D-01,
     -   0.471891376724191505832444201805D-01,
     -   0.751904796544568950188645776282D-01,
     -   0.102663026501639351910471524367D+00,
     -   0.126519220292683774040685450171D+00,
     -   0.142403640599617524011167964331D+00,
     -   0.145946174521533176393768464488D+00,
     -   0.134309015746636306623468461808D+00,
     -   0.107527851227008124291726237890D+00,
     -   0.695719580813518094264071952286D-01,
     -   0.286683944483077621746522457810D-01/
        data xs18 /
     -   0.997098718755707617181911911435D-02,
     -   0.524770724633077640146393552869D-01,
     -   0.128369289721710683268931753017D+00,
     -   0.235459033133174867016384305515D+00,
     -   0.368435910470727875095815360604D+00,
     -   0.517839137927429690522618683351D+00,
     -   0.670126660885701505591483493335D+00,
     -   0.809089677725004329825387371219D+00,
     -   0.918313528002299588270589765234D+00,
     -   0.984054509168588106362777091623D+00/
        data whts18 /
     -   0.255852633284893536040097704000D-01,
     -   0.593757302584501612599554817522D-01,
     -   0.920691152825024714815815113131D-01,
     -   0.121244378280209715588465089565D+00,
     -   0.143132034543318582558106669009D+00,
     -   0.153393026122470918792233301276D+00,
     -   0.148437886857287157121784790064D+00,
     -   0.126702843806542811874246432584D+00,
     -   0.894085177177049915115071696193D-01,
     -   0.406512038030238362081097844182D-01/
        data xs19 /
     -   0.461732218706212294382670226045D-02,
     -   0.368692151431274655801409407409D-01,
     -   0.107574000279418197385090903570D+00,
     -   0.214206421728040746021797121393D+00,
     -   0.349862904356744751473831955155D+00,
     -   0.503544164343300937230140479558D+00,
     -   0.660487789090811610454761659800D+00,
     -   0.803632448869660412765452977661D+00,
     -   0.916024869674320545396064090721D+00,
     -   0.983613401211981751471327256548D+00/
        data whts19 /
     -   0.144370077801060251049583156263D-01,
     -   0.514165218941931827276738339606D-01,
     -   0.894913968289652260698795363892D-01,
     -   0.122630410701586942795661137948D+00,
     -   0.146844966827761183932568884020D+00,
     -   0.158032153590255207418411984095D+00,
     -   0.152969958314126492389297453461D+00,
     -   0.130445786133814451707213968684D+00,
     -   0.919525256130767991413093674620D-01,
     -   0.417792723161144887130255183535D-01/
        data xs20 /
     -   0.156889701587204112794149004524D-01,
     -   0.811706502682805995053452538695D-01,
     -   0.192742211287158159514373036518D+00,
     -   0.339417563007075397110350047426D+00,
     -   0.505328321658836076411803246627D+00,
     -   0.671108982196823469161566721220D+00,
     -   0.816622714502513470938348666487D+00,
     -   0.924813898113130502657895473644D+00,
     -   0.985845573790254769432661340124D+00/
        data whts20 /
     -   0.400941704301966864071629159995D-01,
     -   0.899225556329085365788813941325D-01,
     -   0.131385705553603108636714562915D+00,
     -   0.159297856412149476985455936418D+00,
     -   0.169250023233432195502019777765D+00,
     -   0.158889667049668643789840366496D+00,
     -   0.129219775824413643329530473531D+00,
     -   0.854857548131586769189306414553D-01,
     -   0.364544910504690318514639312874D-01/
        data xs21 /
     -   0.190199675916229066263918701521D-01,
     -   0.977861927260432189017624269018D-01,
     -   0.229667812764747600661635272922D+00,
     -   0.398378724657822076654903896261D+00,
     -   0.582014186228040908266956544376D+00,
     -   0.755631730510426024198509946775D+00,
     -   0.894783467296951025950346296360D+00,
     -   0.979394920911282499657608024587D+00/
        data whts21 /
     -   0.485337120749998721811603565175D-01,
     -   0.107478562812679649628360638312D+00,
     -   0.153545688394080695084385902479D+00,
     -   0.180182069178932420139572718472D+00,
     -   0.182871264528208702749295316900D+00,
     -   0.160222743514375468073051841001D+00,
     -   0.114674607871354822262339382884D+00,
     -   0.524913516253683698818338434337D-01/
        data xs22 /
     -   0.220922298057625552310965966946D-01,
     -   0.112869990283335657920416989496D+00,
     -   0.262212841080976853249841875887D+00,
     -   0.448009797227360684865433696485D+00,
     -   0.642266581695416884038659242118D+00,
     -   0.815251898200548535109270431970D+00,
     -   0.940136796147598816962008127573D+00,
     -   0.997203804126895081110133214774D+00/
        data whts22 /
     -   0.562887244066005285581484763056D-01,
     -   0.123077585304808483097566510366D+00,
     -   0.171869249119174501933630258475D+00,
     -   0.194988726846591661731089879493D+00,
     -   0.188490231455041639659611678632D+00,
     -   0.152934801065166485106883850286D+00,
     -   0.935108456243787951381632183122D-01,
     -   0.188398361782379047749061281301D-01/
        data xs23 /
     -   0.254460438286207377369051579761D-01,
     -   0.129234407200302780068067613360D+00,
     -   0.297077424311301416546696793962D+00,
     -   0.500000000000000000000000000000D+00,
     -   0.702922575688698583453303206038D+00,
     -   0.870765592799697219931932386640D+00,
     -   0.974553956171379262263094842024D+00/
        data whts23 /
     -   0.647424830844348466353132383465D-01,
     -   0.139852695744638333950750134285D+00,
     -   0.190915025252559472475207068916D+00,
     -   0.208979591836734693877575300384D+00,
     -   0.190915025252559472475207068916D+00,
     -   0.139852695744638333950750134285D+00,
     -   0.647424830844348466353132383465D-01/

        ier = 0

        delta = r0/r1

        if( deltas(1,01) .le. delta .AND. delta .le. deltas(2,01)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs01(l)
        whts(l) = whts01(l)
        end do
        goto 1000
        end if


        if( deltas(1,02) .le. delta .AND. delta .le. deltas(2,02)) then
        nquad = 20
        do l=1,20
        xs(l)   = xs02(l)
        whts(l) = whts02(l)
        end do
        goto 1000
        end if


        if( deltas(1,03) .le. delta .AND. delta .le. deltas(2,03)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs03(l)
        whts(l) = whts03(l)
        end do
        goto 1000
        end if


        if( deltas(1,04) .le. delta .AND. delta .le. deltas(2,04)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs04(l)
        whts(l) = whts04(l)
        end do
        goto 1000
        end if


        if( deltas(1,05) .le. delta .AND. delta .le. deltas(2,05)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs05(l)
        whts(l) = whts05(l)
        end do
        goto 1000
        end if


        if( deltas(1,06) .le. delta .AND. delta .le. deltas(2,06)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs06(l)
        whts(l) = whts06(l)
        end do
        goto 1000
        end if


        if( deltas(1,07) .le. delta .AND. delta .le. deltas(2,07)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs07(l)
        whts(l) = whts07(l)
        end do
        goto 1000
        end if


        if( deltas(1,08) .le. delta .AND. delta .le. deltas(2,08)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs08(l)
        whts(l) = whts08(l)
        end do
        goto 1000
        end if


        if( deltas(1,09) .le. delta .AND. delta .le. deltas(2,09)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs09(l)
        whts(l) = whts09(l)
        end do
        goto 1000
        end if


        if( deltas(1,10) .le. delta .AND. delta .le. deltas(2,10)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs10(l)
        whts(l) = whts10(l)
        end do
        goto 1000
        end if


        if( deltas(1,11) .le. delta .AND. delta .le. deltas(2,11)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs11(l)
        whts(l) = whts11(l)
        end do
        goto 1000
        end if


        if( deltas(1,12) .le. delta .AND. delta .le. deltas(2,12)) then
        nquad = 20
        do l=1,20
        xs(l)   = xs12(l)
        whts(l) = whts12(l)
        end do
        goto 1000
        end if


        if( deltas(1,13) .le. delta .AND. delta .le. deltas(2,13)) then
        nquad = 19
        do l=1,19
        xs(l)   = xs13(l)
        whts(l) = whts13(l)
        end do
        goto 1000
        end if


        if( deltas(1,14) .le. delta .AND. delta .le. deltas(2,14)) then
        nquad = 18
        do l=1,18
        xs(l)   = xs14(l)
        whts(l) = whts14(l)
        end do
        goto 1000
        end if


        if( deltas(1,15) .le. delta .AND. delta .le. deltas(2,15)) then
        nquad = 13
        do l=1,13
        xs(l)   = xs15(l)
        whts(l) = whts15(l)
        end do
        goto 1000
        end if


        if( deltas(1,16) .le. delta .AND. delta .le. deltas(2,16)) then
        nquad = 12
        do l=1,12
        xs(l)   = xs16(l)
        whts(l) = whts16(l)
        end do
        goto 1000
        end if


        if( deltas(1,17) .le. delta .AND. delta .le. deltas(2,17)) then
        nquad = 11
        do l=1,11
        xs(l)   = xs17(l)
        whts(l) = whts17(l)
        end do
        goto 1000
        end if


        if( deltas(1,18) .le. delta .AND. delta .le. deltas(2,18)) then
        nquad = 10
        do l=1,10
        xs(l)   = xs18(l)
        whts(l) = whts18(l)
        end do
        goto 1000
        end if


        if( deltas(1,19) .le. delta .AND. delta .le. deltas(2,19)) then
        nquad = 10
        do l=1,10
        xs(l)   = xs19(l)
        whts(l) = whts19(l)
        end do
        goto 1000
        end if


        if( deltas(1,20) .le. delta .AND. delta .le. deltas(2,20)) then
        nquad = 09
        do l=1,09
        xs(l)   = xs20(l)
        whts(l) = whts20(l)
        end do
        goto 1000
        end if


        if( deltas(1,21) .le. delta .AND. delta .le. deltas(2,21)) then
        nquad = 08
        do l=1,08
        xs(l)   = xs21(l)
        whts(l) = whts21(l)
        end do
        goto 1000
        end if


        if( deltas(1,22) .le. delta .AND. delta .le. deltas(2,22)) then
        nquad = 08
        do l=1,08
        xs(l)   = xs22(l)
        whts(l) = whts22(l)
        end do
        goto 1000
        end if


        if( deltas(1,23) .le. delta .AND. delta .le. deltas(2,23)) then
        nquad = 07
        do l=1,07
        xs(l)   = xs23(l)
        whts(l) = whts23(l)
        end do
        goto 1000
        end if

        ier = 4
        return

 1000 continue

        do i=1,nquad
        xs(i)   = (r1-r0)*xs(i)+r0
        whts(i) = (r1-r0)*whts(i)
        end do

        end


        subroutine pvquad16(ier,r0,r1,nquad,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(*),whts(*)
        dimension deltas(2,23)
        dimension xs01(23),whts01(23)
        dimension xs02(21),whts02(21)
        dimension xs03(21),whts03(21)
        dimension xs04(21),whts04(21)
        dimension xs05(21),whts05(21)
        dimension xs06(21),whts06(21)
        dimension xs07(21),whts07(21)
        dimension xs08(21),whts08(21)
        dimension xs09(21),whts09(21)
        dimension xs10(21),whts10(21)
        dimension xs11(22),whts11(22)
        dimension xs12(21),whts12(21)
        dimension xs13(21),whts13(21)
        dimension xs14(20),whts14(20)
        dimension xs15(15),whts15(15)
        dimension xs16(13),whts16(13)
        dimension xs17(13),whts17(13)
        dimension xs18(11),whts18(11)
        dimension xs19(11),whts19(11)
        dimension xs20(10),whts20(10)
        dimension xs21(10),whts21(10)
        dimension xs22(09),whts22(09)
        dimension xs23(09),whts23(09)
        data deltas /
     -   0.1000000000000000D-14, 0.1000000000000000D-13,
     -   0.1000000000000000D-13, 0.1000000000000000D-12,
     -   0.1000000000000000D-12, 0.1000000000000000D-11,
     -   0.1000000000000000D-11, 0.1000000000000000D-10,
     -   0.1000000000000000D-10, 0.1000000000000000D-09,
     -   0.1000000000000000D-09, 0.1000000000000000D-08,
     -   0.1000000000000000D-08, 0.1000000000000000D-07,
     -   0.1000000000000000D-07, 0.1000000000000000D-06,
     -   0.1000000000000000D-06, 0.1000000000000000D-05,
     -   0.1000000000000000D-05, 0.1000000000000000D-04,
     -   0.1000000000000000D-04, 0.1000000000000000D-03,
     -   0.1000000000000000D-03, 0.1000000000000000D-02,
     -   0.1000000000000000D-02, 0.1000000000000000D-01,
     -   0.1000000000000000D-01, 0.1000000000000000D+00,
     -   0.1000000000000000D+00, 0.2000000000000000D+00,
     -   0.2000000000000000D+00, 0.3000000000000000D+00,
     -   0.3000000000000000D+00, 0.4000000000000000D+00,
     -   0.4000000000000000D+00, 0.5000000000000000D+00,
     -   0.5000000000000000D+00, 0.6000000000000000D+00,
     -   0.6000000000000000D+00, 0.7000000000000000D+00,
     -   0.7000000000000000D+00, 0.8000000000000000D+00,
     -   0.8000000000000000D+00, 0.9000000000000000D+00,
     -   0.9000000000000000D+00, 0.1000000000000000D+01/
        data xs01 /
     -   0.269233817274214178188713945923D-16,
     -   0.148060496890395143194289807974D-15,
     -   0.393259057524734971442178764441D-15,
     -   0.818277748772183409708918583834D-15,
     -   0.152904299561166381872536894905D-14,
     -   0.272935691245894513106324175101D-14,
     -   0.484093420591204780702359383600D-14,
     -   0.885341045638309302978650179539D-14,
     -   0.176024026405808838034337545430D-13,
     -   0.426543534978953580671495317876D-13,
     -   0.251354068743036723654290065390D-12,
     -   0.568290148009936760220319204824D-12,
     -   0.335879075726172089765667991354D-11,
     -   0.842798471245011997468631687799D-02,
     -   0.324245430234960591286741451140D-01,
     -   0.942462753143482028220020509242D-01,
     -   0.203242507086456618259106974252D+00,
     -   0.345786243657723355104468475172D+00,
     -   0.505896730681790373333133635617D+00,
     -   0.666084101000342651876614085701D+00,
     -   0.808941277142593012276479861516D+00,
     -   0.918970091991557874705780235791D+00,
     -   0.984265280476467895059485207071D+00/
        data whts01 /
     -   0.697848588242272959537907964528D-16,
     -   0.176763148757403604664253871112D-15,
     -   0.322611070090104880033871640931D-15,
     -   0.544371676499712845053049051740D-15,
     -   0.909402375583856906666531378714D-15,
     -   0.155675690446701679716242042290D-14,
     -   0.281559243634282203098783517061D-14,
     -   0.561327826547393591115246508470D-14,
     -   0.133264124935126121224280366334D-13,
     -   0.455189883376008177390416531201D-13,
     -   0.147041910670562490475068666651D-11,
     -  -0.643999677354532745632309258458D-11,
     -   0.992768997604316699701640899180D-10,
     -   0.181979342929725894856368894554D-01,
     -   0.363241029096334853066170509039D-01,
     -   0.874074747996957778994638459741D-01,
     -   0.128254810574921181998255934002D+00,
     -   0.154159309755275456862733345290D+00,
     -   0.163117957128309263238441364202D+00,
     -   0.154324431165956772643038991875D+00,
     -   0.128778551706106074474470427518D+00,
     -   0.892709867773738584759839961996D-01,
     -   0.401644407953773635601990346517D-01/
        data xs02 /
     -   0.267715999880700989093785327055D-15,
     -   0.147171698764440986586007304793D-14,
     -   0.390622645416216158294165845642D-14,
     -   0.811869913149315915701787744492D-14,
     -   0.151442910975351071357704101035D-13,
     -   0.269598490786696847334900662990D-13,
     -   0.476072602864426309891777876422D-13,
     -   0.863796571231189550125632152458D-13,
     -   0.168871036502593366823095346484D-12,
     -   0.389583849545032575742362622756D-12,
     -   0.137282866810065024092944273310D-11,
     -   0.103081299006986531789360940701D-09,
     -   0.159198823488037064570056655541D-01,
     -   0.819844482981433057836074837266D-01,
     -   0.193314285373295147933395543729D+00,
     -   0.337873289712816520139998173466D+00,
     -   0.500000001068315898578995092182D+00,
     -   0.662126712423815276048598400751D+00,
     -   0.806685716763336644248315240507D+00,
     -   0.918015553838488472219066016972D+00,
     -   0.984080119787827967291044476033D+00/
        data whts02 /
     -   0.693856604489784684026221678254D-15,
     -   0.175624625714766330437967991787D-14,
     -   0.320096833342287495812863980974D-14,
     -   0.538973165327336513170393726146D-14,
     -   0.897469792554909307181046463185D-14,
     -   0.152862712385021312259780393101D-13,
     -   0.274176319457733655288153283489D-13,
     -   0.538245063620222152696066685655D-13,
     -   0.123520689510010413169064050113D-12,
     -   0.380089645515235284570658817321D-12,
     -   0.233179729399505297693504031171D-11,
     -   0.213367947846668687812982387328D-08,
     -   0.406371940939608153115873981992D-01,
     -   0.903240801544394272555511506605D-01,
     -   0.130305347923053187506998899873D+00,
     -   0.156173538186316073935335260235D+00,
     -   0.165119677147829928979112850405D+00,
     -   0.156173538186316071848308270481D+00,
     -   0.130305347923053180727850548365D+00,
     -   0.903240801544393997583222353866D-01,
     -   0.406371940939604846709060290969D-01/
        data xs03 /
     -   0.267704556226923302133413457412D-14,
     -   0.147165002032164257970134208188D-13,
     -   0.390602806279917410463270060822D-13,
     -   0.811821793268453286615221442924D-13,
     -   0.151431972105844154199481270526D-12,
     -   0.269573640002415991341973748314D-12,
     -   0.476013395849062309761465294008D-12,
     -   0.863640105020050163249635071005D-12,
     -   0.168820893074363176678197710743D-11,
     -   0.389350237646886615941131951346D-11,
     -   0.137017429318015704522398734073D-10,
     -   0.910129711544653774387941855473D-09,
     -   0.159198967553825173238067725721D-01,
     -   0.819844617375663786853914146944D-01,
     -   0.193314297182889470874900400137D+00,
     -   0.337873299406118190160895367942D+00,
     -   0.500000008388139715631533774718D+00,
     -   0.662126717370161182129910045139D+00,
     -   0.806685719593389657660317821918D+00,
     -   0.918015555038711887279946123981D+00,
     -   0.984080120020889406405371709244D+00/
        data whts03 /
     -   0.693826511161392612161027161015D-14,
     -   0.175616053712744612062273181792D-13,
     -   0.320077945770930206337729691269D-13,
     -   0.538932748440759248630231122237D-13,
     -   0.897380915194056604589170510689D-13,
     -   0.152841928746848426571157580866D-12,
     -   0.274122503746328714022495241550D-12,
     -   0.538080972094867548518277786686D-12,
     -   0.123455098014190811071903859933D-11,
     -   0.379659755343643106842915511060D-11,
     -   0.232295910886498014583928282012D-10,
     -   0.167468311718621932655560090060D-07,
     -   0.406371934990663903832139822193D-01,
     -   0.903240788321283543005545492601D-01,
     -   0.130305346015429203172246599450D+00,
     -   0.156173535899990630522089562695D+00,
     -   0.165119674730536073670004447024D+00,
     -   0.156173535899990503559522064778D+00,
     -   0.130305346015428790768371224706D+00,
     -   0.903240788321266815298305397922D-01,
     -   0.406371934990462761737341162003D-01/
        data xs04 /
     -   0.267689493313039651888571306488D-13,
     -   0.147156187401809438310848990209D-12,
     -   0.390576693199680483941726463844D-12,
     -   0.811758457242831151780686958832D-12,
     -   0.151417574669131171298722568770D-11,
     -   0.269540933634932490837399887594D-11,
     -   0.475935478833810596012292656353D-11,
     -   0.863434221836000639008207788406D-11,
     -   0.168754931219092660115143814736D-10,
     -   0.389043156554891687767155483353D-10,
     -   0.136669652113133038137573787859D-09,
     -   0.789456610323760749647861885730D-08,
     -   0.159200056397950903046878204301D-01,
     -   0.819845633125315471133221291271D-01,
     -   0.193314386439666843248852045903D+00,
     -   0.337873372667989031631081255517D+00,
     -   0.500000063711293115019699786297D+00,
     -   0.662126754754593854710460922249D+00,
     -   0.806685740982902222277040205510D+00,
     -   0.918015564109988611206182216867D+00,
     -   0.984080121782365469791364004226D+00/
        data whts04 /
     -   0.693786900344057084133085685020D-13,
     -   0.175604770832371349186629434380D-12,
     -   0.320053085619511717894859374440D-12,
     -   0.538879552805999229987686672813D-12,
     -   0.897263942723806813669695599540D-12,
     -   0.152814577109985385595750412269D-11,
     -   0.274051689395283996277272340441D-11,
     -   0.537865093355929947064315249059D-11,
     -   0.123368842203567090764156946369D-10,
     -   0.379095057348058957041915958600D-10,
     -   0.231140616053980210037676919306D-09,
     -   0.127128284346596606446556638455D-06,
     -   0.406371890038309935506793113340D-01,
     -   0.903240688381951076214837233342D-01,
     -   0.130305331597646217295533578317D+00,
     -   0.156173518619972752013482745385D+00,
     -   0.165119656460655886516219958328D+00,
     -   0.156173518619965553299680382594D+00,
     -   0.130305331597622834215483568753D+00,
     -   0.903240688381002625923273517774D-01,
     -   0.406371890026905466487320306938D-01/
        data xs05 /
     -   0.267668772971650399337765867875D-12,
     -   0.147144062232976376662508528742D-11,
     -   0.390540773401198954436896770867D-11,
     -   0.811671337963413438573199866000D-11,
     -   0.151397771701900305111113930570D-10,
     -   0.269495950522801320487126656621D-10,
     -   0.475828325794037513578745296539D-10,
     -   0.863151139357735976393621677136D-10,
     -   0.168664270403055854022657340870D-09,
     -   0.388621514922644511080650857029D-09,
     -   0.136194244977516758523064404997D-08,
     -   0.668804821714595162610836687896D-07,
     -   0.159207919483119544341584202785D-01,
     -   0.819852968502879039474228791186D-01,
     -   0.193315031021670265695818369608D+00,
     -   0.337873901740686904788698219615D+00,
     -   0.500000463236729466906822092153D+00,
     -   0.662127024732599349863733982785D+00,
     -   0.806685895450902239900456818695D+00,
     -   0.918015629619758906311386369425D+00,
     -   0.984080134503165317126216304430D+00/
        data whts05 /
     -   0.693732412357111512741511442159D-12,
     -   0.175589250572212260162769381110D-11,
     -   0.320018890106507190617655118300D-11,
     -   0.538806384728913406720437647445D-11,
     -   0.897103063017462003856309816900D-11,
     -   0.152776962282732708224985281910D-10,
     -   0.273954318408314393026649506908D-10,
     -   0.537568338589622558478169000520D-10,
     -   0.123250339789188556944476584879D-09,
     -   0.378320418651718298734479762521D-09,
     -   0.229566316357522699173580323554D-08,
     -   0.923494393063641490468244514471D-06,
     -   0.406371565903802979356499457480D-01,
     -   0.903239966694365735340352504995D-01,
     -   0.130305227478206455461046349494D+00,
     -   0.156173393829750879763471419684D+00,
     -   0.165119524521752467113863471548D+00,
     -   0.156173393829379115525336766127D+00,
     -   0.130305227476998884226050753228D+00,
     -   0.903239966645385626728827637484D-01,
     -   0.406371565314909078853639381647D-01/
        data xs06 /
     -   0.267638484770013938722166064497D-11,
     -   0.147126338371864148207070980872D-10,
     -   0.390488269311959389587611195875D-10,
     -   0.811544000598357834723956851209D-10,
     -   0.151368828572792956665559563479D-09,
     -   0.269430211294483673801741105139D-09,
     -   0.475671753532085200720506064081D-09,
     -   0.862737608483359166628079438647D-09,
     -   0.168531904467501788223509025741D-08,
     -   0.388006803186993202254488926755D-08,
     -   0.135505545280441574175539634149D-07,
     -   0.548224896459115781072933489924D-06,
     -   0.159261240555753660322517162556D-01,
     -   0.819902717529181829068838606611D-01,
     -   0.193319402698935244458996569157D+00,
     -   0.337877490032276781129225466505D+00,
     -   0.500003172915926799272608902307D+00,
     -   0.662128855791739215481290417739D+00,
     -   0.806686943092685456276024543815D+00,
     -   0.918016073924074107760217279620D+00,
     -   0.984080220778957417031005477823D+00/
        data whts06 /
     -   0.693652764154215056832588861096D-11,
     -   0.175566564314929066687707791677D-10,
     -   0.319968908130403481563730014447D-10,
     -   0.538699445346187789383495767361D-10,
     -   0.896867949825170917187423461052D-10,
     -   0.152721998759462773349910682580D-09,
     -   0.273812069854004435374261793252D-09,
     -   0.537134985089358782083530377840D-09,
     -   0.123077432699162601522008431655D-08,
     -   0.377192598301241397095822732869D-08,
     -   0.227296011658471888894416179097D-07,
     -   0.631396987452536220453626216552D-05,
     -   0.406369389741291528451478630888D-01,
     -   0.903235073870628975723974607729D-01,
     -   0.130304521359439605038549125260D+00,
     -   0.156172547486591427536812431879D+00,
     -   0.165118629684428093355157969536D+00,
     -   0.156172547469717757805254290495D+00,
     -   0.130304521304631645587337089735D+00,
     -   0.903235071647772191644786304899D-01,
     -   0.406369363033303317753944600030D-01/
        data xs07 /
     -   0.267590157757228787602544326088D-10,
     -   0.147098059273937934792255118279D-09,
     -   0.390404500378959445859210085631D-09,
     -   0.811340849806430769884121791455D-09,
     -   0.151322657776468704742743412884D-08,
     -   0.269325357483847425687882217845D-08,
     -   0.475422079411560428507768093197D-08,
     -   0.862078455230501362410867516116D-08,
     -   0.168321099023537347943964410541D-07,
     -   0.387030029145443363668723949723D-07,
     -   0.134421980178181364382030044528D-06,
     -   0.427988312219873790701797193330D-05,
     -   0.159591379003406587616382441943D-01,
     -   0.820210976142619117248082020990D-01,
     -   0.193346493704155775877244029977D+00,
     -   0.337899727195146928026090363974D+00,
     -   0.500019965435706499476632204528D+00,
     -   0.662140203380775113362640564137D+00,
     -   0.806693435650524925543631843953D+00,
     -   0.918018827420651044814050832747D+00,
     -   0.984080755458352721808702586970D+00/
        data whts07 /
     -   0.693525680305812597873660413509D-10,
     -   0.175530368378818744828696166016D-09,
     -   0.319889167382086111529682562108D-09,
     -   0.538528852466128359531853327482D-09,
     -   0.896492944259014673322272346163D-09,
     -   0.152634350870021971527106706351D-08,
     -   0.273585311691389748628885787448D-08,
     -   0.536444609236909455426430574308D-08,
     -   0.122802330544464809510273080198D-07,
     -   0.375404286001350725118416505021D-07,
     -   0.223749124018926652192750675946D-06,
     -   0.395343320395954334312790000460D-04,
     -   0.406356716988623408607987429102D-01,
     -   0.903204819565439075119416059256D-01,
     -   0.130300147023402119876028662898D+00,
     -   0.156167303038263435422676115741D+00,
     -   0.165113084335668455261179463655D+00,
     -   0.156167302402109141509513730678D+00,
     -   0.130300144957404541541042637907D+00,
     -   0.903204735819920380462274966378D-01,
     -   0.406355714774921325291137042038D-01/
        data xs08 /
     -   0.267502426087763062021743074528D-09,
     -   0.147046723930367728048138275091D-08,
     -   0.390252444009758384867445547448D-08,
     -   0.810972132165022642988120348803D-08,
     -   0.151238871516834539638000300168D-07,
     -   0.269135126246674764430853949106D-07,
     -   0.474969291411216638365367830889D-07,
     -   0.860883942141587557275748588344D-07,
     -   0.167939657244333670667230489631D-06,
     -   0.385269652448668694128654848525D-06,
     -   0.132502625272655927564135771135D-05,
     -   0.309459651594429537888817978477D-04,
     -   0.161377040303195268975729753739D-01,
     -   0.821884742770188328145502780850D-01,
     -   0.193493671023278810400786532268D+00,
     -   0.338020555783911178376189016315D+00,
     -   0.500111216628540712554948580163D+00,
     -   0.662201868957474014832029009482D+00,
     -   0.806728718540911630364362775356D+00,
     -   0.918033791084415174602766807677D+00,
     -   0.984083661148598990245542654888D+00/
        data whts08 /
     -   0.693294977271851406880587996022D-09,
     -   0.175464664124189416021653653009D-08,
     -   0.319744435884567331924404903169D-08,
     -   0.538219275262262160930271499281D-08,
     -   0.895812585888929039477825442921D-08,
     -   0.152475394545766644995432899801D-07,
     -   0.273174320508714367099660483445D-07,
     -   0.535194708764436380719778551220D-07,
     -   0.122305406721419760010558701993D-06,
     -   0.372193361900652757207561672571D-06,
     -   0.217543424152024926027251735249D-05,
     -   0.216508554663722157020833192867D-03,
     -   0.406309735408749800631600259852D-01,
     -   0.903042283822516978913059998210D-01,
     -   0.130276422350194762490463065051D+00,
     -   0.156138819458310884145266819597D+00,
     -   0.165082955378957539307239683125D+00,
     -   0.156138801126883543932496002010D+00,
     -   0.130276362864104125487471072309D+00,
     -   0.903039879360244586961382650313D-01,
     -   0.406281544045771727445905713302D-01/
        data xs09 /
     -   0.267312035193972447956081795407D-08,
     -   0.146935326190140585823740361650D-07,
     -   0.389922524793368086013279522098D-07,
     -   0.810172287119081326490266890947D-07,
     -   0.151057175737168151512666007406D-06,
     -   0.268722805778850216830343783237D-06,
     -   0.473988706651058616873845449384D-06,
     -   0.858300956502548967761380589519D-06,
     -   0.167117447323496552744363815974D-05,
     -   0.381506635727743336721619116662D-05,
     -   0.128542368405223929112541649038D-04,
     -   0.197981979697210415834798409364D-03,
     -   0.169348973261434806309467241598D-01,
     -   0.829470475734262292366272452002D-01,
     -   0.194162189043735260585112638533D+00,
     -   0.338569781447305329891262560346D+00,
     -   0.500526126837434759000249601134D+00,
     -   0.662482300652349850354346857241D+00,
     -   0.806889185814318141960046314899D+00,
     -   0.918101849429972898356070783624D+00,
     -   0.984096877240283075596028868311D+00/
        data whts09 /
     -   0.692794324777313107409203634817D-08,
     -   0.175322096766615769908769891126D-07,
     -   0.319430464511411983527002628738D-07,
     -   0.537547927772757883306012683567D-07,
     -   0.894337907452926198634475569518D-07,
     -   0.152131123118707604764353802888D-06,
     -   0.272285321898909606075765905220D-06,
     -   0.532497343748920234578227362342D-06,
     -   0.121238150729568087885329433334D-05,
     -   0.365383187800975816642778180772D-05,
     -   0.205058236589306632889114688913D-04,
     -   0.967896765581231202972197181055D-03,
     -   0.406458900310317334295460099178D-01,
     -   0.902337777966872332108097050449D-01,
     -   0.130169405758500604180230181770D+00,
     -   0.156009589609409366075148393057D+00,
     -   0.164946050477644194706792322849D+00,
     -   0.156009220419723289413704490501D+00,
     -   0.130168211836716776508232198176D+00,
     -   0.902290096099603959037022715605D-01,
     -   0.405944191521292745847680073327D-01/
        data xs10 /
     -   0.266810163668385772642786197223D-07,
     -   0.146641729151392098100696001444D-06,
     -   0.389053271652060177957130680403D-06,
     -   0.808065972040735394206405322696D-06,
     -   0.150579076491505834280541551545D-05,
     -   0.267639221573887372186246280319D-05,
     -   0.471417125747031608597179362792D-05,
     -   0.851553048018355139697011015300D-05,
     -   0.164986623149657655629175427269D-04,
     -   0.371955729480519783128355969583D-04,
     -   0.119298156468135784748999936371D-03,
     -   0.106378759715357008820672650813D-02,
     -   0.197746514761457938915884776813D-01,
     -   0.857621310365499013574343809703D-01,
     -   0.196661037953700902309041488661D+00,
     -   0.340627607733685871016433236562D+00,
     -   0.502082328756126225606418467547D+00,
     -   0.663534679702053254208427906865D+00,
     -   0.807491552531444160645057799344D+00,
     -   0.918357371198875407251172803928D+00,
     -   0.984146500425124074836118983975D+00/
        data whts10 /
     -   0.691474650129262424709930035669D-07,
     -   0.174946418676252357504008287010D-06,
     -   0.318603579206235026895769743723D-06,
     -   0.535781327306848201381059832027D-06,
     -   0.890462259999248198989969330175D-06,
     -   0.151228093959531319123248920125D-05,
     -   0.269960998299849473212907961686D-05,
     -   0.525486287215876344009829426566D-05,
     -   0.118497637017496400387912328870D-04,
     -   0.348434542143528416217958407705D-04,
     -   0.177603532024881202765772241334D-03,
     -   0.331345048436682546388560076111D-02,
     -   0.409990848765333292376246216106D-01,
     -   0.900099339224654265706841822832D-01,
     -   0.129778636523514653993856635748D+00,
     -   0.155528446043871106746324667893D+00,
     -   0.164433686526176631151631365744D+00,
     -   0.155523353660370350322199462876D+00,
     -   0.129762364327526779242460694224D+00,
     -   0.899475373418294791645399700884D-01,
     -   0.404677538485594803409757377938D-01/
        data xs11 /
     -   0.806732320908936258995754118078D-07,
     -   0.785671884747863445565810027137D-06,
     -   0.252805533303108110338887414481D-05,
     -   0.568847238036037415605155789538D-05,
     -   0.109740342006558507146125285250D-04,
     -   0.196830819340660439235623040540D-04,
     -   0.342680747064751901011355837861D-04,
     -   0.596981642120427024258553906465D-04,
     -   0.107191662078099510095575076708D-03,
     -   0.206276355687652643783499400928D-03,
     -   0.453819787222894591330973823759D-03,
     -   0.129148882423476445939517047257D-02,
     -   0.574712766944209152782332275324D-02,
     -   0.302220763858912912471989793268D-01,
     -   0.969498758363987226435518516271D-01,
     -   0.206810294597232648914413804051D+00,
     -   0.349052918681322668088722019510D+00,
     -   0.508477385195084355760612766007D+00,
     -   0.667867697619855628754691937090D+00,
     -   0.809974401216215769015903954281D+00,
     -   0.919411226429360166467601940610D+00,
     -   0.984351223383972708231764340618D+00/
        data whts11 /
     -   0.282204510732279507527071980654D-06,
     -   0.117779313884514802933333631860D-05,
     -   0.236713987770083660556766876084D-05,
     -   0.406764088360005200895845229833D-05,
     -   0.671144692117428812887695004588D-05,
     -   0.110955665372551351618007429943D-04,
     -   0.188505000339979317467170383488D-04,
     -   0.337211992090200648590903363949D-04,
     -   0.656055764845387776379342948721D-04,
     -   0.145921013304166067426882630042D-03,
     -   0.403508850445720505917456294543D-03,
     -   0.160157789449565977633516489029D-02,
     -   0.987916077749004249500995156626D-02,
     -   0.436773194001990707186882228702D-01,
     -   0.895343694297171654158872632062D-01,
     -   0.128314671591402365004073285945D+00,
     -   0.153601998132959853245475455337D+00,
     -   0.162344701536245260118221929482D+00,
     -   0.153529020636597272295255343703D+00,
     -   0.128091535607020987665722805379D+00,
     -   0.887871114231159909445493998144D-01,
     -   0.399452246394095812337510329905D-01/
        data xs12 /
     -   0.845690376195442453732293400098D-06,
     -   0.842330098649461571873564404320D-05,
     -   0.272912679162794524092454988766D-04,
     -   0.618604638092027420800464327202D-04,
     -   0.120493754932278438606893816190D-03,
     -   0.218933320976879797825842557206D-03,
     -   0.388004251713444208795493114786D-03,
     -   0.693254738075982803981638155442D-03,
     -   0.129173763659163058618595476152D-02,
     -   0.262307374544854747892620449458D-02,
     -   0.615251779461003276476065564735D-02,
     -   0.171680612084187843675554317872D-01,
     -   0.495914110408000878322159496325D-01,
     -   0.118845666674467155502552033448D+00,
     -   0.227347065107696205324512084388D+00,
     -   0.366367842168136726336845124708D+00,
     -   0.521723605551230091499372441436D+00,
     -   0.676881679697770739800519105285D+00,
     -   0.815152354846237897317752752292D+00,
     -   0.921612143868334376255838255330D+00,
     -   0.984779074020081675059043221547D+00/
        data whts12 /
     -   0.300092997656108003966718080127D-05,
     -   0.127003854435306322150549739616D-04,
     -   0.257331062907764558267362923883D-04,
     -   0.447486678531752214817342680305D-04,
     -   0.750276711893977465871520949151D-04,
     -   0.126707288197244728994251328969D-03,
     -   0.221578576276716888087041795615D-03,
     -   0.412606378321420256584305803973D-03,
     -   0.848323736198836854837956724023D-03,
     -   0.201994878780193128873929428725D-02,
     -   0.580154441831909622529373300260D-02,
     -   0.186425255962154765262634861703D-01,
     -   0.491473619236842202523223937741D-01,
     -   0.896487732823524428716411897990D-01,
     -   0.125803587544555475442011168922D+00,
     -   0.149827533788309268255749043168D+00,
     -   0.158093019557586864985095251848D+00,
     -   0.149409323943284657553682636180D+00,
     -   0.124616724768915701037251278081D+00,
     -   0.863658672288233696438369507167D-01,
     -   0.388533624204038360534596735879D-01/
        data xs13 /
     -   0.261764463579121931270153500221D-04,
     -   0.143516810568560037315826092148D-03,
     -   0.378940786328330781539611205489D-03,
     -   0.780850388974619280354364799493D-03,
     -   0.143699464322599145708845573534D-02,
     -   0.250393006348323054916953414703D-02,
     -   0.426848823243610240502589018783D-02,
     -   0.728066332162544336967035401174D-02,
     -   0.126465709344630512306069123839D-01,
     -   0.226468456985292436254879158055D-01,
     -   0.417875393037321430482743826586D-01,
     -   0.775423600319177216393584448120D-01,
     -   0.138470305937790040700597998098D+00,
     -   0.229325194490484838995942873453D+00,
     -   0.347756692282941144260843826626D+00,
     -   0.484968918780701437130404031169D+00,
     -   0.628126220533099301795532988192D+00,
     -   0.762775834741686248449411216519D+00,
     -   0.874987632685706751276640148796D+00,
     -   0.953403785797307828846136407882D+00,
     -   0.992497190558614214274943499888D+00/
        data whts13 /
     -   0.678024460066445673433066422711D-04,
     -   0.170711590144012621824013655586D-03,
     -   0.307977717187655295595382941543D-03,
     -   0.510011297424787508139031141371D-03,
     -   0.827428957728678489982395438825D-03,
     -   0.135222599594049081679699150318D-02,
     -   0.226408006659710094881078574759D-02,
     -   0.393484056973087882318130437366D-02,
     -   0.715853841186480094704190093201D-02,
     -   0.135738404720997164744802843220D-01,
     -   0.259801848968287996133383652640D-01,
     -   0.470540536362876958477376490561D-01,
     -   0.756554154477088643995328630302D-01,
     -   0.105642957787768985885612675027D+00,
     -   0.129710143212170857132251091781D+00,
     -   0.142534454716775123278320836152D+00,
     -   0.141336033862334685277810070746D+00,
     -   0.125608970296720079667476347521D+00,
     -   0.968793564153900492303868241302D-01,
     -   0.588917009670732177464022102974D-01,
     -   0.205392712362168754279356702962D-01/
        data xs14 /
     -   0.179286353841196370102218290510D-03,
     -   0.101551488420419225000924838754D-02,
     -   0.275569512672506537153256329154D-02,
     -   0.573875573867052296700201526763D-02,
     -   0.104773176768694621771162494208D-01,
     -   0.177566668997193917331231554487D-01,
     -   0.287561034493108691487153578856D-01,
     -   0.451944997439334445428651191266D-01,
     -   0.694654040657942901382920450817D-01,
     -   0.104651812532939084986733066377D+00,
     -   0.154222495585414360473695403489D+00,
     -   0.221228751238752013606260075546D+00,
     -   0.307085733484447375001890100483D+00,
     -   0.410413721949366653366153450250D+00,
     -   0.526510866333922507871272321221D+00,
     -   0.647687046071216520303235636045D+00,
     -   0.764279774831454855199003898558D+00,
     -   0.866013512568645802439411878108D+00,
     -   0.943407009383635747890306866147D+00,
     -   0.989033809402553422227932005979D+00/
        data whts14 /
     -   0.469107068127068630239063253147D-03,
     -   0.124190416236177184977573960750D-02,
     -   0.229366277279072422200575945065D-02,
     -   0.375625486853859135180767923614D-02,
     -   0.584974714008354051836907999553D-02,
     -   0.890300877268229547469429950081D-02,
     -   0.133798432456044915822360600049D-01,
     -   0.198951215544128674022973700939D-01,
     -   0.291650297868582617249623336651D-01,
     -   0.417994160536781331251693294103D-01,
     -   0.578720862740541236866552902213D-01,
     -   0.763971880483500673581172650585D-01,
     -   0.950961704710795557617286831649D-01,
     -   0.110770037891183607700498477821D+00,
     -   0.120124741904325455996933856507D+00,
     -   0.120593835045214976309412321134D+00,
     -   0.110855105824201994687185966725D+00,
     -   0.910109804119974825772315605452D-01,
     -   0.625208379927925659594131753711D-01,
     -   0.280059207116624240812666892341D-01/
        data xs15 /
     -   0.135477504683912779712588307195D-02,
     -   0.843741080105820133418496432889D-02,
     -   0.238956638708842374036648592469D-01,
     -   0.498646612309145782177368103718D-01,
     -   0.887951023199147398735148121574D-01,
     -   0.143319967502468335861019001569D+00,
     -   0.215508853944249835393186565681D+00,
     -   0.305848193928133789942800972629D+00,
     -   0.412331671477780509120154896943D+00,
     -   0.530057430381574240272632138107D+00,
     -   0.651511336858890588401277796213D+00,
     -   0.767438096430255314774747018833D+00,
     -   0.868052656414148244443572881751D+00,
     -   0.944335820683085236020590012135D+00,
     -   0.989220793588268854629720746738D+00/
        data whts15 /
     -   0.368274452085692680287103259751D-02,
     -   0.109139252340795741767518297245D-01,
     -   0.203416361524414853764441600707D-01,
     -   0.320086191994286660424669279032D-01,
     -   0.463040915389965016147853120926D-01,
     -   0.631085362350147758533721757293D-01,
     -   0.813665352004145645211117425987D-01,
     -   0.989849205414053560579773838685D-01,
     -   0.113160034648364628920226188498D+00,
     -   0.121023957648334974373305388402D+00,
     -   0.120320484891588397441260556344D+00,
     -   0.109880119767530346568372265500D+00,
     -   0.898258245350209043983913258694D-01,
     -   0.615460234159069730036479302200D-01,
     -   0.275325464706159248490157805823D-01/
        data xs16 /
     -   0.520507685623446100623838017673D-02,
     -   0.277726750614504115895206885790D-01,
     -   0.696054541222504588684725488098D-01,
     -   0.132101163838012738434662629846D+00,
     -   0.215811100030263160131539278229D+00,
     -   0.319316000978549655197787517906D+00,
     -   0.438457626994173414834529207181D+00,
     -   0.566238056472764989686424215045D+00,
     -   0.693416011635685154202268997868D+00,
     -   0.809604204146419554451735938633D+00,
     -   0.904595042649759109126779082661D+00,
     -   0.969665420603268184972930300397D+00,
     -   0.998600799181947759767372263985D+00/
        data whts16 /
     -   0.133992808857727073383698625770D-01,
     -   0.319536169935204762912997614360D-01,
     -   0.519627479257941546010892813434D-01,
     -   0.731404964231806057346112655830D-01,
     -   0.940564464006603447982596631042D-01,
     -   0.112262029257513796881561388553D+00,
     -   0.124844425590123596400025652674D+00,
     -   0.129157241634065501797432398743D+00,
     -   0.123444761105411460160250522597D+00,
     -   0.107208131748643552615437128110D+00,
     -   0.813013365185693123801105374763D-01,
     -   0.477999144339560352579562151278D-01,
     -   0.946957108278845574359632267551D-02/
        data xs17 /
     -   0.502265403227479240802913218476D-02,
     -   0.271651793391456352523049751517D-01,
     -   0.689653431130008847542994082563D-01,
     -   0.131761562972504830907699151408D+00,
     -   0.215383943579257301014985478946D+00,
     -   0.317642299525798332953969003725D+00,
     -   0.434012070649203094158200248496D+00,
     -   0.557702331514406508980919712534D+00,
     -   0.680218925571818383882073174592D+00,
     -   0.792362205092731933835480542310D+00,
     -   0.885509616044683020047413504703D+00,
     -   0.952982563638380772282805837663D+00,
     -   0.991069196073705655140507880765D+00/
        data whts17 /
     -   0.129815744763462106481037470357D-01,
     -   0.316558626685117203300243633829D-01,
     -   0.521809625416362406836746312087D-01,
     -   0.733913275045159615659040647880D-01,
     -   0.934943145459710496374369198336D-01,
     -   0.110267448554157901749677773657D+00,
     -   0.121330696088063727740004078571D+00,
     -   0.124617120122408444771249819639D+00,
     -   0.118860998741436475443505016149D+00,
     -   0.103964866211975181356003147337D+00,
     -   0.811974386122744918457551805879D-01,
     -   0.531309459124967677126905236333D-01,
     -   0.229264440202058265159707341766D-01/
        data xs18 /
     -   0.911194497086687847920336695178D-02,
     -   0.477842982729496130971544900734D-01,
     -   0.116238343500286657362585668732D+00,
     -   0.211896322862001885623586545546D+00,
     -   0.329947794006667213992946955432D+00,
     -   0.463012794648014852070821226656D+00,
     -   0.601389837112043529189483615347D+00,
     -   0.733879991654847450413085216556D+00,
     -   0.848993954565832866120242996671D+00,
     -   0.936294925948209347406729746918D+00,
     -   0.987663908812533074495822148771D+00/
        data whts18 /
     -   0.233603993974726134633648696714D-01,
     -   0.538377709849750025704986951931D-01,
     -   0.826523460766963622475392062841D-01,
     -   0.107868213789263558404836688501D+00,
     -   0.127000517968195137492287905510D+00,
     -   0.137498319323060077156255996600D+00,
     -   0.137368540618493033428123162565D+00,
     -   0.125673963758485949712904780746D+00,
     -   0.102793150889826117562836212884D+00,
     -   0.704374104929226832919436954492D-01,
     -   0.315093667006094646694087865954D-01/
        data xs19 /
     -   0.704203313586598907412194545210D-02,
     -   0.398530846626565174959555510779D-01,
     -   0.103560047761084424030081677573D+00,
     -   0.197245458135684267675512412961D+00,
     -   0.315735628322886642168729425190D+00,
     -   0.450827388846398702845160975228D+00,
     -   0.592036882878433119547457952453D+00,
     -   0.727538783936462060831178315301D+00,
     -   0.845373169903625212667405464710D+00,
     -   0.934764355913078164681843670281D+00,
     -   0.987367404970115331604543491160D+00/
        data whts19 /
     -   0.185039891246347776246901358493D-01,
     -   0.479561051458629421646243631282D-01,
     -   0.792633213372176764678535093925D-01,
     -   0.107227874975310755064368921584D+00,
     -   0.128374322937569140974707106983D+00,
     -   0.140046239648492100935693016531D+00,
     -   0.140381311434726772905800737246D+00,
     -   0.128606706343880066534599755086D+00,
     -   0.105245095373798344527017742676D+00,
     -   0.721283532823875835857878284008D-01,
     -   0.322666803961198392148568831238D-01/
        data xs20 /
     -   0.123666998654864799250685201887D-01,
     -   0.641828082759067665249694019060D-01,
     -   0.153397001239567347665132346632D+00,
     -   0.273143528931098041397545418078D+00,
     -   0.413653797595948280844646729042D+00,
     -   0.562842348540594558050316524550D+00,
     -   0.707300529135197411439544726035D+00,
     -   0.833601267464778710082720353312D+00,
     -   0.929721170464558602340513360027D+00,
     -   0.986383000877785863975217815636D+00/
        data whts20 /
     -   0.316258206820542880220571045165D-01,
     -   0.714004130838845693296178435709D-01,
     -   0.105879435305947140194336082652D+00,
     -   0.131976017541708318286571756876D+00,
     -   0.147013430880516252848640952416D+00,
     -   0.149106656155861373877265509531D+00,
     -   0.137547647433894689453831137759D+00,
     -   0.113027294191361836196508231626D+00,
     -   0.776468044957233410015453808695D-01,
     -   0.347764802290481907918731309148D-01/
        data xs21 /
     -   0.136953108398016769853407259633D-01,
     -   0.707838189915762047060774719846D-01,
     -   0.167970586390630161185388726976D+00,
     -   0.296219486642809682381277178234D+00,
     -   0.443393180518695226341973274061D+00,
     -   0.595379118145679379068473972251D+00,
     -   0.737543420886840668511128579555D+00,
     -   0.856467502649638415438857917418D+00,
     -   0.941962815483254569911609233938D+00,
     -   0.989164052747515203536506189375D+00/
        data whts21 /
     -   0.349887606305164027892476247364D-01,
     -   0.783394375694550012316438950859D-01,
     -   0.114514217949533174939976011596D+00,
     -   0.139942016274001283277517513113D+00,
     -   0.152028872216320410541563790649D+00,
     -   0.149475062833532040049555115155D+00,
     -   0.132583671856556036969226381196D+00,
     -   0.103527953889707190709227902117D+00,
     -   0.666356361526470099002332070954D-01,
     -   0.279643706277314495918085592548D-01/
        data xs22 /
     -   0.159198802461869550822118985482D-01,
     -   0.819844463366821028502851059651D-01,
     -   0.193314283649704801345648980329D+00,
     -   0.337873288298095535480730992678D+00,
     -   0.500000000000000000000000000000D+00,
     -   0.662126711701904464519269007322D+00,
     -   0.806685716350295198654351019671D+00,
     -   0.918015553663317897149714894035D+00,
     -   0.984080119753813044917788101452D+00/
        data whts22 /
     -   0.406371941807872059859477384733D-01,
     -   0.903240803474287020292397040014D-01,
     -   0.130305348201467731159376755723D+00,
     -   0.156173538520001420034321580632D+00,
     -   0.165119677500629881582269277298D+00,
     -   0.156173538520001420034321580632D+00,
     -   0.130305348201467731159376755723D+00,
     -   0.903240803474287020292397040014D-01,
     -   0.406371941807872059859477384733D-01/
        data xs23 /
     -   0.159198802461869550822118985482D-01,
     -   0.819844463366821028502851059651D-01,
     -   0.193314283649704801345648980329D+00,
     -   0.337873288298095535480730992678D+00,
     -   0.500000000000000000000000000000D+00,
     -   0.662126711701904464519269007322D+00,
     -   0.806685716350295198654351019671D+00,
     -   0.918015553663317897149714894035D+00,
     -   0.984080119753813044917788101452D+00/
        data whts23 /
     -   0.406371941807872059859460790575D-01,
     -   0.903240803474287020292360156265D-01,
     -   0.130305348201467731159371434717D+00,
     -   0.156173538520001420034315203301D+00,
     -   0.165119677500629881582262534653D+00,
     -   0.156173538520001420034315203301D+00,
     -   0.130305348201467731159371434717D+00,
     -   0.903240803474287020292360156265D-01,
     -   0.406371941807872059859460790575D-01/

        ier = 0

        delta = r0/r1

        if( deltas(1,01) .le. delta .AND. delta .le. deltas(2,01)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs01(l)
        whts(l) = whts01(l)
        end do
        goto 1000
        end if


        if( deltas(1,02) .le. delta .AND. delta .le. deltas(2,02)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs02(l)
        whts(l) = whts02(l)
        end do
        goto 1000
        end if


        if( deltas(1,03) .le. delta .AND. delta .le. deltas(2,03)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs03(l)
        whts(l) = whts03(l)
        end do
        goto 1000
        end if


        if( deltas(1,04) .le. delta .AND. delta .le. deltas(2,04)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs04(l)
        whts(l) = whts04(l)
        end do
        goto 1000
        end if


        if( deltas(1,05) .le. delta .AND. delta .le. deltas(2,05)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs05(l)
        whts(l) = whts05(l)
        end do
        goto 1000
        end if


        if( deltas(1,06) .le. delta .AND. delta .le. deltas(2,06)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs06(l)
        whts(l) = whts06(l)
        end do
        goto 1000
        end if


        if( deltas(1,07) .le. delta .AND. delta .le. deltas(2,07)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs07(l)
        whts(l) = whts07(l)
        end do
        goto 1000
        end if


        if( deltas(1,08) .le. delta .AND. delta .le. deltas(2,08)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs08(l)
        whts(l) = whts08(l)
        end do
        goto 1000
        end if


        if( deltas(1,09) .le. delta .AND. delta .le. deltas(2,09)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs09(l)
        whts(l) = whts09(l)
        end do
        goto 1000
        end if


        if( deltas(1,10) .le. delta .AND. delta .le. deltas(2,10)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs10(l)
        whts(l) = whts10(l)
        end do
        goto 1000
        end if


        if( deltas(1,11) .le. delta .AND. delta .le. deltas(2,11)) then
        nquad = 22
        do l=1,22
        xs(l)   = xs11(l)
        whts(l) = whts11(l)
        end do
        goto 1000
        end if


        if( deltas(1,12) .le. delta .AND. delta .le. deltas(2,12)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs12(l)
        whts(l) = whts12(l)
        end do
        goto 1000
        end if


        if( deltas(1,13) .le. delta .AND. delta .le. deltas(2,13)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs13(l)
        whts(l) = whts13(l)
        end do
        goto 1000
        end if


        if( deltas(1,14) .le. delta .AND. delta .le. deltas(2,14)) then
        nquad = 20
        do l=1,20
        xs(l)   = xs14(l)
        whts(l) = whts14(l)
        end do
        goto 1000
        end if


        if( deltas(1,15) .le. delta .AND. delta .le. deltas(2,15)) then
        nquad = 15
        do l=1,15
        xs(l)   = xs15(l)
        whts(l) = whts15(l)
        end do
        goto 1000
        end if


        if( deltas(1,16) .le. delta .AND. delta .le. deltas(2,16)) then
        nquad = 13
        do l=1,13
        xs(l)   = xs16(l)
        whts(l) = whts16(l)
        end do
        goto 1000
        end if


        if( deltas(1,17) .le. delta .AND. delta .le. deltas(2,17)) then
        nquad = 13
        do l=1,13
        xs(l)   = xs17(l)
        whts(l) = whts17(l)
        end do
        goto 1000
        end if


        if( deltas(1,18) .le. delta .AND. delta .le. deltas(2,18)) then
        nquad = 11
        do l=1,11
        xs(l)   = xs18(l)
        whts(l) = whts18(l)
        end do
        goto 1000
        end if


        if( deltas(1,19) .le. delta .AND. delta .le. deltas(2,19)) then
        nquad = 11
        do l=1,11
        xs(l)   = xs19(l)
        whts(l) = whts19(l)
        end do
        goto 1000
        end if


        if( deltas(1,20) .le. delta .AND. delta .le. deltas(2,20)) then
        nquad = 10
        do l=1,10
        xs(l)   = xs20(l)
        whts(l) = whts20(l)
        end do
        goto 1000
        end if


        if( deltas(1,21) .le. delta .AND. delta .le. deltas(2,21)) then
        nquad = 10
        do l=1,10
        xs(l)   = xs21(l)
        whts(l) = whts21(l)
        end do
        goto 1000
        end if


        if( deltas(1,22) .le. delta .AND. delta .le. deltas(2,22)) then
        nquad = 09
        do l=1,09
        xs(l)   = xs22(l)
        whts(l) = whts22(l)
        end do
        goto 1000
        end if


        if( deltas(1,23) .le. delta .AND. delta .le. deltas(2,23)) then
        nquad = 09
        do l=1,09
        xs(l)   = xs23(l)
        whts(l) = whts23(l)
        end do
        goto 1000
        end if

        ier = 4
        return

 1000 continue

        do i=1,nquad
        xs(i)   = (r1-r0)*xs(i)+r0
        whts(i) = (r1-r0)*whts(i)
        end do

        end



        subroutine pvquad20(ier,r0,r1,nquad,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(*),whts(*)
        dimension deltas(2,23)
        dimension xs01(24),whts01(24)
        dimension xs02(23),whts02(23)
        dimension xs03(23),whts03(23)
        dimension xs04(23),whts04(23)
        dimension xs05(23),whts05(23)
        dimension xs06(23),whts06(23)
        dimension xs07(23),whts07(23)
        dimension xs08(23),whts08(23)
        dimension xs09(23),whts09(23)
        dimension xs10(23),whts10(23)
        dimension xs11(25),whts11(25)
        dimension xs12(23),whts12(23)
        dimension xs13(23),whts13(23)
        dimension xs14(21),whts14(21)
        dimension xs15(16),whts15(16)
        dimension xs16(14),whts16(14)
        dimension xs17(13),whts17(13)
        dimension xs18(13),whts18(13)
        dimension xs19(12),whts19(12)
        dimension xs20(11),whts20(11)
        dimension xs21(11),whts21(11)
        dimension xs22(11),whts22(11)
        dimension xs23(11),whts23(11)
        data deltas /
     -   0.1000000000000000D-14, 0.1000000000000000D-13,
     -   0.1000000000000000D-13, 0.1000000000000000D-12,
     -   0.1000000000000000D-12, 0.1000000000000000D-11,
     -   0.1000000000000000D-11, 0.1000000000000000D-10,
     -   0.1000000000000000D-10, 0.1000000000000000D-09,
     -   0.1000000000000000D-09, 0.1000000000000000D-08,
     -   0.1000000000000000D-08, 0.1000000000000000D-07,
     -   0.1000000000000000D-07, 0.1000000000000000D-06,
     -   0.1000000000000000D-06, 0.1000000000000000D-05,
     -   0.1000000000000000D-05, 0.1000000000000000D-04,
     -   0.1000000000000000D-04, 0.1000000000000000D-03,
     -   0.1000000000000000D-03, 0.1000000000000000D-02,
     -   0.1000000000000000D-02, 0.1000000000000000D-01,
     -   0.1000000000000000D-01, 0.1000000000000000D+00,
     -   0.1000000000000000D+00, 0.2000000000000000D+00,
     -   0.2000000000000000D+00, 0.3000000000000000D+00,
     -   0.3000000000000000D+00, 0.4000000000000000D+00,
     -   0.4000000000000000D+00, 0.5000000000000000D+00,
     -   0.5000000000000000D+00, 0.6000000000000000D+00,
     -   0.6000000000000000D+00, 0.7000000000000000D+00,
     -   0.7000000000000000D+00, 0.8000000000000000D+00,
     -   0.8000000000000000D+00, 0.9000000000000000D+00,
     -   0.9000000000000000D+00, 0.1000000000000000D+01/
        data xs01 /
     -   0.268963679618091426033377149495D-16,
     -   0.147902294238008495686601903427D-15,
     -   0.392789711994945650672809892587D-15,
     -   0.817136809515630553198823800402D-15,
     -   0.152644067750194625371083811448D-14,
     -   0.272341475724709028133735222357D-14,
     -   0.482665887939379219037633571007D-14,
     -   0.881511361520763086746393469840D-14,
     -   0.174756712939488851972037322501D-13,
     -   0.420070214962656824769538452956D-13,
     -   0.249627288750256483217656425718D-12,
     -   0.334963745572518063420510996011D-12,
     -   0.547786400602319400297613179722D-11,
     -   0.108856710607008486844973507349D-01,
     -   0.564687002435188204591262536048D-01,
     -   0.134923997329934565525832858566D+00,
     -   0.240451935499285826256985068980D+00,
     -   0.365228422109649329752663507671D+00,
     -   0.500000000067600550260918773685D+00,
     -   0.634771578025551770766543437795D+00,
     -   0.759548064635915274252467575826D+00,
     -   0.865076002805266534957663669660D+00,
     -   0.943531299891682279938738900144D+00,
     -   0.989114329074500251091304490450D+00/
        data whts01 /
     -   0.697138086375967402815827416191D-16,
     -   0.176560474325354006480065920274D-15,
     -   0.322163394236949750533636923625D-15,
     -   0.543410283322994785506746145748D-15,
     -   0.907277429468799742561899290758D-15,
     -   0.155174946385571323355106679342D-14,
     -   0.280246299912321292148400215363D-14,
     -   0.557232059309049119847968085127D-14,
     -   0.131544649519737684105467164194D-13,
     -   0.442189274081422250131871307244D-13,
     -   0.242791045798278265263994168581D-11,
     -  -0.348109775809756764145416628545D-11,
     -   0.136184966571859540083499489797D-09,
     -   0.278342835543236094433617710987D-01,
     -   0.627901847239630104037898284932D-01,
     -   0.931451054512738049874975034003D-01,
     -   0.116596882280231213170398864618D+00,
     -   0.131402272237357599276444625050D+00,
     -   0.136462543370500429313337669039D+00,
     -   0.131402272237357599270913064679D+00,
     -   0.116596882280231213155167170732D+00,
     -   0.931451054512738049455923972759D-01,
     -   0.627901847239630102390880173211D-01,
     -   0.278342835543236074718581560924D-01/
        data xs02 /
     -   0.267714296318294147460617448776D-15,
     -   0.147170701851079887178331309859D-14,
     -   0.390619692035903770167862707565D-14,
     -   0.811862749662217261364723576738D-14,
     -   0.151441282514454871603043611361D-13,
     -   0.269594791206003189829234527144D-13,
     -   0.476063788364567071546869308746D-13,
     -   0.863773276023145756994659699575D-13,
     -   0.168863570218950292487815545405D-12,
     -   0.389549055672443746076492200888D-12,
     -   0.137243285519011406965253208359D-11,
     -   0.101080059303423979446094373375D-09,
     -   0.108856729612571398693780241813D-01,
     -   0.564687020564886932387964045507D-01,
     -   0.134923998992154723851284388106D+00,
     -   0.240451936958736908233530854074D+00,
     -   0.365228423329345766040558861011D+00,
     -   0.500000001028337041540374643053D+00,
     -   0.634771578727328316442457220792D+00,
     -   0.759548065097937172033083970115D+00,
     -   0.865076003064519350517504354999D+00,
     -   0.943531300000185361672577133036D+00,
     -   0.989114329095416773694951282093D+00/
        data whts02 /
     -   0.693852124636199958715637659703D-15,
     -   0.175623349634453602185535507504D-14,
     -   0.320094021599756744297859699025D-14,
     -   0.538967148499266591798483601094D-14,
     -   0.897456561226562494915304897013D-14,
     -   0.152859618208219024972278097727D-13,
     -   0.274168307276474050991406025288D-13,
     -   0.538220631554233552247638786941D-13,
     -   0.123510921896324520971200194309D-12,
     -   0.380025602047433547076423961800D-12,
     -   0.233047824727150062920678553192D-11,
     -   0.205372302846263306928714243546D-08,
     -   0.278342835008412315451576068150D-01,
     -   0.627901846033134040252943655671D-01,
     -   0.931451052722980107435724643734D-01,
     -   0.116596882056193457374764593526D+00,
     -   0.131402271984871684662961100149D+00,
     -   0.136462543108291339578471963615D+00,
     -   0.131402271984871683406063468742D+00,
     -   0.116596882056193453913773737705D+00,
     -   0.931451052722980012217695702849D-01,
     -   0.627901846033133666012606491838D-01,
     -   0.278342835008407835744240227963D-01/
        data xs03 /
     -   0.267702349551471520115692648622D-14,
     -   0.147163710708564473596764009448D-13,
     -   0.390598980747240677106068026500D-13,
     -   0.811812514527927094866956792165D-13,
     -   0.151429862846539642023506684799D-12,
     -   0.269568848330412917084690552852D-12,
     -   0.476001980133079116116326549970D-12,
     -   0.863609938839647088428254131361D-12,
     -   0.168811226972585713602320339641D-11,
     -   0.389305221588815695455943267902D-11,
     -   0.136966366523763807450190891920D-10,
     -   0.890118657554898937143413390563D-09,
     -   0.108856868209320979114514895250D-01,
     -   0.564687152774527006063464716362D-01,
     -   0.134924011113787198554223313796D+00,
     -   0.240451947601689316661103411185D+00,
     -   0.365228432223903051897856281033D+00,
     -   0.500000008034445732729094767968D+00,
     -   0.634771583844988377575448167789D+00,
     -   0.759548068467201979379857113799D+00,
     -   0.865076004955103742424360112874D+00,
     -   0.943531300791437068992412637857D+00,
     -   0.989114329247949162292384563897D+00/
        data whts03 /
     -   0.693820708281040396664317165552D-14,
     -   0.175614400789672006379935343511D-13,
     -   0.320074303759388672448631542144D-13,
     -   0.538924955155122301770860638865D-13,
     -   0.897363778056395772657498763171D-13,
     -   0.152837921437223407964818595176D-12,
     -   0.274112128105022085209893149404D-12,
     -   0.538049338645004767145490945726D-12,
     -   0.123442456067320815813898019612D-11,
     -   0.379576946609534375379531510687D-11,
     -   0.232126091962106640268482862678D-10,
     -   0.160394534422147537071291035793D-07,
     -   0.278342831108477173651324716553D-01,
     -   0.627901837234858954442634690861D-01,
     -   0.931451039671291016103128903899D-01,
     -   0.116596880422412797923668886924D+00,
     -   0.131402270143634562845777598240D+00,
     -   0.136462541196148547879099863708D+00,
     -   0.131402270143634487177709875976D+00,
     -   0.116596880422412589564234348603D+00,
     -   0.931451039671285283764230863819D-01,
     -   0.627901837234836424343078170283D-01,
     -   0.278342831108207486022906504007D-01/
        data xs04 /
     -   0.267686522719226961500023153083D-13,
     -   0.147154449055675760810941579200D-12,
     -   0.390571543448152128078238193913D-12,
     -   0.811745966945498853682799670835D-12,
     -   0.151414735456613686381259386608D-11,
     -   0.269534484061357790239305387850D-11,
     -   0.475920114689550965969657385753D-11,
     -   0.863393628316889316666572497000D-11,
     -   0.168741928158123006098845577416D-10,
     -   0.388982652220017983068227925464D-10,
     -   0.136601282174951445294323098402D-09,
     -   0.769447934800507235348931590774D-08,
     -   0.108857908772574273319029741804D-01,
     -   0.564688145388156528912360727770D-01,
     -   0.134924102121553575784135892186D+00,
     -   0.240452027507716156115171796388D+00,
     -   0.365228499003191198905092907468D+00,
     -   0.500000060635486342379668652208D+00,
     -   0.634771622267779473159740552205D+00,
     -   0.759548093763247052858068120430D+00,
     -   0.865076019149389773974401784383D+00,
     -   0.943531306732062179807613831971D+00,
     -   0.989114330393144467185694363165D+00/
        data whts04 /
     -   0.693779088607258798818768886520D-13,
     -   0.175602545731113322745530113318D-12,
     -   0.320048183016953532120952354560D-12,
     -   0.538869062479735220153125891012D-12,
     -   0.897240876178376219712398763724D-12,
     -   0.152809183728042005623954744244D-11,
     -   0.274037726829446314554989706504D-11,
     -   0.537822534171303947158174693956D-11,
     -   0.123351842291600544501262808660D-10,
     -   0.378983847088233194793381627536D-10,
     -   0.230913855737315904153042867029D-09,
     -   0.120976494776320731162258544143D-06,
     -   0.278342801841043901121816410976D-01,
     -   0.627901771179511953161117910205D-01,
     -   0.931450941681012501097991853613D-01,
     -   0.116596868156189366130237979511D+00,
     -   0.131402256319846600561798512318D+00,
     -   0.136462526840006696401427875591D+00,
     -   0.131402256319842368321028379723D+00,
     -   0.116596868156177712244365179622D+00,
     -   0.931450941680691882214688294012D-01,
     -   0.627901771178251812021711696572D-01,
     -   0.278342801825960173598241895094D-01/
        data xs05 /
     -   0.267664561062760787506823617034D-12,
     -   0.147141597516936008385678923762D-11,
     -   0.390533471979336155557182972210D-11,
     -   0.811653629569469067538875047346D-11,
     -   0.151393746550076665136524743012D-10,
     -   0.269486807666356701247329490026D-10,
     -   0.475806548430608930700573021355D-10,
     -   0.863093614233199142519520312728D-10,
     -   0.168645852162843181295662084816D-09,
     -   0.388535916481259408530167136437D-09,
     -   0.136098030742159212449705304225D-08,
     -   0.648802398945146819165288657929D-07,
     -   0.108865355798783753450916968760D-01,
     -   0.564695249443752328542204956767D-01,
     -   0.134924753459229956807896385650D+00,
     -   0.240452599391518353107546555806D+00,
     -   0.365228976939768361667934619864D+00,
     -   0.500000437099038957399915632684D+00,
     -   0.634771897258207515724786108417D+00,
     -   0.759548274806079169062895464724D+00,
     -   0.865076120737360769727879341655D+00,
     -   0.943531349248894039136044748773D+00,
     -   0.989114338589264602952642654172D+00/
        data whts05 /
     -   0.693721336377307415434985057063D-12,
     -   0.175586095750571584874858608389D-11,
     -   0.320011939297603695759829039590D-11,
     -   0.538791512569590563985679196020D-11,
     -   0.897070364055981192269096605781D-11,
     -   0.152769317548831928590744860018D-10,
     -   0.273934531155548708972788049111D-10,
     -   0.537508045289962839666434881463D-10,
     -   0.123226272688976208578057928399D-09,
     -   0.378163260990743469618038825431D-09,
     -   0.229248409159990781929324275577D-08,
     -   0.871202482304909957438782314220D-06,
     -   0.278342593018580897268459227022D-01,
     -   0.627901298477639376200581186256D-01,
     -   0.931450240382002577240639324198D-01,
     -   0.116596780367805537971811326711D+00,
     -   0.131402157383744666788764833373D+00,
     -   0.136462424093746041995681666868D+00,
     -   0.131402157383530105289990513328D+00,
     -   0.116596780367214723322540925618D+00,
     -   0.931450240365748320763200734023D-01,
     -   0.627901298413755831095434435449D-01,
     -   0.278342592254007843322451816279D-01/
        data xs06 /
     -   0.267632060425911027063870138371D-11,
     -   0.147122579052593473525901646761D-10,
     -   0.390477133145719848554162912269D-10,
     -   0.811516993029403206815276212087D-10,
     -   0.151362690163968938693204362027D-09,
     -   0.269416269919567023808568015670D-09,
     -   0.475638552731990480261463256002D-09,
     -   0.862649937158265270870497118046D-09,
     -   0.168503853206348061809181277851D-08,
     -   0.387876669188483568594339708027D-08,
     -   0.135360419026749565392374202498D-07,
     -   0.528248640724461243249485030511D-06,
     -   0.108915214669748337840516562669D-01,
     -   0.564742820421866162130015574451D-01,
     -   0.134929115126311619200826978291D+00,
     -   0.240456429027125455211165433919D+00,
     -   0.365232177465962687383174956475D+00,
     -   0.500002958110452646200061390882D+00,
     -   0.634773738750430870274485736596D+00,
     -   0.759549487172538717154948243457D+00,
     -   0.865076801028837146167291117572D+00,
     -   0.943531633966116101668676292486D+00,
     -   0.989114393475217400372960967375D+00/
        data whts06 /
     -   0.693635870239917701480551055205D-11,
     -   0.175561752501101268565853526113D-10,
     -   0.319958307174201011906213352229D-10,
     -   0.538676765057133994034453924431D-10,
     -   0.896818089117468422053018101930D-10,
     -   0.152710343766519137938570967715D-09,
     -   0.273781911027731585602386847677D-09,
     -   0.537043134404453348174822875115D-09,
     -   0.123040806511712873494569111253D-08,
     -   0.376954074284839712178345229211D-08,
     -   0.226819172770770534867235212797D-07,
     -   0.588360782818418614626752041344D-05,
     -   0.278341222609568502570453066942D-01,
     -   0.627898135334765809877420349495D-01,
     -   0.931445544678236133604194902353D-01,
     -   0.116596192508758443516813333339D+00,
     -   0.131401494860651808134212980983D+00,
     -   0.136461736050327640298075502667D+00,
     -   0.131401494851164672460634632908D+00,
     -   0.116596192482635119368262471765D+00,
     -   0.931445543959568686795969816220D-01,
     -   0.627898132510561962566290048918D-01,
     -   0.278341188839246981657489305847D-01/
        data xs07 /
     -   0.267579255379801141219274707555D-10,
     -   0.147091679733585358579884199572D-09,
     -   0.390385603375445301674236389107D-09,
     -   0.811295024342667142303487412930D-09,
     -   0.151312243633504373170096084034D-08,
     -   0.269301709635843739588251221373D-08,
     -   0.475365780175460415859448901630D-08,
     -   0.861929869798342691239964962850D-08,
     -   0.168273610860043952468121929346D-07,
     -   0.386810371638372449919323433197D-07,
     -   0.134180134416298675967719266251D-06,
     -   0.408145193734912454732061770977D-05,
     -   0.109218149085827588127245393708D-01,
     -   0.565032148102677661349305587721D-01,
     -   0.134955646636720777756751668481D+00,
     -   0.240479725295942684336386196189D+00,
     -   0.365251647145465864393978176077D+00,
     -   0.500018294270916680642220647909D+00,
     -   0.634784941232614675009295551903D+00,
     -   0.759556862474979193473923674428D+00,
     -   0.865080939518852774460591801738D+00,
     -   0.943533366019187674052340296102D+00,
     -   0.989114727369480678390787143670D+00/
        data whts07 /
     -   0.693497010817498978318789937967D-10,
     -   0.175522203017038268021651984040D-09,
     -   0.319871179834753987425122137021D-09,
     -   0.538490373977178170436946614227D-09,
     -   0.896408368524595542453418751236D-09,
     -   0.152614586794607040168574781650D-08,
     -   0.273534192897810816620981115794D-08,
     -   0.536289050548686212980663808865D-08,
     -   0.122740404684735594039363509529D-07,
     -   0.375002779822450760095010870762D-07,
     -   0.222961723360860208870332647018D-06,
     -   0.361682234022504770131592794733D-04,
     -   0.278333869210745135068564184003D-01,
     -   0.627878975420433647500722389864D-01,
     -   0.931416999955390614528649282524D-01,
     -   0.116592617117896768990439906367D+00,
     -   0.131397464815002735174379486409D+00,
     -   0.136457550566093438369461626831D+00,
     -   0.131397464470672265143274957717D+00,
     -   0.116592616169843983383521188460D+00,
     -   0.931416973879714934458240529244D-01,
     -   0.627878873023835741296390192599D-01,
     -   0.278332651280146107514521694722D-01/
        data xs08 /
     -   0.267480836452085520364054605909D-09,
     -   0.147034091393862451119596702437D-08,
     -   0.390215028549634899450689609873D-08,
     -   0.810881413297530231365773542572D-08,
     -   0.151218259843884383686041920653D-07,
     -   0.269088339040143348771412310760D-07,
     -   0.474857967533920827096492088989D-07,
     -   0.860590436177809210833295337372D-07,
     -   0.167846050342834214423114167681D-06,
     -   0.384839065587687397029019488254D-06,
     -   0.132039698298693905606283310006D-05,
     -   0.290243756278616246122343386783D-04,
     -   0.110811761819231226565387951181D-01,
     -   0.566561716082703106222147016944D-01,
     -   0.135096006930221033397866461911D+00,
     -   0.240602997404842219805543768512D+00,
     -   0.365354680848568800179210960002D+00,
     -   0.500099457441257686066692336626D+00,
     -   0.634844229559845925322267608954D+00,
     -   0.759595896419071740671486518149D+00,
     -   0.865102842804715654930102944018D+00,
     -   0.943542533107225820827355239202D+00,
     -   0.989116494548764098516776798549D+00/
        data whts08 /
     -   0.693238204660600422881324829821D-09,
     -   0.175448496314544137915709944419D-08,
     -   0.319708825787437095211464777462D-08,
     -   0.538143118054026786833189760792D-08,
     -   0.895645252007055676715135193427D-08,
     -   0.152436312047246558323680363957D-07,
     -   0.273073322920702315142578674561D-07,
     -   0.534887839062861192552467329432D-07,
     -   0.122183634392048535301432166830D-06,
     -   0.371410380417476971980561611824D-06,
     -   0.216061662548947021212164751799D-05,
     -   0.192535279173455567753081783543D-03,
     -   0.278319704659329638014401964782D-01,
     -   0.627779723338580183136195835193D-01,
     -   0.931266477634879805541123911071D-01,
     -   0.116573714887427080138456450651D+00,
     -   0.131376144699289949760114610502D+00,
     -   0.136435402855817680848513123378D+00,
     -   0.131376135291933319890594842546D+00,
     -   0.116573688996906143298390240302D+00,
     -   0.931265766331238422033750234205D-01,
     -   0.627776940439143475747649097863D-01,
     -   0.278287465160523896809021513766D-01/
        data xs09 /
     -   0.267260082339909213007709355279D-08,
     -   0.146904931136334485363382869262D-07,
     -   0.389832519884835961740227011914D-07,
     -   0.809954134686622785914111217070D-07,
     -   0.151007636852211247361219195641D-06,
     -   0.268610446459253990983365055033D-06,
     -   0.473721712818813551889525294967D-06,
     -   0.857598668537218562661303996620D-06,
     -   0.166894543194720268805230883852D-05,
     -   0.380494038120030252802024379271D-05,
     -   0.127509050878168777584780033253D-04,
     -   0.180990125105267091958923078723D-03,
     -   0.117672165951606430841044079917D-01,
     -   0.573267851841418856809980659645D-01,
     -   0.135713093302128748278226920273D+00,
     -   0.241145435705010282042214435850D+00,
     -   0.365808239351547106122443739728D+00,
     -   0.500456813146884845627199829388D+00,
     -   0.635105302250114516348034660216D+00,
     -   0.759767792312295431464669856006D+00,
     -   0.865199303740348660327183161529D+00,
     -   0.943582905552575757294284577154D+00,
     -   0.989124277424315539748566211372D+00/
        data whts09 /
     -   0.692657712006512233430933740782D-08,
     -   0.175283200658981661169703425537D-07,
     -   0.319344827910350647348673419844D-07,
     -   0.537364885370597321087596317598D-07,
     -   0.893936050005209192580944692147D-07,
     -   0.152037379726445987135104663877D-06,
     -   0.272043546301485486540440253216D-06,
     -   0.531765318798631275875531006968D-06,
     -   0.120949763746634987284036786744D-05,
     -   0.363563390560938576709185209969D-05,
     -   0.201872185423833469200319492487D-04,
     -   0.826296233122839453126306394608D-03,
     -   0.278628162114905501343738976658D-01,
     -   0.627379572768000062647870447115D-01,
     -   0.930613264306981187366774865549D-01,
     -   0.116490837023893770890865679482D+00,
     -   0.131282414869586535582713770743D+00,
     -   0.136337941475147330731127531041D+00,
     -   0.131282236836126979373251344627D+00,
     -   0.116490347909009202413448663539D+00,
     -   0.930599889065731720067059275792D-01,
     -   0.627328030593926179450480496760D-01,
     -   0.278088460523550762435602297414D-01/
        data xs10 /
     -   0.266661601230813282161112597242D-07,
     -   0.146554840804492594798503079166D-06,
     -   0.388796137629465810710817227899D-06,
     -   0.807443332931235300030405401013D-06,
     -   0.150437889419854884984565360256D-05,
     -   0.267319702353941671838798849192D-05,
     -   0.470660584989431766443727055574D-05,
     -   0.849575716355908533061479560843D-05,
     -   0.164367133353139603171560202041D-04,
     -   0.369233312915882025539350468178D-04,
     -   0.116854771599272572679920896905D-03,
     -   0.942065164963600019690107727105D-03,
     -   0.141272350936745103655033213646D-01,
     -   0.597432204291798826946708572716D-01,
     -   0.137955831999423853963228994525D+00,
     -   0.243122524407815796096051509737D+00,
     -   0.367463485278628007410336488854D+00,
     -   0.501761841503866098438592512697D+00,
     -   0.636059082531671832688525677667D+00,
     -   0.760395930030773451030082171423D+00,
     -   0.865551840612753299864447514635D+00,
     -   0.943730468152978793686353905293D+00,
     -   0.989152725358661802157247212346D+00/
        data whts10 /
     -   0.691084025946395910485673617827D-07,
     -   0.174835268984774680866411930892D-06,
     -   0.318359123084204903918110087543D-06,
     -   0.535259626715994855014291109263D-06,
     -   0.889319445938652268431246534214D-06,
     -   0.150962393743254414560393174203D-05,
     -   0.269279416469810919571102512167D-05,
     -   0.523442379027785595819074177038D-05,
     -   0.117708020547867719299927605747D-04,
     -   0.343692926137501105744283693956D-04,
     -   0.170736988491111940110769945496D-03,
     -   0.270236803258543239643662511665D-02,
     -   0.282398353099354025432201999630D-01,
     -   0.626318877636348735354819607562D-01,
     -   0.928338882909146240191049784488D-01,
     -   0.116192317394997192806347489396D+00,
     -   0.130941822035656819346179268237D+00,
     -   0.135982668494310822936112493178D+00,
     -   0.130939479666647255107187173682D+00,
     -   0.116185920839911635726987130981D+00,
     -   0.928166700925224927517767355855D-01,
     -   0.625687345711484331434489415322D-01,
     -   0.277361067008156400895030277231D-01/
        data xs11 /
     -   0.165773595586393687688305586940D-06,
     -   0.987007527401057263390961503627D-06,
     -   0.278749857686556487273803143730D-05,
     -   0.599950530734901068499509018736D-05,
     -   0.113577355340006294292148262592D-04,
     -   0.202005926495574255900987495781D-04,
     -   0.350898121758786185852872291356D-04,
     -   0.613406672304223573376604015077D-04,
     -   0.111387693283958586319485501586D-03,
     -   0.219707163938073517834836511045D-03,
     -   0.507668311994946343637415777374D-03,
     -   0.155489296497452395800015000678D-02,
     -   0.650028838676379176984196585595D-02,
     -   0.248077863338179274959918321891D-01,
     -   0.663792818131063305436032051967D-01,
     -   0.134055318351603466317525899725D+00,
     -   0.225400835977478360938350252275D+00,
     -   0.335183774514175457032840912251D+00,
     -   0.456429065751093509734581903783D+00,
     -   0.581103394159723278185450663201D+00,
     -   0.700813510090293886684822208028D+00,
     -   0.807560286063491301305906786545D+00,
     -   0.894521913485484787409057446140D+00,
     -   0.956789652720423366388206360461D+00,
     -   0.991793155537284842212933938414D+00/
        data whts11 /
     -   0.441264708587354076126869212982D-06,
     -   0.125314313922452308132252840103D-05,
     -   0.241655162741216475129376831387D-05,
     -   0.412538882674355248893978098316D-05,
     -   0.680494552139394996536721574169D-05,
     -   0.112856858802988037205038131768D-04,
     -   0.193185974559892558000553325604D-04,
     -   0.350641628391287207498593934287D-04,
     -   0.700266883457546640952588763148D-04,
     -   0.163022789232142809048181126556D-03,
     -   0.484589011774134350142053335780D-03,
     -   0.202275770306322243018267732984D-02,
     -   0.957992813314665778259078443459D-02,
     -   0.288411890581389732256238403326D-01,
     -   0.547076409412421589419050467352D-01,
     -   0.801989582366871629957366067817D-01,
     -   0.101605494876393672008624313916D+00,
     -   0.116787072712481014567168955204D+00,
     -   0.124353002944099472754387369289D+00,
     -   0.123584143184007407727934415676D+00,
     -   0.114489722715172313774393150447D+00,
     -   0.978540676816897358868048972929D-01,
     -   0.752425125663864023205199094135D-01,
     -   0.488701418604193093014716693113D-01,
     -   0.210650191577216861347374017963D-01/
        data xs12 /
     -   0.260960129465153694668711910407D-05,
     -   0.143122684072595061305179762571D-04,
     -   0.378196658468606655414872321949D-04,
     -   0.780599434293019020424346363043D-04,
     -   0.144122814250838302780279485980D-03,
     -   0.252714796286855347073150721581D-03,
     -   0.436011856293258357479271207780D-03,
     -   0.761155997674179080161289839097D-03,
     -   0.138442845090607347156451661966D-02,
     -   0.272180277437471479199405561781D-02,
     -   0.603869594720797436752118909977D-02,
     -   0.152845002378902915295903590858D-01,
     -   0.395381751329528968741085597046D-01,
     -   0.888380208347826173652149322331D-01,
     -   0.166335005810439836729054593788D+00,
     -   0.268696623506074655417730844016D+00,
     -   0.389112943296235999285271515101D+00,
     -   0.518936538695611521282851429617D+00,
     -   0.648657791156365170097309893283D+00,
     -   0.768712365346300902904191586640D+00,
     -   0.870226176825445766947075591664D+00,
     -   0.945688745366183385373098539886D+00,
     -   0.989530422766363869465582116013D+00/
        data whts12 /
     -   0.675986818646091969990803823281D-05,
     -   0.170315756629347270023351140391D-04,
     -   0.307797952809939207750478482916D-04,
     -   0.511632391263327351389750369360D-04,
     -   0.836300709985647714462048090868D-04,
     -   0.138634376734875980256980360976D-03,
     -   0.238393101818338986475808434635D-03,
     -   0.435648266794447799117710221665D-03,
     -   0.872391125360902825575693416475D-03,
     -   0.198560441356809481609007837675D-02,
     -   0.524257238431768152844482122566D-02,
     -   0.148264522888801639671918316258D-01,
     -   0.355451226183075025216550050323D-01,
     -   0.635351273568980666128392996735D-01,
     -   0.908491563724856637927908581943D-01,
     -   0.112714832610891317924617853722D+00,
     -   0.126661267439901076615948829961D+00,
     -   0.131387682280474042063822103121D+00,
     -   0.126448302192242007514541300932D+00,
     -   0.112170524402461011289832167277D+00,
     -   0.895957832332046553773597857151D-01,
     -   0.603926806201014143813821380063D-01,
     -   0.267704603663034489279952638569D-01/
        data xs13 /
     -   0.255184284168003991470618043153D-04,
     -   0.139702387695868648240015558749D-03,
     -   0.367787214470717550109998507783D-03,
     -   0.754135373274810312303528365478D-03,
     -   0.137698667994145213476606929494D-02,
     -   0.237035152084715060648796150629D-02,
     -   0.396633047356742464992768932844D-02,
     -   0.657823497087154578548247546015D-02,
     -   0.109635693636378456963335077909D-01,
     -   0.185267573572395335197290013081D-01,
     -   0.317893229932433394219713891260D-01,
     -   0.548027888975051800073892064966D-01,
     -   0.927621722948761005752618634701D-01,
     -   0.150224593084596748332392857949D+00,
     -   0.228960510899335066877004330793D+00,
     -   0.327109204376686262452007985967D+00,
     -   0.439702528810363015951610188309D+00,
     -   0.559672062301687161273139268231D+00,
     -   0.678834987263728125832016155994D+00,
     -   0.788762859372694171432103877949D+00,
     -   0.881544098560213488628730588978D+00,
     -   0.950442892436817881859089654226D+00,
     -   0.990448627297968080575160521448D+00/
        data whts13 /
     -   0.660761559794395653460387181862D-04,
     -   0.165874915829758235856726977307D-03,
     -   0.297507922350347025329232843994D-03,
     -   0.487885231950121123621440682166D-03,
     -   0.779492051117246406824791648668D-03,
     -   0.124462527464087687346420755025D-02,
     -   0.201365536400831884147972271388D-02,
     -   0.333150466415599414963820960708D-02,
     -   0.566504854805033781113694636149D-02,
     -   0.987277572128480471019908196034D-02,
     -   0.173316096764730412194068826477D-01,
     -   0.295937384080642247244946096726D-01,
     -   0.471258877515523884824816787883D-01,
     -   0.681075337056833509953823817151D-01,
     -   0.890427018310177480184671037962D-01,
     -   0.106408805179134336524200434643D+00,
     -   0.117583252524164449046472220084D+00,
     -   0.120977590185036052227241674130D+00,
     -   0.115930814673935911598605727458D+00,
     -   0.102595314217166673825187376137D+00,
     -   0.818354420280113817661481936115D-01,
     -   0.551192062498095710985880237229D-01,
     -   0.244236577205836257304272945291D-01/
        data xs14 /
     -   0.255620118286698236593623243707D-03,
     -   0.139529815173999277362604891093D-02,
     -   0.365117643857364172245052441454D-02,
     -   0.740789208402439896041673791295D-02,
     -   0.132924753570324232001855699689D-01,
     -   0.222520965750031001463918446238D-01,
     -   0.356419901326391687523419863309D-01,
     -   0.552880334040421002669845054888D-01,
     -   0.834528839374245293455070701521D-01,
     -   0.122611411148390814700485075280D+00,
     -   0.174981552629270045448068775415D+00,
     -   0.241886149381019235661855516525D+00,
     -   0.323158757145366236802927707224D+00,
     -   0.416816795924166714866947223698D+00,
     -   0.519095046421213331384873647907D+00,
     -   0.624786438056543680865135919467D+00,
     -   0.727771978150004968382920198445D+00,
     -   0.821629005471921763485910289821D+00,
     -   0.900237992686259218054179268073D+00,
     -   0.958334484420391012169849137700D+00,
     -   0.991976844522788369627739830890D+00/
        data whts14 /
     -   0.661460163853639282285754410062D-03,
     -   0.165066844620440621583106622338D-02,
     -   0.292424961211204245632909914117D-02,
     -   0.469237052857572175143233064178D-02,
     -   0.723311364138835119287276977376D-02,
     -   0.109107870131761110458292697881D-01,
     -   0.161736269683573622532714564623D-01,
     -   0.234991383631077678365169701286D-01,
     -   0.332510295063762612773746804267D-01,
     -   0.454468870542399720350907864669D-01,
     -   0.595243742379457025186539892980D-01,
     -   0.742627490854230932327045399807D-01,
     -   0.879513094972176949552737819238D-01,
     -   0.987312394302742140900606516263D-01,
     -   0.104948313524108008835772391594D+00,
     -   0.105404313986115926771361097012D+00,
     -   0.994860857607223310260098949164D-01,
     -   0.871987663624295826569173723217D-01,
     -   0.691340460686524471960186899166D-01,
     -   0.463953172998705895793437234323D-01,
     -   0.205201534498487737910496845150D-01/
        data xs15 /
     -   0.224878755802716701158898892336D-02,
     -   0.120978441451046350231217671129D-01,
     -   0.308010595288933550618497213370D-01,
     -   0.599007211841317118347310619055D-01,
     -   0.101238700758448386489421898375D+00,
     -   0.156472845604625144895887745341D+00,
     -   0.226469953776845163020037600698D+00,
     -   0.310745922381071052403874785567D+00,
     -   0.407138865165486828616759733871D+00,
     -   0.511808476950210078921105444312D+00,
     -   0.619532519225550057842970342845D+00,
     -   0.724201777562086340710448946892D+00,
     -   0.819408185576637775373713348731D+00,
     -   0.899043266341426000136481319347D+00,
     -   0.957848511044927462239314254803D+00,
     -   0.991884610125741879406251886235D+00/
        data whts15 /
     -   0.579981455141374776030977420111D-02,
     -   0.140616241959305547578683669699D-01,
     -   0.236064728488424861764670548228D-01,
     -   0.349080664994175133689398803908D-01,
     -   0.480547600935841734972117482382D-01,
     -   0.625678230629187072715698701715D-01,
     -   0.773513908037960630392214546386D-01,
     -   0.908401449009560445347792959706D-01,
     -   0.101300485377155168056557937309D+00,
     -   0.107159854779374366123576055100D+00,
     -   0.107260028997342519538822442960D+00,
     -   0.101000366190171636898247903878D+00,
     -   0.883851058915990086514884370775D-01,
     -   0.700008559625108744958384922932D-01,
     -   0.469463534679938482307195329563D-01,
     -   0.207568523769932875983817530223D-01/
        data xs16 /
     -   0.458849318721402125801204605826D-02,
     -   0.243824898995194101441372878371D-01,
     -   0.607030321190286710240555584291D-01,
     -   0.114287024238690622102246875135D+00,
     -   0.185268852214541964722038987301D+00,
     -   0.272558413166421274952991402464D+00,
     -   0.373456967272527844350424089652D+00,
     -   0.483618123659061296715748101316D+00,
     -   0.597325895505307761757764267100D+00,
     -   0.707988920509439451172599826265D+00,
     -   0.808741228229900913864900441481D+00,
     -   0.893061329285272811648240079303D+00,
     -   0.955346062604670722675744052617D+00,
     -   0.991402329456371986670375045434D+00/
        data whts16 /
     -   0.118003104023517973910349882650D-01,
     -   0.279179096171351173393330264437D-01,
     -   0.448569113929324728417964582458D-01,
     -   0.623444612648363615874905577212D-01,
     -   0.794441657665294476712848268693D-01,
     -   0.946857191970948830773180132808D-01,
     -   0.106382087827940363439416553076D+00,
     -   0.112979051715145575784047173658D+00,
     -   0.113327709140967997487895723878D+00,
     -   0.106843545232139922936112549910D+00,
     -   0.935644897539231560178351917563D-01,
     -   0.741336945635054014592323288119D-01,
     -   0.497298407091702578638931877252D-01,
     -   0.219901034163272451033094203592D-01/
        data xs17 /
     -   0.636686804373680661692959691881D-02,
     -   0.335389468182532697676009618757D-01,
     -   0.822565040403119421130085395349D-01,
     -   0.151772204336725861831406221530D+00,
     -   0.240199439899564910406030864436D+00,
     -   0.344143951776306330555142848942D+00,
     -   0.458640212965061843423712636537D+00,
     -   0.577404439344222559648064287204D+00,
     -   0.693317583068545243516958190313D+00,
     -   0.799029641357120636056253623013D+00,
     -   0.887592492147654876964130334999D+00,
     -   0.953052189744624013874024303991D+00,
     -   0.990959627067690458415587096025D+00/
        data whts17 /
     -   0.163399073935820539878321623664D-01,
     -   0.379995130838314092843639240668D-01,
     -   0.593198540744717530564190368123D-01,
     -   0.794037139335967621445449176309D-01,
     -   0.968860947164711007457054390048D-01,
     -   0.110170885968530591996674490108D+00,
     -   0.117766016421497247733649456948D+00,
     -   0.118565634285020287747364798855D+00,
     -   0.112024942620742142010371774159D+00,
     -   0.982311696218582656820687396083D-01,
     -   0.778937323971055574596908753293D-01,
     -   0.522767532896572832139953157464D-01,
     -   0.231217821936355450487039646273D-01/
        data xs18 /
     -   0.542463833175336435357423617512D-02,
     -   0.297460765201093519581047030093D-01,
     -   0.759496873372568959207175963510D-01,
     -   0.144306061435524010204963639011D+00,
     -   0.232806065105593800908487503522D+00,
     -   0.337605308083923075124117055606D+00,
     -   0.453324505155217740534635788505D+00,
     -   0.573393063871370009950573815542D+00,
     -   0.690520828875537830879257614391D+00,
     -   0.797266491425120807075842336693D+00,
     -   0.886638078590384985700357901137D+00,
     -   0.952663036764617210031135426210D+00,
     -   0.990885717184502181805668524724D+00/
        data whts18 /
     -   0.140900107050361978359013576230D-01,
     -   0.350033236829359129515199850767D-01,
     -   0.574427833593893518517421315009D-01,
     -   0.789272359512492998811609264637D-01,
     -   0.974277422045807687042996534461D-01,
     -   0.111268619345890655734933450420D+00,
     -   0.119066705482167043491045893711D+00,
     -   0.119845747329461976383102753883D+00,
     -   0.113159404104991890268372110461D+00,
     -   0.991569734643610367907930540647D-01,
     -   0.785825132788112287563500124571D-01,
     -   0.527175067883675225617246712246D-01,
     -   0.233114343027571147890539996670D-01/
        data xs19 /
     -   0.754490802642326661397879954159D-02,
     -   0.401709054985051799770663753398D-01,
     -   0.994300314981247707748489199152D-01,
     -   0.183848016720589490919373740433D+00,
     -   0.289336499333354591530486538411D+00,
     -   0.409644336296503701114177295372D+00,
     -   0.537017446804686363511342483124D+00,
     -   0.662885286561024055236459659997D+00,
     -   0.778562898774300337348116758755D+00,
     -   0.875942564294336725220958287492D+00,
     -   0.948132272620363553445778183833D+00,
     -   0.990006707054016128055540695770D+00/
        data whts19 /
     -   0.194253789907764398316504870084D-01,
     -   0.459695812646557596593129067333D-01,
     -   0.723089934411016624606484471392D-01,
     -   0.958361268047063154442418387979D-01,
     -   0.114087815437017924368818433069D+00,
     -   0.125224085745238765152674755616D+00,
     -   0.128081940396601387080838759133D+00,
     -   0.122194609382496598897982563975D+00,
     -   0.107798584480839372152599976557D+00,
     -   0.858029374971536432048800500511D-01,
     -   0.577142300584027660270561493626D-01,
     -   0.255557165010093657294232085307D-01/
        data xs20 /
     -   0.108856709269715035980309994386D-01,
     -   0.564687001159523504624211153480D-01,
     -   0.134923997212975337953291873984D+00,
     -   0.240451935396594092037137165271D+00,
     -   0.365228422023827513834234007300D+00,
     -   0.500000000000000000000000000000D+00,
     -   0.634771577976172486165765992700D+00,
     -   0.759548064603405907962862834729D+00,
     -   0.865076002787024662046708126016D+00,
     -   0.943531299884047649537578884652D+00,
     -   0.989114329073028496401969000561D+00/
        data whts20 /
     -   0.278342835580868332421930390504D-01,
     -   0.627901847324523123191883329622D-01,
     -   0.931451054638671257157800950076D-01,
     -   0.116596882295995239962680798222D+00,
     -   0.131402272255123331094197516074D+00,
     -   0.136462543388950315361243226561D+00,
     -   0.131402272255123331094197516074D+00,
     -   0.116596882295995239962680798222D+00,
     -   0.931451054638671257157800950076D-01,
     -   0.627901847324523123191883329622D-01,
     -   0.278342835580868332421930390504D-01/
        data xs21 /
     -   0.108856709269715035980309994386D-01,
     -   0.564687001159523504624211153480D-01,
     -   0.134923997212975337953291873984D+00,
     -   0.240451935396594092037137165271D+00,
     -   0.365228422023827513834234007300D+00,
     -   0.500000000000000000000000000000D+00,
     -   0.634771577976172486165765992700D+00,
     -   0.759548064603405907962862834729D+00,
     -   0.865076002787024662046708126016D+00,
     -   0.943531299884047649537578884652D+00,
     -   0.989114329073028496401969000561D+00/
        data whts21 /
     -   0.278342835580868332413771817394D-01,
     -   0.627901847324523123173478749112D-01,
     -   0.931451054638671257130498966495D-01,
     -   0.116596882295995239959263199250D+00,
     -   0.131402272255123331090345952793D+00,
     -   0.136462543388950315357243340468D+00,
     -   0.131402272255123331090345952793D+00,
     -   0.116596882295995239959263199250D+00,
     -   0.931451054638671257130498966495D-01,
     -   0.627901847324523123173478749112D-01,
     -   0.278342835580868332413771817394D-01/
        data xs22 /
     -   0.108856709269715035980309994386D-01,
     -   0.564687001159523504624211153479D-01,
     -   0.134923997212975337953291873984D+00,
     -   0.240451935396594092037137165270D+00,
     -   0.365228422023827513834234007299D+00,
     -   0.500000000000000000000000000000D+00,
     -   0.634771577976172486165765992700D+00,
     -   0.759548064603405907962862834729D+00,
     -   0.865076002787024662046708126015D+00,
     -   0.943531299884047649537578884652D+00,
     -   0.989114329073028496401969000561D+00/
        data whts22 /
     -   0.278342835580868332413768602322D-01,
     -   0.627901847324523123173471496367D-01,
     -   0.931451054638671257130488207525D-01,
     -   0.116596882295995239959261852468D+00,
     -   0.131402272255123331090344434997D+00,
     -   0.136462543388950315357241764222D+00,
     -   0.131402272255123331090344434997D+00,
     -   0.116596882295995239959261852468D+00,
     -   0.931451054638671257130488207528D-01,
     -   0.627901847324523123173471496369D-01,
     -   0.278342835580868332413768602323D-01/
        data xs23 /
     -   0.108856709269715035980309994386D-01,
     -   0.564687001159523504624211153479D-01,
     -   0.134923997212975337953291873984D+00,
     -   0.240451935396594092037137165270D+00,
     -   0.365228422023827513834234007299D+00,
     -   0.500000000000000000000000000000D+00,
     -   0.634771577976172486165765992700D+00,
     -   0.759548064603405907962862834729D+00,
     -   0.865076002787024662046708126015D+00,
     -   0.943531299884047649537578884652D+00,
     -   0.989114329073028496401969000561D+00/
        data whts23 /
     -   0.278342835580868332413768602212D-01,
     -   0.627901847324523123173471496119D-01,
     -   0.931451054638671257130488207157D-01,
     -   0.116596882295995239959261852421D+00,
     -   0.131402272255123331090344434945D+00,
     -   0.136462543388950315357241764168D+00,
     -   0.131402272255123331090344434945D+00,
     -   0.116596882295995239959261852422D+00,
     -   0.931451054638671257130488207160D-01,
     -   0.627901847324523123173471496121D-01,
     -   0.278342835580868332413768602213D-01/

        ier = 0

        delta = r0/r1

        if( deltas(1,01) .le. delta .AND. delta .le. deltas(2,01)) then
        nquad = 24
        do l=1,24
        xs(l)   = xs01(l)
        whts(l) = whts01(l)
        end do
        goto 1000
        end if


        if( deltas(1,02) .le. delta .AND. delta .le. deltas(2,02)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs02(l)
        whts(l) = whts02(l)
        end do
        goto 1000
        end if


        if( deltas(1,03) .le. delta .AND. delta .le. deltas(2,03)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs03(l)
        whts(l) = whts03(l)
        end do
        goto 1000
        end if


        if( deltas(1,04) .le. delta .AND. delta .le. deltas(2,04)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs04(l)
        whts(l) = whts04(l)
        end do
        goto 1000
        end if


        if( deltas(1,05) .le. delta .AND. delta .le. deltas(2,05)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs05(l)
        whts(l) = whts05(l)
        end do
        goto 1000
        end if


        if( deltas(1,06) .le. delta .AND. delta .le. deltas(2,06)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs06(l)
        whts(l) = whts06(l)
        end do
        goto 1000
        end if


        if( deltas(1,07) .le. delta .AND. delta .le. deltas(2,07)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs07(l)
        whts(l) = whts07(l)
        end do
        goto 1000
        end if


        if( deltas(1,08) .le. delta .AND. delta .le. deltas(2,08)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs08(l)
        whts(l) = whts08(l)
        end do
        goto 1000
        end if


        if( deltas(1,09) .le. delta .AND. delta .le. deltas(2,09)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs09(l)
        whts(l) = whts09(l)
        end do
        goto 1000
        end if


        if( deltas(1,10) .le. delta .AND. delta .le. deltas(2,10)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs10(l)
        whts(l) = whts10(l)
        end do
        goto 1000
        end if


        if( deltas(1,11) .le. delta .AND. delta .le. deltas(2,11)) then
        nquad = 25
        do l=1,25
        xs(l)   = xs11(l)
        whts(l) = whts11(l)
        end do
        goto 1000
        end if


        if( deltas(1,12) .le. delta .AND. delta .le. deltas(2,12)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs12(l)
        whts(l) = whts12(l)
        end do
        goto 1000
        end if


        if( deltas(1,13) .le. delta .AND. delta .le. deltas(2,13)) then
        nquad = 23
        do l=1,23
        xs(l)   = xs13(l)
        whts(l) = whts13(l)
        end do
        goto 1000
        end if


        if( deltas(1,14) .le. delta .AND. delta .le. deltas(2,14)) then
        nquad = 21
        do l=1,21
        xs(l)   = xs14(l)
        whts(l) = whts14(l)
        end do
        goto 1000
        end if


        if( deltas(1,15) .le. delta .AND. delta .le. deltas(2,15)) then
        nquad = 16
        do l=1,16
        xs(l)   = xs15(l)
        whts(l) = whts15(l)
        end do
        goto 1000
        end if


        if( deltas(1,16) .le. delta .AND. delta .le. deltas(2,16)) then
        nquad = 14
        do l=1,14
        xs(l)   = xs16(l)
        whts(l) = whts16(l)
        end do
        goto 1000
        end if


        if( deltas(1,17) .le. delta .AND. delta .le. deltas(2,17)) then
        nquad = 13
        do l=1,13
        xs(l)   = xs17(l)
        whts(l) = whts17(l)
        end do
        goto 1000
        end if


        if( deltas(1,18) .le. delta .AND. delta .le. deltas(2,18)) then
        nquad = 13
        do l=1,13
        xs(l)   = xs18(l)
        whts(l) = whts18(l)
        end do
        goto 1000
        end if


        if( deltas(1,19) .le. delta .AND. delta .le. deltas(2,19)) then
        nquad = 12
        do l=1,12
        xs(l)   = xs19(l)
        whts(l) = whts19(l)
        end do
        goto 1000
        end if


        if( deltas(1,20) .le. delta .AND. delta .le. deltas(2,20)) then
        nquad = 11
        do l=1,11
        xs(l)   = xs20(l)
        whts(l) = whts20(l)
        end do
        goto 1000
        end if


        if( deltas(1,21) .le. delta .AND. delta .le. deltas(2,21)) then
        nquad = 11
        do l=1,11
        xs(l)   = xs21(l)
        whts(l) = whts21(l)
        end do
        goto 1000
        end if


        if( deltas(1,22) .le. delta .AND. delta .le. deltas(2,22)) then
        nquad = 11
        do l=1,11
        xs(l)   = xs22(l)
        whts(l) = whts22(l)
        end do
        goto 1000
        end if


        if( deltas(1,23) .le. delta .AND. delta .le. deltas(2,23)) then
        nquad = 11
        do l=1,11
        xs(l)   = xs23(l)
        whts(l) = whts23(l)
        end do
        goto 1000
        end if

        ier = 4
        return

 1000 continue

        do i=1,nquad
        xs(i)   = (r1-r0)*xs(i)+r0
        whts(i) = (r1-r0)*whts(i)
        end do

        end
