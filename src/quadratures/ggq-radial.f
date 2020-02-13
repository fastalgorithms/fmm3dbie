C
C (c) Jim Bremer, all rights reserved
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the 
c       quadrature code.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains subroutines for generating quadrature rules
c       for integrating a class of radially singular functions over
c       more or less arbitrary triangles containing the origin.
c
c       The quadrature rules are designed to integrate functions
c       which admit an expansion of the form
c
c                                N    (j-1)
c       f(r\cos(t),r\sin(t)) = \sum  r      p  (t) ,                       (1)
c                               j=0          3j+2
c
c       where p_{3j+2} denotes a trigonometric polynomial of order 3j+2,
c       over more or less arbitrary triangles containing the origin.
c       We say that the order of the quadrature rule (1) is N.
c
c       These rules are constructed on the fly using precomputed tables 
c       stored on the disk.  The possible orders of the quadrature 
c       rules are fixed at the time of precomputation.
c
c       Note also that the routines will fail for sufficiently 
c       degenerate triangles; the precise tolerances depend on the
c       limits of the parameters for the precomputed quadrature rules.
c
c       The quadrature rules are accurate to roughly 25 digits.
c
c       The following subroutines are user-callable:
c
c   radial_init - read a precomputed table of quadraturs from a 
c       text file on the disk.
c
c   raddiag - return a quadrature for integrating a function of the 
c        form (1) over a triangle containing the origin.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



        subroutine radial_init_mem(norder,lkeep)
        implicit real *8 (a-h,o-z)

        character*54 filename
        dimension rs(2,1000),ts(2,1000)

c
 0050 format ("../../src/quadratures/ggq-self-quads/radquads",I2.2,
     1   ".txt")
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
        open(iw,FILE=filename,STATUS='old',ERR=1000)



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
cccc        call prina("in radial_init, filename = *",filename,23)
cccccall prin2("in radial_init, rs = *",rs,2*nr)
cccc        call prin2("in radial_init, ts = *",ts,2*nt)
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
 1000 continue

        return
        end


        subroutine radial_init(ier,norder,rad,lrad,lkeep)
        implicit double precision (a-h,o-z)
        dimension rad(1)
        dimension rs(2,1000),ts(2,1000)
        character*2 str
        character*54 filename
c
c       Read a table of precomputed quadrature rules from a file
c       on the disk into the user-supplied array.  The table is
c       stored in the file ``radquads??.txt'' where ?? is the order
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
c                  attempting to read the input file
c
c   rad - on return, this user-supplied array will contain a quadrature
c       table and a structure header describing the quadratures
c
c
 0050 format ("../../src/quadratures/ggq-self-quads/radquads",I2.2,
     1   ".txt")
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
        open(iw,FILE=filename,STATUS='old',ERR=1000)


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
cccc        call prina("in radial_init, filename = *",filename,23)
cccccall prin2("in radial_init, rs = *",rs,2*nr)
cccc        call prin2("in radial_init, ts = *",ts,2*nt)
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
c       Construct the Legendre quadrature.
c
cccc        call legequad(nlege,rad(ixslege),rad(iwhtslege))
        itype111 = 1
        call legerts(itype111, nlege, rad(ixslege), rad(iwhtslege))

c
c       Copy the as and bs into the array.
c
        call radmove(nr*2,rs,rad(irs))
        call radmove(nt*2,ts,rad(its))
c
c       Call the auxillary routine to read the quadrature rules.
c
        call radinit0(ier,iw,max,nr,nt,rad(ins),rad(ixs),rad(iwhts),
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
        subroutine radinit0(ier,iw,max,nr,nt,ns,xs,whts,rs,ts,norder)
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
        return
      
        end



        subroutine raddiag(ier,rad,verts,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension rad(1),xs(1),ys(1),whts(1),verts(2,3)
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
c       Build quadratures for each subtriangle.
c
        call raddiag0(ier,rad,x1,y1,x2,y2,nquad1,xs(1),ys(1),whts(1),
     1    w,lw)
        if (ier .ne. 0) return
c
        call raddiag0(ier,rad,x1,y1,x3,y3,nquad2,xs(nquad1+1),
     1    ys(nquad1+1),whts(nquad1+1),w,lw)
        if (ier .ne. 0) return
c
        call raddiag0(ier,rad,x2,y2,x3,y3,nquad3,xs(nquad1+nquad2+1),
     1    ys(nquad1+nquad2+1),whts(nquad1+nquad2+1),w,lw)
        if (ier .ne. 0) return
c
        nquad=nquad1+nquad2+nquad3
c

        end
c
c
c
        subroutine raddiag0(ier,rad,x1_,y1_,x2_,y2_,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension rad(1),xs(1),ys(1),whts(1),verts(2,3)
        dimension amatr(2,2),ts(1000),twhts(1000)
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
c        
c
c       Fetch parameters from the rad structure.
c
        nlege     = rad(10) 
        ixslege   = rad(11) 
        iwhtslege = rad(12) 
c
c       Compute the radii of the vertices.
c
        rr1 = sqrt(x1_**2+y1_**2)
        rr2 = sqrt(x2_**2+y2_**2)
c      
c       Perform a sanity check ... make sure neither of the specified
c       vertices is 0.
c      
        if (rr1 .eq. 0 .OR. rr2 .eq. 0) then
        ier = 128
        return
        endif
c
c       Swap the two vertices in order to ensure (x2,y2) is the most
c       distant from the origin.
c      
        if (rr1 .gt. rr2) then
        x2 = x1_
        y2 = y1_
        x1 = x2_
        y1 = y2_
        else
        x1 = x1_
        y1 = y1_
        x2 = x2_
        y2 = y2_
        endif
c
c       Build a transformation matrix taking the vertex of largest
c       radius to (1,0).
c
        rr = (x2**2+y2**2)
c
        amatr(1,1) = x2/rr
        amatr(2,2) = x2/rr
        amatr(1,2) = y2/rr
        amatr(2,1) = -y2/rr
c
        u = amatr(1,1)*x1+amatr(1,2)*y1
        v = amatr(2,1)*x1+amatr(2,2)*y1
c
        if (v .lt. 0) then
        amatr(2,1)=-amatr(2,1)
        amatr(2,2)=-amatr(2,2)
        v=-v
        endif
c
c       Fetch the theta quadrature.
c
        r = sqrt(u**2+v**2)
        t = atan2(v,u)
c
        lw = 10 000 000
c
        call radfetch(ier,rad,r,t,nt,its,itwhts)
c
        if (ier .ne. 0) then
        call prin2("radfetch failed with r = *",r,1)
        call prin2("and t = *",t,1)
c
        print *,"r = ",r
        print *,"t = ",t
        stop
c
        ier = 64
        return
        endif
c
c       Call an auxillary routine to build the product quadrature.
c
        call raddiag1(nt,rad(its),rad(itwhts),
     1    nlege,rad(ixslege),rad(iwhtslege),
     2    amatr,nquad,xs,ys,whts,r,t)
c
        end

c
        subroutine raddiag1(nt,ts,twhts,nlege,xslege,whtslege,amatr,
     1     nquad,xs,ys,whts,a,b)
        implicit double precision (a-h,o-z)
        dimension ts(1),twhts(1),xslege(1),whtslege(1)
        dimension amatr(2,2),ainv(2,2)
        dimension xs(1),ys(1),whts(1)
c
        nquad=0
c
c        print *,"raddiag: ",a,b
c        print *,""
c
c       Build the product quadrature.
c
        do 1000 j=1,nt
        t    = ts(j)*b
        twht = twhts(j)*b
c
        r1 = 0
        r2 = a*sin(b)/(a*sin(b-t)+sin(t))
c
        alpha = (r2-r1)/2
        beta  = (r2+r1)/2
c
        do 1100 i=1,nlege
        r = xslege(i)*alpha+beta
        rwht = whtslege(i)*alpha
c
        x = r*cos(t)
        y = r*sin(t)
        wht = rwht*twht*r
c
        nquad=nquad+1
        xs(nquad)=x
        ys(nquad)=y
        whts(nquad)=wht
c
 1100 continue
 1000 continue
c
c       Transform the quadrature.
c
        call rad2x2inv(amatr,ainv,det)
        det=abs(1/det)
c
        do 2000 j=1,nquad
        x = xs(j)
        y = ys(j)
        wht = whts(j)
c
        u = ainv(1,1)*x+ainv(1,2)*y
        v = ainv(2,1)*x+ainv(2,2)*y
        wht = wht*det
c
        xs(j)=u
        ys(j)=v
        whts(j)=wht
 2000 continue
c
        end



        subroutine radfetch(ier,rad,r,t,n,ixs0,iwhts0)
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
        call radfetch0(ier,max,nr,nt,rad(irs),rad(its),rad(ins),
     1     rad(ixs),rad(iwhts),n,ii,r,t)
c
        ixs0   = ixs+ii-1
        iwhts0 = iwhts+ii-1
c
        end
c
c
c
        subroutine radfetch0(ier,max,nr,nt,rs,ts,ns,xs,xwhts,
     1    n,ii,r,t)
        implicit double precision (a-h,o-z)
        dimension rs(2,nr),ts(2,nt)
        dimension ns(nt,nr),xs(max,nt,nr),whts(max,nt,nr)
c
        data eps / 1.0d-12 /
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


        subroutine rayintersect(ier,x1,y1,x2,y2,th,x,y)
        implicit double precision (a-h,o-z)
c
c       Find the point of intersection of the ray \theta=th
c       and the line segment between two points --- presuming
c       there is such an intersection.  
c
c                          Input Parameters:
c
c   (x1,y1) - coordinates of one point forming the line segment
c   (x2,y2) - coordinates of the second point forming the line segment
c   th - real number in the interval [-pi,pi] specifiying the angle
c       of the ray
c               
c                         Output Parameters:
c
c   ier - an error return code;
c       ier = 0   indicates successful execution
c       ier = 16  means that there were no intersections
c       ier = 64  means that there are an infinite number of 
c                 intersections
c
c  (x,y) - coordinates of the intersection
c
        ier = 0
c        eps0 = 1.0d-15
c
c       Compute some initial parameters.
c
        ct = cos(th)
        st = sin(th)
        dx = x2-x1
        dy = y2-y1
c
c       If the line segments are parallel, there are either no
c       intersections or an infinite number, depending on the 
c       direction of the ray.
c
        det = -dx*st+dy*ct
c
        if (abs(det) .le. 0) then
        ier=16
        th1 = atan2(y1,x1)        
        if (th1*th .gt. 0) ier = 64
        return
        endif
c
c       If not, then solve the linear system to check for 
c       intersection.
c
        t = (st*x1-ct*y1)/det
        s = (dy*x1-dx*y1)/det
c
        ier = 0
        if (0 .le. t .AND. t. le. 1 .AND. s .gt. 0) then
        x=(1-t)*x1+t*x2
        y=(1-t)*y1+t*y2
        return
        endif
c
        ier = 16
c
        end
c
c
c
        subroutine rad2x2inv(amatr,ainv,det)
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
        subroutine radinsort2(k,a,b)
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
        subroutine radmove(k,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
        do 1000 j=1,k
        b(j)=a(j)
 1000 continue
        end
