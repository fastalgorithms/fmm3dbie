!
!
!  TODO: fix speed of determining clashes
!
!
      subroutine get_npols(iptype,norder,npols)
!
!
!  This subroutine determines the number of polynomials 
!  for a given order corresponding to patch type
!
!  Input arguments
!    - iptype: integer *8
!        type of patch
!        * iptype = 1, triangular patch with RV nodes
!        * iptype = 11, quad patch with GL nodes and full order
!        * iptype = 12, quad patch with cheb nodes and full order
!
!    - norder: integer *8
!        order of discretization
!  Output arguments
!    - npols: integer *8
!        number of nodes/basis functions in the discretization
!       
      implicit real *8(a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 iptype,norder,npols

      if(iptype.eq.1) then
        npols = (norder+1)*(norder+2)/2
      elseif (iptype.eq.11.or.iptype.eq.12) then
        npols = (norder+1)*(norder+1)
      else
        print *, "Invalid order in get_npols"
        print *, "returning with npols=0"
        npols = 0
      endif


      return
      end
!
!
!
!
      

      subroutine get_disc_nodes(norder,npols,iptype,uvs)
!
!  Given a patch type and order, this subroutine gets the corresponding
!  nodes and weights
!
!  Input arguments:
!    - norder: integer *8
!        order of discretization
!    - npols: integer *8
!        number of nodes/basis functions in the discretization
!    - iptype: integer *8
!        type of patch
!        * iptype = 1, triangular patch with RV nodes
!        * iptype = 11, quad patch with GL nodes and full order
!        * iptype = 12, quad patch with cheb nodes and full order
!  Output arguments:
!    - uvs:  real *8(2,npols)
!        discretization nodes
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: norder,npols,iptype
      real *8, intent(out) :: uvs(2,npols)
!
!   temporary variable
!
      integer *8 ipoly,itype
      character * 1, ttype
      real *8 umatr,vmatr
      integer *8 ldu,ldv
      real *8 wts(1)
      

      itype = 0
      ldu = 1
      ldv = 1

      if(iptype.eq.1) then
        call get_vioreanu_nodes(norder,npols,uvs)
      endif

      if(iptype.eq.11) ipoly = 0
      if(iptype.eq.12) ipoly = 1

      ttype = "F"

      if(iptype.eq.11.or.iptype.eq.12) then
        call polytens_exps_2d(ipoly,itype,norder+1,ttype,uvs,umatr,ldu, &
          vmatr,ldv,wts)
      endif

     
      return
      end
!
!
!
!
!
!
!

      subroutine get_disc_wts(norder,npols,iptype,wts)
!
!  Given a patch type and order, this subroutine gets the corresponding
!  nodes and weights
!
!  Input arguments:
!    - norder: integer *8
!        order of discretization
!    - npols: integer *8
!        number of nodes/basis functions in the discretization
!    - iptype: integer *8
!        type of patch
!        * iptype = 1, triangular patch with RV nodes
!        * iptype = 11, quad patch with GL nodes and full order
!        * iptype = 12, quad patch with cheb nodes and full order
!  Output arguments:
!    - wts:  real *8(npols)
!        discretization weights
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: norder,npols,iptype
      real *8, intent(out) :: wts(2,npols)
!
!   temporary variable
!
      integer *8 ipoly,itype
      character * 1, ttype
      real *8 umatr,vmatr
      integer *8 ldu,ldv
      real *8 uvs(2,npols)

      

      itype = 1
      ldu = 1
      ldv = 1

      if(iptype.eq.1) then
        call get_vioreanu_wts(norder,npols,wts)
      endif

      if(iptype.eq.11) ipoly = 0
      if(iptype.eq.12) ipoly = 1

      ttype = "F"

      if(iptype.eq.11.or.iptype.eq.12) then
        call polytens_exps_2d(ipoly,itype,norder+1,ttype,uvs,umatr,ldu, &
          vmatr,ldv,wts)
      endif

     
      return
      end
!
!
!
!

      subroutine get_disc_nodes_wts(norder,npols,iptype,uvs,wts)
!
!  Given a patch type and order, this subroutine gets the corresponding
!  nodes and weights
!
!  Input arguments:
!    - norder: integer *8
!        order of discretization
!    - npols: integer *8
!        number of nodes/basis functions in the discretization
!    - iptype: integer *8
!        type of patch
!        * iptype = 1, triangular patch with RV nodes
!        * iptype = 11, quad patch with GL nodes and full order
!        * iptype = 12, quad patch with cheb nodes and full order
!  Output arguments:
!    - uvs:  real *8(2,npols)
!        discretization nodes
!    - wts: real *8(npols)
!        the corresponding quadrature weights
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: norder,npols,iptype
      real *8, intent(out) :: uvs(2,npols),wts(npols)
!
!   temporary variable
!
      integer *8 ipoly,itype
      character * 1, ttype
      real *8 umatr,vmatr
      integer *8 ldu,ldv
      

      itype = 1
      ldu = 1
      ldv = 1

      if(iptype.eq.1) then
        call get_vioreanu_nodes_wts(norder,npols,uvs,wts)
      endif

      if(iptype.eq.11) ipoly = 0
      if(iptype.eq.12) ipoly = 1

      ttype = "F"

      if(iptype.eq.11.or.iptype.eq.12) then
        call polytens_exps_2d(ipoly,itype,norder+1,ttype,uvs,umatr,ldu, &
          vmatr,ldv,wts)
      endif

     
      return
      end
!
!
!
!
!
!
!
!

      subroutine get_disc_exps(norder,npols,iptype,uvs,umatr,vmatr,wts)
!
!  Given a patch type and order, this subroutine gets the corresponding
!  nodes and weights
!
!  Input arguments:
!    - norder: integer *8
!        order of discretization
!    - npols: integer *8
!        number of nodes/basis functions in the discretization
!    - iptype: integer *8
!        type of patch
!        * iptype = 1, triangular patch with RV nodes
!        * iptype = 11, quad patch with GL nodes and full order
!        * iptype = 12, quad patch with cheb nodes and full order
!  Output arguments:
!    - uvs:  real *8(2,npols)
!        discretization nodes
!    - umatr: real *8(npols,npols)
!        values to coefs matrix
!    - vmatr: real *8(npols,npols)
!        coefs to values matrix
!    - wts: real *8(npols)
!        the corresponding quadrature weights
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: norder,npols,iptype
      real *8, intent(out) :: uvs(2,npols),wts(npols)
      real *8, intent(out) :: umatr(npols,npols),vmatr(npols,npols)
!
!   temporary variable
!
      integer *8 ipoly,itype
      character * 1, ttype
      integer *8 ldu,ldv
      

      itype = 2
      ldu = npols
      ldv = npols

      if(iptype.eq.1) then
        call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
      endif

      if(iptype.eq.11) ipoly = 0
      if(iptype.eq.12) ipoly = 1

      ttype = "F"


      if(iptype.eq.11.or.iptype.eq.12) then
        call polytens_exps_2d(ipoly,itype,norder+1,ttype,uvs,umatr,ldu, &
          vmatr,ldv,wts)
      endif

     
      return
      end
!
!
!
!
      subroutine get_basis_pols(uvs,norder,npols,iptype,pols)
!
!  This subroutine evaluates the basis functions corresponding
!  to patch type at a particular point on the patch
!
!  Input arguments:
!    - uvs: real *8(2)
!        location on patch where basis functions are desired
!    - norder: integer *8
!        order of discretization
!    - npols: integer *8
!        number of nodes/basis functions in the discretization
!    - iptype: integer *8
!        type of patch
!        * iptype = 1, triangular patch with RV nodes
!        * iptype = 11, quad patch with GL nodes and full order
!        * iptype = 12, quad patch with cheb nodes and full order
!  Output arguments:
!    - pols:  real *8(npols)
!        the evaluated polynomials
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: norder,npols,iptype
      real *8, intent(in) :: uvs(2)
      real *8, intent(out) :: pols(npols)

      integer *8 ipoly
      character *1 ttype

      if(iptype.eq.11) ipoly = 0
      if(iptype.eq.12) ipoly = 1

      ttype = "F"
      
      if(iptype.eq.1) then
        call koorn_pols(uvs,norder,npols,pols)
      elseif(iptype.eq.11.or.iptype.eq.12) then
        call polytens_pols_2d(ipoly,uvs,norder,ttype,pols) 
      endif

      return
      end
!
!
!
!
!
!
      subroutine get_basis_ders(uvs,norder,npols,iptype,pols,ders)
!
!  This subroutine evaluates the basis functions corresponding
!  to patch type at a particular point on the patch
!
!  Input arguments:
!    - uvs: real *8(2)
!        location on patch where basis functions are desired
!    - norder: integer *8
!        order of discretization
!    - npols: integer *8
!        number of nodes/basis functions in the discretization
!    - iptype: integer *8
!        type of patch
!        * iptype = 1, triangular patch with RV nodes
!        * iptype = 11, quad patch with GL nodes and full order
!        * iptype = 12, quad patch with cheb nodes and full order
!  Output arguments:
!    - pols:  real *8(npols)
!        the evaluated polynomials
!    - pols:  real *8(2,npols)
!        the evaluated derivatives
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: norder,npols,iptype
      real *8, intent(in) :: uvs(2)
      real *8, intent(out) :: pols(npols),ders(2,npols)

      integer *8 ipoly
      character *1 ttype

      if(iptype.eq.11) ipoly = 0
      if(iptype.eq.12) ipoly = 1

      ttype = "F"
      
      if(iptype.eq.1) then
        call koorn_ders(uvs,norder,npols,pols,ders)
      elseif(iptype.eq.11.or.iptype.eq.12) then
        call polytens_ders_2d(ipoly,uvs,norder,ttype,pols,ders) 
      endif

      return
      end
!
!
!
!
!
      subroutine get_tail_coefs(norder,npols,iptype,itailcoefs,ntailcoefs)
!----------------
!  This subroutine estimates the set of coefficients that define
!  the tail of the expansion
!
!  Input arguments:
!    - norder: integer *8
!        order of discretization
!    - npols: integer *8
!        number of nodes/basis functions in the discretization
!    - iptype: integer *8
!        type of patch
!        * iptype = 1, triangular patch with RV nodes
!        * iptype = 11, quad patch with GL nodes and full order
!        * iptype = 12, quad patch with cheb nodes and full order
!
!  Output arguments:
!    - itailcoefs: integer *8(ntailcoefs)
!        set of coefs in 1:npols which correspond to the tail
!        of the expansion
!    - ntailcoefs: integer *8
!        number of tail coefs
!----------------------        
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: norder,npols
      integer *8, intent(out) :: itailcoefs(npols),ntailcoefs

      integer *8 iind2p(2,npols)

      if(iptype.eq.1) then
        norderhead = norder-1
        i1 = (norderhead+1)*(norderhead+2)/2
        ntailcoefs = 0
        do i=i1+1,npols
          ntailcoefs = ntailcoefs +1
          itailcoefs(ntailcoefs) = i
        enddo
      endif

      if(iptype.eq.11.or.iptype.eq.12) then
        ntailcoefs = 0
        call polytens_ind2pow_2d(norder,"F",iind2p)
        do i=1,npols
          if(iind2p(1,i).eq.norder+1.or.iind2p(2,i).eq.norder+1) then
            ntailcoefs = ntailcoefs+1
            itailcoefs(ntailcoefs) = i
          endif
        enddo
      endif




      

      return
      end
!
!
!
!
!
      subroutine get_boundary_vertices(iptype,uv,nv)
!
!  This subroutine returns the number and collection of boundary
!  vertices for a given patch type
!
!  Input arguments:
!    - iptype: integer *8
!        type of patch
!        * iptype = 1, triangular patch with RV nodes
!        * iptype = 11, quad patch with GL nodes and full order
!        * iptype = 12, quad patch with cheb nodes and full order
!
!  Output arguments:
!    - uv: real *8 (2,nv)
!        set of boundary vertices
!    - nv: integer *8
!        number of boundary vertices
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 iptype,nv
      real *8 uv(2,*)

      if(iptype.eq.1) then
        nv = 3
        uv(1,1) = 0
        uv(2,1) = 0 

        uv(1,2) = 1
        uv(2,2) = 0

        uv(1,3) = 0
        uv(2,3) = 1
      endif

      if(iptype.eq.11.or.iptype.eq.12) then
        nv = 4
        uv(1,1) = -1
        uv(2,1) = -1

        uv(1,2) = 1
        uv(2,2) = -1

        uv(1,3) = 1
        uv(2,3) = 1

        uv(1,4) = -1
        uv(2,4) = 1
      endif

      return
      end
!
!
!
!
!
!
      subroutine get_clashing_indices(norder,npols,iptype,nfar,nfar_pols, &
        iclash,iclashfar,nclash)
!
!  This subroutine estimates the clashing indicies between two discretizations
!  on a patch, in particular given two sets of nodes u1(i),u2(i), and v1(j),v2(j)
!  this subroutine identifices all pairs of indices (i,j) such that
!  abs(u1(i)-v1(j))**2 + abs(u2(i)-v2(j))**2 < 1e-24
!  
!
!
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: norder,npols,iptype,nfar
      integer *8 iclash(npols),iclashfar(npols)
      integer *8, intent(out) :: nclash
      real *8 us(2,npols),ufar(2,nfar_pols)

      tol = 1.0d-12

      if(iptype.eq.1) then
        n1 = mod(norder,3)
        n2 = mod(nfar,3)
        if(n1.eq.0.and.n2.eq.0) then
          nclash = 1
          iclash(nclash) = 1
          iclashfar(nclash) = 1
          if(norder.eq.15) iclash(nclash) = 76
          if(norder.eq.18) iclash(nclash) = 106

          if(nfar.eq.15) iclashfar(nclash) = 76
          if(nfar.eq.18) iclashfar(nclash) = 106
        else
          nclash = 0
        endif
      elseif(iptype.eq.11.or.iptype.eq.12) then
        call get_disc_nodes(norder,npols,iptype,us) 
        call get_disc_nodes(nfar,nfar_pols,iptype,ufar)

        nclash = 0
        do i=1,npols
          do j=1,nfar_pols
            dd = (us(1,i)-ufar(1,j))**2 + (us(2,i)-ufar(2,j))**2
            if(dd.le.tol**2) then
              nclash = nclash + 1
              iclash(nclash) = i
              iclashfar(nclash) = j
            endif
          enddo
        enddo
      endif
      
      return
      end

