      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcvals_tot(:,:),srccoefs_tot(:,:)
      real *8, allocatable :: coefs_pols(:,:,:),vals_pols(:,:,:)
      real *8, allocatable :: xyzvals(:,:,:), srccoefs(:,:,:)
      real *8, allocatable :: dxyzvals(:,:,:)
      real *8, allocatable :: cintvals(:,:)
      real *8, allocatable :: cintex(:), cinttest(:), cinttest2(:)
      real *8, allocatable :: umatr(:,:), vmatr(:,:)
      real *8 uvs(2,1000), wts(1000)
      real *8 xyztarg(3,2)
      real *8 uv_targ(2), uv(2), pols(1000)
      real *8 v1(3), v2(3), v3(3), v4(3)
      real *8 verts(3,3,2)
      real *8 dpars
      complex *16 zpars
      integer ipars

      integer itargptr(3), ntargptr(3)
      external lslp


      call prini(6,13)
      norder = 8
      npols = (norder+1)*(norder+2)/2
      allocate(srccoefs(9,npols,2), xyzvals(3,npols,2))
      allocate(dxyzvals(6,npols,2))
      allocate(umatr(npols,npols), vmatr(npols,npols))

      allocate(coefs_pols(npols,npols,2), vals_pols(npols,npols,2))

      call vioreanu_simplex_quad(norder, npols, uvs, umatr, vmatr, 
     1       wts)

      uv_targ(1) = 0.3d0
      uv_targ(2) = 0.0d0
      iedge = 1

      v1(1:3) = 0

      v2(1:3) = 0
      v2(1) = 1

      v3(1:3) = 0
      v3(2) = 1

      v4(1) = uv_targ(1)
      v4(2) = uv_targ(2)
      v4(3) = 0

      if (iedge.eq.1) then
c
c   vertices for triangle 1
c
         verts(1:3,1,1) = v1
         verts(1:3,2,1) = v4 
         verts(1:3,3,1) = v3
c
c  vertices for triangle 2
c
c
         verts(1:3,1,2) = v4
         verts(1:3,2,2) = v2 
         verts(1:3,3,2) = v3
      elseif (iedge.eq.2) then
         verts(1:3,1,1) = v1
         verts(1:3,2,1) = v2
         verts(1:3,3,1) = v4

         verts(1:3,1,2) = v1
         verts(1:3,2,2) = v4
         verts(1:3,3,2) = v3
      else
         verts(1:3,1,1) = v1
         verts(1:3,2,1) = v2
         verts(1:3,3,1) = v4

         verts(1:3,1,2) = v2
         verts(1:3,2,2) = v3
         verts(1:3,3,2) = v4
      endif

      npatches = 2
      do i=1,npatches
        do j=1,npols
           xyzvals(1:3,j,i) = verts(1:3,1,i) + 
     1        uvs(1,j)*(verts(1:3,2,i) - verts(1:3,1,i)) + 
     2        uvs(2,j)*(verts(1:3,3,i) - verts(1:3,1,i))
           uv(1) = xyzvals(1,j,i)
           uv(2) = xyzvals(2,j,i)
           call koorn_pols(uv, norder, npols, pols)
           do l=1,npols
             vals_pols(j,l,i) = pols(l)
           enddo

           dxyzvals(1:3,j,i) = verts(1:3,2,i) - verts(1:3,1,i)
           dxyzvals(4:6,j,i) = verts(1:3,3,i) - verts(1:3,1,i)
        enddo
      enddo

      npts = npatches*npols
      do ipatch =1,npatches
        do i=1,npols
          do j=1,3
            srccoefs(j,i,ipatch) = 0
            do l=1,npols
              srccoefs(j,i,ipatch) = srccoefs(j,i,ipatch)+umatr(i,l)*
     1                                 xyzvals(j,l,ipatch)
            enddo
          enddo

          do j=1,6
            srccoefs(j+3,i,ipatch) = 0
            do l=1,npols
              srccoefs(j+3,i,ipatch)=srccoefs(j+3,i,ipatch)+umatr(i,l)*
     1                                   dxyzvals(j,l,ipatch)
            enddo
          enddo

          do j=1,npols
            coefs_pols(i,j,ipatch) = 0
            do l=1,npols
              coefs_pols(i,j,ipatch) = coefs_pols(i,j,ipatch) +
     1                          umatr(i,l)*vals_pols(l,j,ipatch)      
            enddo
          enddo
        enddo
      enddo
      call prin2('coefs_pols 1=*',coefs_pols(1,3,1),npols)
      call prin2('coefs_pols 2=*',coefs_pols(1,7,2),npols)

      ntarg = 2
      xyztarg(1:3,1) = v4 
      xyztarg(1:3,2) = v4 
      
      itargptr(1) = 1
      itargptr(2) = 2
     
      ntargptr(1) = 1
      ntargptr(2) = 1

      ntrimax = 3000
      istrat = 2
      intype = 1
      ifp = 0
      ifmetric = 0

      nqorder = 16

      eps = 1.0d-9

      allocate(cintvals(npols,ntarg), cintex(npols))
      allocate(cinttest(npols), cinttest2(npols))

      ndd = 0 
      ndz = 0
      ndi = 0
      
      call dtriaints(eps, istrat, intype, npatches, norder, npols,
     1  srccoefs, 3, ntarg, xyztarg, ifp, tmp, itargptr, ntargptr,
     2  norder, npols, lslp, ndd, dpars, ndz, zpars, ndi, ipars,
     3  nqorder, ntrimax, rfac, cintvals, ifmetric, rn1, n2)
      
      do i=1,npols
        cinttest(i) = 0
        do j=1,npols
          cinttest(i) = cinttest(i) + 
     1       coefs_pols(j,i,1)*cintvals(j,1)
          cinttest(i) = cinttest(i) + 
     1       coefs_pols(j,i,2)*cintvals(j,2)
        enddo
      enddo

      call prin2('cinttest=*',cinttest,npols)


      allocate(srccoefs_tot(9,npols),srcvals_tot(9,npols))
      do i=1,npols
        srcvals_tot(1,i) = uvs(1,i)
        srcvals_tot(2,i) = uvs(2,i)
        srcvals_tot(3,i) = 0

        srcvals_tot(4,i) = 1
        srcvals_tot(5,i) = 0
        srcvals_tot(6,i) = 0

        srcvals_tot(7,i) = 0
        srcvals_tot(8,i) = 1
        srcvals_tot(9,i) = 0
      enddo

      do i=1,npols
        do j=1,9
          srccoefs_tot(j,i) = 0
          do l=1,npols
            srccoefs_tot(j,i) = srccoefs_tot(j,i) + 
     1         umatr(i,l)*srcvals_tot(j,l)
          enddo
        enddo
      enddo


      
      ipv = 0

      call dget_ggq_self_quad_pt(ipv, norder, npols, uv_targ, 1,
     1  umatr, srccoefs_tot, 3, xyztarg, 2, lslp, ndd, dpars, 
     2  ndz, zpars, ndi, ipars, cintex)
      do i=1,npols 
        cinttest2(i) = 0
      enddo

      do j=1,npols
        call koorn_pols(uvs(1,j), norder, npols, pols)
        do i=1,npols
          cinttest2(i) = cinttest2(i) + pols(i)*cintex(j)
        enddo
      enddo

      call prin2('cinttest2=*',cinttest2,npols)
      


      erra = 0
      ra = 0
      do i=1,npols
        erra = erra + abs(cinttest(i) - cinttest2(i))**2
        ra = ra + abs(cinttest(i))**2
      enddo
      erra = sqrt(erra/ra)
      call prin2('relative error=*',erra,1)

      stop
      end



      subroutine lslp(x,ndt,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(*),y(ndt),dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 f
      
      rr = (x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2
      f = 1/sqrt(rr)

      return
      end
      

      
      
