     implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:), wts(:)
      integer *8 ipars(2)
      real *8 dpars(3), dpars_solver(2)
      
      integer *8, allocatable :: norders(:), ixyzs(:), iptype(:)

      real *8, allocatable :: dintvals(:,:,:,:)
      real *8, allocatable :: dintcoefs(:,:,:,:)
      
      real *8, allocatable :: targuse(:,:)
      real *8, allocatable :: amat(:,:), xs(:), ys(:), ws(:)
      real *8, allocatable :: srctmp(:,:), qwts(:), sigvals(:,:)
      real *8, allocatable :: fval(:,:)
      real *8, allocatable :: rhs(:,:), soln(:,:), soln2(:,:), errs(:)
      real *8, allocatable :: uvs(:,:), umat(:,:), vmat(:,:), wtmp(:)
      real *8, allocatable :: rtmp(:,:)
      real *8 strengths(3), pot(3,1000), pot_ex(3,1000)
      real *8, allocatable :: bmat(:,:), rhsb(:), solnb(:), resb(:)
      real *8, allocatable :: dmul1(:,:), dmul2(:,:)
      real *8, allocatable :: dhs(:), dpars_quad(:)

      procedure (), pointer :: fker


      real *8 c0(3), abc(3), verts(2,4), xyz_out(12), xyz_in(3,1000)
      integer *8 nabc(3)
      
      real *8, allocatable :: cms(:,:), rads(:), rad_near(:)
      integer *8, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      integer *8, allocatable :: col_ptr(:), row_ind(:), iper(:)
      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8, allocatable :: wnear(:,:,:), wnear2(:,:,:)

      complex *16 zpars
      real *8 uv(2)
      character *1 transa, transb
      external l3d_slp, st3d_slp_vec, st3d_comb_vec
      external el3d_elastlet_string_mindlin_normalstress_vec
      external el3d_elastlet_mindlin_normalstress_vec

      
      
      
      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4
      
      npatches = 0
      npts = 0

      abc(1) = 2.1d0
      abc(2) = 1.0d0
      abc(3) = 4.0d0

      abc(1) = 1.0d0
      abc(2) = 1.0d0
      abc(3) = 1.0d0

      abc(1:3) = abc(1:3)

      nabc(1) = 5
      nabc(2) = 3
      nabc(3) = 10

      na = 2
      huse = 0.5d0

      nabc(1) = na 
      nabc(2) = na
      nabc(3) = na
      
      c0(1) = 0
      c0(2) = 0
      c0(3) = 0

      xyz_out(1) = 3.1d0
      xyz_out(2) = 3.19d0
      xyz_out(3) = -3.19d0

      xyz_out(4:9) = 0
      xyz_out(10) = 1.0d0/sqrt(3.0d0)
      xyz_out(11) = 1.0d0/sqrt(3.0d0)
      xyz_out(12) = -1.0d0/sqrt(3.0d0)


      nin = 100
      do i=1,nin
        rr = hkrand(0)*0.25d0
        thet = hkrand(0)*pi
        phi = hkrand(0)*2*pi
        xyz_in(1,i) = rr*sin(thet)*cos(phi) 
        xyz_in(2,i) = rr*sin(thet)*sin(phi) 
        xyz_in(3,i) = rr*cos(thet)
      enddo

      allocate(bmat(3*nin,6), rhsb(3*nin), solnb(6), resb(3*nin))

      norder = 9
      iptype0 = 1

      

      call get_ellipsoid_npat_mem(abc, nabc, c0, norder, iptype0, &
        npatches, npts)

      allocate(srcvals(12,npts), srccoefs(9,npts))
      allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))
      allocate(wts(npts))
      
      call get_ellipsoid_npat(abc, nabc, c0, norder, iptype0, &
        npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)
 
      call prin2('srcvals=*',srcvals(1:3,1:10),30)
      call prin2('srcnorms=*',srcvals(10:12,1:10),30)
      
      
      call get_qwts(npatches, norders, ixyzs, iptype, npts, & 
        srcvals, wts)

      npols = ixyzs(2) - ixyzs(1)
      allocate(uvs(2,npols), umat(npols,npols),vmat(npols,npols))
      allocate(wtmp(npols))
      
      call get_disc_exps(norder, npols, iptype0, uvs, umat, vmat, wtmp)
      
      nv = 0
      call get_boundary_vertices(iptype0, verts, nv)
      call prin2('verts=*',verts,2*nv)

      istrat = 1
      intype = 1
      npatches0 = 1
      eps = 0.51d-12

      isd = 0
      ndsc = 9
      ndtarg = 12

      ifp = 0
      
      fker => el3d_elastlet_string_mindlin_normalstress_vec 
      ndim = 3
      nd = ndim*ndim
      

      ntarguse = npts - npols
      allocate(amat(ndim*npts, ndim*npts))
      allocate(targuse(12,ntarguse))
      allocate(fval(ndim,ndim))

      dlam = 1.0d0
      dmu = 1.0d0
      dpars(1) = dlam
      dpars(2) = dmu

      itargptr = 1

      ntrimax = 3000

      ndd = 3
      ndi = 0
      ndz = 0

      nqorder = 8
      call get_quadparams_adap(eps, iptype0, nqorder, eps_adap, nlev, &
        nqorder_f)
      call prin2('eps_adap=*',eps_adap, 1)
      call prinf('nlev=*', nlev, 1)
      if(nqorder.gt.20) nqorder = 20
      if(nqorder_f.gt.20) nqorder_f = 20
      print *, "nqorder=",nqorder
      rfac = 2.0d0
      
      ifmetric = 0
      
      nmax = 50000
      allocate(xs(nmax), ys(nmax), ws(nmax))
      allocate(srctmp(12,nmax), qwts(nmax))
      allocate(sigvals(npols,nmax))

      ipv = 0

      allocate(rhs(ndim,npts), soln(ndim,npts), soln2(ndim,npts))
      allocate(rtmp(ndim,ndim))

      strengths(1) = 1.1d0
      strengths(2) = 2.1d0
      strengths(3) = 0.3d0

      hout = huse 
      dpars(3) = hout

      ndim3 = 3
      ndd = 3
      ndz = 0
      ndi = 0
      
      do i=1,npts
        call el3d_elastlet_string_mindlin_normalstress_vec(nd, xyz_out, ndim3, &
          srcvals(1,i), ndd, dpars, ndz, zpars, ndi, ipars, rtmp)
        rhs(1,i) = rtmp(1,1)*strengths(1) + rtmp(1,2)*strengths(2) + &
                   rtmp(1,3)*strengths(3)
        rhs(2,i) = rtmp(2,1)*strengths(1) + rtmp(2,2)*strengths(2) + &
                   rtmp(2,3)*strengths(3)
        rhs(3,i) = rtmp(3,1)*strengths(1) + rtmp(3,2)*strengths(2) + &
                   rtmp(3,3)*strengths(3)
        soln(1:ndim,i) = 0
      enddo
!
!  find near
!
      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))
      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, cms, rads)
      
      rfac = 1.25d0
!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO      
      call prin2('rad_near=*', rad_near, npatches)

      allocate(dhs(npatches))
      allocate(dpars_quad(npatches+2))
      dpars_quad(1) = dlam
      dpars_quad(2) = dmu
      do i=1,npatches
        dhs(i) = rads(i)*2 
        dpars_quad(2+i) = dhs(i)
      enddo
      
!
!    find near quadrature correction interactions
!
      print *, "entering find near mem"
      ndtarg = 12
      call findnearmem(cms, npatches, rad_near, ndtarg, srcvals, npts, &
         nnz)
      

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms, npatches, rad_near, ndtarg, srcvals, npts, &
        row_ptr, col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind, &
         iquad)
      nquad = iquad(nnz+1) - 1
      
      allocate(row_ind(nnz), col_ptr(npatches+1), iper(nnz))
      call rsc_to_csc(npatches, npts, nnz, row_ptr, col_ind, col_ptr, &
        row_ind, iper)
      
      allocate(ipatch_id(npts), uvs_targ(2,npts))
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ) 
      allocate(wnear(ndim,ndim,nquad))
      print *, "Starting quadrature gen"
      rfac_tri = 2.0d0
      allocate(wnear2(ndim,ndim,nquad))
      iquadtype = 1
      call getnearquad_el_mstring_trans(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, dpars_quad, &
        iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear2)
 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ipatch, i,istart,j,targuse) &
!$OMP PRIVATE(h, dpars, rn1, n2, dintcoefs, ier, ipv, ns, xs, ys) &
!$OMP PRIVATE(ws, l, uv, sigvals, transa, transb, alpha, beta) &
!$OMP PRIVATE(lda, ldb, ldc, srctmp, rr, qwts, fval, k, isrc) &
!$OMP PRIVATE(dintvals, jtarg, jpatch, ntt, ii, iqstart, iind)
      do ipatch = 1,npatches
        istart = ixyzs(ipatch)

!  Extract non-self targets
!
        ntt = 0
        do i = col_ptr(ipatch), col_ptr(ipatch+1)-1
          jtarg = row_ind(i)
          jpatch = ipatch_id(jtarg)
          if(jpatch.ne.ipatch) then
            ntt = ntt + 1
            targuse(1:12,ntt) = srcvals(1:12,jtarg)
          endif
        enddo

!
!  For now set h to be constant
!
        dpars(1) = dlam 
        dpars(2) = dmu 
        dpars(3) = dhs(ipatch)
        rn1 = 0
        n2 = 0

        allocate(dintcoefs(ndim,ndim,npols,ntt+npols))
        allocate(dintvals(ndim,ndim,npols,ntt+npols))
        dintcoefs = 0

        call dtriaints_vec(eps, istrat, intype, npatches0, &
          norders(ipatch), npols, isd, ndsc, srccoefs(1,istart), ndtarg, &
          ntt, targuse, ifp, targuse, itargptr, ntt, &
          norders(ipatch), npols, fker, nd, ndd, dpars, ndz, zpars, ndi, &
          ipars, nqorder, ntrimax, rfac_tri, dintcoefs, ifmetric, rn1, n2)
!
!  Now start getting stuff for self
!
        ier = 0
        ipv = 0
        do j=1,npols
          ier = 0
          ns = 0
          call self_quadrature(norder, ipv, verts, nv, uvs(1,j), & 
            uvs(2,j), srcvals(4,istart+j-1), ns, xs, ys, ws, ier) 
          do l=1,ns
            uv(1) = xs(l)
            uv(2) = ys(l)
            call get_basis_pols(uv, norder, npols, iptype(ipatch), &
               sigvals(1,l))
          enddo
          transa = 'n'
          transb = 'n'
          alpha = 1.0d0
          beta = 0.0d0
          lda = ndsc
          ldb = npols
          ldc = 12
          call dgemm_guru(transa, transb, ndsc, ns, npols, alpha, &
             srccoefs(1,istart), lda, sigvals, ldb, beta, srctmp, ldc)
          do l=1,ns
            call cross_prod3d(srctmp(4,l), srctmp(7,l), srctmp(10,l))
            rr = sqrt(srctmp(10,l)**2 + srctmp(11,l)**2 + &
              srctmp(12,l)**2)
            qwts(l) = rr*ws(l)
            srctmp(10:12,l) = srctmp(10:12,l)/rr
            fval(1:ndim,1:ndim) = 0
             call fker(nd, srctmp(1,l), ndtarg, srcvals(1,istart+j-1), &
               ndd, dpars, ndz, zpars, ndi, ipars, fval)
            do k=1,npols
              dintcoefs(1:ndim,1:ndim,k,ntt+j) = & 
                dintcoefs(1:ndim,1:ndim,k,ntt+j) + & 
                fval(1:ndim,1:ndim)*sigvals(k,l)*qwts(l)
            enddo
          enddo
        enddo

!
!  convert coefs to vals
!
        do k=1,ntt+npols
          do j=1,npols
            dintvals(1:ndim,1:ndim,j,k) = 0
            do l=1,npols
              dintvals(1:ndim,1:ndim,j,k) = & 
                dintvals(1:ndim,1:ndim,j,k) + & 
                umat(l,j)*dintcoefs(1:ndim,1:ndim,l,k)
            enddo
          enddo
        enddo
!
!  start filling out relevant sections of wnear
!
        ii = 1
        do i=col_ptr(ipatch),col_ptr(ipatch+1)-1
          jtarg = row_ind(i)
          jpatch = ipatch_id(jtarg)
          iqstart = iquad(iper(i))-1
          if(ipatch.ne.jpatch) then
            do l=1,npols
              wnear(1:ndim,1:ndim,iqstart+l) = dintvals(1:ndim,1:ndim,l,ii)
            enddo
            ii = ii + 1
          else
            iind = 0
            rr = 0
            do l=1,npols
              rr = (uvs_targ(1,jtarg) - uvs(1,l))**2 + &
                   (uvs_targ(2,jtarg) - uvs(2,l))**2
              if(rr.le.1.0d-16) iind = l
            enddo
            do l=1,npols
              wnear(1:ndim,1:ndim,iqstart+l) = dintvals(1:ndim,1:ndim,l,iind+ntt)
            enddo
          endif
        enddo
      enddo 
!$OMP END PARALLEL DO
      print *, "done quadrature gen"


      errq = 0
      rq = 0
      do i=1,nquad
        do j=1,3
          do k=1,3
            errq = errq + abs(wnear(k,j,i) - wnear2(k,j,i))**2
            rq = rq + abs(wnear(k,j,i))**2
          enddo
        enddo
      enddo

      errq = sqrt(errq/rq)
      call prin2('wnear=*',wnear,60)
      call prin2('wnear2=*',wnear2,60)
      call prin2('errq=*',errq,1)

      do i=1,npts
        ipatch = ipatch_id(i)
        dpars(1) = dlam 
        dpars(2) = dmu 
        dpars(3) = dhs(ipatch)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j, rtmp)     
        do j=1,npts
          call fker(nd, srcvals(1,i), ndtarg, srcvals(1,j), &
            ndd, dpars, ndz, zpars, ndi, ipars, rtmp)
          amat(3*j-2,3*i-2) = rtmp(1,1)*wts(i)
          amat(3*j-1,3*i-2) = rtmp(2,1)*wts(i)
          amat(3*j-0,3*i-2) = rtmp(3,1)*wts(i)

          amat(3*j-2,3*i-1) = rtmp(1,2)*wts(i)
          amat(3*j-1,3*i-1) = rtmp(2,2)*wts(i)
          amat(3*j-0,3*i-1) = rtmp(3,2)*wts(i)

          amat(3*j-2,3*i-0) = rtmp(1,3)*wts(i)
          amat(3*j-1,3*i-0) = rtmp(2,3)*wts(i)
          amat(3*j-0,3*i-0) = rtmp(3,3)*wts(i)
        enddo
!$OMP END PARALLEL DO        
      enddo

!
!  Fix quadrature corrections
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,isrc)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            isrc = jstart+l-1
            amat(3*i-2,3*isrc-2) = wnear(1,1,jquadstart+l-1)
            amat(3*i-1,3*isrc-2) = wnear(2,1,jquadstart+l-1)
            amat(3*i-0,3*isrc-2) = wnear(3,1,jquadstart+l-1)

            amat(3*i-2,3*isrc-1) = wnear(1,2,jquadstart+l-1)
            amat(3*i-1,3*isrc-1) = wnear(2,2,jquadstart+l-1)
            amat(3*i-0,3*isrc-1) = wnear(3,2,jquadstart+l-1)

            amat(3*i-2,3*isrc-0) = wnear(1,3,jquadstart+l-1)
            amat(3*i-1,3*isrc-0) = wnear(2,3,jquadstart+l-1)
            amat(3*i-0,3*isrc-0) = wnear(3,3,jquadstart+l-1)
          enddo

        enddo
      enddo
!$OMP END PARALLEL DO

      if(1.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(isrc,j)
      do isrc = 1,npts
        do j=1,npts
!
!  fix null space issue, forces
!
            amat(3*j-2,3*isrc-2) = amat(3*j-2,3*isrc-2) + wts(isrc)
            amat(3*j-1,3*isrc-1) = amat(3*j-1,3*isrc-1) + wts(isrc)
            amat(3*j-0,3*isrc-0) = amat(3*j-0,3*isrc-0) + wts(isrc)
!
!  fix null space issue, torques
!
            amat(3*j-2,3*isrc-2) = amat(3*j-2,3*isrc-2) + &
                        (srcvals(2,j)*srcvals(2,isrc) + &
                        srcvals(3,j)*srcvals(3,isrc))*wts(isrc) 
            amat(3*j-2,3*isrc-1) = amat(3*j-2,3*isrc-1) - &
                        (srcvals(2,j)*srcvals(1,isrc))*wts(isrc) 
            amat(3*j-2,3*isrc-0) = amat(3*j-2,3*isrc-0) - &
                        (srcvals(3,j)*srcvals(1,isrc))*wts(isrc) 
            
            amat(3*j-1,3*isrc-2) = amat(3*j-1,3*isrc-2) - &
                        (srcvals(1,j)*srcvals(2,isrc))*wts(isrc) 
            amat(3*j-1,3*isrc-1) = amat(3*j-1,3*isrc-1) + &
                        (srcvals(1,j)*srcvals(1,isrc) + &
                        srcvals(3,j)*srcvals(3,isrc))*wts(isrc) 
            amat(3*j-1,3*isrc-0) = amat(3*j-1,3*isrc-0) - &
                        (srcvals(3,j)*srcvals(2,isrc))*wts(isrc) 

            amat(3*j-0,3*isrc-2) = amat(3*j-0,3*isrc-2) - &
                        (srcvals(1,j)*srcvals(3,isrc))*wts(isrc) 
            amat(3*j-0,3*isrc-1) = amat(3*j-0,3*isrc-1) - &
                        (srcvals(2,j)*srcvals(3,isrc))*wts(isrc) 
            amat(3*j-0,3*isrc-0) = amat(3*j-0,3*isrc-0) + &
                        (srcvals(1,j)*srcvals(1,isrc) + &
                        srcvals(2,j)*srcvals(2,isrc))*wts(isrc)
        enddo
      enddo
!$OMP END PARALLEL DO
      endif

      call prin2('amat=*', amat, 24)
      call prinf('ndim=*', ndim, 1)
      do i=1,ndim*npts
         amat(i,i) = amat(i,i) + 1.0d0
      enddo
!
!  test matmul
!
      allocate(dmul1(3,npts), dmul2(3,npts))
      call dmatvec(3*npts, 3*npts, amat, rhs, dmul1)
      print *, "done with matvec1"
      
      call matmul(npts, srcvals, wts, ipatch_id, &
        uvs_targ, dhs, dlam, dmu, ixyzs, &
        row_ptr, col_ind, iquad, wnear, rhs, dmul2)
      call prin2('dmul1=*',dmul1,24)
      call prin2('dmul2=*',dmul2,24)
      erra = 0
      ra = 0
      do i=1,npts
        erra = erra + abs(dmul1(1,i) - dmul2(1,i))**2*wts(i)
        erra = erra + abs(dmul1(2,i) - dmul2(2,i))**2*wts(i)
        erra = erra + abs(dmul1(3,i) - dmul2(3,i))**2*wts(i)

        ra = ra + abs(dmul1(1,i))**2*wts(i)
        ra = ra + abs(dmul1(2,i))**2*wts(i)
        ra = ra + abs(dmul1(3,i))**2*wts(i)
      enddo
      errmv = sqrt(erra/ra)
      read *, i

      info = 0
      call dgausselim(npts*ndim, amat, rhs, info, soln, dcond)

      print *, "dcond=", dcond
      print *, "info=", info

      call prin2('dpars=*',dpars,3)
      call prin2('xyz_out=*',xyz_out,12)

      do j=1,nin
        dpars(3) = hout
        rtmp(1:3,1:3) = 0
        call el3d_elastlet_string_mindlin_vec(nd, xyz_out, ndim3, &
            xyz_in(1,j), ndd, dpars, ndz, zpars, ndi, ipars, rtmp)
      
        pot_ex(1,j) = rtmp(1,1)*strengths(1) + rtmp(1,2)*strengths(2) + &
                  rtmp(1,3)*strengths(3)
        pot_ex(2,j) = rtmp(2,1)*strengths(1) + rtmp(2,2)*strengths(2) + &
                  rtmp(2,3)*strengths(3)
        pot_ex(3,j) = rtmp(3,1)*strengths(1) + rtmp(3,2)*strengths(2) + &
                  rtmp(3,3)*strengths(3)

        pot(1:3,j) = 0
        do i=1,npts
          rtmp(1:3,1:3) = 0 
          ipatch = ipatch_id(i)
          dpars(3) = dhs(ipatch)
          call el3d_elastlet_string_mindlin_vec(nd, srcvals(1,i), ndim3, &
            xyz_in(1,j), ndd, dpars, ndz, zpars, ndi, ipars, rtmp)
         
          pot(1,j) = pot(1,j) + (rtmp(1,1)*soln(1,i) + rtmp(1,2)*soln(2,i) + &
                    rtmp(1,3)*soln(3,i))*wts(i)
          pot(2,j) = pot(2,j) + (rtmp(2,1)*soln(1,i) + rtmp(2,2)*soln(2,i) + &
                    rtmp(2,3)*soln(3,i))*wts(i)
          pot(3,j) = pot(3,j) + (rtmp(3,1)*soln(1,i) + rtmp(3,2)*soln(2,i) + &
                    rtmp(3,3)*soln(3,i))*wts(i)
        enddo
      enddo
      pot(1:3,1:nin) = pot(1:3,1:nin)

      call prin2('pot=*',pot,12)
      call prin2('pot_ex=*',pot_ex,12)

      rarhs = 0
      rasoln = 0
      do i=1,npts
        rarhs = rarhs + abs(rhs(1,i))**2*wts(i)
        rarhs = rarhs + abs(rhs(2,i))**2*wts(i)
        rarhs = rarhs + abs(rhs(3,i))**2*wts(i)

        rasoln = rasoln + abs(soln(1,i))**2*wts(i)
        rasoln = rasoln + abs(soln(2,i))**2*wts(i)
        rasoln = rasoln + abs(soln(3,i))**2*wts(i)
      enddo

      rarhs = sqrt(rarhs)
      rasoln = sqrt(rasoln)

      call prin2('rarhs=*',rarhs,1)
      call prin2('rasoln=*',rasoln,1)



      erra = 0.0d0
      bmat(1:3*nin,1:6) = 0

      do i=1,nin
        bmat(3*i-2,1) = 1
        bmat(3*i-1,2) = 1
        bmat(3*i-0,3) = 1

        bmat(3*i-2,4) = 0
        bmat(3*i-2,5) = xyz_in(3,i)
        bmat(3*i-2,6) = -xyz_in(2,i)

        bmat(3*i-1,4) = -xyz_in(3,i)
        bmat(3*i-1,5) = 0 
        bmat(3*i-1,6) = xyz_in(1,i)

        bmat(3*i-0,4) = xyz_in(2,i)
        bmat(3*i-0,5) = -xyz_in(1,i)
        bmat(3*i-0,6) = 0

        rhsb(3*i-2) = pot_ex(1,i) - pot(1,i)
        rhsb(3*i-1) = pot_ex(2,i) - pot(2,i)
        rhsb(3*i-0) = pot_ex(3,i) - pot(3,i)
      enddo

      m = 3*nin
      n = 6
      nrhs = 1
      eps = 1.0d-12
      info = 0
      irank = 0
      call dleastsq(m, n, bmat, nrhs, rhsb, eps, info, solnb, irank)
      call prinf('irank=*', irank, 1)
      call prinf('info=*', info, 1)
      call prinf('m=*',m,1)
      call prinf('n=*',n,1)
      
      erra = 0
      do i=1,3*nin
        resb(i) = -rhsb(i)
        do j = 1,6
          resb(i) = resb(i) + bmat(i,j)*solnb(j)
        enddo

        erra = erra + resb(i)**2
      enddo
      erra = sqrt(erra)
      call prin2('resb=*',resb,12)
      call prin2('soln=*',solnb,6)
      print *, "Error in solution=", erra
      print *, "error in quadrature=",errq
      print *, "error in matvec=",errmv

      stop
      end



      subroutine getnearquad_el_mstring_trans(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, dpars, &
        iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature for the
!  representation:
!
!  u = K[\sigma]
!
!  where K is the mindlin Green's function with strings, and returns 
!  the corrections corresponding to normal stress associated with
!  this representation.
!
!  The quadrature is computed using adaptive integration
!
!  NOTES:
!    - wnear must be of size (3,3,nquad) 
!    - the parameters dpars must be dpars(2+npatches)
!        dpars(1) = \lambda Lame parameter
!        dpars(2) = \mu Lame parameter
!        dpars(3:end) = length of string on patch i
!
!  Input arguments:
!    - npatches: integer *8
!        number of patches
!    - norders: integer *8(npatches)
!        order of discretization on each patch
!    - ixyzs: integer *8(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer *8(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer *8
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch.
!        For each point
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - eps: real *8
!        precision requested
!    - dpars: complex *16 (2+npatches)
!        kernel parameters (See notes above)
!        dpars(1) = lambda Lame parameter
!        dpars(2) = mu Lame parameter 
!        dpars(3:end) = Length of strings on patch i 
!    - iquadtype: integer *8
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer *8
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer *8(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer *8 (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer *8(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: integer *8
!        radius parameter for near field (currently unused)
!    - nquad: integer *8
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is (3,3,nquad) since there are 9 kernels
!        per source target pair
!
!  Output arguments
!    - wnear: real *8(3,3,nquad)
!        The desired near field quadrature
!

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: npatches, norders(npatches), npts
      integer *8, intent(in) :: iptype(npatches), ixyzs(npatches+1)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, dpars(npatches+2)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(npatches+1), col_ind(nnz) 
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      real *8, intent(out) :: wnear(3,3,nquad)

      real *8 dpars_use(3)
      real *8, allocatable :: xs(:), ys(:), ws(:)
      real *8, allocatable :: srctmp(:,:), qwts(:), sigvals(:,:)

      integer *8, allocatable :: col_ptr(:), row_ind(:), iper(:)
      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)

      integer *8 ntt, i, j, k, l, isrc, ipatch, jpatch
      real *8, allocatable :: dintvals(:,:,:,:), dintcoefs(:,:,:,:)
      real *8, allocatable :: targuse(:,:)
      real *8 fval(3,3), uv(2), verts(2,4)
      character *1 transa, transb
      real *8, allocatable :: umat(:,:),vmat(:,:),uvs_r(:,:),wts_r(:)
      procedure (), pointer :: fker
      integer *8 ipars
      complex *16 zpars

      external el3d_elastlet_string_mindlin_normalstress_vec 
      
!
!  Get adaptive integration parameters
!
      iptype0 = iptype(1)
      nqorder = 8
      call get_quadparams_adap(eps, iptype0, nqorder, eps_adap, nlev, &
        nqorder_f)
      if(nqorder.gt.20) nqorder = 20
      if(nqorder_f.gt.20) nqorder_f = 20
      rfac = 2.0d0
      ntrimax = 3000
      
      ifmetric = 0
      
      nmax = 50000
      allocate(xs(nmax), ys(nmax), ws(nmax))
      allocate(srctmp(12,nmax), qwts(nmax))

      ipv = 0
!
!  transpose row_ptr and col_ind to get col_ptr and row_ind
!
      allocate(row_ind(nnz), col_ptr(npatches+1), iper(nnz))
      call rsc_to_csc(npatches, npts, nnz, row_ptr, col_ind, col_ptr, &
        row_ind, iper)
!
!  Get patch id and uvs targ
!
      allocate(ipatch_id(npts), uvs_targ(2,npts))
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)
      allocate(targuse(12,npts))

      istrat = 1
      intype = 1
      npatches0 = 1

      isd = 0
      ndsc = 9
      ndtarg = 12
      itargptr = 1

      ifp = 0
      
      fker => el3d_elastlet_string_mindlin_normalstress_vec 
      ndim = 3
      nd = ndim*ndim
      call prin2('dpars=*',dpars,24)
      rfac_tri = 2.0d0
      ndd = 3
      ndz = 0
      ndi = 0
      
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ipatch, i,istart,j,targuse) &
!$OMP PRIVATE(h, dpars_use, rn1, n2, dintcoefs, ier, ipv, ns, xs, ys) &
!$OMP PRIVATE(ws, l, uv, sigvals, transa, transb, alpha, beta) &
!$OMP PRIVATE(lda, ldb, ldc, srctmp, rr, qwts, fval, k, isrc) &
!$OMP PRIVATE(dintvals, jtarg, jpatch, ntt, ii, iqstart, iind) &
!$OMP PRIVATE(umat, vmat, uvs_r, wts_r, npols, verts, nv)
      do ipatch = 1,npatches
        istart = ixyzs(ipatch)
        npols = ixyzs(ipatch+1)-ixyzs(ipatch)
        allocate(umat(npols,npols), vmat(npols,npols))
        allocate(uvs_r(2,npols),wts_r(npols))
        allocate(sigvals(npols,nmax))
        nv = 0
        call get_boundary_vertices(iptype(ipatch), verts, nv)
        call get_disc_exps(norders(ipatch), npols, iptype(ipatch), &
               uvs_r, umat, vmat, wts_r)

!  Extract non-self targets
!
        ntt = 0
        do i = col_ptr(ipatch), col_ptr(ipatch+1)-1
          jtarg = row_ind(i)
          jpatch = ipatch_id(jtarg)
          if(jpatch.ne.ipatch) then
            ntt = ntt + 1
            targuse(1:12,ntt) = srcvals(1:12,jtarg)
          endif
        enddo

        dpars_use(1) = dpars(1)
        dpars_use(2) = dpars(2) 
        dpars_use(3) = dpars(ipatch+2)
        rn1 = 0
        n2 = 0

        allocate(dintcoefs(3,3,npols,ntt+npols))
        allocate(dintvals(3,3,npols,ntt+npols))
        dintcoefs = 0

        call dtriaints_vec(eps, istrat, intype, npatches0, &
          norders(ipatch), npols, isd, ndsc, srccoefs(1,istart), ndtarg, &
          ntt, targuse, ifp, targuse, itargptr, ntt, &
          norders(ipatch), npols, fker, nd, ndd, dpars_use, ndz, zpars, ndi, &
          ipars, nqorder, ntrimax, rfac_tri, dintcoefs, ifmetric, rn1, n2)
!
!  Now start getting stuff for self
!
        ier = 0
        ipv = 0
        do j=1,npols
          ier = 0
          ns = 0
          call self_quadrature(norders(ipatch), ipv, verts, nv, uvs_r(1,j), &
            uvs_r(2,j), srcvals(4,istart+j-1), ns, xs, ys, ws, ier) 
          do l=1,ns
            uv(1) = xs(l)
            uv(2) = ys(l)
            call get_basis_pols(uv, norders(ipatch), npols, iptype(ipatch), &
               sigvals(1,l))
          enddo
          transa = 'n'
          transb = 'n'
          alpha = 1.0d0
          beta = 0.0d0
          lda = ndsc
          ldb = npols
          ldc = 12
          call dgemm_guru(transa, transb, ndsc, ns, npols, alpha, &
             srccoefs(1,istart), lda, sigvals, ldb, beta, srctmp, ldc)
          do l=1,ns
            call cross_prod3d(srctmp(4,l), srctmp(7,l), srctmp(10,l))
            rr = sqrt(srctmp(10,l)**2 + srctmp(11,l)**2 + &
              srctmp(12,l)**2)
            qwts(l) = rr*ws(l)
            srctmp(10:12,l) = srctmp(10:12,l)/rr
            fval(1:3,1:3) = 0
             call fker(nd, srctmp(1,l), ndtarg, srcvals(1,istart+j-1), &
               ndd, dpars_use, ndz, zpars, ndi, ipars, fval)
            do k=1,npols
              dintcoefs(1:3,1:3,k,ntt+j) = & 
                dintcoefs(1:3,1:3,k,ntt+j) + & 
                fval(1:3,1:3)*sigvals(k,l)*qwts(l)
            enddo
          enddo
        enddo

!
!  convert coefs to vals
!
        do k=1,ntt+npols
          do j=1,npols
            dintvals(1:3,1:3,j,k) = 0
            do l=1,npols
              dintvals(1:3,1:3,j,k) = & 
                dintvals(1:3,1:3,j,k) + & 
                umat(l,j)*dintcoefs(1:3,1:3,l,k)
            enddo
          enddo
        enddo
!
!  start filling out relevant sections of wnear
!
        ii = 1
        do i=col_ptr(ipatch),col_ptr(ipatch+1)-1
          jtarg = row_ind(i)
          jpatch = ipatch_id(jtarg)
          iqstart = iquad(iper(i))-1
          if(ipatch.ne.jpatch) then
            do l=1,npols
              wnear(1:3,1:3,iqstart+l) = dintvals(1:3,1:3,l,ii)
            enddo
            ii = ii + 1
          else
            iind = 0
            rr = 0
            do l=1,npols
              rr = (uvs_targ(1,jtarg) - uvs_r(1,l))**2 + &
                   (uvs_targ(2,jtarg) - uvs_r(2,l))**2
              if(rr.le.1.0d-16) iind = l
            enddo
            do l=1,npols
              wnear(1:3,1:3,iqstart+l) = dintvals(1:3,1:3,l,iind+ntt)
            enddo
          endif
        enddo
        deallocate(dintvals, dintcoefs, umat, vmat, uvs_r, wts_r, sigvals)
      enddo 
!$OMP END PARALLEL DO




      return
      end
      
      

      subroutine matmul(npts, srcvals, wts, ipatch_id, &
        uvs_targ, dhs, dla, dmu, ixyzs, &
        row_ptr, col_ind, iquad, wnear, x, y)
      implicit real *8 (a-h,o-z)
      implicit integer *8(i-n)
      integer *8, intent(in) :: npts
      real *8, intent(in) :: srcvals(12,npts), wts(npts), dhs(*)
      integer *8, intent(in) :: ipatch_id(npts)
      real *8, intent(in) :: uvs_targ(2,npts)
      real *8, intent(in) :: dmu, dla
      integer *8, intent(in) :: ixyzs(*), row_ptr(*), col_ind(*)
      integer *8, intent(in) :: iquad(*)
      real *8, intent(in) :: wnear(3,3,*)
      real *8, intent(in) :: x(3,npts)
      real *8, intent(out) :: y(3,npts)
      
      real *8 rtmp(3,3), dpars(3)
      real *8 rf(3), rtorque(3), dtmp(3)
      integer *8 ipars
      complex *16 zpars

      procedure (), pointer :: fker
      external el3d_elastlet_string_mindlin_normalstress_vec

      fker => el3d_elastlet_string_mindlin_normalstress_vec 
!
!  the identity term
!
      dpars(1) = dla
      dpars(2) = dmu
      y(1:3,1:npts) = x(1:3,1:npts)
      nd = 9
      ndd = 3
      ndi = 0
      ndz = 0
      ndtarg = 12
      do i=1,npts
        ipatch = ipatch_id(i)
        dpars(3) = dhs(ipatch)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j, rtmp)        
        do j=1,npts
          if(j.ne.i) then
            call fker(nd, srcvals(1,i), ndtarg, srcvals(1,j), &
              ndd, dpars, ndz, zpars, ndi, ipars, rtmp)
            y(1:3,j) = y(1:3,j) + rtmp(1:3,1)*x(1,i)*wts(i)
            y(1:3,j) = y(1:3,j) + rtmp(1:3,2)*x(2,i)*wts(i)
            y(1:3,j) = y(1:3,j) + rtmp(1:3,3)*x(3,i)*wts(i)
          endif
        enddo
!$OMP END PARALLEL DO
      enddo

! Add and subrtract
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,isrc,rtmp,dpars)
      do i=1,npts
        dpars(1) = dla 
        dpars(2) = dmu
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          dpars(3) = dhs(jpatch)
          do l=1,npols
            isrc = jstart+l-1
            rtmp(1:3,1:3) = 0
            if(i.ne.isrc) then
              call fker(nd, srcvals(1,isrc), ndtarg, srcvals(1,i), &
              ndd, dpars, ndz, zpars, ndi, ipars, rtmp)
            endif
            y(1:3,i) = y(1:3,i) + (wnear(1:3,1,jquadstart+l-1) - &
                 rtmp(1:3,1)*wts(isrc))*x(1,isrc)
            y(1:3,i) = y(1:3,i) + (wnear(1:3,2,jquadstart+l-1) - &
                 rtmp(1:3,2)*wts(isrc))*x(2,isrc)
            y(1:3,i) = y(1:3,i) + (wnear(1:3,3,jquadstart+l-1) - &
                 rtmp(1:3,3)*wts(isrc))*x(3,isrc)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!
!  Ones matrix
!
      rf(1:3) = 0
      rtorque(1:3) = 0
      do i=1,npts
        rf(1:3) = rf(1:3) + x(1:3,i)*wts(i)
        call cross_prod3d(srcvals(1,i),x(1,i),dtmp)
        rtorque(1:3) = rtorque(1:3) + dtmp(1:3)*wts(i)
      enddo

      do i=1,npts
         y(1:3,i) = y(1:3,i) + rf(1:3)
         call cross_prod3d(srcvals(1,i), rtorque, dtmp)
         y(1:3,i) = y(1:3,i) - dtmp(1:3)
      enddo

      return
      end
      

