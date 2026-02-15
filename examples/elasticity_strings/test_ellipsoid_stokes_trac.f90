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

      procedure (), pointer :: fker


      real *8 c0(3), abc(3), verts(2,4), xyz_out(12), xyz_in(3,1000)
      integer *8 nabc(3)

      complex *16 zpars
      real *8 uv(2)
      character *1 transa, transb
      external l3d_slp, st3d_slp_vec, st3d_comb_vec
      external el3d_elastlet_string_mindlin_normalstress_vec
      external el3d_elastlet_mindlin_normalstress_vec
      external st3d_strac_vec
      
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

      na = 3
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


      nin = 10
      do i=1,nin
        rr = hkrand(0)*0.25d0
        thet = hkrand(0)*pi
        phi = hkrand(0)*2*pi
        xyz_in(1,i) = rr*sin(thet)*cos(phi) 
        xyz_in(2,i) = rr*sin(thet)*sin(phi) 
        xyz_in(3,i) = rr*cos(thet)
      enddo

      allocate(bmat(3*nin,6), rhsb(3*nin), solnb(6), resb(3*nin))

      norder = 7
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
     
      fker => st3d_strac_vec
      ndim = 3
      nd = ndim*ndim
      

      ntarguse = npts - npols
      allocate(dintcoefs(ndim,ndim,npols,npts))
      allocate(dintvals(ndim,ndim,npols,npts))
      allocate(amat(ndim*npts, ndim*npts))
      allocate(targuse(12,ntarguse))
      allocate(fval(ndim,ndim))

      dlam = 1.0d0
      dmu = 1.0d0
      dpars(1) = dlam
      dpars(2) = dmu

      itargptr = 1
      ntargptr = ntarguse

      ntrimax = 3000

      ndd = 3
      ndi = 0
      ndz = 0

      nqorder = 8
      call get_quadparams_adap(eps, iptype0, nqorder, eps_adap, nlev, &
        nqorder_f)
      call prin2('eps_adap=*',eps_adap, 1)
      call prinf('nlev=*', nlev, 2)
      if(nqorder.gt.20) nqorder = 20
      if(nqorder_f.gt.20) nqorder_f = 20
      print *, "nqorder=",nqorder
      rfac = 2.0d0
      
      ifmetric = 0
      
      nmax = 50000
      allocate(xs(nmax), ys(nmax), ws(nmax))
      allocate(srctmp(12,nmax), qwts(nmax))
      allocate(sigvals(npols,nmax))

      ipv = 1

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
        call st3d_slp_vec(nd, xyz_out, ndim3, &
          srcvals(1,i), ndd, dpars, ndz, zpars, ndi, ipars, rtmp)
        rhs(1,i) = rtmp(1,1)*strengths(1) + rtmp(1,2)*strengths(2) + &
                   rtmp(1,3)*strengths(3)
        rhs(2,i) = rtmp(2,1)*strengths(1) + rtmp(2,2)*strengths(2) + &
                   rtmp(2,3)*strengths(3)
        rhs(3,i) = rtmp(3,1)*strengths(1) + rtmp(3,2)*strengths(2) + &
                   rtmp(3,3)*strengths(3)
        soln(1:ndim,i) = 0
      enddo

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,istart,j,targuse) &
!$OMP PRIVATE(h, dpars, rn1, n2, dintcoefs, ier, ipv, ns, xs, ys) &
!$OMP PRIVATE(ws, l, uv, sigvals, transa, transb, alpha, beta) &
!$OMP PRIVATE(lda, ldb, ldc, srctmp, rr, qwts, fval, k, isrc) &
!$OMP PRIVATE(dintvals)
      do i = 1,npatches
        istart = ixyzs(i)

!  Extract non-self targets
!
        do j=1,ixyzs(i)-1
          targuse(1:12,j) = srcvals(1:12,j)
        enddo

        do j=ixyzs(i+1),npts
          targuse(1:12,j-npols) = srcvals(1:12,j)
        enddo

!        call prin2('targuse=*',targuse,12*ntarguse)
!
!  For now set h to be constant
!
        h = huse
        dpars(1) = dlam 
        dpars(2) = dmu 
        dpars(3) = h
        rn1 = 0
        n2 = 0
!
!        call dtriaints(eps, istrat, intype, npatches0, &
!          norders(i), npols, isd, ndsc, srccoefs(1,istart), ndtarg, &
!          ntarguse, targuse, ifp, targuse, itargptr, ntargptr, &
!          norders(i), npols, fker, ndd, dpars, ndz, zpars, ndi, &
!          ipars, nqorder, ntrimax, rfac, dintcoefs, ifmetric, rn1, n2)
        call dtriaints_vec(eps, istrat, intype, npatches0, &
          norders(i), npols, isd, ndsc, srccoefs(1,istart), ndtarg, &
          ntarguse, targuse, ifp, targuse, itargptr, ntargptr, &
          norders(i), npols, fker, nd, ndd, dpars, ndz, zpars, ndi, &
          ipars, nqorder, ntrimax, rfac, dintcoefs, ifmetric, rn1, n2)
        

!
!  Move dintcoefs to make room for self
!

        do j=npts,ixyzs(i+1),-1
          dintcoefs(1:ndim,1:ndim,1:npols,j) = &
            dintcoefs(1:ndim,1:ndim,1:npols,j-npols)
        enddo
!
!  Reset values for the self part
!
        do j=ixyzs(i),ixyzs(i+1)-1
          dintcoefs(1:ndim,1:ndim,1:npols,j) = 0
        enddo
!
!  Now start getting stuff for self
!
        ier = 0
        ipv = 1
        do j=1,npols
          ier = 0
          ns = 0
          call self_quadrature(norder, ipv, verts, nv, uvs(1,j), & 
            uvs(2,j), srcvals(4,ixyzs(i)+j-1), ns, xs, ys, ws, ier) 
          do l=1,ns
            uv(1) = xs(l)
            uv(2) = ys(l)
            call get_basis_pols(uv, norder, npols, iptype(i), &
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
             srccoefs(1,ixyzs(i)), lda, sigvals, ldb, beta, srctmp, ldc)
          do l=1,ns
            call cross_prod3d(srctmp(4,l), srctmp(7,l), srctmp(10,l))
            rr = sqrt(srctmp(10,l)**2 + srctmp(11,l)**2 + &
              srctmp(12,l)**2)
            qwts(l) = rr*ws(l)
            srctmp(10:12,l) = srctmp(10:12,l)/rr
            fval(1:ndim,1:ndim) = 0
!            call fker(srctmp(1,l), ndtarg, srcvals(1,ixyzs(i)+j-1), &
!              ndd, dpars, ndz, zpars, ndi, ipars, fval)
             call fker(nd, srctmp(1,l), ndtarg, srcvals(1,ixyzs(i)+j-1), &
               ndd, dpars, ndz, zpars, ndi, ipars, fval)
            do k=1,npols
              dintcoefs(1:ndim,1:ndim,k,ixyzs(i)+j-1) = & 
                dintcoefs(1:ndim,1:ndim,k,ixyzs(i)+j-1) + & 
                fval(1:ndim,1:ndim)*sigvals(k,l)*qwts(l)
            enddo
          enddo
        enddo

!
!  convert coefs to vals
!
        do k=1,npts
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
!  start filling in matrix entries
!
        do j=1,npts
          do l=1,npols
            isrc = ixyzs(i)+l-1
            amat(3*j-2,3*isrc-2) = dintvals(1,1,l,j)
            amat(3*j-2,3*isrc-1) = dintvals(1,2,l,j)
            amat(3*j-2,3*isrc-0) = dintvals(1,3,l,j)

            amat(3*j-1,3*isrc-2) = dintvals(2,1,l,j)
            amat(3*j-1,3*isrc-1) = dintvals(2,2,l,j)
            amat(3*j-1,3*isrc-0) = dintvals(2,3,l,j)

            amat(3*j-0,3*isrc-2) = dintvals(3,1,l,j)
            amat(3*j-0,3*isrc-1) = dintvals(3,2,l,j)
            amat(3*j-0,3*isrc-0) = dintvals(3,3,l,j)

            if(1.eq.1) then

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
            endif
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      call prin2('amat=*', amat, 24)
      call prinf('ndim=*', ndim, 1)
      do i=1,ndim*npts
         amat(i,i) = amat(i,i) + 0.5d0
      enddo

      ipt = 15
      jpt = 455

      i1 = 1
      j1 = 3
      
      i = 3*(ipt-1) + i1
      j = 3*(jpt-1) + j1
      h = huse 
      dpars(1) = dlam
      dpars(2) = dmu
      dpars(3) = h
      call fker(nd, srcvals(1,jpt), ndtarg, srcvals(1,ipt), &
         ndd, dpars, ndz, zpars, ndi, ipars, rtmp)
      ftmp = rtmp(i1,j1) 
      print *, "ftmp=", ftmp*wts(jpt)
      print *, "amat=",amat(i,j)
      info = 0
      call dgausselim(npts*ndim, amat, rhs, info, soln, dcond)

      print *, "dcond=", dcond
      print *, "info=", info

      call prin2('dpars=*',dpars,3)
      call prin2('xyz_out=*',xyz_out,12)

      do j=1,nin
        dpars(3) = hout
        rtmp(1:3,1:3) = 0
        call st3d_slp_vec(nd, xyz_out, ndim3, &
            xyz_in(1,j), ndd, dpars, ndz, zpars, ndi, ipars, rtmp)
      
        pot_ex(1,j) = rtmp(1,1)*strengths(1) + rtmp(1,2)*strengths(2) + &
                  rtmp(1,3)*strengths(3)
        pot_ex(2,j) = rtmp(2,1)*strengths(1) + rtmp(2,2)*strengths(2) + &
                  rtmp(2,3)*strengths(3)
        pot_ex(3,j) = rtmp(3,1)*strengths(1) + rtmp(3,2)*strengths(2) + &
                  rtmp(3,3)*strengths(3)

        pot(1:3,j) = 0
        h = huse 
        dpars(3) = h
        do i=1,npts
          rtmp(1:3,1:3) = 0 
          call st3d_slp_vec(nd, srcvals(1,i), ndim3, &
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

      do i=1,nin
        bmat(3*i-2,1) = 1.0d0
        bmat(3*i-1,2) = 1.0d0
        bmat(3*i-0,3) = 1.0d0

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
      call prin2('bmat=*',bmat,3*nin*6)
      call prin2('rhs=*',rhs,3*nin)

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
      call prin2('solnb=*',solnb,6)
      
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
      

      stop
      end
      
      
