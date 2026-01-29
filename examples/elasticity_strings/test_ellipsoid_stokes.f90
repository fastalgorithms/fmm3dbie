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
      real *8 strengths(3)

      procedure (), pointer :: fker


      real *8 c0(3), abc(3), verts(2,4), xyz_out(3), xyz_in(3)
      integer *8 nabc(3)

      complex *16 zpars
      real *8 uv(2)
      character *1 transa, transb
      external l3d_slp, st3d_slp_vec, st3d_comb_vec
      
      call prini(6,13)
      
      npatches = 0
      npts = 0

      abc(1) = 2.1d0
      abc(2) = 1.0d0
      abc(3) = 4.0d0

      abc(1) = 1.0d0
      abc(2) = 1.0d0
      abc(3) = 1.0d0

      abc(1:3) = abc(1:3)*1.1d0

      nabc(1) = 5
      nabc(2) = 3
      nabc(3) = 10

      nabc(1) = 2
      nabc(2) = 2
      nabc(3) = 2
      
      c0(1) = 0
      c0(2) = 0
      c0(3) = 0

      xyz_out(1) = 3.1d0
      xyz_out(2) = 1.19d0
      xyz_out(3) = -4.5d0

      xyz_in(1) = 0.31d0
      xyz_in(2) = -0.07d0
      xyz_in(3) = 0.11d0

      norder = 8
      iptype0 = 1

      

      call get_ellipsoid_npat_mem(abc, nabc, c0, norder, iptype0, &
        npatches, npts)

      allocate(srcvals(12,npts), srccoefs(9,npts))
      allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))
      allocate(wts(npts))
      
      call get_ellipsoid_npat(abc, nabc, c0, norder, iptype0, &
        npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)
      
      
      call get_qwts(npatches, norders, ixyzs, iptype, npts, & 
        srcvals, wts)

      npols = ixyzs(2) - ixyzs(1)
      allocate(uvs(2,npols), umat(npols,npols),vmat(npols,npols))
      allocate(wtmp(npols))
      
      call get_disc_exps(norder, npols, iptype0, uvs, umat, vmat, wtmp)
      
      nv = 0
      call get_boundary_vertices(iptype0, verts, nv)

      istrat = 2
      intype = 1
      npatches0 = 1
      eps = 0.51d-7

      isd = 0
      ndsc = 9
      ndtarg = 12

      ifp = 0
      
      fker => st3d_comb_vec
      ndim = 3
      nd = ndim*ndim
      

      ntarguse = npts - npols
      allocate(dintcoefs(ndim,ndim,npols,npts))
      allocate(dintvals(ndim,ndim,npols,npts))
      allocate(amat(ndim*npts, ndim*npts))
      allocate(targuse(12,ntarguse))
      allocate(fval(ndim,ndim))

      dlam = 1.1d0
      dmu = 3.21d0
      dpars(1) = dlam
      dpars(2) = dmu

      dpars(1) = 1.1d0
      dpars(2) = 1.0d0

      itargptr = 1
      ntargptr = ntarguse

      ntrimax = 3000

      ndd = 3
      ndi = 0
      ndz = 0

      nqorder = 8
      call get_quadparams_adap(eps, 1, nqorder, eps_adap, nlev, &
        nqorder_f)
      rfac = 2.0d0
      
      ifmetric = 0
      
      nmax = 20000
      allocate(xs(nmax), ys(nmax), ws(nmax))
      allocate(srctmp(12,nmax), qwts(nmax))
      allocate(sigvals(npols,nmax))

      ipv = 0

      allocate(rhs(ndim,npts), soln(ndim,npts), soln2(ndim,npts))
      allocate(rtmp(ndim,ndim))

      strengths(1) = 1.1d0
      strengths(2) = 2.1d0
      strengths(3) = 0.3d0

      
      do i=1,npts
!        call l3d_slp(xyz_out, 3, srcvals(1,i), 0, dpars, 1, zpars, &
!          0, ipars, rhs(1,i))
        call st3d_slp_vec(nd, xyz_out, 3, srcvals(1,i), 0, dpars, 1, &
          zpars, 0, ipars, rtmp)
        rhs(1,i) = rtmp(1,1)*strengths(1) + rtmp(1,2)*strengths(2) + &
                   rtmp(1,3)*strengths(3)
        rhs(2,i) = rtmp(2,1)*strengths(1) + rtmp(2,2)*strengths(2) + &
                   rtmp(2,3)*strengths(3)
        rhs(3,i) = rtmp(3,1)*strengths(1) + rtmp(3,2)*strengths(2) + &
                   rtmp(3,3)*strengths(3)
        soln(1:ndim,i) = 0
        soln2(1:ndim,i) = 0
      enddo

      numit = 200
      ifinout = 0
      niter = 0
      allocate(errs(numit+1))
      eps_gmres = 0.5d-9
      dpars_solver(1) = 1.0d0
      dpars_solver(2) = 0.0d0
      print *, "npts=",npts

      call stok_comb_vel_solver(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, eps, dpars, numit, ifinout, rhs, &
        eps_gmres, niter, errs, rres, soln)
     
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
        h = 0.1d0
        dpars(1) = 1.1d0
        dpars(2) = 1.0d0
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
        ipv = 0
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
!            print *, "isrc=",isrc
!            print *, "ndim=",ndim
!            amat((ndim*(j-1)+1):ndim*j,(ndim*(isrc-1)+1):ndim*isrc) = & 
!              dintvals(1:ndim,1:ndim,l,j)
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
            amat(3*j-2,3*isrc-2) = amat(3*j-2,3*isrc-2) + &
              srcvals(10,j)*srcvals(10,isrc)*wts(isrc)
            amat(3*j-2,3*isrc-1) = amat(3*j-2,3*isrc-1) + &
              srcvals(10,j)*srcvals(11,isrc)*wts(isrc)
            amat(3*j-2,3*isrc-0) = amat(3*j-2,3*isrc-0) + &
              srcvals(10,j)*srcvals(12,isrc)*wts(isrc)

            amat(3*j-1,3*isrc-2) = amat(3*j-1,3*isrc-2) + &
              srcvals(11,j)*srcvals(10,isrc)*wts(isrc)
            amat(3*j-1,3*isrc-1) = amat(3*j-1,3*isrc-1) + &
              srcvals(11,j)*srcvals(11,isrc)*wts(isrc)
            amat(3*j-1,3*isrc-0) = amat(3*j-1,3*isrc-0) + &
              srcvals(11,j)*srcvals(12,isrc)*wts(isrc)

            amat(3*j-0,3*isrc-2) = amat(3*j-0,3*isrc-2) + &
              srcvals(12,j)*srcvals(10,isrc)*wts(isrc)
            amat(3*j-0,3*isrc-1) = amat(3*j-0,3*isrc-1) + &
              srcvals(12,j)*srcvals(11,isrc)*wts(isrc)
            amat(3*j-0,3*isrc-0) = amat(3*j-0,3*isrc-0) + &
              srcvals(12,j)*srcvals(12,isrc)*wts(isrc)
            endif
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      call prin2('amat=*', amat, 24)
      call prinf('ndim=*', ndim, 1)
      do i=1,ndim*npts
         amat(i,i) = amat(i,i) - 0.5d0
      enddo
      ipt = 15
      jpt = 455

      i1 = 1
      j1 = 3
      
      i = 3*(ipt-1) + i1
      j = 3*(jpt-1) + j1
      call fker(nd, srcvals(1,jpt), ndtarg, srcvals(1,ipt), &
         ndd, dpars, ndz, zpars, ndi, ipars, rtmp)
      ftmp = rtmp(i1,j1) + srcvals(9+i1,ipt)*srcvals(9+j1,jpt)
      print *, "ftmp=", ftmp*wts(jpt)
      print *, "amat=",amat(i,j)
      info = 0
      call dgausselim(npts*ndim, amat, rhs, info, soln2, dcond)

      print *, "dcond=", dcond
      print *, "info=", info

      erra = 0
      ra = 0
      do i=1,npts
       do j=1,ndim
         erra = erra + (soln(j,i) - soln2(j,i))**2*wts(i)
         ra = ra + soln(j,i)**2*wts(i)
       enddo
      enddo
      call prin2('soln=*', soln, 24)
      call prin2('soln2=*', soln2, 24)

      erra = sqrt(erra/ra)
      print *, "Error in solution=", erra
      

      stop
      end
      
      
