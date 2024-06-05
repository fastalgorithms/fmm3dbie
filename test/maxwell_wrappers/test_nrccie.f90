      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)
      integer, allocatable :: norders(:), ixyzs(:), iptype(:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: ru(:,:), rv(:,:)
      integer ipars(2)


      real *8 xyz_out(3), xyz_in(3)

      real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 
      integer nover, npolso, nptso
      integer nnz, nquad
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      complex *16, allocatable :: wnear(:,:)

      real *8, allocatable :: srcover(:,:), wover(:)
      integer, allocatable :: ixyzso(:), novers(:)

      complex *16, allocatable :: einc(:,:), hinc(:,:)
      complex *16, allocatable :: zjvec(:,:), rho(:)
      complex *16, allocatable :: rhs(:,:), pot(:,:), sigma(:,:)
      complex *16, allocatable :: soln(:,:), zynm(:)
      complex *16, allocatable :: unm(:,:), xnm(:,:)
      complex *16, allocatable :: unm_uv(:,:)

      complex *16, allocatable :: w(:,:)
      complex *16 vf(3)
      complex *16 e_ex(3), h_ex(3)
      complex *16 e_comp(3), h_comp(3)
      complex *16 zvec(3), zrhsvec(3), zsolvec(3)
      complex *16 zju, zjx, zrhoy


      real *8, allocatable :: errs(:)
      real *8 thet,phi
      complex * 16 zpars(3)
      integer numit,niter
      character *200 title, fname, fname1, fname2

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0, isout1

      complex *16 potex, ztmp, ima, zk
      complex *16 alpha_rhs
      real *8 ptinfo_out(12)
      real *8 c0(3)

      complex *16 fjvals(0:100), fhvals(0:100), fjder(0:100)
      complex *16 fhder(0:100)

      complex *16 zrho3, zrho12, zjx3, zju3, zju12, zjx12
      complex *16 zmat(3,3), zmatuse(3,3), zmatinv(3,3)
      complex *16 zhu, zhx, zeu, zex, zey

      integer count1

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      eps = 1.0d-6

      zk = 1.1d0 
      zpars(1) = zk 
      zpars(2) = 1.0d0 

      norder = 4 
!
!  patch type, iptype0 = 1, for triangles with RV nodes
!                      = 11, for quadrangles with GL nodes
!                      = 12, for quadrangles with Cheb nodes      
!      
      iptype0 = 1

      a = 1
      na = 4
      c0(1) = 0
      c0(2) = 0
      c0(3) = 0

      call get_sphere_npat_mem(a, na, c0, norder, iptype0, &
       npatches, npts)

      allocate(srcvals(12,npts), srccoefs(9,npts))
      allocate(ixyzs(npatches+1), iptype(npatches))
      allocate(norders(npatches))
      
      allocate(ru(3,npts), rv(3,npts))

      call get_sphere_npat(a, na, c0, norder, iptype0, &
        npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

      xyz_out(1) = 3.17d0
      xyz_out(2) = -0.03d0
      xyz_out(3) = 3.15d0

      xyz_in(1) = 0.17d0
      xyz_in(2) = 0.23d0
      xyz_in(3) = -0.11d0



      allocate(wts(npts))
      call get_qwts(npatches, norders, ixyzs, iptype, npts, srcvals, &
        wts)

      call orthonormalize_all(srcvals(4:6,:), srcvals(10:12,:), ru, &
         rv, npts)
      iquadtype = 1

      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg, iptype_avg, rfac, rfac0)


      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))

      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, & 
        srccoefs, cms, rads)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO      

!
!    find near quadrature correction interactions
!
      print *, "entering find near mem"
      ndtarg = 12
      call findnearmem(cms, npatches, rad_near, ndtarg, srcvals, npts, &
        nnz)

      allocate(row_ptr(npts+1), col_ind(nnz))
      
      call findnear(cms, npatches, rad_near, ndtarg, srcvals, npts, &
        row_ptr, col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches, ixyzs, npts, nnz, row_ptr, col_ind, &
        iquad)


      nquad = iquad(nnz+1)-1
      allocate(wnear(9,nquad))
      print *, "beginning quadrature generation"
      call getnearquad_em_nrccie_pec(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, &
        nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)
      
      print *, "nquad=", nquad

      ikerorder = 0

!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches), ixyzso(npatches+1))

      print *, "beginning far order estimation"
      ndtarg = 12
      call get_far_order(eps, npatches, norders, ixyzs, iptype, cms, &
        rads, npts, srccoefs, ndtarg, npts, srcvals, ikerorder, &
        zpars(1), nnz, row_ptr, col_ind, rfac, novers, ixyzso)
      call prinf('novers=*',novers,10)

      npts_over = ixyzso(npatches+1)-1
      

      allocate(srcover(12,npts_over), wover(npts_over))

      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, npts_over, srcover)

      call get_qwts(npatches, novers, ixyzso, iptype, npts_over, &
        srcover, wover)
      print *, "done oversampling geometry"


      allocate(zjvec(3,npts), rho(npts), sigma(3,npts))
      allocate(pot(3,npts), rhs(3,npts), zynm(npts))
      allocate(unm(3,npts), xnm(3,npts))
      allocate(unm_uv(2,npts))
!
!
!  Get spherical harmonic data
!
      
      nn = 2
      mm = 1
      nmax = nn
      allocate(w(0:nmax,0:nmax))


      call l3getsph(nmax, mm, nn, ndtarg, srcvals, zynm, npts, w)
      call get_surf_grad(2, npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, zynm, unm_uv)

!
!  Set strengths of ynm, unm, and xnm terms
!
      zju = 1.0d0
      zjx = 0.0d0 + ima*0.3d0
      zrhoy = 1.2d0

      zrhsvec(1) = zju
      zrhsvec(2) = zjx
      zrhsvec(3) = zrhoy


      do i=1,npts
        unm(1:3,i) = unm_uv(1,i)*srcvals(4:6,i) + &
          unm_uv(2,i)*srcvals(7:9,i) 
        call dzcross_prod3d(srcvals(10,i), unm(1,i), xnm(1,i))
        zjvec(1:3,i) = zju*unm(1:3,i) + zjx*xnm(1:3,i)
        rho(i) = zrhoy*zynm(i)
      enddo

      do i=1,npts
        rhs(1,i) = ru(1,i)*zjvec(1,i) + ru(2,i)*zjvec(2,i) + &
          ru(3,i)*zjvec(3,i)
        
        rhs(2,i) = rv(1,i)*zjvec(1,i) + rv(2,i)*zjvec(2,i) + &
          rv(3,i)*zjvec(3,i)
        
        rhs(3,i) = rho(i)
      enddo


!
!  Apply the nrrcie integral operator 
!
      ndd = 0
      ndz = 2
      ndi = 0
      lwork = 0
      ndim = 3
      nker = 9
      call lpcomp_em_nrccie_pec_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, npts_over, ixyzso, srcover, wover, lwork, work, &
        ndim, rhs, pot)

      njh = nn + 5
      ifder = 1
      rscale = 1.0d0
      call besseljs3d(njh, zk, rscale, fjvals, ifder, fjder)
      call h3dall(njh, zk, rscale, fhvals, ifder, fhder)

      zrho3 = ima*zk*zk*(fjder(nn)*fhvals(nn) + fhder(nn)*fjvals(nn))/2
      zrho3 = zrho3 - ima*zk*(ima*zk)*zpars(2)*fhvals(nn)*fjvals(nn)

      zrho12 = ima*zk*zpars(2)*fjvals(nn)*fhvals(nn)

      zjx3 = 0
      zju3 = zk*(fjvals(nn)*fhvals(nn-1) - fjvals(nn+1)*fhvals(nn))
      zju3 = zju3 - zpars(2)*ima*zk*fjvals(nn)*fhvals(nn)
      zju3 = zju3*nn*(nn+1.0d0)

      zju12 = 0.5d0 - ima*(zk*fhvals(nn))*(fjvals(nn) + zk*fjder(nn))
      zju12 = zju12 + & 
        zpars(2)*ima*zk*(-ima*((nn+1.0d0)*fjvals(nn)*fhvals(nn-1) + & 
                                nn*fjvals(nn+1)*fhvals(nn) - &
                                zk*fjvals(nn+1)*fhvals(nn-1)))
      zjx12 = 0.5d0 + ima*(zk*fjvals(nn))*(fhvals(nn) + zk*fhder(nn)) 
      zjx12 = zjx12 + (zk*fjvals(nn))*(zk*fhvals(nn))

      zmat(1,1) = zju12
      zmat(1,2) = 0
      zmat(1,3) = zrho12
      
      zmat(2,1) = 0
      zmat(2,2)= zjx12
      zmat(2,3) = 0
      
      zmat(3,1) = zju3
      zmat(3,2) = zjx3
      zmat(3,3) = zrho3


      do i=1,3
        zsolvec(i) = 0
        do j=1,3
          zsolvec(i) = zsolvec(i) + zmat(i,j)*zrhsvec(j)
        enddo
      enddo
      
      
      
      erra = 0
      ra = 0

      err_vec = 0
      r_vec = 0
      do i=1,npts
        ra = ra + abs(zynm(i))**2*wts(i)
        erra = erra + abs(pot(3,i) - zsolvec(3)*zynm(i))**2*wts(i)

        zvec(1:3) = pot(1,i)*ru(1:3,i) + pot(2,i)*rv(1:3,i)

        err_vec = err_vec + abs(zvec(1) - zsolvec(1)*unm(1,i) - &
                                  zsolvec(2)*xnm(1,i))**2*wts(i)
        err_vec = err_vec + abs(zvec(2) - zsolvec(1)*unm(2,i) - &
                                  zsolvec(2)*xnm(2,i))**2*wts(i)
        err_vec = err_vec + abs(zvec(3) - zsolvec(1)*unm(3,i) - &
                                  zsolvec(2)*xnm(3,i))**2*wts(i)
      enddo
      erra = sqrt(erra/ra)
      err_vec = sqrt(err_vec/ra)
      call prin2('error in rho contrib to scalar part=*',erra,1)
      call prin2('error in rho contrib to vector part=*',err_vec,1)

!
!
!  Now test the solver
!
      allocate(errs(numit+1)) 
      numit = 100
      niter = 0
      rres = 0
      eps_gmres = eps
      allocate(hinc(3,npts), einc(3,npts))
      zhu = 0.0d0
      zhx = 0.0d0
      
      zeu = 0.0d0
      zex = 1.0d0
      zey = 0.0d0
      do i=1,npts
        hinc(1:3,i) = zhu*unm(1:3,i) + zhx*xnm(1:3,i)
        einc(1:3,i) = zeu*unm(1:3,i) + zex*xnm(1:3,i) + &
          zey*zynm(i)*srcvals(10:12,i)
      enddo

      call em_nrccie_pec_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, &
        einc, hinc, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, npts_over, ixyzso, srcover, wover, eps_gmres, niter, &
        errs, rres, zjvec, rho)
      
      do i=1,3
         do j=1,3
           zmatuse(j,i) = zmat(j,i)
         enddo
         zmatuse(i,i) = zmatuse(i,i) + 0.5d0
      enddo
      zrhsvec(1) = -zhx + zpars(2)*zeu
      zrhsvec(2) = zhu + zpars(2)*zex
      zrhsvec(3) = zey
      
      info = 0 
      call zinverse(3, zmatuse, info, zmatinv)

      do i=1,3
        zsolvec(i) = 0
        do j=1,3
          zsolvec(i) = zsolvec(i) + zmatinv(i,j)*zrhsvec(j)
        enddo
      enddo

!
!  Now test accuracy of solution
!
      errj = 0
      err_rho = 0
      do i=1,npts
        zvec(1:3) = zsolvec(1)*unm(1:3,i) + zsolvec(2)*xnm(1:3,i)
        errj = errj + abs(zvec(1) - zjvec(1,i))**2*wts(i)
        errj = errj + abs(zvec(2) - zjvec(2,i))**2*wts(i)
        errj = errj + abs(zvec(3) - zjvec(3,i))**2*wts(i)

        err_rho = err_rho + abs(zsolvec(3)*zynm(i) - rho(i))**2*wts(i)
      enddo

      errj = sqrt(errj/ra)
      err_rho = sqrt(err_rho/ra)

      call prin2('error in zjvec=*', errj, 1)
      call prin2('error in rho=*', err_rho, 1)

      stop
      end
!
!
!
!
!
      subroutine l3getsph(nmax, mm, nn, ndx, xyzs, ynms, npts, ynm)
      implicit real *8 (a-h,o-z)
      real *8 :: xyzs(ndx,npts)
      complex *16 ynms(npts), ima
      real *8 rat1(10000), rat2(10000)
      real *8 ynm(0:nmax,0:nmax)
      data ima/(0.0d0,1.0d0)/
  
      call ylgndrini(nmax, rat1, rat2)
  
      do i=1,npts
        x=xyzs(1,i)
        y=xyzs(2,i)
        z=xyzs(3,i)
        r=sqrt(x**2 + y**2 + z**2)
        call cart2polar(xyzs(1,i),r,theta,phi)
        ctheta = cos(theta)
        call ylgndrf(nmax, ctheta, ynm, rat1, rat2)
        ynms(i) = ynm(nn,abs(mm))*exp(ima*mm*phi)        
      enddo
       
      return
      end




