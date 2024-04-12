!
!
!  This file contains various maxwell kernels that
!  are to be used in the representation.
!
!  Unlike the helmholtz case, it is beneficial
!  to construct the kernels on a representation by representation
!  basis
!
!  The subroutine names will take the following form
!
!  em3d_<rep>_kers
!  em3d_<rep>_kers_postproc
!  em3d_<rep>_kers_comp  (temporarily used for until vectorized
!                         versions of adaptive integration is
!                         available)
!  em3d_<rep>_kers_postproc_comp  (temporarily used for until vectorized
!                         versions of adaptive integration is
!                         available)
!  For consistency of notation, we will stick to the following convention
!  electic current: j
!  magnetic current: k
!  electric charge: r
!  magnetic charge: q
!  
!  for the physical variables
!
!  The routines rely on srcinfo and targinfo arrays, standardized
!  in the following way
!
!  *info(1:3) = xyz
!  *info(4:6) = dxyz/du
!  *info(7:9) = dxyz/dv
!  *info(10:12) = normal vector
!
!  The u and v components above correspond to a local orthonormal 
!  basis consutrcted from the above quantities as follows
!
!  uvec = dxyz/du/|dxyz/du|
!  vvec = (n \times dxyz/du) = (n \times dxyz/du)/|n\times dxyz/du|
!
!  so by ju and jv above, we mean their projections onto the
!  uvec and vvec directions
!
!
!--------------------------
!  NRCCIE routine
subroutine em3d_nrccie_kers(nd, src, ndt, targ, ndd, dpars, ndz, &
   zpars, ndi, ipars, val)
!
!  returns all the relevant kernels for nrccie
!  Input arguments:
!    - nd: integer
!        number of kernels, should be 9
!    - src: real *8 (12)
!        source info
!          * src(1:3)   - xyz
!          * src(4:6)   - dxyz/du
!          * src(7:9)   - dxyz/dv
!          * src(10:12) - normals
!    - ndt: integer
!        size of targ array, must be atleast 12, and only
!        first 12 components used
!    - targ: real *8 (ndt)
!        target info
!          * targ(1:3)   - xyz
!          * targ(4:6)   - dxyz/du
!          * targ(7:9)   - dxyz/dv
!          * targ(10:12) - normals
!    - ndd: integer
!        number of double precision parameters (unused for this routine)
!    - dpars: real *8 (ndd)
!        real precision parameters (unused for this routine)
!    - ndz: integer
!        number of complex parameters, only first 2 components used 
!    - zpars: complex *16 (ndz)
!        zpars(1) - Maxwell wavenumber
!        zpars(2) - scaling parameter alpha above
!    - ndi: integer
!        number of integer parameters, unused for this routine
!    - ipars: integer (ndi)
!        integer parameters, unused for this routine
!
!  Output arguments:
!    - val: complex *16(9)
!        value of the kernel above
!        * val(1) corresponds to
!          \ru \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{u} \ru]
!        * val(2) corresponds to
!          \ru \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{v} \rv]
!        * val(3) corresponds to
!            \ru \cdot (-\alpha n \times n \times \nabla S_{k}[\rho])        
!        * val(4) corresponds to
!          \rv \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{u} \ru]
!        * val(5) corresponds to
!          \rv \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{v} \rv]
!        * val(6) corresponds to
!          \rv \cdot (-\alpha n \times n \times \nabla S_{k}[\rho])
!        * val(7) corresponds to
!         (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{u} \ru]
!        * val(8) corresponds to
!         (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{v} \rv]
!        * val(9) corresponds to S_{k}'[\rho]        

!
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nd
      real *8 :: src(12), targ(ndt), dpars(ndd)
      real *8 :: over4pi, ru_s(3), rv_s(3), rr, ru_t(3), rv_t(3)
      integer ipars(ndi)
      complex *16 :: zpars(ndz), val(nd), zk, alpha, zexp, zgrad(3)
      complex *16 :: zfac, zjvec_u(3), zjvec_v(3), ztmp
      complex *16 :: zcurljvec_u(3), zcurljvec_v(3)
      complex *16 :: zncurljvec_u(3), zncurljvec_v(3)
      complex *16 :: znnjvec_u(3), znnjvec_v(3)
      complex *16 :: zvec(3), zn, znnrho(3)
      complex *16 :: ima

      data ima/(0.0d0,1.0d0)/
      data over4pi/0.07957747154594767d0/

      zk = zpars(1)
      alpha = zpars(2)

      dx = targ(1) - src(1)
      dy = targ(2) - src(2)
      dz = targ(3) - src(3)

      rr = sqrt(dx**2 + dy**2 + dz**2)
      rinv = 1.0d0/rr
      rinv2 = rinv*rinv
!   zexp also serves as a proxy for S \rho
      zexp = exp(ima*zk*rr)*rinv*over4pi
      zfac = (ima*zk*r - 1.0d0)*zexp*rinv2

!   zgrad also  serves as a proxy for \nabla S \rho      
      zgrad(1) = zfac*dx
      zgrad(2) = zfac*dy
      zgrad(3) = zfac*dz

      zn = zgrad(1)*targ(10) + zgrad(2)*targ(11) + zgrad(3)*targ(12)
      znnrho(1:3) = zn*targ(10:12) - zgrad(1:3)


      call orthonormalize(src(4:6), src(10:12), ru_s, rv_s)
      call orthonormalize(targ(4:6), targ(10:12), ru_t, rv_t)

!  Compute S_{k}[\ru] and store it in zjvec_u      

      zjvec_u(1) = zexp*ru_s(1)
      zjvec_u(2) = zexp*ru_s(2)
      zjvec_u(3) = zexp*ru_s(3)

      ztmp = zjvec_u(1)*targ(10) + zjvec_u(2)*targ(11) + &
        zjvec_u(3)*targ(12)
      znnjvec_u(1:3) = ztmp*targ(10:12) - zjvec_u(1:3)
      
!  Compute S_{k}[\rv] and store it in zjvec_v      

      zjvec_v(1) = zexp*rv_s(1)
      zjvec_v(2) = zexp*rv_s(2)
      zjvec_v(3) = zexp*rv_s(3)

      ztmp = zjvec_v(1)*targ(10) + zjvec_v(2)*targ(11) + &
        zjvec_v(3)*targ(12)
      znnjvec_v(1:3) = ztmp*targ(10:12) - zjvec_v(1:3)

!  Compute \nabla \times S_{k}[\ru] and store it in zcurljvec_u
      
      zcurljvec_u(1) = zgrad(2)*ru_s(3) - zgrad(3)*ru_s(2)
      zcurljvec_u(2) = zgrad(3)*ru_s(1) - zgrad(1)*ru_s(3)
      zcurljvec_u(3) = zgrad(1)*ru_s(2) - zgrad(2)*ru_s(1)

      call dzcross_prod3d(targ(10:12), zcurljvec_u, zncurljvec_u)

!  Compute \nabla \times S_{k}[\rv] and store it in zcurljvec_v

      zcurljvec_v(1) = zgrad(2)*rv_s(3) - zgrad(3)*rv_s(2)
      zcurljvec_v(2) = zgrad(3)*rv_s(1) - zgrad(1)*rv_s(3)
      zcurljvec_v(3) = zgrad(1)*rv_s(2) - zgrad(2)*rv_s(1)

      call dzcross_prod3d(targ(10:12), zcurljvec_v, zncurljvec_v)

!
!   Now start combining the kernels to compute vals 
      zvec = -zncurljvec_u + alpha*znnjvec_u
      
      val(1) = ru_t(1)*zvec(1) + ru_t(2)*zvec(2) + ru_t(3)*zvec(3)
      val(4) = rv_t(1)*zvec(1) + rv_t(2)*zvec(2) + rv_t(3)*zvec(3)

      zvec = -zncurljvec_v + alpha*znnjvec_v
      val(2) = ru_t(1)*zvec(1) + ru_t(2)*zvec(2) + ru_t(3)*zvec(3)
      val(5) = rv_t(1)*zvec(1) + rv_t(2)*zvec(2) + rv_t(3)*zvec(3)

      zvec = -alpha*znnrho
      val(3) = ru_t(1)*zvec(1) + ru_t(2)*zvec(2) + ru_t(3)*zvec(3)
      val(6) = rv_t(1)*zvec(1) + rv_t(2)*zvec(2) + rv_t(3)*zvec(3)
      

      val(9) = zn

!        * val(7) corresponds to
!         (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{u} \ru]
!        * val(8) corresponds to
!         (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{v} \rv]


! 
!

return
end
