c
c
c   This file contains subroutines for computing the subtraction
c   step corresponding to the computed quadrature correction
c
c   Using this feature would increase the memory requirements of
c   the code but imporve the CPU time performance of the solver
c
c
c
c
c     This file has the following user callable routines
c
c        ?getnearquadsub_guru - guru interface for 
c         computing subraction corresponding to near field quadrature 
c        
c
c         z - complex
c         d - double precision
c
c
c
c
c
c

      subroutine zgetnearquadsub_guru(npatches,norders,
     1   ixyzso,iptype,npts_over,srcover,wover,ndtarg,ntarg,targvals,
     2   ipatch_id,uvs_targ,thresh,ipv,fker,ndd,dpars,ndz,zpars,
     3   ndi,ipars,nnz,row_ptr,col_ind,iquad,nquad,wnear)
c
c------------------------
c  This subroutine generates the fmm subtraction assocaited
c  with the near field quadrature 
c  for a given kernel which is assumed to be
c  a compact/principal value integral operator
c  where the near field is specified by the user 
c  in row sparse compressed format.
c
c
c  Input arguments:
c
c    - npatches: integer(8)
c        number of patches
c    - norders: integer(8)(npatches)
c        order of discretization on each patch
c    - ixyzso: integer(8)(npatches+1)
c        ixyzso(i) denotes the starting location in srcover,
c    - iptype: integer(8)(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c    - npts_over: integer(8)
c        total number of oversampled points on the boundary
c    - srcover: double precision (12,npts_over)
c        xyz(u,v) and derivative info sampled at the 
c        oversampled nodes on the surface
c        * srcover(1:3,i) - xyz info
c        * srcoverls(4:6,i) - dxyz/du info
c        * srcoverls(7:9,i) - dxyz/dv info
c        * srcoverls(10:12,i) - normals info
c    - wover: double precision (npts_over)
c        quadrature weights for smooth functions
c        discretized at oversampled nodes
c    - ndtarg: integer(8)
c        leading dimension of target array
c    - ntarg: integer(8)
c        number of targets
c    - targvals: double precision (ndtarg,ntarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - ipatch_id: integer(8)(ntarg)
c        patch id of target, patch_id = -1, if target off-surface
c    - uvs_targ: double precision (2,ntarg)
c        local uv coordinates on patch if target on surface
c    - thresh: double precision
c        threshold for ignoring close interactions
c    - ipv: integer(8)
c        Flag for choosing type of self-quadrature
c        * ipv = 0, for compact/weakly singular operators
c        * ipv = 1, for singular operators
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence 
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c    - ndd: integer(8)
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer(8)
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer(8)
c        number of integer(8) parameters
c    - ipars: integer(8)(ndi)
c        integer(8) parameters
c    - nnz: integer(8)
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(8)(ntarg+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer(8) (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(8)(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - nquad: integer(8)
c        number of entries in wnear
c
c  Output parameters:
c
c    - wnear: double complex(nquad)
c        near field quadrature subtractions 
c----------------------------------------------------               
c
      implicit real *8 (a-h,o-z)
      implicit integer(8) (i-n)
      integer(8), intent(in) :: ndi,ndd,ndz,ipv
      integer(8), intent(in) :: ipars(ndi)
      integer(8), intent(in) :: ndtarg
      integer(8), intent(in) :: npatches,norders(npatches),npts_over
      integer(8), intent(in) :: ixyzso(npatches+1),iptype(npatches)
      real *8, intent(in) :: srcover(12,npts_over),thresh
      real *8, intent(in) :: wover(npts_over)
      integer(8), intent(in) :: ntarg,ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: targvals(ndtarg,ntarg)
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer(8), intent(in) :: nnz
      integer(8), intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer(8), intent(in) :: iquad(nnz+1)
      complex *16, intent(out) :: wnear(nquad)

      integer(8) i,j,ipatch,istart,npols,iquadstart,l
      real *8 rr,threshsq
      complex *16 val

      external fker

      threshsq = thresh**2


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,ipatch)
C$OMP$PRIVATE(istart,npols,iquadstart,l,val,rr)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          ipatch = col_ind(j)
          istart = ixyzso(ipatch)
          npols = ixyzso(ipatch+1)-ixyzso(ipatch)
          iquadstart = iquad(j)
          do l=1,npols
            rr = (srcover(1,l+istart-1) - targvals(1,i))**2 + 
     1          (srcover(2,l+istart-1) - targvals(2,i))**2 + 
     2          (srcover(3,l+istart-1) - targvals(3,i))**2 
            if(rr.gt.threshsq) then
              call fker(srcover(1,l+istart-1),ndtarg,targvals(1,i),ndd,
     1          dpars,ndz,zpars,ndi,ipars,val)
              wnear(iquadstart+l-1) = val*wover(l+istart-1)
            endif
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
