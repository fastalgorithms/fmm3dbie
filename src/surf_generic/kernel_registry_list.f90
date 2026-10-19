!
!  kernel_registry_list.f90
!
!  mwrap-callable queries of the runtime kernel registry, so MATLAB asks
!  the registry what exists instead of mirroring the list by hand.
!
!  Names cross the mwrap boundary as integer vectors of ASCII character
!  codes, never as Fortran strings (mwrap has no portable path for passing
!  a character string to Fortran).
!
!-----------------------------------------------------------------------
      subroutine kernel_registry_count(n)
!
!  Number of registered kernels.
!
      use kernel_registry, only: kernel_count
      implicit none
      integer *8, intent(out) :: n

      call kernel_count(n)

      return
      end subroutine kernel_registry_count

!-----------------------------------------------------------------------
      subroutine kernel_registry_names(nmax, nchar, codes, lens)
!
!  ASCII codes of the first nmax registered names, one column per kernel:
!  codes(1:lens(i), i) is name i, zero-padded to nchar rows. lens(i) = 0
!  marks an out-of-range slot.
!
      use kernel_registry, only: kernel_count, kernel_name
      implicit none
      integer *8, intent(in)  :: nmax, nchar
      integer *8, intent(out) :: codes(nchar, nmax)
      integer *8, intent(out) :: lens(nmax)

      character(len=64) :: name
      integer *8 :: i, j, n, ln

      call kernel_count(n)

      do i = 1, nmax
        do j = 1, nchar
          codes(j, i) = 0
        enddo
        lens(i) = 0
        if (i .gt. n) cycle
        call kernel_name(i, name, ln)
        ln = min(ln, nchar)
        do j = 1, ln
          codes(j, i) = iachar(name(j:j))
        enddo
        lens(i) = ln
      enddo

      return
      end subroutine kernel_registry_names

!-----------------------------------------------------------------------
      subroutine kernel_registry_info(nchar, kername, ipv, ifcomplex, ier)
!
!  Look up one kernel by name.
!
!  Input arguments:
!    - nchar: integer *8
!        length of kername
!    - kername: integer *8(nchar)
!        ASCII codes of the kernel name
!
!  Output arguments:
!    - ipv: integer *8
!        self-quadrature type of the kernel
!    - ifcomplex: integer *8
!        1 if the kernel is complex valued, 0 if real valued
!    - ier: integer *8
!        0 if the kernel is registered, 1 otherwise
!
      use kernel_registry, only: lookup_kernel
      implicit none
      integer *8, intent(in)  :: nchar
      integer *8, intent(in)  :: kername(max(nchar,1_8))
      integer *8, intent(out) :: ipv, ifcomplex, ier

      character(len=64) :: name
      integer *8 :: ln
      procedure(), pointer :: fker => null()

      call kername_to_fstr(kername, nchar, name, ln)
      call lookup_kernel(trim(name), fker, ipv, ifcomplex, ier)

      return
      end subroutine kernel_registry_info

!-----------------------------------------------------------------------
      subroutine kername_to_fstr(kername, nchar, name, ln)
!
!  Rebuild a Fortran string from a vector of ASCII codes. At most 64
!  characters are kept; a code <= 0 ends the name.
!
      implicit none
      integer *8,        intent(in)  :: kername(*)
      integer *8,        intent(in)  :: nchar
      character(len=64), intent(out) :: name
      integer *8,        intent(out) :: ln
      integer *8 :: i, nmax, ic

      name = ''
      nmax = min(nchar, 64_8)
      ln = 0
      do i = 1, nmax
        ic = kername(i)
        if (ic .le. 0) exit
        name(i:i) = achar(int(ic))
        ln = i
      enddo

      return
      end subroutine kername_to_fstr
