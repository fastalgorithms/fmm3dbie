!
!  kernel_registry.f90
!
!  Runtime name -> (fker, ipv, ifcomplex) registry. Entries are generated
!  from !@kernel annotations by gen_registry_init.sh.
!
      module kernel_registry
        implicit none

        integer, parameter :: KR_NAME_LEN = 64
        integer, parameter :: KR_MAX      = 1024

        type :: kr_entry
          character(len=KR_NAME_LEN) :: name = ''
          procedure(), pointer, nopass :: fker => null()
          integer *8 :: ipv = 0
          integer *8 :: ifcomplex = 0
        end type kr_entry

        type(kr_entry), save :: kr_table(KR_MAX)
        integer *8,     save :: kr_count = 0
        logical,        save :: kr_initialised = .false.

      contains

!-----------------------------------------------------------------------
        subroutine register_kernel(name, fker, ipv, ifcomplex)
!
!  Register one kernel by name; a repeated name overwrites. ipv is the
!  self-quadrature type (0 weakly singular, 1 pv, 2 hypersingular),
!  ifcomplex is 1 if fker returns complex *16.
!
          implicit none
          character(len=*), intent(in) :: name
          external                     :: fker
          integer *8,       intent(in) :: ipv, ifcomplex

          integer *8 :: i, islot

          islot = 0
          do i = 1, kr_count
            if (trim(kr_table(i)%name) .eq. trim(name)) then
              islot = i
              exit
            endif
          enddo

          if (islot .eq. 0) then
            if (kr_count .ge. KR_MAX) then
              call prinf('register_kernel: KR_MAX exceeded*', kr_count, 1)
              return
            endif
            kr_count = kr_count + 1
            islot = kr_count
          endif

          kr_table(islot)%name = trim(name)
          kr_table(islot)%fker => fker
          kr_table(islot)%ipv  = ipv
          kr_table(islot)%ifcomplex = ifcomplex

          return
        end subroutine register_kernel

!-----------------------------------------------------------------------
        subroutine lookup_kernel(name, fker, ipv, ifcomplex, ier)
!
!  Resolve a name to (fker, ipv, ifcomplex). ier = 0 on success, 1 if the
!  name is not registered. Initializes the registry on first use.
!
          implicit none
          character(len=*), intent(in)  :: name
          procedure(), pointer, intent(out) :: fker
          integer *8,       intent(out) :: ipv, ifcomplex, ier

          integer *8 :: i

          if (.not. kr_initialised) call init_kernel_registry()

          ier = 1
          fker => null()
          ipv = 0
          ifcomplex = 0
          do i = 1, kr_count
            if (trim(kr_table(i)%name) .eq. trim(name)) then
              fker => kr_table(i)%fker
              ipv  = kr_table(i)%ipv
              ifcomplex = kr_table(i)%ifcomplex
              ier  = 0
              return
            endif
          enddo

          return
        end subroutine lookup_kernel

!-----------------------------------------------------------------------
        subroutine kernel_count(n)
!
!  Number of registered kernels. Initializes the registry on first use.
!
          implicit none
          integer *8, intent(out) :: n
          if (.not. kr_initialised) call init_kernel_registry()
          n = kr_count
          return
        end subroutine kernel_count

!-----------------------------------------------------------------------
        subroutine kernel_name(i, name, ln)
!
!  Name of the i-th registered kernel (1-based); ln = 0 if i is out of
!  range. Initializes the registry on first use.
!
          implicit none
          integer *8,       intent(in)  :: i
          character(len=*), intent(out) :: name
          integer *8,       intent(out) :: ln

          if (.not. kr_initialised) call init_kernel_registry()

          name = ''
          ln   = 0
          if (i .lt. 1 .or. i .gt. kr_count) return
          name = kr_table(i)%name
          ln   = len_trim(kr_table(i)%name)

          return
        end subroutine kernel_name

      end module kernel_registry
