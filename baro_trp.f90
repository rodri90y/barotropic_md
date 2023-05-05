program baro_trp
! compile command:  gfortran -I/usr/local/netcdf/include module_constants.f90 module_solve.f90 baro_trp.f90 -o barotrp -fbounds-check
!

    use module_configure
    use module_solve

    implicit none
    integer             :: error
    integer :: i, j
    real(dp), allocatable :: psi(:,:)

    call init()

    allocate(psi(e_ns, e_we),stat=error)

    if (error /= 0) stop "*** Not enough memory ***"


    call psi_field(psi)

    call psi_prog(psi, dx, dy, dt)


!
!write(*,*) 'JAC'


!
!!    !saida jacobiano
!    do i=1, e_ns
!    write(*,*)
!      do j=1, e_we
!        write(*,'(E10.1)', advance='no') advf(i,j)
!      end do
!    end do
!write(*,*)

CONTAINS


end program baro_trp
