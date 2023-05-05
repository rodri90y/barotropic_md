!
! Model contants
!
MODULE module_configure

  implicit none
  integer, parameter          :: dp=kind(0.d0)                        ! double precision
  real(dp), parameter         :: earthr        = 6.37122e06, &
                                 pi            = 3.141592653589793_dp, &
                                 omega         = 7.292e-05, &
                                 rad           = pi/180.

  real                        :: dx, &
                                 dy, &
                                 dt, &
                                 dphi, &
                                 dtheta
  real(dp)                    :: beta, &
                                 sigma, &
                                 alpha, &
                                 lat, &
                                 lon, &
                                 A0
  integer(kind=4)             :: e_ns, e_we, nint
  character(len=255)          :: io_out
  integer                     :: jactyp, io_dt




CONTAINS
  subroutine init( )
!<descricao>
!
!<descricao>
    implicit none

    call initial_config( )

    beta = 2.0*omega*cos(lat*rad)/earthr
    dphi = dy/(111.1*1.E3)
    dtheta = dx/(111.1*1.E3)*cos(lat*rad)

  end subroutine init


  subroutine initial_config( )
!<descricao>
!
!
! tip: unexpected characters, such as DOS CR-LF characters cause reading error
!<descricao>
    implicit none
    integer                     :: io_status
    integer, parameter          :: nml_read_unit = 10
    namelist/init/ sigma, A0, alpha, jactyp
    namelist/domains/ dx, dy, dt, e_ns, e_we, lat, lon, nint
    namelist/io/ io_out, io_dt

    open(UNIT   = nml_read_unit ,      &
         FILE   = "namelist.dat",      &
         IOSTAT = io_status)

    if ( io_status .ne. 0 ) then
        print *, 'ERROR OPENING namelist.dat'
    end if

    read(UNIT=nml_read_unit,NML=init)
    read(UNIT=nml_read_unit,NML=domains)
    read(UNIT=nml_read_unit,NML=io)

    close(UNIT=nml_read_unit)

  end subroutine initial_config

END MODULE module_configure
