MODULE module_io

  use netcdf
  use module_configure
  implicit none

  integer, private, parameter   :: NDIMS = 4,&
                                   NLVLS = 1
  integer(kind=4), private      :: ncid
  integer, private              :: psi_varid, vort_varid, u_varid, v_varid
  ! The start and count arrays will tell the netCDF library where to
  ! write our data.
  integer, private               :: start(NDIMS), count(NDIMS)



CONTAINS

  subroutine open_ncdf()
  !<descricao>
  !
  !<descricao>
    implicit none

    integer(kind=4)               :: i, j


    character(len=*), parameter   :: LVL_NAME = "level", &
                                     LAT_NAME = "lat", &
                                     LON_NAME = "lon", &
                                     REC_NAME = "time"

    ! These program variables hold the latitudes and longitudes.
    real                          :: lats(e_ns), lons(e_we)
    integer                       :: lvl_dimid, lon_dimid, lat_dimid, rec_dimid
    integer                       :: lon_varid, lat_varid

    ! We will create two netCDF variables, one each for temperature and
    ! pressure fields.
    character(len=*), parameter   :: PSI_NAME = "streamfunction", &
                                     VORT_NAME = "vorticity", &
                                     U_NAME = "u", &
                                     V_NAME = "v"

    integer                       :: dimids(NDIMS)

    ! We recommend that each variable carry a "units" attribute.
    character(len=*), parameter   :: UNITS = "units", &
                                     PSI_UNITS = "m^2_s^-1", &
                                     VORT_UNITS = "s^-1", &
                                     U_UNITS = "m_s^-1", &
                                     V_UNITS = "m_s^-1", &
                                     LAT_UNITS = "degrees", &
                                     LON_UNITS = "degrees"

    ! Create pretend data. If this wasn't an example program, we would
    ! have some real data to write, for example, model output.
    do i = 1, e_ns
       lats(i) = lat + (i - e_ns/2.)*dphi
    end do
    do i = 1, e_we
       lons(i) = lon + (i - e_we/2.)*dtheta
    end do

    ! Create the file.
    call check( nf90_create(trim(io_out), nf90_clobber, ncid) )


    ! Define the dimensions. The record dimension is defined to have
    ! unlimited length - it can grow as needed. In this example it is
    ! the time dimension.
    call check( nf90_def_dim(ncid, LVL_NAME, NLVLS, lvl_dimid) )
    call check( nf90_def_dim(ncid, LAT_NAME, e_ns, lat_dimid) )
    call check( nf90_def_dim(ncid, LON_NAME, e_we, lon_dimid) )
    call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )


    ! Define the coordinate variables. We will only define coordinate
    ! variables for lat and lon.  Ordinarily we would need to provide
    ! an array of dimension IDs for each variable's dimensions, but
    ! since coordinate variables only have one dimension, we can
    ! simply provide the address of that dimension ID (lat_dimid) and
    ! similarly for (lon_dimid).
    call check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
    call check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )


    ! Assign units attributes to coordinate variables.
    call check( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
    call check( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )

    ! The dimids array is used to pass the dimids of the dimensions of
    ! the netCDF variables. Both of the netCDF variables we are creating
    ! share the same four dimensions. In Fortran, the unlimited
    ! dimension must come last on the list of dimids.
    dimids = (/ lon_dimid, lat_dimid, lvl_dimid, rec_dimid /)

    ! Define the netCDF variables for the variable data.
    call check( nf90_def_var(ncid, PSI_NAME, NF90_REAL, dimids, psi_varid) )
    call check( nf90_def_var(ncid, VORT_NAME, NF90_REAL, dimids, vort_varid) )
    call check( nf90_def_var(ncid, U_NAME, NF90_REAL, dimids, u_varid) )
    call check( nf90_def_var(ncid, V_NAME, NF90_REAL, dimids, v_varid) )

    ! Assign units attributes to the netCDF variables.
    call check( nf90_put_att(ncid, psi_varid, UNITS, PSI_UNITS) )
    call check( nf90_put_att(ncid, vort_varid, UNITS, VORT_UNITS) )
    call check( nf90_put_att(ncid, u_varid, UNITS, U_UNITS) )
    call check( nf90_put_att(ncid, v_varid, UNITS, V_UNITS) )

    ! End define mode.
    call check( nf90_enddef(ncid) )


    ! Write the coordinate variable data. This will put the latitudes
    ! and longitudes of our data grid into the netCDF file.
    call check( nf90_put_var(ncid, lat_varid, lats) )
    call check( nf90_put_var(ncid, lon_varid, lons) )

  ! These settings tell netcdf to write one timestep of data. (The
    ! setting of start(4) inside the loop below tells netCDF which
    ! timestep to write.)
    count = (/ e_we, e_ns, NLVLS, 1 /)
    start = (/ 1, 1, 1, 1 /)

  end subroutine open_ncdf




  subroutine write_ncdf(n, psi, zeta, u, v)

    implicit none
    integer, intent(in)     :: n
    real(dp), intent(in), allocatable  :: psi(:,:), &
                                          zeta(:,:), &
                                          u(:,:), &
                                          v(:,:)

    integer(kind=4)                       :: i, j


    start(NDIMS) = n
    call check( nf90_put_var(ncid, psi_varid, transpose(psi), start = start,&
                              count = count) )
    call check( nf90_put_var(ncid, vort_varid, transpose(zeta), start = start, &
                              count = count) )
    call check( nf90_put_var(ncid, u_varid, transpose(u), start = start,&
                              count = count) )
    call check( nf90_put_var(ncid, v_varid, transpose(v), start = start, &
                              count = count) )


!    write(*,*) start
!    !saida psi
!    do i=1, e_ns
!    write(*,*)
!      do j=1, e_we
!        write(*,'(F7.4)', advance='no') psi(i,j)
!      end do
!    end do
!    write(*,*)

  end subroutine write_ncdf



  subroutine close_ncdf()

    implicit none

    ! Close the file. This causes netCDF to flush all buffers and make
    ! sure your data are really written to disk.
    call check( nf90_close(ncid) )

  end subroutine close_ncdf



  subroutine check(status)

    implicit none
    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check
END MODULE module_io
