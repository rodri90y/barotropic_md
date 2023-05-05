MODULE module_solve

  use module_configure
  use module_io
  implicit none




CONTAINS


  subroutine psi_field(psi)
  !<descricao>
  ! Calcula campo de psi inicial: gaussiana psi = exp(-1*((x-x0)^2+(x-x0)^2)/sigma**)
  !<descricao>
    implicit none
    real(dp), intent(inout), allocatable :: psi(:,:)

    integer(kind=4)                      :: i, j

    do i=1, e_ns
      do j=1, e_we
        psi(i,j) = A0*exp(-1.*((j-e_we/2.)**2 + (i-e_ns/2.)**2)/sigma**2)
      end do
    end do

  end subroutine psi_field



  subroutine lap(psi, zeta, dx, dy)
  !<descricao>
  ! Calcula campo de zeta
  !<descricao>
    implicit none
    real(dp), intent(in), allocatable     :: psi(:,:)
    real(dp), intent(inout), allocatable  :: zeta(:,:)
    real, intent(in)                      :: dx, dy
    integer(kind=4)                       ::  i, j, im, ip, jm ,jp, im_i1, im_in

    do i=1, e_ns

      im_i1 = 1.
      im_in = 1.

      im = i-1
      ip = i+1

      if(i == 1) then
        im_i1 = 0
        im = i
      else if(i == e_ns) then
        im_in = 0
        ip = i
      end if

      do j=1, e_we

        jm = j-1
        jp = j+1

        if(j == 1) then
          jm = e_we-1
        else if(j == e_we) then
          jp = 2
        end if

        zeta(i,j) = 1./dy**2 * (psi(i,jp) + psi(i,jm) - 2*psi(i,j)) &
                  + 1./dx**2 * (psi(ip,j)*im_in + psi(im,j)*im_i1 - 2*psi(i,j))

      end do
    end do

  end subroutine lap


  subroutine jcb_arak1(psi, zeta, jac, dx, dy)
  !<descricao>
  ! Calcula Jacobiano de Arakawa
  !<descricao>
    implicit none
    real(dp), intent(in), allocatable     :: psi(:,:)
    real(dp), intent(inout), allocatable  :: zeta(:,:),&
                                             jac(:,:)
    real, intent(in)                      :: dx, dy
    real(dp)                              :: dm, J1, J2, J3
    integer(kind=4)                       :: i, j, im, ip, jm ,jp, im_i1, im_in

    dm = 4.*dy*dx

    do i=1, e_ns

      im_i1 = 1.
      im_in = 1.

      im = i-1
      ip = i+1

      if(i == 1) then
        im_i1 = 0       !cond. contorno nulo fora das fronteiras
        im = i
      else if(i == e_ns) then
        im_in = 0       !cond. contorno nulo fora das fronteiras
        ip = i
      end if

      do j=1, e_we

        jm = j-1
        jp = j+1

        if(j == 1) then
          jm = e_we-1
        else if(j == e_we) then
          jp = 2
        end if


        J1  = 1./dm * ( ((zeta(i,jp)-zeta(i,jm))*(psi(ip,j)*im_in-psi(im,j)*im_i1)) &
                  - ((zeta(ip,j)*im_in-zeta(im,j)*im_i1)*(psi(i,jp)-psi(i,jm))) )


        J2  = 1./dm * ( psi(i,jm)*(zeta(ip,jm)*im_in-zeta(im,jm)*im_i1) &
                  - psi(i,jp)*(zeta(ip,jp)*im_in-zeta(im,jp)*im_i1)*im_in &
                  + psi(ip,j)*(zeta(ip,jp)-zeta(ip,jm)*im_in) &
                  - psi(im,j)*(zeta(im,jp)-zeta(im,jm)*im_i1) )

        J3  = 1./dm * ( zeta(i,jp)*(psi(ip,jp)*im_in-psi(im,jp)*im_i1) &
                  - zeta(i,jm)*(psi(ip,jm)*im_in-psi(im,jm)*im_i1) &
                  - zeta(ip,j)*(psi(ip,jp)-psi(ip,jm))*im_in &
                  + zeta(im,j)*(psi(im,jp)-psi(im,jm))*im_i1 )

        jac(i,j) = (J1 + J2 + J3)/3.

      end do
    end do

  end subroutine jcb_arak1

  subroutine jcb_arak2(psi, zeta, jac, dx, dy)
  !<descricao>
  ! Calcula Jacobiano de Arakawa nao conservativo
  !<descricao>
    implicit none
    real(dp), intent(in), allocatable     :: psi(:,:)
    real(dp), intent(inout), allocatable  :: zeta(:,:),&
                                             jac(:,:)
    real, intent(in)                      :: dx, dy
    real(dp)                              :: dm, J1
    integer(kind=4)                       :: i, j, im, ip, jm ,jp, im_i1, im_in

    dm = 4.*dy*dx

    do i=1, e_ns

      im_i1 = 1.
      im_in = 1.

      im = i-1
      ip = i+1

      if(i == 1) then
        im_i1 = 0       !cond. contorno nulo fora das fronteiras
        im = i
      else if(i == e_ns) then
        im_in = 0       !cond. contorno nulo fora das fronteiras
        ip = i
      end if

      do j=1, e_we

        jm = j-1
        jp = j+1

        if(j == 1) then
          jm = e_we-1
        else if(j == e_we) then
          jp = 2
        end if


        jac(i,j)  = 1./dm * ( ((zeta(i,jp)-zeta(i,jm))*(psi(ip,j)*im_in-psi(im,j)*im_i1)) &
                  - ((zeta(ip,j)*im_in-zeta(im,j)*im_i1)*(psi(i,jp)-psi(i,jm))) )

      end do
    end do

  end subroutine jcb_arak2

  subroutine relax(psi, zeta, dx, dy)
  !<descricao>
  ! Calcula gaus-sidel relaxation
  !<descricao>
    implicit none
    real(dp), intent(inout), allocatable  :: psi(:,:)
    real(dp), intent(in), allocatable     :: zeta(:,:)

    real(dp), allocatable                 :: psi_n(:,:)
    real(dp)                              :: err_max, err
    real, intent(in)                      :: dx, dy
    integer(kind=4)                       :: i, j, im, ip, jm ,jp, im_i1, im_in
    integer(kind=4)                       :: int_count = 0, int_max = 5000

    allocate(psi_n(e_ns, e_we))

    err = 0
    int_count = 1

    err_max = 1.E-8                ! set it between 10-5 to 10-7 (similar about 10-4 order) best 10-6

    do
      do i=1, e_ns

        im_i1 = 1.
        im_in = 1.

        im = i-1
        ip = i+1

        if(i == 1) then
          im_i1 = 0
          im = i
        else if(i == e_ns) then
          im_in = 0
          ip = i
        end if

        do j=1, e_we

          jm = j-1
          jp = j+1

          if(j == 1) then
            jm = e_we-1
          else if(j == e_we) then
            jp = 2
          end if

          psi_n(i,j) = psi(i,j)*(1.-alpha) + alpha/4. * (dx**2 * zeta(i,j) &         ! Randall
                     + psi_n(i,jm)+psi(i,jp)+psi_n(im,j)*im_i1+psi(ip,j)*im_in)

        end do
      end do

      err = maxval(abs(psi_n-psi))/maxval(abs(psi))

      if(err < err_max .or. int_count>=int_max) exit

      int_count = int_count + 1
      psi = psi_n

    end do

    psi = -1.*psi
    !write(*,*) int_count, err

  end subroutine relax

  subroutine adv_f(advf, psi, dx)
  !<descricao>
  ! Calcula adveccao: Beta * d(psi)/dx
  !<descricao>
    implicit none
    real(dp), intent(in), allocatable     :: psi(:,:)
    real(dp), intent(inout), allocatable  :: advf(:,:)
    real, intent(in)                      :: dx
    integer(kind=4)                       :: i, j, jm, jp

    do i=1, e_ns
      do j=1, e_we

        jm = j-1
        jp = j+1

        if(j == 1) then
          jm = e_we-1
        else if(j == e_we) then
          jp = 2
        end if

        advf(i,j) = beta/(2.*dx) * (psi(i,jp) - psi(i,jm))  !tem inversao no sinal apos l/2 (x)

      end do
    end do

    end subroutine adv_f


    subroutine psi2uv(psi, u, v, dx, dy)
  !<descricao>
  ! Calcula campo de u e a partir de psi
  !<descricao>
    implicit none
    real(dp), intent(in), allocatable     :: psi(:,:)
    real(dp), intent(inout), allocatable  :: u(:,:), &
                                             v(:,:)
    real, intent(in)                      :: dx, dy
    integer(kind=4)                       :: i, j, im, ip, jm ,jp, im_i1, im_in

    do i=1, e_ns

      im_i1 = 1.
      im_in = 1.

      im = i-1
      ip = i+1

      if(i == 1) then
        im_i1 = 0
        im = i
      else if(i == e_ns) then
        im_in = 0
        ip = i
      end if

      do j=1, e_we

        jm = j-1
        jp = j+1

        if(j == 1) then
          jm = e_we-1
        else if(j == e_we) then
          jp = 2
        end if

        u(i,j) = -1.*( psi(i,jp) - psi(i,jm) )/(2.*dy)
        v(i,j) = ( psi(ip,j)*im_in - psi(im,j)*im_i1 )/(2.*dx)

      end do
    end do

  end subroutine psi2uv


    subroutine psi_prog(psi, dx, dy, dt)
  !<descricao>
  ! Calcula conservação da vorticidade absoluta
  !<descricao>
    implicit none
    real(dp), intent(inout), allocatable  :: psi(:,:)
    real(dp), allocatable                 :: zeta_3d(:,:,:), &
                                             zeta_2d(:,:), &
                                             jac(:,:), &
                                             advf(:,:), &
                                             u(:,:), &
                                             v(:,:)
    real, intent(in)                      :: dt, dx, dy
    integer(kind=4)                       :: i, j, t
    integer(kind=4)                       :: error


    allocate(zeta_3d(e_ns, e_we,3),zeta_2d(e_ns, e_we),stat=error)
    allocate(jac(e_ns, e_we),stat=error)
    allocate(advf(e_ns, e_we),stat=error)

    allocate(u(e_ns, e_we),v(e_ns, e_we),stat=error)

    if (error /= 0) stop "*** Not enough memory ***"

    call open_ncdf()

    do t=0, nint


      call lap(psi, zeta_2d, dx, dy)                       ! initialize vorticity field (zeta) from psi fiel
                                                           ! calculate jacobian - advection of relative vorticity
      select case (jactyp)
      case (1)
       call jcb_arak1(psi, zeta_2d, jac, dx, dy)
      case (2)
        call jcb_arak2(psi, zeta_2d, jac, dx, dy)
      end select


      call adv_f(advf, psi, dx)                             ! calculate advection of planetary vortivit


      !time integration
      if(t==0) then
        zeta_3d(:,:,1) = zeta_2d

        zeta_3d(:,:,2) = zeta_3d(:,:,1) - 1.*dt*(jac + advf)
        zeta_2d = zeta_3d(:,:,2)
      else
        zeta_3d(:,:,3) = zeta_3d(:,:,1) - 2.*dt*(jac + advf)

        zeta_2d = zeta_3d(:,:,3)

        zeta_3d = CSHIFT(zeta_3d, shift=1, dim=3)
      end if

      call relax(psi, zeta_2d, dx, dy)                    ! calculate new stream function (psi) field from zeta field

      call psi2uv(psi, u, v, dx, dy)

      if (mod(t,io_dt) == 0) then
        print*,t+1
        call write_ncdf((1+t/io_dt), psi, zeta_2d, u, v)
      end if

    end do

    call close_ncdf()

  end subroutine psi_prog

END MODULE module_solve
