program main
    use iso_fortran_env, only: real64, int64
    use iso_c_binding
    use mod_rhpix, only: pi, deg2rad, rad2deg, wrap_longitude, wrap_latitude
    implicit none

    real(kind=real64) :: wrapped

    real(8),  parameter :: pi_8  = 4 * atan (1.0_8)
    real(16), parameter :: pi_16 = 4 * atan (1.0_16)
    real(real64), parameter :: pi_64 = 4 * atan (1.0_real64)

    print *, 'fortran_quick2 library demo'
    write (*,*) pi_8
    write (*,*) pi_16
    write (*,*) pi
    write (*,*) pi_64

    print *, "Hello, rhpix!"

    ! >>> wrap_longitude(-185, radians=False)
    ! 175.0
    ! wrapped = wrap_longitude(-185.0_real64, .false.)
    wrapped = -185.0_real64
    call wrap_longitude(wrapped, 0)
    print *, 'wrap_longitude', -185.0_real64, wrapped
    ! >>> wrap_longitude(-180, radians=False)
    ! -180.0
    ! wrapped = wrap_longitude(-180.0_real64, .false.)
    !
    ! print *, 'wrap_longitude', -180.0_real64, wrapped
    ! ! >>> wrap_longitude(185, radians=False)
    ! ! -175.0
    ! wrapped = wrap_longitude(185.0_real64, .false.)
    !
    ! print *, 'wrap_longitude', 185.0_real64, wrapped
    !
    ! wrapped = wrap_longitude(180.0_real64, .false.)
    !
    ! print *, 'wrap_longitude', 180.0_real64, wrapped

    ! >>> wrap_latitude(-45, radians=False)
    ! -45.0
    wrapped = wrap_latitude(-45.0_real64, 0)
    print *, 'wrap_latitude', -45.0_real64, wrapped
    ! >>> wrap_latitude(90, radians=False)
    ! 90.0
    wrapped = wrap_latitude(90.0_real64, 0)
    print *, 'wrap_latitude', 90.0_real64, wrapped
    ! >>> wrap_latitude(-90, radians=False)
    ! -90.0
    wrapped = wrap_latitude(-90.0_real64, 0)
    print *, 'wrap_latitude', -90.0_real64, wrapped
    ! >>> wrap_latitude(135, radians=False)
    ! -45.0
    wrapped = wrap_latitude(135.0_real64, 0)
    print *, 'wrap_latitude', 135.0_real64, wrapped
    ! >>> wrap_latitude(-135, radians=False)
    ! 45.0
    wrapped = wrap_latitude(-135.0_real64, 0)
    print *, 'wrap_latitude', -135.0_real64, wrapped

end program main
