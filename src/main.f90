program main
    use iso_fortran_env, only: real64, int64
    use iso_c_binding
    use mod_rhpix_utils, only: pi, deg2rad, rad2deg, wrap_longitude, wrap_latitude, auth_lat, auth_rad
    use mod_rhpix_ellipsoids, only: ellipsoid, init_ellipsoid, WGS84_A, WGS84_F, WGS84_B, WGS84_E, WGS84_R_A, R_EM
    implicit none

    real(kind=real64) :: wrapped
    real(kind=real64) :: beta
    integer :: i

    real(8),  parameter :: pi_8  = 4 * atan (1.0_8)
    real(16), parameter :: pi_16 = 4 * atan (1.0_16)
    real(real64), parameter :: pi_64 = 4 * atan (1.0_real64)

    ! testing using ellipsoid
    real(kind=real64), parameter :: oo = 0.0_real64
    type(ellipsoid), parameter :: param = ellipsoid(0, 0, 0.0_real64)
    type(ellipsoid) :: wgs84_base

    type(ellipsoid) :: WGS84_ELLIPSOID
    type(ellipsoid) :: WGS84_ELLIPSOID_RADIANS
    type(ellipsoid) :: WGS84_ASPHERE
    type(ellipsoid) :: WGS84_ASPHERE_RADIANS
    type(ellipsoid) :: EMR_SPHERE
    type(ellipsoid) :: EMR_SPHERE_RADIANS
    type(ellipsoid) :: UNIT_SPHERE
    type(ellipsoid) :: UNIT_SPHERE_RADIANS

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

    ! >>> beta = auth_lat(pi/4, 0.5, radians=True)
    !     >>> print(my_round(beta, 15))
    !     0.68951821243544
    beta = auth_lat(pi_64/4, 0.5_real64, 0, 1)
    print *, 'auth_lat', pi_64/4, beta

    beta = auth_lat(pi_64/4, 0.0_real64, 0, 1)
    print *, 'auth_lat e=0.0', beta

    ! >>> auth_rad(1, 0)
    !     1
    !     >>> for i in range(2, 11):
    !     ...     e = 1.0/i**2
    !     ...     print(my_round((e, auth_rad(1, 1.0/i**2)), 15))
    beta = auth_rad(1.0_real64, 0.0_real64, 0)
    print *, 'auth_rad 1, 0: ', beta

    do i = 2, 10, 1
        beta = 1.0_real64/i**2
        print *, 'auth_rad i', i, beta, auth_rad(1.0_real64, beta, 0)
    end do

    print *, 'ellipsoid', param

    wgs84_base = param
    call init_ellipsoid(wgs84_base)
    print *, 'wgs84_base', wgs84_base

    ! `is_sphere` 1 `R` radius 2 `a` major radius 3 `b` minor radius of the ellipsoid in meters
    ! 4 `e` eccentricity 5 `f` Flattening 6 `R_A` authalic radius meters
    ! 7 `lon_0` Central meridian 8 `lat_0` Latitude of origin 10 `is_radians` 9 `phi_0` latitude separating equatorial and north polar region in (r)HEALPix
    ! # Define some common ellipsoids.
    ! WGS84_ELLIPSOID = Ellipsoid(a=WGS84_A, f=WGS84_F)
    WGS84_ELLIPSOID%is_sphere = 0
    WGS84_ELLIPSOID%is_radians = 0
    WGS84_ELLIPSOID%el_params = oo
    WGS84_ELLIPSOID%el_params(2) = WGS84_A
    WGS84_ELLIPSOID%el_params(5) = WGS84_F
    call init_ellipsoid(WGS84_ELLIPSOID)
    print *, 'WGS84_ELLIPSOID', WGS84_ELLIPSOID

    ! WGS84_ELLIPSOID_RADIANS = Ellipsoid(a=WGS84_A, f=WGS84_F, radians=True)
    WGS84_ELLIPSOID_RADIANS%is_sphere = 0
    WGS84_ELLIPSOID_RADIANS%is_radians = 1
    WGS84_ELLIPSOID_RADIANS%el_params = oo
    WGS84_ELLIPSOID_RADIANS%el_params(2) = WGS84_A
    WGS84_ELLIPSOID_RADIANS%el_params(5) = WGS84_F
    call init_ellipsoid(WGS84_ELLIPSOID_RADIANS)
    print *, 'WGS84_ELLIPSOID_RADIANS', WGS84_ELLIPSOID_RADIANS

    ! WGS84_ASPHERE = Ellipsoid(R=WGS84_R_A)
    WGS84_ASPHERE%is_sphere = 1
    WGS84_ASPHERE%is_radians = 0
    WGS84_ASPHERE%el_params = oo
    WGS84_ASPHERE%el_params(1) = WGS84_R_A
    call init_ellipsoid(WGS84_ASPHERE)
    print *, 'WGS84_ASPHERE', WGS84_ASPHERE

    ! WGS84_ASPHERE_RADIANS = Ellipsoid(R=WGS84_R_A, radians=True)
    WGS84_ASPHERE_RADIANS%is_sphere = 1
    WGS84_ASPHERE_RADIANS%is_radians = 1
    WGS84_ASPHERE_RADIANS%el_params = oo
    WGS84_ASPHERE_RADIANS%el_params(1) = WGS84_R_A
    call init_ellipsoid(WGS84_ASPHERE_RADIANS)
    print *, 'WGS84_ASPHERE_RADIANS', WGS84_ASPHERE_RADIANS

    ! EMR_SPHERE = Ellipsoid(R=R_EM)
    EMR_SPHERE%is_sphere = 1
    EMR_SPHERE%is_radians = 0
    EMR_SPHERE%el_params = oo
    EMR_SPHERE%el_params(1) = R_EM
    call init_ellipsoid(EMR_SPHERE)
    print *, 'EMR_SPHERE', EMR_SPHERE

    ! EMR_SPHERE_RADIANS = Ellipsoid(R=R_EM, radians=True)
    EMR_SPHERE_RADIANS%is_sphere = 1
    EMR_SPHERE_RADIANS%is_radians = 1
    EMR_SPHERE_RADIANS%el_params = oo
    EMR_SPHERE_RADIANS%el_params(1) = R_EM
    call init_ellipsoid(EMR_SPHERE_RADIANS)
    print *, 'EMR_SPHERE_RADIANS', EMR_SPHERE_RADIANS

    ! UNIT_SPHERE = Ellipsoid(R=1)
    UNIT_SPHERE%is_sphere = 1
    UNIT_SPHERE%is_radians = 0
    UNIT_SPHERE%el_params = oo
    UNIT_SPHERE%el_params(1) = 1.0_real64
    call init_ellipsoid(UNIT_SPHERE)
    print *, 'UNIT_SPHERE', UNIT_SPHERE

    ! UNIT_SPHERE_RADIANS = Ellipsoid(R=1, radians=True)
    UNIT_SPHERE_RADIANS%is_sphere = 1
    UNIT_SPHERE_RADIANS%is_radians = 1
    UNIT_SPHERE_RADIANS%el_params = oo
    UNIT_SPHERE_RADIANS%el_params(1) = 1.0_real64
    call init_ellipsoid(UNIT_SPHERE_RADIANS)
    print *, 'UNIT_SPHERE_RADIANS', UNIT_SPHERE_RADIANS


end program main
