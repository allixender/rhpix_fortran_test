module mod_rhpix_ellipsoids
    use iso_fortran_env, only: real64, int64
    use iso_c_binding
    use mod_rhpix_utils, only: pi, deg2rad, rad2deg, auth_lat, auth_rad
    implicit none

    ! Parameters of some common ellipsoids.
    real(real64), parameter :: WGS84_A = 6378137.0_real64
    ! # ORIGINAL: GRS80 from EPSG:42310 298.257222101, based on WGS84+GRS80
    real(real64), parameter :: WGS84_F = 1.0_real64 / 298.257222101_real64
    ! WGS84_F = 1 / 298.257222100882711  # GRS80 from https://en.wikipedia.org/wiki/World_Geodetic_System
    ! WGS84_F = 1 / 298.257223563  # new value from EPSG:7030
    real(real64), parameter :: WGS84_B = WGS84_A * (1.0_real64 - WGS84_F)
    real(real64), parameter :: WGS84_E = sqrt(WGS84_F * (1.0_real64 - WGS84_F))
    real(real64), parameter :: WGS84_R_A = sqrt(WGS84_A ** 2 / 2.0_real64 + WGS84_B ** 2 / 2.0_real64 * (atanh(WGS84_E) / WGS84_E))
    ! # Earth's mean radius
    real(real64), parameter :: R_EM = 6371000.0_real64

    ! setup declarations, treat the ellipsoid class as a fixed length array of values

    ! Represents an ellipsoid of revolution (possibly a sphere) with a
    ! geodetic longitude-latitude coordinate frame.
    !
    ! INSTANCE ATTRIBUTES:
    !
    ! is_sphere - True (1) if the ellipsoid is a sphere, and False (0) otherwise.
    ! 1 - `R` - The radius of the ellipsoid in meters, implying that it is a
    !   sphere.
    ! 2 - `a` - Major radius of the ellipsoid in meters.
    ! 3 - `b` - Minor radius of the ellipsoid in meters.
    ! 4 - `e` - Eccentricity of the ellipsoid.
    ! 5 - `f` - Flattening of the ellipsoid.
    ! 6 - `R_A` - Authalic radius of the ellipsoid in meters.
    ! 7 - `lon_0` - Central meridian.
    ! 8 - `lat_0` - Latitude of origin.
    ! is_radians - If True (1), use angles measured in radians for all calculations.
    !   else false (0) Use degrees otherwise.
    ! 9 - `phi_0` - The latitude separating the equatorial region and
    !   the north polar region in the context of the (r)HEALPix projection.
    !
    ! Except for phi_0, these attribute names match the names of the
    ! `PROJ.4 ellipsoid parameters
    ! Fortran array index start with 1

    type ellipsoid
        integer(c_int) :: is_sphere
        integer(c_int) :: is_radians
        real(real64), dimension(9) :: el_params
    end type ellipsoid

    ! # Define some common ellipsoids.
    ! WGS84_ELLIPSOID = Ellipsoid(a=WGS84_A, f=WGS84_F)
    ! WGS84_ELLIPSOID_RADIANS = Ellipsoid(a=WGS84_A, f=WGS84_F, is_radians=True)
    ! WGS84_ASPHERE = Ellipsoid(R=WGS84_R_A)
    ! WGS84_ASPHERE_RADIANS = Ellipsoid(R=WGS84_R_A, is_radians=True)
    ! EMR_SPHERE = Ellipsoid(R=R_EM)
    ! EMR_SPHERE_RADIANS = Ellipsoid(R=R_EM, is_radians=True)
    ! UNIT_SPHERE = Ellipsoid(R=1)
    ! UNIT_SPHERE_RADIANS = Ellipsoid(R=1, is_radians=True)

    private

    public :: WGS84_A, WGS84_F, WGS84_B, WGS84_E, WGS84_R_A, R_EM, ellipsoid, init_ellipsoid

contains

    ! class methods implement as procedure over the same type ellipsoid

    ! is_sphere` 1 `R` radius 2 `a` major radius 3 `b` minor radius of the ellipsoid in meters
    ! 4 `e` eccentricity 5 `f` Flattening 6 `R_A` authalic radius meters
    ! 7 `lon_0` Central meridian 8 `lat_0` Latitude of origin `is_radians` 9 `phi_0` latitude separating equatorial and north polar region in (r)HEALPix
    subroutine init_ellipsoid(params)
        type(ellipsoid), intent(in out) :: params

        params%el_params(1) = 1.0_real64

        if (params%is_sphere == 1) then
            ! The ellipsoid is a sphere.
            ! Override the other geometric parameters.
            params%el_params(2) = params%el_params(1)
            params%el_params(3) = params%el_params(1)
            params%el_params(4) = 0
            params%el_params(5) = 0
            params%el_params(6) = params%el_params(1)
        else

            if (params%el_params(3) > 0) then
                ! # Derive the other geometric parameters from a and b.

                params%el_params(4) = sqrt(1 - (params%el_params(3) / params%el_params(2)) ** 2)
                params%el_params(5) = (params%el_params(2) - params%el_params(3)) / params%el_params(2)
            else if (params%el_params(4) > 0) then
                ! # Derive the other geometric parameters from a and e.

                params%el_params(3) = params%el_params(2) * sqrt(1 - params%el_params(4) ** 2)
                params%el_params(5) = 1 - sqrt(1 - params%el_params(4) ** 2)
            else
                params%el_params(3) = params%el_params(4) * (1 - params%el_params(5))
                params%el_params(4) = sqrt(params%el_params(5) * (1 - params%el_params(5)))
            end if
            params%el_params(6)= auth_rad(params%el_params(2), params%el_params(4), 0)
        end if
        params%el_params(9) = auth_lat(asin(2.0_real64 / 3.0_real64), params%el_params(4), 1, 1)
        if (params%is_sphere <= 0) then
            ! Convert to degrees.
            params%el_params(9) = rad2deg(params%el_params(9))
        end if
    end subroutine




end module mod_rhpix_ellipsoids
