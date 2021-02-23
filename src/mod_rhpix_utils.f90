module mod_rhpix_utils
    use iso_fortran_env, only: real64, int64
    use iso_c_binding
    implicit none
    real(kind=real64), parameter :: pi    = 3.1415926535897932384626433832795_real64
    real(real64), parameter :: pi_64 = 4 * atan (1.0_real64)

    private

    public :: pi, deg2rad, rad2deg, wrap_longitude, wrap_latitude, auth_lat, auth_rad

contains

    real(real64) pure function deg2rad(degrees)
        real(real64), intent(in) :: degrees
        deg2rad = degrees * (pi_64/180.)
    end function deg2rad

    real(real64) pure function rad2deg(radians)
        real(real64), intent(in) :: radians
        rad2deg = radians * (180./pi_64)
    end function

    subroutine wrap_longitude(lam, is_radians) bind(c, name='c_wrap_longitude')
        real(C_DOUBLE), intent(in out) :: lam
        integer(c_int), intent(in) :: is_radians
        real(C_DOUBLE) :: res
        real(C_DOUBLE) :: t_lam

        ! Given a point p on the unit circle at angle `lam` from the positive
        ! x-axis, return its angle theta in the range -pi <= theta < pi.
        ! If `is_radians` = True, then `lam` and the output are given in radians.
        ! Otherwise, they are given in degrees.
        !
        ! EXAMPLES::
        !
        !     >>> wrap_longitude(2*pi + pi, is_radians=True)
        !     -3.141592653589793
        !
        ! NOTES:: .. Issue #1 was ..
        !     -3.1415926535897931

        if (is_radians <= 0) then
            t_lam = deg2rad(lam)
        else
            t_lam = lam
        end if

        if ((t_lam < -pi_64) .or. (t_lam >= pi_64)) then
            res = t_lam - 2 * pi_64 * floor(t_lam / (2 * pi_64))
            if (res >= pi_64) then
                res = res - 2 * pi_64
            end if
        else
            res = t_lam
        end if

        if (is_radians <= 0) then
            ! Convert to degrees.
            res = rad2deg(res)
        end if
        lam = res
    end subroutine


    function wrap_latitude(phi, is_radians) bind(c, name='c_wrap_latitude') result(res)

        real(C_DOUBLE), intent(in) :: phi
        integer(c_int), intent(in) :: is_radians
        real(C_DOUBLE) :: res
        real(C_DOUBLE) :: t_phi

        ! Given a point p on the unit circle at angle `phi` from the positive x-axis,
        ! if p lies in the right half of the circle, then return its angle that lies
        ! in the interval [-pi/2, pi/2].
        ! If p lies in the left half of the circle, then reflect it through the
        ! origin, and return the angle of the reflected point that lies in the
        ! interval [-pi/2, pi/2].
        ! If `is_radians` = True, then `phi` and the output are given in radians.
        ! Otherwise, they are given in degrees.
        !
        ! EXAMPLES::
        !
        !     >>> wrap_latitude(45, is_radians=False)
        !     45.0

        if (is_radians <= 0) then
            ! # Convert to radians.
            t_phi = deg2rad(phi)
        else
            t_phi = phi
        end if

        !# Put phi in range -pi <= phi < pi.
        ! t_phi =
        call wrap_longitude(t_phi, 1)
        if ( abs(t_phi) <= pi_64 / 2) then
            res = t_phi
        else
            res = t_phi - sign(1.0_real64, t_phi) * pi_64
        end if

        if (is_radians <= 0) then
            !# Convert to degrees.
            res = rad2deg(res)
        end if
    end function


    function auth_lat(phi, e, inverse, is_radians) bind(c, name='c_auth_lat') result(res)

        real(C_DOUBLE), intent(in) :: phi
        real(C_DOUBLE), intent(in) :: e
        integer(c_int), intent(in) :: is_radians
        integer(c_int), intent(in) :: inverse
        real(C_DOUBLE) :: res
        real(C_DOUBLE) :: t_phi
        real(C_DOUBLE) :: q
        real(C_DOUBLE) :: qp
        real(C_DOUBLE) :: ratio

        ! Given a point of geographic latitude `phi` on an ellipse of
        ! eccentricity `e`, return the authalic latitude of the point.
        ! If `inverse` =True, then compute its inverse approximately.
        !
        ! EXAMPLES::
        !
        !     >>> beta = auth_lat(pi/4, 0.5, is_radians=True)
        !     >>> print(my_round(beta, 15))
        !     0.68951821243544

        ! The power series approximation used for the inverse is
        ! standard in cartography (PROJ.4 uses it, for instance)
        ! and accurate for small eccentricities.

        if (e == 0.0_real64) then
            res = phi
            return
        end if

        if (is_radians <= 0.0_real64) then
            ! Convert to radians to do calculations below.
            t_phi = deg2rad(phi)
        else
            t_phi = phi
        end if

        if (inverse <= 0.0_real64) then
            ! Compute authalic latitude from latitude phi.
            q = ((1.0_real64 - e ** 2) * sin(t_phi)) / (1.0_real64 - (e * sin(t_phi)) ** 2) - (1.0_real64 - e ** 2) / &
                ( 2.0 * e ) * log((1.0_real64 - e * sin(t_phi)) / (1.0_real64 + e * sin(t_phi)))
            qp = 1.0_real64 - (1.0_real64 - e ** 2) / (2.0_real64 * e) * log((1.0_real64 - e) / (1.0_real64 + e))
            ratio = q / qp
            ! Avoid rounding errors.
            if (abs(ratio) > 1.0_real64) then
                ! Make abs(ratio) = 1
                ratio = sign(1.0_real64, ratio)
            end if
            res = asin(ratio)
        else
            ! Compute an approximation of latitude from authalic latitude phi.
            res = ( t_phi + (e ** 2 / 3.0 + 31 * e ** 4 / 180.0_real64 + 517 * e ** 6 / 5040.0_real64) * &
                sin(2 * t_phi) + (23 * e ** 4 / 360.0_real64 + 251 * e ** 6 / 3780.0_real64) * sin(4 * t_phi) + &
                (761 * e ** 6 / 45360.0_real64) * sin(6 * t_phi) )
        end if
        if (is_radians <= 0) then
            !# Convert to degrees.
            res = rad2deg(res)
        end if
    end function


    function auth_rad(a, e, inverse) bind(c, name='c_auth_rad') result(res)

        real(C_DOUBLE), intent(in) :: a
        real(C_DOUBLE), intent(in) :: e
        integer(c_int), intent(in) :: inverse
        real(C_DOUBLE) :: res
        real(C_DOUBLE) :: k

        ! Return the radius of the authalic sphere of the ellipsoid with major
        ! radius `a` and eccentricity `e`.
        ! If `inverse` = True, then return the major radius of the ellipsoid
        ! with authalic radius `a` and eccentricity `e`.
        !
        ! EXAMPLES::
        !
        !     >>> auth_rad(1, 0)
        !     1
        !     >>> for i in range(2, 11):
        !     ...     e = 1.0/i**2
        !     ...     print(my_round((e, auth_rad(1, 1.0/i**2)), 15))
        !     (0.25, 0.989393259670095)
        !     (0.111111111111111, 0.997935147429943)
        !     (0.0625, 0.999348236455825)
        !     (0.04, 0.99973321235361)
        !     (0.027777777777778, 0.99987137105188)
        !     (0.020408163265306, 0.999930576285614)
        !     (0.015625, 0.999959307080847)
        !     (0.012345679012346, 0.999974596271211)
        !     (0.01, 0.999983332861089)

        if (e == 0.0_real64) then
            res = a
            return
        end if

        k = sqrt(0.5_real64 * (1.0_real64 - (1.0_real64 - e ** 2) / (2.0_real64 * e) * log((1.0_real64 - e) / (1.0_real64 + e))))

        if (inverse <= 0.0_real64) then
            ! The expression below is undefined when e=0 (sphere),
            ! but its limit as e tends to 0 is a, as expected.
            res = a * k
        else
            ! Then a is the authalic radius and output major radius of ellipsoid.
            res = a / k
        end if
    end function

end module mod_rhpix_utils
