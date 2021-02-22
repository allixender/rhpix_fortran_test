module mod_rhpix
    use iso_fortran_env, only: real64, int64
    use iso_c_binding
    implicit none
    real(kind=real64), parameter :: pi    = 3.1415926535897932384626433832795_real64
    real(real64), parameter :: pi_64 = 4 * atan (1.0_real64)

    private

    public :: pi, deg2rad, rad2deg, wrap_longitude, wrap_latitude

contains

    real(real64) pure function deg2rad(degrees)
        real(real64), intent(in) :: degrees
        deg2rad = degrees * (pi_64/180.)
    end function deg2rad

    real(real64) pure function rad2deg(radians)
        real(real64), intent(in) :: radians
        rad2deg = radians * (180./pi_64)
    end function

    subroutine wrap_longitude(lam, radians) bind(c, name='c_wrap_longitude')
        real(C_DOUBLE), intent(in out) :: lam
        integer(c_int), intent(in) :: radians
        real(C_DOUBLE) :: res
        real(C_DOUBLE) :: t_lam

        ! Given a point p on the unit circle at angle `lam` from the positive
        ! x-axis, return its angle theta in the range -pi <= theta < pi.
        ! If `radians` = True, then `lam` and the output are given in radians.
        ! Otherwise, they are given in degrees.
        !
        ! EXAMPLES::
        !
        !     >>> wrap_longitude(2*pi + pi, radians=True)
        !     -3.141592653589793
        !
        ! NOTES:: .. Issue #1 was ..
        !     -3.1415926535897931

        if (radians <= 0) then
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

        if (radians <= 0) then
            ! Convert to degrees.
            res = rad2deg(res)
        end if
        lam = res
    end subroutine


    function wrap_latitude(phi, radians) bind(c, name='c_wrap_latitude') result(res)

        real(C_DOUBLE), intent(in) :: phi
        integer(c_int), intent(in) :: radians
        real(C_DOUBLE) :: res
        real(C_DOUBLE) :: t_phi

        ! Given a point p on the unit circle at angle `phi` from the positive x-axis,
        ! if p lies in the right half of the circle, then return its angle that lies
        ! in the interval [-pi/2, pi/2].
        ! If p lies in the left half of the circle, then reflect it through the
        ! origin, and return the angle of the reflected point that lies in the
        ! interval [-pi/2, pi/2].
        ! If `radians` = True, then `phi` and the output are given in radians.
        ! Otherwise, they are given in degrees.
        !
        ! EXAMPLES::
        !
        !     >>> wrap_latitude(45, radians=False)
        !     45.0

        if (radians <= 0) then
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

        if (radians <= 0) then
            !# Convert to degrees.
            res = rad2deg(res)
        end if
    end function

end module mod_rhpix
