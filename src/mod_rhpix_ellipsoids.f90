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

    subroutine pi_ellipsoid(params)
      type(ellipsoid), intent(in out) :: params
        ! """
        ! Return pi if `self.radians` = True and 180 otherwise.
        ! """
        ! if self.radians:
        !     return pi
        ! else:
        !     return 180.0
        end subroutine

    ! self, lam_min=None, lam_max=None, phi_min=None, phi_max=None
    subroutine random_point_ellipsoid(params)
      type(ellipsoid), intent(in out) :: params
        ! """
        ! Return a point (given in geodetic coordinates) sampled uniformly at
        ! random from the section of this ellipsoid with longitude in the range
        ! `lam_min <= lam < lam_max` and latitude in the range
        ! `phi_min <= phi < phi_max`.
        ! But avoid the poles.
        !
        ! EXAMPLES::
        !
        !    >>> E = UNIT_SPHERE
        !    >>> print(E.random_point()) # doctest: +SKIP
        !    (-1.0999574573422948, 0.21029104897701129)
        !
        ! """
        ! PI = self.pi()
        ! if lam_min is None:
        !     lam_min = -PI
        ! if lam_max is None:
        !     lam_max = PI
        ! if phi_min is None:
        !     phi_min = -PI / 2
        ! if phi_max is None:
        !     phi_max = PI / 2
        ! if not self.radians:
        !     # Convert to radians.
        !     lam_min, lam_max, phi_min, phi_max = deg2rad(
        !         [lam_min, lam_max, phi_min, phi_max]
        !     )
        ! # Pick a longitude.
        ! while True:
        !     u = uniform(0, 1)
        !     lam = (lam_max - lam_min) * u + lam_min
        !     # Don't include lam_max.
        !     if lam < lam_max:
        !         # Success.
        !         break
        ! # Pick a latitude.
        ! delta = pi / 360
        ! while True:
        !     v = uniform(0, 1)
        !     if self.sphere:
        !         phi = arcsin((sin(phi_max) - sin(phi_min)) * v + sin(phi_min))
        !     else:
        !         # Sample from the authalic sphere.
        !         # The map from the ellipsoid to the authalic sphere is
        !         # an equiareal diffeomorphism.
        !         # So a uniform distribution on the authalic sphere gives
        !         # rise to a uniform distribution on the ellipsoid.
        !         beta0 = auth_lat(phi_min, e=self.e, radians=True)
        !         beta1 = auth_lat(phi_max, e=self.e, radians=True)
        !         beta = arcsin((sin(beta1) - sin(beta0)) * v + sin(beta0))
        !         phi = auth_lat(beta, e=self.e, radians=True, inverse=True)
        !     # Avoid the poles.
        !     if abs(phi) <= pi / 2 - delta:
        !         # Success.
        !         break
        ! if not self.radians:
        !     # Convert back to degrees.
        !     lam, phi = rad2deg([lam, phi])
        ! return lam, phi
        end subroutine

    ! (self, n=90
    subroutine lattice_ellipsoid(params)
      type(ellipsoid), intent(in out) :: params
        ! """
        ! Return a 2n x n square lattice of longitude-latitude points.
        !
        ! EXAMPLES::
        !
        !     >>> E = UNIT_SPHERE
        !     >>> for p in E.lattice(n=3):
        !     ...     print(p)
        !     (-150.0, -60.0)
        !     (-150.0, 0.0)
        !     (-150.0, 60.0)
        !     (-90.0, -60.0)
        !     (-90.0, 0.0)
        !     (-90.0, 60.0)
        !     (-30.0, -60.0)
        !     (-30.0, 0.0)
        !     (-30.0, 60.0)
        !     (30.0, -60.0)
        !     (30.0, 0.0)
        !     (30.0, 60.0)
        !     (90.0, -60.0)
        !     (90.0, 0.0)
        !     (90.0, 60.0)
        !     (150.0, -60.0)
        !     (150.0, 0.0)
        !     (150.0, 60.0)
        !
        ! """
        ! PI = self.pi()
        ! # Longitudinal and latitudinal spacing between points.
        ! delta = PI / n
        ! return [
        !     (-PI + delta * (0.5 + i), -PI / 2 + delta * (0.5 + j))
        !     for i in range(2 * n)
        !     for j in range(n)
        ! ]
        end subroutine

    ! (self, lam, n=200
    subroutine meridian_ellipsoid(params)
      type(ellipsoid), intent(in out) :: params
        ! """
        ! Return a list of `n` equispaced longitude-latitude
        ! points lying along the meridian of longitude `lam`.
        ! Avoid the poles.
        ! """
        ! PI = self.pi()
        ! delta = PI / n
        ! return [(lam, -PI / 2 + delta * (0.5 + i)) for i in range(n)]
        end subroutine

    ! (self, phi, n=200
    subroutine parallel_ellipsoid(params)
      type(ellipsoid), intent(in out) :: params
        ! """
        ! Return a list of `2*n` equispaced longitude-latitude
        ! points lying along the parallel of latitude `phi`.
        ! """
        ! PI = self.pi()
        ! delta = PI / n
        ! return [(-PI + delta * (0.5 + i), phi) for i in range(2 * n)]
        end subroutine

    ! (self, n=400, spacing=None
    subroutine graticule_ellipsoid(params)
      type(ellipsoid), intent(in out) :: params
        ! """
        ! Return a list of longitude-latitude points sampled from a
        ! longitude-latitude graticule on this ellipsoid with the given
        ! spacing between meridians and between parallels.
        ! The number of points on longitude and latitude per pi radians is `n`.
        ! The spacing should be specified in the angle units used for this
        ! ellipsoid.
        ! If `spacing=None`, then a default spacing of pi/16 radians will be set.
        !
        ! EXAMPLES::
        !
        !     >>> E = UNIT_SPHERE
        !     >>> print(len(E.graticule(n=400)))
        !     25600
        !
        ! """
        ! PI = self.pi()
        ! result = []
        ! # delta = PI/n
        ! # Set default spacing.
        ! if spacing is None:
        !     spacing = PI / 16
        ! # Longitude lines.
        ! lam = -PI
        ! while lam < PI:
        !     # result.extend([(lam, -PI/2 + delta*(0.5 + i)) for i in range(n)])
        !     result.extend(self.meridian(lam, n))
        !     lam += spacing
        ! # Latitude lines. Avoid the poles.
        ! eps = PI / 360
        ! phi = -PI / 2 + eps
        ! while phi < PI / 2:
        !     # result.extend([(-PI + delta*(0.5 + i), phi) for i in range(2*n)])
        !     result.extend(self.parallel(phi, n))
        !     phi += spacing
        ! return result
        end subroutine

    ! (self, filename
    subroutine get_points_ellipsoid(params)
      type(ellipsoid), intent(in out) :: params
        ! """
        ! Return a list of longitude-latitude points contained in
        ! the file with filename `filename`.
        ! Assume the file is a text file containing at most one
        ! longitude-latitude point per line with the coordinates separated by
        ! whitespace and angles given in degrees.
        ! """
        ! result = []
        ! for line in open(filename, "rb"):
        !     if line[0] not in ["-", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
        !         # Ignore line.
        !         continue
        !     else:
        !         # Split coordinate pair on whitespace.
        !         p = [float(x) for x in line.split()]
        !         result.append(p)
        ! if self.radians:
        !     # Convert to radians.
        !     result = [deg2rad(p) for p in result]
        ! return result
        end subroutine

    ! (self, lam, phi
    subroutine xyz_ellipsoid(params)
      type(ellipsoid), intent(in out) :: params
        ! """
        ! Given a point on this ellipsoid with longitude-latitude coordinates
        ! `(lam, phi)`, return the point's 3D rectangular coordinates.
        !
        ! EXAMPLES::
        !
        !     >>> E = UNIT_SPHERE
        !     >>> print(my_round(E.xyz(0, 45), 15))
        !     (0.707106781186548, 0.0, 0.707106781186548)
        !
        ! NOTES:: .. Issue #1 was ..
        !     (0.70710678118654802, 0.0, 0.70710678118654802)
        !
        ! """
        end subroutine


end module mod_rhpix_ellipsoids
