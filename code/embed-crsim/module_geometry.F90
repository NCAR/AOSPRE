module module_geometry
  use module_configuration, only : RKIND
  use module_llxy, only : RAD_PER_DEG
  use module_llxy, only : DEG_PER_RAD
  
  implicit none

contains

  subroutine coordinate_transformation ( coordinate_system_type, tilt_degrees, rotation_degrees, heading_degrees, roll_degrees, pitch_degrees, &
       &                                 bwdeg1, ref_angle, azimuth, elevation, x, y, z, bwuse_h, bwuse_v )

    !
    ! PURPOSE:
    !
    !   Convert beam angles from aircraft-relative tilt and rotation, to ground-relative azimuth and elevation.
    !
    !   This routine implements equations 15, 16, 17 from Lee et al. (1994), omitting the
    !   gate range term ("r" from Lee94) which drops out when we compute azimuth and elevation.
    !
    ! INPUT:
    !
    !    coordinate_system_type:  As in CfRadial documentation on sensor type: X, Y, Y-Prime, Z
    !
    !    ref_angle:        Each panel has it's own orientation relative to the aircraft; this is used
    !                      to compute the beamwidth in the correct context to the orientation of each
    !                      panel.  Units of degrees.
    !   
    !    tilt_degrees
    !                      Angle between radar beam (when it is in a plane containing the 
    !                      longitudinal axis of the platform) and a line perpendicular to
    !                      the longitudinal axis.  Zero is perpendicular to the longitudinal
    !                      axis, positive is towards the front of the platform.
    !                      Units of degrees.
    !
    !    rotation_degrees
    !                      Angle between the radar beam and the vertical axis of the platform.
    !                      Zero is along the vertical axis, positive is clockwise looking 
    !                      forward from behind the platform.  Units of degrees.
    !  
    !    heading_degrees
    !                      Heading of the platform relative to true north, looking down from above.
    !                      Units of degrees.
    !
    !    roll_degrees
    !                      Roll about the longitudinal axis of the platform.  Postive is
    !                      left side up, looking forward.  Units of degrees.
    !
    !    pitch_degrees
    !                      Pitch about the lateral axis of the platform.  Positive is up 
    !                      at the front.  Units of degrees.
    !
    !    bwdeg1            Initial broadside beamwidth; the same for h and v.
    !
    !
    ! OUTPUT:
    !
    !    azimuth
    !                      Azimuth of antenna, relative to true north.  The azimuth should refer
    !                      to the center of the dwell.  Units of degrees.
    !
    !    elevation
    !                      Elevation of antenna, relative to the horizontal plane.  The elevation
    !                      should refer to the center of the dwell.  Units of degrees.
    !
    !    x
    !                      West-east component of the unit vector defined by azimuth and elevation.
    !
    !    y
    !                      South-north component of the unit vector defined by azimuth and elevation.
    !
    !    z
    !                      Vertical component of the unit vector defined by azimuth and elevation.
    !
    !    bwuse_h           Updated beamwidth - horizontal polarization
    !
    !    bwuse_v           Updated beamwidth - vertical polarization
    !
    !
    ! LOCAL:
    !
    !    tauA
    !                      Tilt term in radians
    !
    !    thetaA
    !                      Rotation term in radians
    !
    !    H
    !                      Heading term in radians
    !
    !    R
    !                      Roll term in radians
    !
    !    P
    !                      Pitch term in radians
    !
    !    RFA
    !                      Reference angle in radians
    !
    !
    ! INHERITED FROM MODULE:
    !
    !    RAD_PER_DEG
    !                      Multiplication factor to convert from degrees to radians.
    !
    !    DEG_PER_RAD
    !                      Multiplication factor to convert from radians to degrees.
    !

    implicit none

    ! Input:
    character(len=*), intent(in)  :: coordinate_system_type
    real(kind=RKIND), intent(in)  :: ref_angle
    real(kind=RKIND), intent(in)  :: tilt_degrees
    real(kind=RKIND), intent(in)  :: rotation_degrees
    real(kind=RKIND), intent(in)  :: heading_degrees
    real(kind=RKIND), intent(in)  :: roll_degrees
    real(kind=RKIND), intent(in)  :: pitch_degrees
    real(kind=RKIND), intent(in)  :: bwdeg1

    ! Output:
    real(kind=RKIND), intent(out) :: azimuth
    real(kind=RKIND), intent(out) :: elevation
    real(kind=RKIND), intent(out) :: x
    real(kind=RKIND), intent(out) :: y
    real(kind=RKIND), intent(out) :: z
    real(kind=RKIND), intent(out) :: bwuse_h
    real(kind=RKIND), intent(out) :: bwuse_v

    ! Local:
    real(kind=RKIND)              :: tauA
    real(kind=RKIND)              :: thetaA
    real(kind=RKIND)              :: H
    real(kind=RKIND)              :: R
    real(kind=RKIND)              :: P
    real(kind=RKIND)              :: RFA
    
    real(kind=RKIND), dimension(3,3) :: M      ! Transformation Matrix
    real(kind=RKIND), dimension(3)   :: Xa     ! Aircraft-relative x/y/z coordinates
    real(kind=RKIND), dimension(3)   :: Xearth ! Earth-relative x/y/z coordinates

    tauA   = tilt_degrees     * RAD_PER_DEG
    thetaA = rotation_degrees * RAD_PER_DEG
    H      = heading_degrees  * RAD_PER_DEG
    R      = roll_degrees     * RAD_PER_DEG
    P      = pitch_degrees    * RAD_PER_DEG
    RFA    = ref_angle        * RAD_PER_DEG

#ifdef _OLD_TO_BE_REMOVED_
    x = -cos(thetaA+R) * sin(H)*cos(tauA)*sin(P) + cos(H)*sin(thetaA+R)*cos(tauA) + sin(H)*cos(P)*sin(tauA)
    y = -cos(thetaA+R) * cos(H)*cos(tauA)*sin(P) - sin(H)*sin(thetaA+R)*cos(tauA) + cos(P)*cos(H)*sin(tauA)
    z = cos(P)*cos(tauA)*cos(thetaA+R) + sin(P)*sin(tauA)

    azimuth = atan2(x,y) * DEG_PER_RAD
    elevation = asin(z)  * DEG_PER_RAD
#endif

    !
    ! Transformation Matrix M:
    !

    M(1,1) = cos(H)*cos(R) + sin(H)*sin(P)*sin(R)
    M(1,2) = sin(H)*cos(P)
    M(1,3) = cos(H)*sin(R) - sin(H)*sin(P)*cos(R)

    M(2,1) = -sin(H)*cos(R)+cos(H)*sin(P)*sin(R)
    M(2,2) = cos(H)*cos(P)
    M(2,3) = -sin(H)*sin(R)-cos(H)*sin(P)*cos(R)

    M(3,1) = -cos(P)*sin(R)
    M(3,2) = sin(P)
    M(3,3) = cos(P)*cos(R)

    select case (coordinate_system_type)
    case default
        write(*,'("coordinate_system_type = ''", A,"''")') trim(coordinate_system_type)
        stop "UNRECOGNIZED COORDINATE_SYSTEM_TYPE"
    case ("X","x")
        Xa = (/ sin(tauA) , sin(thetaA+R)*cos(tauA) , cos(thetaA+R)*cos(tauA) /)
        bwuse_v = bwdeg1 + abs(sin(abs(thetaA)))
        bwuse_h = bwdeg1 + abs(sin(abs(tauA-RFA)))
    case ("Y","y")
        Xa = (/ cos(thetaA+R)*cos(tauA) , sin(tauA) , sin(thetaA+R)*cos(tauA) /)
        bwuse_h = bwdeg1 + abs(sin(abs(tauA - RFA - (90.0 * RAD_PER_DEG))))
        bwuse_v = bwdeg1 + abs(sin(abs(thetaA - (90.0 * RAD_PER_DEG))))
    case ("Y-Prime", "Y-prime", "y-Prime", "y-prime", "Y Prime", "Y prime", "y Prime", "y prime")
        Xa = (/ sin(thetaA+R)*cos(tauA) , sin(tauA) , cos(thetaA+R)*cos(tauA) /)
        bwuse_h = bwdeg1 + abs(sin(abs(tauA - RFA)))
        bwuse_v = bwdeg1 + abs(sin(abs(thetaA)))
    case ("Z","z")
        Xa = (/ sin(thetaA+R)*cos(tauA) , cos(thetaA+R)*cos(tauA) , sin(tauA) /)
        bwuse_v = bwdeg1 + abs(sin(abs(tauA)))
        bwuse_h = bwdeg1 + abs(sin(abs(thetaA - RFA)))
    end select

    Xearth = matmul(M,Xa)
    x = Xearth(1)
    y = Xearth(2)
    z = Xearth(3)

    ! added by B.Klotz (3/19/2021)
    select case (coordinate_system_type)
    case default
        write(*,'("coordinate_system_type = ''", A,"''")') trim(coordinate_system_type)
        stop "UNRECOGNIZED COORDINATE_SYSTEM_TYPE"
!    case ("X","x")
!       azimuth = asin(y) * DEG_PER_RAD
!       elevation = atan2(z,x)  * DEG_PER_RAD
    case ("Z","z","X","x","Y-Prime","Y-prime")
       azimuth = atan2(x,y) * DEG_PER_RAD
       elevation = asin(z)  * DEG_PER_RAD
    case ("Y","y")
       elevation = atan2(z,y) * DEG_PER_RAD
       azimuth = asin(x)  * DEG_PER_RAD
      ! z=-z
!    case ("Y-Prime", "Y-prime", "y-Prime", "y-prime", "Y Prime", "Y prime", "y Prime", "y prime")
!       elevation = asin(z) * DEG_PER_RAD
!       azimuth = atan2(x,y)  * DEG_PER_RAD
    end select
    ! end addition
  end subroutine coordinate_transformation

end module module_geometry
