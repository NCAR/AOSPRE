module module_flightpath

   ! Modifications by B. Klotz (6/1/2021)
   ! When the aircraft track is < 0, the value is converted to 
   ! a positive number (i.e., -90 becomes 270). However, computing
   ! sine and cosine with track > 180 reduces the precision and
   ! causes problems with generating a flight track. In order to
   ! solve this issue, a set of variables with tmp_track are added
   ! for the purpose of computing sine and cosine reliant portions
   ! of the code. The positive track values are the only ones
   ! sent to flightpath.ascii. A search through the document will 
   ! help find all instances of this addition.
      
   use module_configuration, only : RKIND
  private
  integer, parameter      :: MAXPATH = 36000
  real(kind=RKIND), public, parameter :: path_dt = 0.1

  type, public :: flightpath_class
      integer :: pathcount
      real(kind=RKIND), dimension(0:MAXPATH)   :: time_history
      real(kind=RKIND), dimension(0:MAXPATH)   :: track_history
      real(kind=RKIND), dimension(0:MAXPATH)   :: tmp_track_history        ! Keep history of tmp_track (BWK, 6/1/2021)
      real(kind=RKIND), dimension(3,0:MAXPATH) :: location_history         ! Vector of x,y,z for each pathcount
      real(kind=RKIND), dimension(3,0:MAXPATH) :: ground_velocity_history  ! Vector of dx/dt,dy/dt,dz/dt for each pathcount
      real(kind=RKIND), dimension(3,0:MAXPATH) :: wind_history             ! Vector of u,v,w at aircraft location for each pathcount

      real(kind=RKIND)                         :: time
      real(kind=RKIND), dimension(3)           :: location                 ! Instantaneous x, y, z.
      real(kind=RKIND), dimension(3)           :: wind                     ! Instantaneous local u, v, w.
      real(kind=RKIND)                         :: track
      real(kind=RKIND)                         :: tmp_track                ! Needed for sine/cosine (added by BWK, 6/1/2021

      integer                      :: waypoint_index
    contains
      procedure :: update => flightpath_update
      procedure :: update_location => flightpath_update_location
      procedure :: getloc => flightpath_getloc
      procedure :: dump => flightpath_dump
      procedure :: compute_at_time_and_location => flightpath_compute_at_time_and_location
      procedure :: precompute => flightpath_precompute
  end type flightpath_class

  interface flightpath_class
      module procedure flightpath_initialize
  end interface flightpath_class

  type, public :: aircraft_metadata_type
      real(kind=RKIND) :: x  = -1.E36 ! horizontal distance (meters) in i direction from grid i=1.
      real(kind=RKIND) :: y  = -1.E36 ! horizontal distance (meters) in j direction from grid j=1.
      real(kind=RKIND) :: z  = -1.E36 ! vertical distance (meters) from surface.

      real(kind=RKIND) :: dxdt = -1.E36 ! (m/s)
      real(kind=RKIND) :: dydt = -1.E36 ! (m/s
      real(kind=RKIND) :: dzdt = -1.E36 ! (m/s)

      real(kind=RKIND) :: xgrid = -1.E36 ! (fractional grid index)
      real(kind=RKIND) :: ygrid = -1.E36 ! (fractional grid index)
      real(kind=RKIND) :: zgrid = -1.E36 ! (fractional grid index)

      real(kind=RKIND) :: u = -1.E36 ! x-component of 3d wind at aircraft location
      real(kind=RKIND) :: v = -1.E36 ! y-component of 3d wind at aircraft location
      real(kind=RKIND) :: w = -1.E36 ! vertical component of 3d wind at aircraft location

      real(kind=RKIND) :: lat = -1.E36
      real(kind=RKIND) :: lon = -1.E36

      real(kind=RKIND) :: air_speed    = -1.E36
      real(kind=RKIND) :: ground_speed = -1.E36

      ! <heading> is the direction (degrees N), in the horizontal, that the aircraft is pointed, 
      !           which is not necessarily the direction the aircraft is moving.
      real(kind=RKIND) :: heading = -1.E36

      ! <track> is the direction (degrees N), in the horizontal, that the aircraft is moving.
      !         <track> = <heading> + <drift>
      real(kind=RKIND) :: track = -1.E36
      real(kind=RKIND) :: tmp_track = -1.E36

      ! <drift> is the difference between heading and track over the ground. 
      !         Positive drift implies track is clockwise from heading looking from above.
      real(kind=RKIND) :: drift = -1.E36

      ! <roll> is the angle about the longitudinal access of the aircraft (i.e., side-to-side roll).
      !        Positive is left side up, looking forward.
      real(kind=RKIND) :: roll = -1.E36

      real(kind=RKIND) :: pitch = -1.E36

  end type aircraft_metadata_type


contains

  !
  !-------------------------------------------------------------------------------
  !

  function flightpath_initialize(waypoint, t, gridspace_dx)
    use module_access_wrf, only : waypoint_type
    use module_llxy, only : DEG_PER_RAD
    implicit none
    type(flightpath_class) :: flightpath_initialize
    type(waypoint_type), dimension(0:), intent(in) :: waypoint
    real(kind=RKIND),                               intent(in) :: t
    real(kind=RKIND),                               intent(in) :: gridspace_dx

    real(kind=RKIND) :: xinc
    real(kind=RKIND) :: yinc
    real(kind=RKIND) :: track
    real(kind=RKIND) :: x_meters
    real(kind=RKIND) :: y_meters
    real(kind=RKIND) :: tmp_track

    x_meters = (waypoint(0)%x_grid - 1) * gridspace_dx
    y_meters = (waypoint(0)%y_grid - 1) * gridspace_dx

    if (waypoint(0)%z_meters > -1.E25) then
        write(*, '("FLIGHTPATH initialize at x_grid,y_grid,z_meters,t = ",4F20.6)') waypoint(0)%x_grid, waypoint(0)%y_grid, waypoint(0)%z_meters, t
    else
        write(*, '("FLIGHTPATH initialize at x_grid,y_grid,p_Pa,t = ",4F20.6)') waypoint(0)%x_grid, waypoint(0)%y_grid, waypoint(0)%p_Pa, t
    endif

    flightpath_initialize%pathcount = 0
    flightpath_initialize%time_history = -1.E36
    flightpath_initialize%location_history = -1.E36
    flightpath_initialize%ground_velocity_history = -1.E36
    flightpath_initialize%wind_history = -1.E36
    flightpath_initialize%track_history = -1.E36
    flightpath_initialize%tmp_track_history = -1.E36

    flightpath_initialize%time_history(0) = t
    flightpath_initialize%location_history(:,0) = (/x_meters, y_meters, waypoint(0)%z_meters/)
    flightpath_initialize%time = t
    flightpath_initialize%location = (/x_meters, y_meters, waypoint(0)%z_meters/)
    flightpath_initialize%wind = -1.E36

    !
    ! Compute the track that will take the aircraft from waypoint(0) to waypoint(1).
    ! Track is the direction the aircraft is moving, projected onto the horizontal plane. 
    ! Track is measured clockwise from North, range [0,360).
    !

    xinc = waypoint(1)%x_grid - waypoint(0)%x_grid
    yinc = waypoint(1)%y_grid - waypoint(0)%y_grid
    track = atan2(xinc,yinc)*DEG_PER_RAD
    tmp_track = track
    print*, 'xinc, yinc, track: ', xinc, yinc, track
    ! print*, 'size(waypoint) = ', size(waypoint)
    if (track<0) track = track + 360.0
    print*, 'Check the track change: ',track

    ! round track to nearest 1e-06 degree
    ! tmp_f = 10.0**(5.0)
    ! tmp_track = NINT(tmp_f * track) / tmp_f
    ! print*, 'Temporary track: ', tmp_track
    ! track = tmp_track

    flightpath_initialize%track = track
    flightpath_initialize%tmp_track = tmp_track
    flightpath_initialize%track_history(0) = track
    flightpath_initialize%tmp_track_history(0) = tmp_track

    flightpath_initialize%waypoint_index = 1 ! I.e., we're in the leg from waypoint(0) to waypoint(1)


  end function flightpath_initialize

  !
  !-------------------------------------------------------------------------------
  !

  subroutine flightpath_compute_at_time_and_location(self, aircraft, helicopter)

    !
    ! From:
    !     self%pathcount
    !     self%wind
    !     aircraft%air_speed (or aircraft%dxdt,dydt,dzdt)
    !     aircraft%track
    ! Update:
    !     self%wind_history
    !     self%ground_velocity_history
    !

    use module_external_attitude, only : use_external_attitudes
    use module_llxy, only : RAD_PER_DEG
    implicit none
    class (flightpath_class), intent(inout) :: self
    type (aircraft_metadata_type), intent(in) :: aircraft
    logical, intent(in) :: helicopter

    real(kind=RKIND), dimension(3) :: ground_velocity ! in x, y, z.
    real(kind=RKIND), dimension(3) :: flightpath_air_velocity ! in x, y, z.

    self%wind_history(:,self%pathcount) = self%wind

    if (helicopter) then
        self%ground_velocity_history(:,self%pathcount) = (/0.0,0.0,0.0/)
        return
    endif

    if ( use_external_attitudes ) then
        if ( aircraft%dxdt < -1.E25 ) stop "dxdt not initialized"
        ground_velocity = (/ aircraft%dxdt, aircraft%dydt, aircraft%dzdt /)
    else

        ! aircraft_air_velocity = aircraft%air_speed * (/ sin(aircraft%track * RAD_PER_DEG ) , &  ! x-component
        !      &                                         cos(aircraft%track * RAD_PER_DEG ) , &  ! y-component
        !      &                                          0.0 /) ! Z-component

        flightpath_air_velocity = aircraft%air_speed * (/ sin(self%tmp_track * RAD_PER_DEG ) , &  ! x-component
             &                                         cos(self%tmp_track * RAD_PER_DEG ) , &  ! y-component
             &                                          0.0_RKIND /) ! Z-component
        
        !print*, 'Sind = ',sin(self%tmp_track * RAD_PER_DEG)
        !print*, 'Cosd = ',cos(self%tmp_track * RAD_PER_DEG)
        !print*, 'Flightpath air vel., wind: ', flightpath_air_velocity,self%wind
        if (self%wind(1) < -9998) then
            ground_velocity = flightpath_air_velocity
        else
            ground_velocity = flightpath_air_velocity * ( 1.0 +  dot_product(flightpath_air_velocity,self%wind) / ( aircraft%air_speed * aircraft%air_speed ) )
        endif
        !print*, 'Ground vel. = ', ground_velocity
    endif
    self%ground_velocity_history(:,self%pathcount) = ground_velocity
    ! print*, 'pathcount = ', self%pathcount, ground_velocity


  end subroutine flightpath_compute_at_time_and_location

  !
  !-------------------------------------------------------------------------------
  !

  subroutine flightpath_update_location(self,aircraft,waypoint,gridspace_dx)

    !
    ! From:
    !     self%pathcount
    !     self%time
    !     self%ground_velocity_history
    !     self%location
    !
    ! Update:
    !     self%pathcount
    !     self%time
    !     self%location
    !     self%location_history
    !     self%time_history
    !

    use module_llxy, only : RAD_PER_DEG
    use module_llxy, only : DEG_PER_RAD
    use module_access_wrf, only : waypoint_type
    use module_external_attitude, only : use_external_attitudes
    implicit none
    class (flightpath_class), intent(inout) :: self
    type (aircraft_metadata_type), intent(in) :: aircraft
    type (waypoint_type), dimension(0:), intent(in) :: waypoint
    real(kind=RKIND), intent(in) :: gridspace_dx
    real(kind=RKIND), dimension(3) :: gvel
    real(kind=RKIND), dimension(3) :: wind
    real(kind=RKIND), dimension(3) :: lm

    real(kind=RKIND) :: x0, xa, xm, x1, x2, xn, y0, ya, ym, y1, y2, yn, p1, p2, dtA, dtB, segment_length, new_track, new_tmp_track
    real(kind=RKIND), dimension(3) :: aircraft_air_velocity, new_gvel ! in x, y, z.
    real(kind=RKIND), dimension(3) :: recomputed_location
    real(kind=RKIND) :: ground_speed

    self%pathcount = self%pathcount + 1

    ! Quick sanity check.
    if (self%location(1) /= self%location_history(1,self%pathcount-1)) stop "LOCATION"
    if (self%location(2) /= self%location_history(2,self%pathcount-1)) stop "LOCATION"
    if (self%location(3) /= self%location_history(3,self%pathcount-1)) stop "LOCATION"
    if (self%time /= self%time_history(self%pathcount-1)) stop "TIME"

    gvel = self%ground_velocity_history(:,self%pathcount-1)

    ! write(*,'("Targeted time = ", F12.4)') self%time + path_dt
    ! write(*,'("Recompute from time/location = ", F12.4, 3F12.6, "  gvel=", 3F12.6)') self%time_history(self%pathcount-1), self%location, gvel
    recomputed_location = self%location + path_dt * gvel

    !  But if the new location is beyond the line from waypoint(n) to waypoint(n+1), then
    !  we need to reach a corner and update the track, heading, etc.

    !
    !  Check if we're between waypoint(waypoint_index-1) and waypoint(waypoint_index)
    !

    ! (xa,ya) is our newly-computed aircraft location.

    xa = recomputed_location(1)
    ya = recomputed_location(2)

    ! print*, 'lbound/ubound (waypoint,1) = ', lbound(waypoint,1), ubound(waypoint,1)
    ! print*, 'size(waypoint) = ', size(waypoint)
    ! print*, 'waypoint_index = ', self%waypoint_index
    if (self%waypoint_index > ubound(waypoint,1)) stop "(1) OUT OF WAYPOINTS"

    x0 = (waypoint(self%waypoint_index-1)%x_grid - 1) * gridspace_dx
    x1 = (waypoint(self%waypoint_index  )%x_grid - 1) * gridspace_dx

    y0 = (waypoint(self%waypoint_index-1)%y_grid - 1) * gridspace_dx
    y1 = (waypoint(self%waypoint_index  )%y_grid - 1) * gridspace_dx

    ! print*, x0,x1,xa
    ! print*, y0,y1,ya
    ! print*, pl_dist(x0,y0,x1,y1,xa,ya)
#ifdef _NOPE_
    segment_length = hypot(x1-x0, y1-y0)
    p1 = hypot(xa-x0,ya-y0)
    p2 = hypot(x1-xa,y1-ya)

    ! print*, 'Seg length, p1, p2: ', segment_length, p1, p2
    ! print*, 'pldist: ', pl_dist(x0,y0,x1,y1,xa,ya)
#endif

    ! print*, 'xy(0) = ', x0, y0
    ! print*, 'xy(1) = ', x1, y1
    ! print*, 'xy(a) = ', xa, ya
    ! write(*,'("Segments:", 4F20.12, " |", 3F20.12, " |", 3F20.12)') p1+p2, segment_length, p1, p2, gvel
    ! print*, (y1-y0)/(x1-x0), (ya-y0)/(xa-x0)

    if (pl_dist(x0,y0,x1,y1,xa,ya) > 5.E-8) then ! edited by B.Klotz (3/30/2021)
        print*, 'out of leg!', pl_dist(x0,y0,x1,y1,xa,ya)
        print*, waypoint(self%waypoint_index), self%waypoint_index, ubound(waypoint,1)

        if (self%waypoint_index+1 > ubound(waypoint,1)) stop "(2) OUT OF WAYPOINTS"

        !
        !  So let's assume we hit waypoint(self%waypoint_index) at (x1,y1), turn a corner,
        !  and head directly toward waypoint(self%waypoint_index+1) at (x2,y2).
        !

        x2 = (waypoint(self%waypoint_index+1)%x_grid - 1) * gridspace_dx
        y2 = (waypoint(self%waypoint_index+1)%y_grid - 1) * gridspace_dx

        ! (xm,ym) is our previous aircraft location.
        lm = self%location_history(:,self%pathcount-1)
        xm = lm(1)
        ym = lm(2)

        !
        !  It takes some fraction of the time step path_dt to go from (xm,ym) to (x1,y1).
        !  The complement of that fraction will take us some distance towards (x2,y2) (at a different ground-relative velocity).
        !

        write(*,'("Length of segment from (xm,ym) to (x1,y1) = ", F20.12)') hypot(x1-xm,y1-ym)
        dtA = hypot(x1-xm,y1-ym)/hypot(gvel(1),gvel(2))  ! meters divided by m/s gives us seconds.
        dtB = path_dt-dtA
        print*, hypot(gvel(1),gvel(2)), 'dtA = ', dtA, dtB, (/xm,ym/) + (/gvel(1),gvel(2)/)*dtA

        !
        !  Compute our new track direction, from (x1,y1) to (x2,y2)
        !

        print*,'x2-x1 = ', x2-x1, x2, x1
        print*,'y2-y1 = ', y2-y1, y2, y1
        new_track = atan2(x2-x1,y2-y1)*DEG_PER_RAD
        new_tmp_track = new_track
        if (new_track<0) new_track = new_track + 360.0
        print*, 'new track = ', new_track

        !
        !  Knowing the new track, compute a new headwind from 3d wind at the prior aircraft location.
        !
        !  We need to know the wind in the vicinity to compute a new ground-relative velocity.
        !  Probably the wind at the most recent location is close enough.
        !

        wind = self%wind_history(:,self%pathcount-1)

        aircraft_air_velocity = aircraft%air_speed * (/ sin(new_tmp_track * RAD_PER_DEG ) , &  ! x-component
             &                                          cos(new_tmp_track * RAD_PER_DEG ) , &  ! y-component
             &                                          0.0_RKIND /) ! Z-component

        if (self%wind(1) < -9998) then
            new_gvel = aircraft_air_velocity
        else
            new_gvel = aircraft_air_velocity * ( 1.0 +  dot_product(aircraft_air_velocity,wind) / ( aircraft%air_speed * aircraft%air_speed ) )
        endif

        write(*,'("aircraft_air_velocity = ", 3F12.6)') aircraft_air_velocity
        write(*,'("aircraft%air_speed    = ", F12.6)') aircraft%air_speed
        write(*,'("new_gvel              = ", 3F12.6)') new_gvel
        write(*,'("wind                  = ", 3F12.6)') wind

        ground_speed = sqrt(dot_product(new_gvel,new_gvel))

        !
        !  Move at this new ground velocity from (x1,y1) towards point (x2,y2)
        !

        xn = x1 + new_gvel(1) * dtB
        yn = y1 + new_gvel(2) * dtB

        !
        !  Insert the exact waypoint location into our flightpath history
        !  At the waypoint, the track will be the direction to the next waypoint.
        !  At the waypoint, the ground-relative velocity will reflect the new track to the next waypoint.
        !

        self%time_history(self%pathcount) = self%time + dtA
        self%track_history(self%pathcount) = new_track
        self%tmp_track_history(self%pathcount) = new_tmp_track
        self%location_history(:,self%pathcount) = (/x1,y1,self%location(3)/)
        self%ground_velocity_history(:,self%pathcount) = new_gvel
        self%wind_history(:,self%pathcount) = wind ! This is actually the wind at the previous point in the path.

        ! print*, 'insert ', self%time_history(self%pathcount)
        ! print*, 'insert ', self%location_history(:,self%pathcount)

        self%pathcount = self%pathcount + 1

        !
        !  And update everything relevant.
        !

        self%waypoint_index = self%waypoint_index + 1
        recomputed_location = (/xn,yn,self%location(3)/)
        self%track = new_track
        self%tmp_track = new_tmp_track
        self%ground_velocity_history(:,self%pathcount) = new_gvel

    endif

    self%time = self%time + path_dt
    self%time_history(self%pathcount) = self%time

    self%track_history(self%pathcount) = self%track
    self%tmp_track_history(self%pathcount) = self%tmp_track

    self%location = recomputed_location
    self%location_history(:,self%pathcount) = self%location

    ! print*, 'insert ', self%time_history(self%pathcount)
    ! print*, 'insert ', self%location_history(:,self%pathcount)


  end subroutine flightpath_update_location

  !
  !-------------------------------------------------------------------------------
  !

  subroutine flightpath_update(self, aircraft)
    ! From:
    !     self%wind
    !     self%pathcount
    !     aircraft%track
    !     aircraft%air_speed
    !
    ! Update:
    !     self%wind_history
    !     self%ground_velocity_history
    !     self%location
    !     self%time
    !     self%pathcount
    !     self%location_history
    !     self%time_history
    use module_llxy, only : RAD_PER_DEG
    use module_external_attitude, only : use_external_attitudes
    implicit none
    class (flightpath_class), intent(inout) :: self
    type (aircraft_metadata_type), intent(in) :: aircraft

    real(kind=RKIND), dimension(3) :: ground_velocity ! in x, y, z.
    real(kind=RKIND), dimension(3) :: aircraft_air_velocity ! in x, y, z.

    ! A simple forward step.
    ! The location at time n is computed from the location and velocity at time n-1.
    ! self%time = self%time + path_dt

    self%wind_history(:,self%pathcount) = self%wind

    if ( use_external_attitudes ) then
        if ( aircraft%dxdt < -1.E25 ) stop "dxdt not initialized"
        ground_velocity = (/ aircraft%dxdt, aircraft%dydt, aircraft%dzdt /)
    else

        ! aircraft_air_velocity = aircraft%air_speed * (/ sin(aircraft%track * RAD_PER_DEG ) , &  ! x-component
        !      &                                         cos(aircraft%track * RAD_PER_DEG ) , &  ! y-component
        !      &                                          0.0 /) ! Z-component

        aircraft_air_velocity = aircraft%air_speed * (/ sin(self%tmp_track * RAD_PER_DEG ) , &  ! x-component
             &                                         cos(self%tmp_track * RAD_PER_DEG ) , &  ! y-component
             &                                          0.0_RKIND /) ! Z-component

        ! The ground velocity is also at time n
        if (self%wind(1) < -9998) then
            ground_velocity = aircraft_air_velocity
        else
            ground_velocity = aircraft_air_velocity * ( 1.0 +  dot_product(aircraft_air_velocity,self%wind) / ( aircraft%air_speed * aircraft%air_speed ) )
        endif
    endif
    self%ground_velocity_history(:,self%pathcount) = ground_velocity
    print*, 'pathcount = ', self%pathcount, ground_velocity

    self%location = self%location + path_dt * ground_velocity

    self%pathcount = self%pathcount + 1
    if ( self%pathcount > MAXPATH ) then
        stop "module_flightpath MAXPATH exceeded"
    endif

    ! The updated location is at time n+1
    self%location_history(:,self%pathcount) = self%location
    ! self%time_history(self%pathcount) = self%time
    ! write(*, '("FLIGHTPATH update at x,y,z,t = ",4F20.6)') self%location(1), self%location(2), self%location(3), self%time

    ! But I don't have the wind at time n+1 yet.

  end subroutine flightpath_update

  !
  !-------------------------------------------------------------------------------
  !

  subroutine flightpath_getloc(self, time_request, proj, aircraft, error_flag)
    ! Input:
    !
    !     from the self of class(flightpath_class):
    !        self%pathcount
    !        self%time_history
    !        self%location_history
    !        self%ground_velocity_history
    !
    !     time_request
    !     proj
    !     aircraft%heading
    !     aircraft%drift
    !
    ! Updates:
    !         aircraft%x     , aircraft%y     , aircraft%z    , 
    !         aircraft%u     , aircraft%v     , aircraft%w    , 
    !         aircraft%dxdt  , aircraft%dydt  , aircraft%dzdt , 
    !         aircraft%xgrid , aircraft%ygrid
    !         aircraft%ground_speed
    !         aircraft%roll  , aircraft%pitch , aircraft%drift ,
    !         aircraft%track
    !         aircraft%lat   , aircraft%lon
    !
    use module_llxy, only : proj_info
    use module_llxy, only : ij_to_latlon
    use module_external_attitude, only : use_external_attitudes
    use module_external_attitude, only : attitude_details_at_time
    implicit none
    class(flightpath_class),       intent(in)    :: self
    real(kind=RKIND),                          intent(in)    :: time_request
    type(proj_info),               intent(in)    :: proj
    type (aircraft_metadata_type), intent(inout) :: aircraft
    integer                      , intent(out)   :: error_flag

    real(kind=RKIND), dimension(3) :: loc
    real(kind=RKIND), dimension(3) :: gvel
    real(kind=RKIND), dimension(3) :: wind
    real(kind=RKIND) :: track
    real(kind=RKIND) :: tmp_track
    integer :: k
    real(kind=RKIND) :: fraction
    error_flag = 0

    k = self%pathcount
    if ( abs(time_request-self%time_history(k)) < 1.E-6 ) then
        loc = self%location_history(:,k)
        gvel = self%ground_velocity_history(:,k)
        wind = self%wind_history(:,k)
        track = self%track_history(k)
        tmp_track = self%tmp_track_history(k)
        
    else
        KLOOP : do k = self%pathcount, 0, -1
            if ( abs(time_request-self%time_history(k)) < 1.E-6 ) then
                loc = self%location_history(:,k)
                gvel = self%ground_velocity_history(:,k)
                wind = self%wind_history(:,k)
                track = self%track_history(k)
                tmp_track = self%tmp_track_history(k)
                exit KLOOP
            endif
            if ( ( self%time_history(k) < time_request ) .and. ( time_request < self%time_history(k+1) ) ) then
                fraction = (time_request - self%time_history(k)) / ( self%time_history(k+1) - self%time_history(k))
                loc = self%location_history(:,k)*(1.0-fraction) + self%location_history(:,k+1)*(fraction)
                gvel = self%ground_velocity_history(:,k)*(1.0-fraction) + self%ground_velocity_history(:,k+1)*(fraction)

                ! Not interpolating track, take the closest one (in time)
                if (fraction>0.5) then
                    track = self%track_history(k+1)
                    
                    tmp_track = self%tmp_track_history(k+1)
                else
                    track = self%track_history(k)
                    tmp_track = self%tmp_track_history(k)
                endif
                    
                if (abs(gvel(1)) > 50000) then
                    print*, 'k = ', k
                    print*, 'pathcount = ', self%pathcount
                    print*, 'times: ', self%time_history(k), time_request, self%time_history(k+1)
                    print*, 'fraction = ', fraction
                    print*, "ground_velocity_history: ", self%ground_velocity_history(:,k)
                    print*, "ground_velocity_history: ", self%ground_velocity_history(:,k+1)
                    print*, gvel
                    stop "GVEL"
                endif
                wind = self%wind_history(:,k)*(1.0-fraction) + self%wind_history(:,k+1)*(fraction)
                exit KLOOP
            endif
            if (k == 0) then
                write (*,'("Problem module_flightpath:flightpath_getloc")')
                write(*,'("time_request = ", F18.8)') time_request
                stop "Problem module_flightpath:flightpath_getloc"
            endif

        enddo KLOOP
    endif

    aircraft%x = loc(1)
    aircraft%y = loc(2)
    aircraft%z = loc(3)

    aircraft%u = wind(1)
    aircraft%v = wind(2)
    aircraft%w = wind(3)

    aircraft%xgrid = 1.0 + ( loc(1) / proj%dx )
    aircraft%ygrid = 1.0 + ( loc(2) / proj%dx )

    if ( use_external_attitudes ) then
        call attitude_details_at_time ( time_request, aircraft%roll, aircraft%pitch, aircraft%drift, aircraft%heading, &
             &                                        aircraft%dxdt, aircraft%dydt, aircraft%dzdt, error_flag )
        if ( error_flag > 0 ) return
        gvel = (/ aircraft%dxdt, aircraft%dydt, aircraft%dzdt /)
        aircraft%track = aircraft%heading + aircraft%drift
    else
        aircraft%roll = 0.0
        aircraft%pitch = 0.0
        aircraft%drift = 0.0
        ! aircraft%heading unchanged

        aircraft%dxdt = gvel(1)
        aircraft%dydt = gvel(2)
        aircraft%dzdt = gvel(3)
        ! Assuming the heading and track are the same:
        aircraft%track = track !aircraft%heading + aircraft%drift
        aircraft%heading = aircraft%track
    endif
    aircraft%ground_speed = sqrt(dot_product(gvel,gvel))

    call ij_to_latlon ( proj, aircraft%xgrid, aircraft%ygrid, aircraft%lat, aircraft%lon )

    if ( abs(aircraft%dxdt) > 5000 ) then
        print*, 'Problem:  aircraft%dxdt = ', aircraft%dxdt
        stop "dxdt"
    endif

!KWM    write(*,'(21F12.3)') aircraft%x, aircraft%y, aircraft%z, &
!KWM         aircraft%dxdt, aircraft%dydt, aircraft%dzdt, &
!KWM         aircraft%xgrid, aircraft%ygrid, & ! aircraft%zgrid, &
!KWM         aircraft%u, aircraft%v, aircraft%w, &
!KWM         aircraft%lat, aircraft%lon, & ! aircraft%air_speed, &
!KWM         aircraft%ground_speed, aircraft%heading, aircraft%track! , &
!KWM         ! aircraft%drift, aircraft%roll, aircraft%pitch


  end subroutine flightpath_getloc

  !
  !-------------------------------------------------------------------------------
  !

  subroutine flightpath_dump(self, outfile, rank)
    implicit none
    class (flightpath_class), intent(in) :: self
    character(len=*), intent(in) :: outfile
    integer, intent(in) :: rank
    integer :: ierr
    integer :: k
    integer :: dumpunit

    open(newunit=dumpunit, file=trim(outfile), form='formatted', action='write', iostat=ierr)
    if ( ierr /= 0 ) then
        print*, 'Problem opening '//trim(outfile)
        return
    endif

    if (rank == 0) then
        write(dumpunit,'(A10,A16,3A20,3A16,3A16,A16)') "Index", "   Time   ", &
             &                                    "       X        ", "       Y        ", "       Z        ", &
             &                                    "   ground dX/dT ", "   ground dY/dT ", "   ground dZ/dT ", &
             &                                    "      U wind    ", "      V wind    ", "      W wind    ", &
             &                                    "     track      "
    endif

    do k = 0, self%pathcount
        write(dumpunit,'(I10, F16.4, 3F20.8, 3F16.8, 3F16.8, F16.8)') k, self%time_history(k), &
             self%location_history(:,k), self%ground_velocity_history(:,k), self%wind_history(:,k), &
             self%track_history(k)
    enddo

    close(dumpunit)

  end subroutine flightpath_dump

  !
  !-------------------------------------------------------------------------------
  !

  real(kind=RKIND) function pl_dist(x1, y1, x2, y2, x3, y3)
    !  Distance of point (x3,y3) from line segment defined
    !  by endpoints (x1,y1) and (x2,y2)
    implicit none
    real(kind=RKIND), intent(in) :: x1, y1, x2, y2, x3, y3
    real(kind=RKIND) :: px, py, norm, u, x, y
    px = x2-x1
    py = y2-y1
    norm = px**2 + py**2
    u =  ((x3 - x1) * px + (y3 - y1) * py) / norm
    if (u > 1) then
        u = 1.0
    else if (u < 0) then
        u = 0.0
    endif
    x = x1 + (u * px)
    y = y1 + (u * py)
   ! print*, 'px,py,norm: ',px, py, norm
   ! print*, 'u,x,y: ', u, x, y

    pl_dist = hypot(x-x3, y-y3)
  end function pl_dist

  !
  !-------------------------------------------------------------------------------
  !

  subroutine flightpath_precompute(self, options, gridspace_dx, wrf_xtimes, wrf_times, wrf_files, wrf_file_index, wrf_time_index)
    use module_access_wrf,       only : options_type
    use module_access_wrf,       only : wrf_metadata_type
    use module_access_wrf,       only : MAX_WRF_FILES
    use module_access_wrf,       only : MAX_WRF_TIMES
    use module_timing_utilities, only : timer_start
    use module_timing_utilities, only : timer_pause
    use module_timing_utilities, only : timer_resume
    implicit none
    class(flightpath_class), intent(inout) :: self
    type(options_type),      intent(in)  :: options
    real(kind=RKIND),        intent(in)  :: gridspace_dx
    character(len=1024), dimension(MAX_WRF_FILES), intent(in) :: wrf_files
    real(kind=RKIND),    dimension(MAX_WRF_FILES*MAX_WRF_TIMES), intent(in) :: wrf_xtimes
    character(len=19),   dimension(MAX_WRF_FILES*MAX_WRF_TIMES), intent(in) :: wrf_times
    integer,             dimension(MAX_WRF_FILES*MAX_WRF_TIMES), intent(in) :: wrf_file_index
    integer,             dimension(MAX_WRF_FILES*MAX_WRF_TIMES), intent(in) :: wrf_time_index

    type (aircraft_metadata_type) :: aircraft
    integer                       :: milliseconds
    real(KIND=RKIND)              :: leg_current_time
    integer                       :: WRFindexA
    integer                       :: WRFindexB
    integer                       :: i
    integer                       :: j
    real(kind=RKIND)              :: time_interpolation_factor
    type (wrf_metadata_type)      :: metaA
    type (wrf_metadata_type)      :: metaB
    real(KIND=RKIND)              :: xgrid
    real(KIND=RKIND)              :: ygrid
    real(KIND=RKIND)              :: zgrid
    real(kind=RKIND)              :: Ua
    real(kind=RKIND)              :: Va
    real(kind=RKIND)              :: Wa
    real(kind=RKIND)              :: Ub
    real(kind=RKIND)              :: Vb
    real(kind=RKIND)              :: Wb

    type(wrf_metadata_type), pointer :: Aptr
    type(wrf_metadata_type), pointer :: Bptr
    real(kind=RKIND),  parameter :: maxgspeed = 200.0
    real(kind=RKIND) :: xlocation, ylocation
    real(kind=RKIND) :: wrf_output_interval_seconds
    real(kind=RKIND) :: maxwindow
    real(kind=RKIND) :: xminimum
    real(kind=RKIND) :: xmaximum
    real(kind=RKIND) :: yminimum
    real(kind=RKIND) :: ymaximum

    aircraft%track = self%track
    aircraft%air_speed = options%air_speed
    aircraft%drift = 0.0
    aircraft%heading = aircraft%track - aircraft%drift

#ifdef _LATER_
    if ( use_external_attitudes ) then
        call attitude_details_at_time ( options%leg_initial_time, aircraft%roll, aircraft%pitch, aircraft%drift, aircraft%heading, &
             &                                                    aircraft%dxdt, aircraft%dydt, aircraft%dzdt, error_flag )
        if ( error_flag  > 0 ) stop "Attitude details?"
        
        call attitude_details_at_time ( options%leg_initial_time, aircraft%roll, aircraft%pitch, aircraft%drift, aircraft%heading, &
             &                                                    aircraft%dxdt, aircraft%dydt, aircraft%dzdt, error_flag )
        if ( error_flag  > 0 ) stop "Attitude details?"
    endif
#endif

    xgrid = 0.0
    ygrid = 0.0
    zgrid = 0.0

    PLEG: do milliseconds = 0, int(options%leg_time_seconds*1000), 100

        leg_current_time = options%leg_initial_time + (milliseconds*1.E-3)

#ifdef _NOPE_
        if (aircraft%x < -1.E25) then
            ! x and y in model grid space [origin at (1.0,1.0)]
            write(*, '("leg initial time, current_time, x, y = " 2F18.8, 2F18.4 )') &
                 &      options%leg_initial_time, leg_current_time, &
                 &      1.0+(self%location_history(1,0)/gridspace_dx), 1.0+(self%location_history(2,0)/gridspace_dx)
        else
            ! x and y in model grid space [origin at (1.0,1.0)]
            write(*, '("leg initial time, current_time, x, y = " 2F18.8, 2F18.4 )') &
                 &      options%leg_initial_time, leg_current_time, &
                 &      1.0+(aircraft%x/gridspace_dx), 1.0+(aircraft%y/gridspace_dx)  
        endif
#endif
        ! write(*, '("leg initial time, current_time, x, y = " 2F18.8, 2F18.4 )') &
        !      &      options%leg_initial_time, leg_current_time, &
        !      &      1.0+(self%location_history(1,self%pathcount)/gridspace_dx), 1.0+(self%location_history(2,self%pathcount)/gridspace_dx)


        !  Nearest point of current aircraft location, used as the center point of a subset of the full WRF field to read from the NetCDF files
        xlocation = nint(self%location_history(1,self%pathcount)/gridspace_dx)
        ylocation = nint(self%location_history(2,self%pathcount)/gridspace_dx)

        ! Find WRF data at (or bracketing) the current time.

        ! call cpu_time(fstart_wrfsearch)

        if (options%time_evolution) then
            WRF_TIME_SEARCH : do i = 1, MAX_WRF_FILES*MAX_WRF_TIMES
                if (wrf_file_index(i) < 0) exit WRF_TIME_SEARCH
                if (abs(wrf_xtimes(i) - leg_current_time) < -1.E3) then
                    ! Time match
                    if (wrf_xtimes(i+1) > -1.E25) then
                        wrf_output_interval_seconds = wrf_xtimes(i+1) - wrf_xtimes(i)
                    endif
                    time_interpolation_factor = 1.0
                    WRFindexA = i
                    WRFindexB = -999999
                    exit WRF_TIME_SEARCH
                else if (wrf_xtimes(i) > leg_current_time) then
                    WRFindexB = i
                    wrf_output_interval_seconds = wrf_xtimes(WRFindexB) - wrf_xtimes(WRFindexA)
                    time_interpolation_factor = ( leg_current_time - wrf_xtimes(WRFindexA) ) / &
                         ( wrf_xtimes(WRFindexB) - wrf_xtimes(WRFindexA) )
                    time_interpolation_factor = 1.0-time_interpolation_factor
                    exit WRF_TIME_SEARCH
                endif
                WRFindexA = i
            enddo WRF_TIME_SEARCH
        else
            time_interpolation_factor = 1.0
            WRFindexA = 1
            WRFindexB = -999999
        endif

        ! call cpu_time(fend_wrfsearch)
        ! actime_wrfsearch = actime_wrfsearch + (fend_wrfsearch-fstart_wrfsearch)
        
        ! call cpu_time(fstart_wrfopen)
        if ((wrf_files(wrf_file_index(WRFindexA)) /= metaA%flnm ) .or. (wrf_time_index(WRFindexA) /= metaA%time_frame)) then

            !
            ! Close an old file if necessary.
            !

            call metaA%close(silent=.TRUE.)

            !
            ! Open a new file
            !

            ! In defining a subset of the full WRF grid to read, assume a maximum ground speed of <maxgspeed> m/s over 
            ! the time period of <wrf_output_interval_seconds>.  Define a square +/- distance of <maxwindow> about the point (xlocation,ylocation).
            maxwindow = 1+(maxgspeed * wrf_output_interval_seconds) / gridspace_dx
            xminimum = xlocation - maxwindow
            xmaximum = xlocation + maxwindow
            yminimum = ylocation - maxwindow
            ymaximum = ylocation + maxwindow

            ! write(*, '("A wrf_file = ", A)') trim(wrf_files(wrf_file_index_A))
            call timer_resume("OpenA")
            call metaA%open_precompute(options%flight_level_coordinate, trim(wrf_files(wrf_file_index(WRFindexA))), wrf_time_index(WRFindexA), &
                 &                     xminimum, xmaximum, yminimum, ymaximum)
            call timer_pause("OpenA")
        endif

        if (time_interpolation_factor < 1) then
            if ( (wrf_files(wrf_file_index(WRFindexB)) /= metaB%flnm) .or. (wrf_time_index(WRFindexB) /= metaB%time_frame)) then
                maxwindow = 1+(maxgspeed * wrf_output_interval_seconds) / gridspace_dx
                xminimum = xlocation - maxwindow
                xmaximum = xlocation + maxwindow
                yminimum = ylocation - maxwindow
                ymaximum = ylocation + maxwindow
                ! write(*, '("B metafile = ", A)') trim(metaB%flnm)
                call metaB%close(silent=.TRUE.)
                ! write(*, '("B wrf_file = ", A)') trim(wrf_files(wrf_file_index_B))
                call timer_resume("OpenB")
                call metaB%open_precompute(options%flight_level_coordinate, trim(wrf_files(wrf_file_index(WRFindexB))), wrf_time_index(WRFindexB), &
                     &                     xminimum, xmaximum, yminimum, ymaximum)
                call timer_pause("OpenB")
            endif
        endif
        ! write(*,'("before read: ", A)') trim(metaA%flnm)

        ! call cpu_time(fend_wrfopen)
        ! actime_wrfopen = actime_wrfopen + (fend_wrfopen-fstart_wrfopen)
        ! call cpu_time(fstart_presweep)

        !
        !  Find the wind velocity at the aircraft location along the flight path.
        !

        if (self%pathcount == 0) then

            !
            !  At the first time, if we're following a pressure surface, we can now 
            !  compute our starting Z, having the model data on hand.
            !


            !
            !  Don't update location yet, we have to compute wind and ground_velocity at the initial point.
            !

            xgrid = 1 + (self%location(1) / gridspace_dx)
            ygrid = 1 + (self%location(2) / gridspace_dx)
            if (options%flight_level_coordinate == "P") then
                self%location(3) = metaA%get_z_at_p(xgrid, ygrid, options%waypoint(0)%p_Pa, z_index=zgrid)
                zgrid = metaA%find_z_index(xgrid, ygrid, self%location(3))
                self%location_history(3,self%pathcount) = self%location(3) ! Correct the history, too.
            else
                zgrid = metaA%find_z_index(xgrid, ygrid, self%location(3))
            endif

            if ( time_interpolation_factor < 1.0 ) then
                ! print*, 'time_interpolation_factor = ', time_interpolation_factor

                call metaA%value ( xgrid, ygrid, zgrid, Ua, Va, Wa )
                call metaB%value ( xgrid, ygrid, zgrid, Ub, Vb, Wb )
                self%wind(1) = Ua*(time_interpolation_factor) + Ub*(1.0-time_interpolation_factor)
                self%wind(2) = Va*(time_interpolation_factor) + Vb*(1.0-time_interpolation_factor)
                self%wind(3) = Wa*(time_interpolation_factor) + Wb*(1.0-time_interpolation_factor)

            else
                call metaA%value ( xgrid, ygrid, zgrid, &
                     &             self%wind(1), self%wind(2), self%wind(3) )
            endif

            !
            ! Now knowing the wind at time zero, we can now update
            ! self%wind_history, self%ground_velocity_history
            ! for time==0.0, self%pathcount==0
            !

            call self%compute_at_time_and_location(aircraft, options%helicopter) 

        endif

        !
        ! The call to self%update_location will update:
        !      self%pathcount  (increment by one)
        !      self%time       (increment by self%path_dt)
        !      self%location   (considering last known location, ground velocity at that place and time, and the delta-time)
        !      self%location_history (add newly-computed location at updated index self%pathcount)
        !      self%time_history     (add newly-computed time at updated index self%pathcount)
        !

        call self%update_location(aircraft,options%waypoint,gridspace_dx)

        xgrid = 1 + (self%location(1) / gridspace_dx)
        ygrid = 1 + (self%location(2) / gridspace_dx)
        ! But if we're following a pressure surface, recompute the new Z location at the desired P level.
        if (options%flight_level_coordinate == "P") then
            self%location(3) = metaA%get_z_at_p(xgrid, ygrid, options%waypoint(0)%p_Pa, z_index=zgrid)
            self%location_history(3,self%pathcount) = self%location(3) ! Correct the history, too.
        else
            zgrid = metaA%find_z_index(xgrid, ygrid, self%location(3))
        endif

        !
        ! Knowing new aircraft location, now we update our knowledge of the wind at the new aircraft location.
        !

        if ( time_interpolation_factor < 1.0 ) then
            call metaA%value ( xgrid, ygrid, zgrid, Ua, Va, Wa)
            call metaB%value ( xgrid, ygrid, zgrid, Ub, Vb, Wb)
            self%wind(1) = Ua*(time_interpolation_factor) + Ub*(1.0-time_interpolation_factor)
            self%wind(2) = Va*(time_interpolation_factor) + Vb*(1.0-time_interpolation_factor)
            self%wind(3) = Wa*(time_interpolation_factor) + Wb*(1.0-time_interpolation_factor)
        else
            call metaA%value ( xgrid, ygrid, zgrid, &
                 &             self%wind(1), self%wind(2), self%wind(3) )
        endif

        !
        ! Now knowing the wind at time self%time_history(self%pathcount), we can now update
        ! self%wind_history, self%ground_velocity_history
        ! for time==self%time_history(self%pathcount)
        !

        !  Update self%wind_history, self%ground_velocity_history:
        call self%compute_at_time_and_location(aircraft, options%helicopter) 

    enddo PLEG
  end subroutine flightpath_precompute

  !
  !-------------------------------------------------------------------------------
  !

end module module_flightpath

