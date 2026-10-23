module module_external_attitude
  use module_configuration, only : RKIND
  implicit none
  private
  public :: read_attitude_file
  public :: attitude_details_at_time
  logical,             public     :: use_external_attitudes
  character(len=2048), public     :: attitude_file
  real(kind=RKIND),                public     :: attitude_orientation_rotate_degrees
  real(kind=RKIND), allocatable, dimension(:) :: time
  real(kind=RKIND), allocatable, dimension(:) :: roll
  real(kind=RKIND), allocatable, dimension(:) :: pitch
  real(kind=RKIND), allocatable, dimension(:) :: drift
  real(kind=RKIND), allocatable, dimension(:) :: heading
  real(kind=RKIND), allocatable, dimension(:) :: altitude
  real(kind=RKIND), allocatable, dimension(:) :: eastward_velocity
  real(kind=RKIND), allocatable, dimension(:) :: northward_velocity
  real(kind=RKIND), allocatable, dimension(:) :: vertical_velocity

  real(kind=RKIND), public :: start_with_time = -1.E36
  real(kind=RKIND), public :: stop_with_time = -1.E36

  integer :: search_start_index = 1
  integer :: replication_factor = 2

contains

  !
  !===============================================================================================================
  !

  subroutine read_attitude_file ( leg_initial_time, start_at, end_at )
    use netcdf, only : NF90_NOWRITE
    use netcdf, only : nf90_open
    use netcdf, only : nf90_close
    use netcdf, only : nf90_inq_dimid
    use netcdf, only : nf90_inquire_dimension
    use netcdf, only : nf90_inq_varid
    use netcdf, only : nf90_get_var
    use module_llxy, only : rad_per_deg
    implicit none
    real(kind=RKIND),             intent(in) :: leg_initial_time
    real(kind=RKIND), optional,   intent(in) :: start_at    ! in a zero-relative <time> array
    real(kind=RKIND), optional,   intent(in) :: end_at      ! in a zero-relative <time> array
    integer :: ncid
    integer :: dimid
    integer :: varid
    integer :: dimlen_time
    integer :: k
    real(kind=RKIND) :: dxdt, dydt
    real(kind=RKIND) :: rotation_angle_radians
    call error_handler ( nf90_open(trim(attitude_file), NF90_NOWRITE, ncid),"Problem open file '"//trim(adjustl(attitude_file))//"'" )
    call error_handler ( nf90_inq_dimid(ncid,"time",dimid),"Problem inq_dimid" )
    call error_handler ( nf90_inquire_dimension(ncid,dimid,len=dimlen_time),"Problem inquire_dimension" )
    allocate (time(dimlen_time*replication_factor))
    allocate (roll(dimlen_time*replication_factor))
    allocate (pitch(dimlen_time*replication_factor))
    allocate (drift(dimlen_time*replication_factor))
    allocate (heading(dimlen_time*replication_factor))
    allocate (altitude(dimlen_time*replication_factor))
    allocate (eastward_velocity(dimlen_time*replication_factor))
    allocate (northward_velocity(dimlen_time*replication_factor))
    allocate (vertical_velocity(dimlen_time*replication_factor))

    call error_handler ( nf90_inq_varid(ncid,"time",varid), "Problem inq_varid" )
    call error_handler ( nf90_get_var(ncid,varid,time(1:dimlen_time)), "Problem get_var" )

    call error_handler ( nf90_inq_varid(ncid,"roll",varid), "Problem inq_varid" )
    call error_handler ( nf90_get_var(ncid,varid,roll(1:dimlen_time)), "Problem get_var" )

    call error_handler ( nf90_inq_varid(ncid,"pitch",varid), "Problem inq_varid" )
    call error_handler ( nf90_get_var(ncid,varid,pitch(1:dimlen_time)), "Problem get_var" )

    call error_handler ( nf90_inq_varid(ncid,"drift",varid), "Problem inq_varid" )
    call error_handler ( nf90_get_var(ncid,varid,drift(1:dimlen_time)), "Problem get_var" )

    call error_handler ( nf90_inq_varid(ncid,"heading",varid), "Problem inq_varid" )
    call error_handler ( nf90_get_var(ncid,varid,heading(1:dimlen_time)), "Problem get_var" )

    call error_handler ( nf90_inq_varid(ncid,"altitude",varid), "Problem inq_varid" )
    call error_handler ( nf90_get_var(ncid,varid,altitude(1:dimlen_time)), "Problem get_var" )

    call error_handler ( nf90_inq_varid(ncid,"eastward_velocity",varid), "Problem inq_varid" )
    call error_handler ( nf90_get_var(ncid,varid,eastward_velocity(1:dimlen_time)), "Problem get_var" )

    call error_handler ( nf90_inq_varid(ncid,"northward_velocity",varid), "Problem inq_varid" )
    call error_handler ( nf90_get_var(ncid,varid,northward_velocity(1:dimlen_time)), "Problem get_var" )

    call error_handler ( nf90_inq_varid(ncid,"vertical_velocity",varid), "Problem inq_varid" )
    call error_handler ( nf90_get_var(ncid,varid,vertical_velocity(1:dimlen_time)), "Problem get_var" )

    call error_handler ( nf90_close(ncid), "Problem close" )

    !
    ! Make the values of the time array relative to the leg start time.
    !

    time = time + leg_initial_time
    if ( present ( start_at ) ) then
        time = time - start_at
        start_with_time = ( start_at + leg_initial_time )
    else
        start_with_time = leg_initial_time
    endif

    !
    ! Rotate everything by <rotation_angle> degrees
    ! attitude_orientation_rotate_degrees = 182.0
    !

    heading = heading + attitude_orientation_rotate_degrees
    where ( heading >= 360 )
        heading = heading - 360
    endwhere
    
    rotation_angle_radians = attitude_orientation_rotate_degrees * RAD_PER_DEG
    do k = 1, dimlen_time
        dxdt =  eastward_velocity(k) 
        dydt = northward_velocity(k)
        ! eastward_velocity(k)  = dxdt * cos(rotation_angle_radians) - dydt*sin(rotation_angle_radians)
        ! northward_velocity(k) = dxdt * sin(rotation_angle_radians) + dydt*cos(rotation_angle_radians)
        ! print '(4F10.4)', dxdt, dydt, eastward_velocity(k), northward_velocity(k)
        eastward_velocity(k)  =   dxdt * cos(rotation_angle_radians) + dydt*sin(rotation_angle_radians)
        northward_velocity(k) = - dxdt * sin(rotation_angle_radians) + dydt*cos(rotation_angle_radians)
    enddo

    if (replication_factor > 1) then
        do k = 2, replication_factor
            roll              ( (dimlen_time*(k-1))+1 : dimlen_time*k) = roll              (1:dimlen_time)
            pitch             ( (dimlen_time*(k-1))+1 : dimlen_time*k) = pitch             (1:dimlen_time)
            drift             ( (dimlen_time*(k-1))+1 : dimlen_time*k) = drift             (1:dimlen_time)
            heading           ( (dimlen_time*(k-1))+1 : dimlen_time*k) = heading           (1:dimlen_time)
            altitude          ( (dimlen_time*(k-1))+1 : dimlen_time*k) = altitude          (1:dimlen_time)
            eastward_velocity ( (dimlen_time*(k-1))+1 : dimlen_time*k) = eastward_velocity (1:dimlen_time)
            northward_velocity( (dimlen_time*(k-1))+1 : dimlen_time*k) = northward_velocity(1:dimlen_time)
            vertical_velocity ( (dimlen_time*(k-1))+1 : dimlen_time*k) = vertical_velocity (1:dimlen_time)
            
            ! Time is a bit trickier:
            time((dimlen_time*(k-1))+1:dimlen_time*k)=(time(1:dimlen_time)-time(1)) + time(dimlen_time*(k-1))
        enddo

    endif

    if ( present ( end_at ) ) then
        stop_with_time = (end_at - start_at) + leg_initial_time
    else
        stop_with_time = time(size(time))
    endif
    
  end subroutine read_attitude_file

  !
  !===============================================================================================================
  !

  subroutine attitude_details_at_time ( request_time , &
       &   roll_return, pitch_return, drift_return, heading_return, &
       &   dxdt_return, dydt_return, dzdt_return, ierr )
    implicit none
    real(kind=RKIND),                         intent(in)  :: request_time
    real(kind=RKIND),                         intent(out) :: roll_return
    real(kind=RKIND),                         intent(out) :: pitch_return
    real(kind=RKIND),                         intent(out) :: drift_return
    real(kind=RKIND),                         intent(out) :: heading_return
    real(kind=RKIND),                         intent(out) :: dxdt_return
    real(kind=RKIND),                         intent(out) :: dydt_return
    real(kind=RKIND),                         intent(out) :: dzdt_return
    integer,                      intent(out) :: ierr

    real(kind=RKIND)    :: fr
    integer :: i
    ierr = 0

    ! print*, 'request_time = ', request_time
    ! print*, 'time(1) = ', time(1)
    ! print*, 'time(<end>) = ', time(size(time))
    ! print*, 'start_with_time = ', start_with_time
    ! print*, 'stop_with_time = ', stop_with_time

    if ( request_time > stop_with_time ) then
        write(*,'( "Beyond desired range" )')
        write(*,*) "request_time = ", request_time
        write(*,*) "stop_with_time = ", stop_with_time
        ierr = 1
        return 
    endif
    if ( request_time > time(size(time)) ) then
        write(*,'("Beyond Time")')
        ierr = 1
        return
    endif
    if ( request_time < time(search_start_index) ) then
        print*, 'request_time = ', request_time
        print*, 'time(1) = ', time(1)
        print*, 'start_with_time = ', start_with_time
        write(*,'("Confusion: before time")')
        ierr = 1
        return
    endif

    ! Interpolate to the given time <t> from the array <time>

    do i = search_start_index, size(time)-1
        if ( request_time >= time(i) ) then
            fr = ( request_time - time(i) ) / ( time(i+1) - time(i) )
            roll_return    = (1.0-fr) *    roll(i) + fr *    roll(i+1)
            pitch_return   = (1.0-fr) *   pitch(i) + fr *   pitch(i+1)
            drift_return   = (1.0-fr) *   drift(i) + fr *   drift(i+1)
            heading_return = (1.0-fr) * heading(i) + fr * heading(i+1)
            dxdt_return = (1.0-fr) *  eastward_velocity(i) + fr *  eastward_velocity(i+1)
            dydt_return = (1.0-fr) * northward_velocity(i) + fr * northward_velocity(i+1)
            dzdt_return = (1.0-fr) *  vertical_velocity(i) + fr *  vertical_velocity(i+1)
        endif
    enddo

  end subroutine attitude_details_at_time
    
  !
  !===============================================================================================================
  !

  subroutine error_handler(ierr,message)
    use netcdf, only : NF90_NOERR
    use netcdf, only : nf90_strerror
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: message
    if ( ierr == NF90_NOERR ) return
    write(*,'("NETCDF>> ",A)') nf90_strerror(ierr)
    write(*,'(A)') trim(adjustl(message))
    stop
  end subroutine error_handler

  !
  !===============================================================================================================
  !

end module module_external_attitude
