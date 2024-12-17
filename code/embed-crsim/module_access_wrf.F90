module module_access_wrf
  use module_configuration, only : RKIND
  use module_llxy, only : PI          ! Added for the variable BW (B. Klotz, 12/18/2023)
  
  implicit none

  private

  integer, public, parameter :: MAX_WRF_FILES = 300
  integer, public, parameter :: MAX_WRF_TIMES = 300
  real(kind=RKIND),parameter :: PI2 = PI*PI
  public :: find_all_wrf_times

  type, public :: waypoint_type
      real(kind=RKIND) :: x_grid   ! x in gridpoint coordinate system, (reference 1)
      real(kind=RKIND) :: y_grid   ! y in gridpoint coordinate system, (reference 1)
      real(kind=RKIND) :: z_meters ! z in meters (AGL)
      real(kind=RKIND) :: p_Pa     ! pressure in Pa
  end type waypoint_type

  type, public :: options_type
      character(len=1024) :: wrf_glob_pattern
      character(len=1024) :: output_filename_format_string
      character(len=1) :: flight_level_coordinate
      type(waypoint_type), allocatable, dimension(:) :: waypoint  ! When allocated with size n, the range should be [0:n-1]
      real(kind=RKIND) :: air_speed
      real(kind=RKIND) :: leg_initial_time
      real(kind=RKIND) :: leg_time_seconds
      integer          :: bwtype  ! Constant or variable beamwdith;
                                  !     0 = does not apply the beamwidth to the calculation
                                  !     1 = applies the constant beamwidth to the calculation
                                  !     2 = applies the variable beamwidth to the calculation
      integer :: conv_minute  ! conversion of XTIME to seconds as nessary
      real(kind=RKIND) :: ref_angle  ! Reference angle for the aircraft fuselage and panel frame of reference; Sides are either 270 or 90, top is 0 and bottom is 180
      logical :: time_evolution
      logical :: helicopter
      logical :: herky_jerky
  end type options_type

  type, public :: wrf_metadata_type
      character(len=1024) :: flnm = " "
      integer :: time_frame
      integer :: ncid
      integer :: ni
      integer :: nj
      integer :: nk
      integer :: ntime
      integer :: varid_u
      integer :: varid_v
      integer :: varid_w
      real(kind=RKIND), allocatable, dimension(:,:,:) :: u
      real(kind=RKIND), allocatable, dimension(:,:,:) :: v
      real(kind=RKIND), allocatable, dimension(:,:,:) :: w
      real(kind=RKIND), allocatable, dimension(:,:,:) :: p
      real(kind=RKIND), allocatable, dimension(:,:,:) :: zf
      real(kind=RKIND), allocatable, dimension(:,:,:) :: zh
      real(kind=RKIND), allocatable, dimension(:) :: xh
      real(kind=RKIND), allocatable, dimension(:) :: yh
      real(kind=RKIND) :: zmin
      real(kind=RKIND) :: zmax
      real(kind=RKIND) :: ymin
      real(kind=RKIND) :: ymax
      real(kind=RKIND) :: xmin
      real(kind=RKIND) :: xmax
      real(kind=RKIND) :: rdx
      real(kind=RKIND) :: rdy
      integer :: nfields
      character(len=256), dimension(1000) :: list_of_fields
      logical :: debug_flag = .false.
    contains
      procedure :: open  => wrf_open
      procedure :: open_precompute  => wrf_open_precompute
      procedure :: close => wrf_close
      procedure :: print => wrf_print
      procedure :: value => wrf_value
      procedure :: interp => wrf_interp
      procedure :: Ku => wrf_Ku
      procedure :: Kv => wrf_Kv
      procedure :: Kw => wrf_Kw
      procedure :: RHO_D => wrf_RHO_D
      procedure :: interp_varbw => wrf_interp_varbw
      procedure :: Ku_varbw => wrf_Ku_varbw
      procedure :: Kv_varbw => wrf_Kv_varbw
      procedure :: Kw_varbw => wrf_Kw_varbw
      procedure :: RHO_D_varbw => wrf_RHO_D_varbw
      procedure :: debug_on => wrf_debug_on
      procedure :: debug_off => wrf_debug_off
      procedure :: find_z_index => wrf_find_z_index
      procedure :: get_z_at_p => wrf_get_z_at_p
  end type wrf_metadata_type
  ! type(wrf_metadata_type), public :: metaA
  ! type(wrf_metadata_type), public :: metaB
contains

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine find_all_wrf_times(wrf_glob_pattern, conv_minute, wrf_reference_time, wrf_files, wrf_times, wrf_xtimes, wrf_file_index, wrf_time_index, dx, idim, jdim)
    use netcdf, only : NF90_NOWRITE
    use netcdf, only : NF90_MPIIO
    use netcdf, only : NF90_NETCDF4
    use netcdf, only : NF90_GLOBAL
    use netcdf, only : NF90_MAX_NAME
    use netcdf, only : NF90_MAX_VAR_DIMS
    use netcdf, only : nf90_open
    use netcdf, only : nf90_inq_varid
    use netcdf, only : nf90_inquire_dimension
    use netcdf, only : nf90_inquire_variable
    use netcdf, only : nf90_get_var
    use netcdf, only : nf90_get_att
    use netcdf, only : nf90_close
    use kwm_date_utilities, only : geth_newdate

    implicit none
    character(len=*),                              intent(in)  :: wrf_glob_pattern
    integer,                                       intent(in)  :: conv_minute
    character(len=19),                             intent(out) :: wrf_reference_time
    character(len=1024), dimension(MAX_WRF_FILES), intent(out) :: wrf_files
    real(kind=RKIND),    dimension(MAX_WRF_FILES*MAX_WRF_TIMES), intent(out) :: wrf_xtimes
    character(len=19),   dimension(MAX_WRF_FILES*MAX_WRF_TIMES), intent(out) :: wrf_times
    integer,             dimension(MAX_WRF_FILES*MAX_WRF_TIMES), intent(out) :: wrf_file_index
    integer,             dimension(MAX_WRF_FILES*MAX_WRF_TIMES), intent(out) :: wrf_time_index
    real(kind=RKIND),                                            intent(out) :: dx
    integer,                                                     intent(out) :: idim
    integer,                                                     intent(out) :: jdim
    character(len=1024) :: string
    integer :: ierr
    integer :: ncid
    integer :: varid
    integer :: i
    integer :: j
    integer :: wrf_ntimes
    integer :: index_count
    integer :: ndims
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
    character(len=NF90_MAX_NAME) :: dimname
    character(len=1024) :: xtime_units_string
    integer :: test_index
    integer :: local_conv_minute

    wrf_files = " "
    wrf_xtimes = -1.E36
    wrf_times = " "
    wrf_file_index = -999999
    wrf_time_index = -999999
    wrf_reference_time = " "
    index_count = 0

    call execute_command_line("ls "//trim(wrf_glob_pattern)//" > .thing")

    open(10,file=".thing", status='old', form='formatted', action='read')
    do i = 1, MAX_WRF_FILES
        read(10,'(A)',end=999) string
        wrf_files(i) = trim(string)

        ierr = nf90_open(trim(string), NF90_NOWRITE, ncid)
        call error_handler ( ierr , "Problem opening file "//trim(string) )

        ierr = nf90_inq_varid ( ncid, "XTIME", varid )
        call error_handler ( ierr, "Problem inquiring varid XTIME" )

        local_conv_minute = conv_minute
        if (local_conv_minute < 0) then
           ierr = nf90_get_att(ncid, varid, 'units', xtime_units_string)
           call error_handler ( ierr, "Problem getting variabl attribute 'units' for XTIME" )

           test_index = index(xtime_units_string, "seconds")
           if (test_index > 0) then
              local_conv_minute = 1
           else
              test_index = index(xtime_units_string, "minutes")
              if (test_index > 0) then
                 local_conv_minute = 60
              endif
           endif
           if (local_conv_minute < 0) then
              stop "CANNOT DETERMINE UNITS FOR XTIME.  ADD CONV_MINUTE TO NAMELIST"
           endif
        endif

        if ( wrf_reference_time == " " ) then
            ierr = nf90_get_att ( ncid, NF90_GLOBAL, "SIMULATION_START_DATE", string)
            call error_handler ( ierr, "Problem get var att SIMULATION_START_DATE" )
            wrf_reference_time = trim(string)
            ! Well, if the first digit of the year is zero, assume it's an idealized case,
            ! and flip that first digit to "2"
            wrf_reference_time(1:1) = "2"
        endif

        ierr = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
        call error_handler ( ierr, "Problem inquire variable XTIME" )

        ! Times should be a 1-d array, only a single dimension
        
        ierr = nf90_inquire_dimension(ncid, dimids(1), name=dimname, len=wrf_ntimes)
        call error_handler ( ierr, "Problem inquire dimension for XIME" )

        ! write(*,'("Dimension ''",A,"'' :",I6)') trim(dimname), wrf_ntimes

        ierr = nf90_get_var ( ncid, varid, wrf_xtimes(index_count+1:index_count+wrf_ntimes) )
        call error_handler ( ierr, "Problem get var XTIME" )
        
        ! Convert wrf_xtimes from minutes to seconds.
        ! This now accounts for files already reporting in seconds
        wrf_xtimes(index_count+1:index_count+wrf_ntimes) = wrf_xtimes(index_count+1:index_count+wrf_ntimes) * local_conv_minute

        do j = 1, wrf_ntimes
            index_count = index_count + 1
            call geth_newdate ( wrf_times(index_count), wrf_reference_time, int(wrf_xtimes(index_count)) )
            wrf_file_index(index_count) = i
            wrf_time_index(index_count) = j
        enddo

        ! While we're here, get some details about the grid.
        call error_handler ( nf90_get_att ( ncid, NF90_GLOBAL, "DX", dx ), "Problem get attribute DX" )
        call error_handler ( nf90_get_att ( ncid, NF90_GLOBAL, "WEST-EAST_GRID_DIMENSION", idim ), &
             "Problem get attribute WEST-EAST_GRID_DIMENSION" )
        call error_handler ( nf90_get_att ( ncid, NF90_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION", jdim ), &
             "Problem get attribute SOUTH-NORTH_GRID_DIMENSION" )

        ierr = nf90_close(ncid)
        call error_handler ( ierr, "Problem close")

    enddo
999 continue
    close(10, status='delete')

  end subroutine find_all_wrf_times

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_open(self, wrf_filename, time_frame)
    use netcdf, only : NF90_NOWRITE
    use netcdf, only : NF90_NETCDF4
    use netcdf, only : NF90_SHARE
    use netcdf, only : NF90_FLOAT
    use netcdf, only : NF90_MAX_NAME
    use netcdf, only : NF90_MAX_DIMS
    use netcdf, only : NF90_MAX_VAR_DIMS
    use netcdf, only : nf90_open
    use netcdf, only : nf90_inq_dimid
    use netcdf, only : nf90_inquire
    use netcdf, only : nf90_inquire_dimension
    use netcdf, only : nf90_inquire_variable
    use netcdf, only : nf90_inq_varid
    use netcdf, only : nf90_get_var
    implicit none

    class(wrf_metadata_type), intent(out) :: self
    character(len=*),         intent(in)  :: wrf_filename
    integer,                  intent(in)  :: time_frame

    integer :: ierr
    integer :: west_east_dimid, south_north_dimid, bottom_top_dimid, time_dimid
    integer :: south_north_stag_dimid, west_east_stag_dimid, bottom_top_stag_dimid
    integer :: varid
    integer :: k
    integer :: nVariables
    integer :: nDims
    integer :: xType
    real(kind=RKIND), allocatable, dimension(:,:,:) :: ph
    real(kind=RKIND), allocatable, dimension(:,:,:) :: pb
    character(len=NF90_MAX_NAME) :: varname
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimid
    integer :: myrank, ierrmpi, comm, prc
   
    self%flnm = wrf_filename
    self%time_frame = time_frame

    ierr = nf90_open(wrf_filename, NF90_NOWRITE, self%ncid)
    call error_handler ( ierr , "Problem opening file "//trim(wrf_filename) )

    write(*,'("Opening WRF file ",A, " and time frame ", I3, " as ncid ", I10)') trim(self%flnm), self%time_frame, self%ncid

    !
    ! Read dimensions <ni>, <nj>, <nk>, <ntime>
    !

    ierr = nf90_inq_dimid ( self%ncid, "west_east", west_east_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'west_east'" )
    ierr = nf90_inquire_dimension ( self%ncid, west_east_dimid, len=self%ni )
    call error_handler ( ierr, "Problem inquiring dimension west_east." )

    ierr = nf90_inq_dimid ( self%ncid, "west_east_stag", west_east_stag_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'west_east_stag'" )

    self%xmin = 1
    self%xmax = self%ni

    ierr = nf90_inq_dimid ( self%ncid, "south_north", south_north_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'south_north'" )
    ierr = nf90_inquire_dimension ( self%ncid, south_north_dimid, len=self%nj )
    call error_handler ( ierr, "Problem inquiring dimension south_north." )

    ierr = nf90_inq_dimid ( self%ncid, "south_north_stag", south_north_stag_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'south_north_stag'" )

    self%ymin = 1
    self%ymax = self%nj

    ierr = nf90_inq_dimid ( self%ncid, "bottom_top", bottom_top_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'bottom_top_dimid'" )
    ierr = nf90_inquire_dimension ( self%ncid, bottom_top_dimid, len=self%nk )
    call error_handler ( ierr, "Problem inquiring dimension bottom_top_dimid." )

    ierr = nf90_inq_dimid ( self%ncid, "bottom_top_stag", bottom_top_stag_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'bottom_top_stag_dimid'" )

    self%zmin = 1
    self%zmax = self%nk

    ierr = nf90_inq_dimid ( self%ncid, "Time", time_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'Time'" )
    ierr = nf90_inquire_dimension ( self%ncid, time_dimid, len=self%ntime )
    call error_handler ( ierr, "Problem inquiring dimension time." )

    call error_handler(nf90_inq_varid(self%ncid, "RDX", varid), "Problem get att RDX varid")
    call error_handler(nf90_get_var(self%ncid, varid, self%rdx), "Problem get att RDX var")

    call error_handler(nf90_inq_varid(self%ncid, "RDY", varid), "Problem get att RDY varid")
    call error_handler(nf90_get_var(self%ncid, varid, self%rdy), "Problem get att RDY var")

    !
    ! Build geopotential height field on full levels.
    !

    allocate ( self%zf ( self%ni, self%nj, self%nk+1 ) )
    ierr = nf90_inq_varid ( self%ncid, "PHB", varid )
    call error_handler ( ierr, "Problem inquiring varid PHB" )
    ierr = nf90_get_var ( self%ncid, varid, self%zf, start=(/1,1,1,time_frame/)  )
    call error_handler ( ierr, "Problem get var PHB" )

    allocate ( ph ( self%ni, self%nj, self%nk+1 ) )
    ierr = nf90_inq_varid ( self%ncid, "PH", varid )
    call error_handler ( ierr, "Problem inquiring varid PH" )
    ierr = nf90_get_var ( self%ncid, varid, ph, start=(/1,1,1,time_frame/)  )
    call error_handler ( ierr, "Problem get var PH" )

    self%zf = ( self%zf + ph ) / 9.81
    deallocate(ph)

    !
    !  Interpolate to half levels.
    !

    allocate(self%zh(self%ni,self%nj,self%nk))
    do k = 1, self%nk
        self%zh(:,:,k) = 0.5*(self%zf(:,:,k) + self%zf(:,:,k+1))
    enddo

    !
    !  Build pressure field on half levels.
    !

    allocate ( self%p ( self%ni, self%nj, self%nk ) )
    ierr = nf90_inq_varid ( self%ncid, "P", varid )
    call error_handler ( ierr, "Problem inquiring varid P" )
    ierr = nf90_get_var ( self%ncid, varid, self%p, start=(/1,1,1,time_frame/)  )
    call error_handler ( ierr, "Problem get var P" )

    allocate (  pb( self%ni, self%nj, self%nk ) )
    ierr = nf90_inq_varid ( self%ncid, "PB", varid )
    call error_handler ( ierr, "Problem inquiring varid PB" )
    ierr = nf90_get_var ( self%ncid, varid, pb, start=(/1,1,1,time_frame/)  )
    call error_handler ( ierr, "Problem get var PB" )

    self%p = ( self%p + pb )
    deallocate(pb)
    
    !
    !  Get list of 3-d fields ( defined with 4d arrays (including time) in wrfout netcdf)
    !

    ierr = nf90_inquire ( self%ncid, nVariables=nVariables )
    call error_handler ( ierr, "Problem inquiring nVariables" )

    self%nfields = 0
    do varid = 1, nVariables
        ierr = nf90_inquire_variable ( self%ncid, varid, xtype=xType, name=varname, ndims=nDims, dimids=dimid )
        call error_handler ( ierr, "Problem inquire variable" )
        if ( xType /= NF90_FLOAT ) cycle  ! Only process variables of type NF90_FLOAT
        if ( ndims /= 4          ) cycle  ! Only process 3d variables (plus time dimension)
        if ( ( dimid(1) /= west_east_dimid   ) .and. ( dimid(1) /= west_east_stag_dimid   ) ) cycle ! Make sure x dimension is as expected
        if ( ( dimid(2) /= south_north_dimid ) .and. ( dimid(2) /= south_north_stag_dimid ) ) cycle ! Make sure y dimension is as expected
        if ( ( dimid(3) /= bottom_top_dimid  ) .and. ( dimid(3) /= bottom_top_stag_dimid  ) ) cycle ! Make sure z dimension is as expected
        if (   dimid(4) /= time_dimid  ) cycle ! Make sure time dimension is as expected
        self%nfields = self%nfields + 1
        self%list_of_fields(self%nfields) = trim(varname)
    enddo

    ierr = nf90_inq_varid ( self%ncid, "U", self%varid_u )
    call error_handler ( ierr, "Problem inquiring varid u" )

    ierr = nf90_inq_varid ( self%ncid, "V", self%varid_v )
    call error_handler ( ierr, "Problem inquiring varid v" )

    ierr = nf90_inq_varid ( self%ncid, "W", self%varid_w )
    call error_handler ( ierr, "Problem inquiring varid w" )

    allocate ( self%u ( self%ni+1, self%nj, self%nk ) )
    ierr = nf90_get_var ( self%ncid, self%varid_u, self%u, start=(/1,1,1,time_frame/) )
    call error_handler ( ierr, "Problem get_var u" )

    allocate ( self%v ( self%ni, self%nj+1, self%nk ) )
    ierr = nf90_get_var ( self%ncid, self%varid_v, self%v, start=(/1,1,1,time_frame/) )
    call error_handler ( ierr, "Problem get_var v" )

    allocate ( self%w ( self%ni, self%nj, self%nk+1 ) )
    ierr = nf90_get_var ( self%ncid, self%varid_w, self%w, start=(/1,1,1,time_frame/) )
    call error_handler ( ierr, "Problem get_var w" )

  end subroutine wrf_open

  subroutine wrf_open_precompute(self, flight_level_coordinate, wrf_filename, time_frame, x1, x2, y1, y2)
    use netcdf, only : NF90_NOWRITE
    use netcdf, only : NF90_MPIIO
    use netcdf, only : NF90_NETCDF4
    use netcdf, only : NF90_SHARE
    use netcdf, only : NF90_FLOAT
    use netcdf, only : NF90_MAX_NAME
    use netcdf, only : NF90_MAX_DIMS
    use netcdf, only : NF90_MAX_VAR_DIMS
    use netcdf, only : nf90_open
    use netcdf, only : nf90_inq_dimid
    use netcdf, only : nf90_inquire
    use netcdf, only : nf90_inquire_dimension
    use netcdf, only : nf90_inquire_variable
    use netcdf, only : nf90_inq_varid
    use netcdf, only : nf90_get_var
    use netcdf, only : nf90_open_par
    implicit none

    class(wrf_metadata_type), intent(out) :: self
    character(len=1),         intent(in)  :: flight_level_coordinate
    character(len=*),         intent(in)  :: wrf_filename
    integer,                  intent(in)  :: time_frame
    real(kind=RKIND),         intent(in)  :: x1
    real(kind=RKIND),         intent(in)  :: y1
    real(kind=RKIND),         intent(in)  :: x2
    real(kind=RKIND),         intent(in)  :: y2

    integer :: ierr
    integer :: west_east_dimid, south_north_dimid, bottom_top_dimid, time_dimid
    integer :: south_north_stag_dimid, west_east_stag_dimid, bottom_top_stag_dimid
    integer :: varid
    integer :: k
    integer :: nVariables
    integer :: nDims
    integer :: xType
    real(kind=RKIND), allocatable, dimension(:,:,:) :: ph
    real(kind=RKIND), allocatable, dimension(:,:,:) :: pb
    character(len=NF90_MAX_NAME) :: varname
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimid
    integer :: myrank, ierrmpi, comm, prc
    integer :: i1, i2, j1, j2
   
    self%flnm = wrf_filename
    self%time_frame = time_frame

    ierr = nf90_open(wrf_filename, NF90_NOWRITE, self%ncid)
    call error_handler ( ierr , "Problem opening file "//trim(wrf_filename) )

    ! write(*,'("Opening WRF file ",A, I10)') trim(self%flnm), self%ncid

    !
    ! Read dimensions <ni>, <nj>, <nk>, <ntime>
    !

    ierr = nf90_inq_dimid ( self%ncid, "west_east", west_east_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'west_east'" )
    ierr = nf90_inquire_dimension ( self%ncid, west_east_dimid, len=self%ni )
    call error_handler ( ierr, "Problem inquiring dimension west_east." )

    ierr = nf90_inq_dimid ( self%ncid, "west_east_stag", west_east_stag_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'west_east_stag'" )

    self%xmin = 1
    self%xmax = self%ni

    i1 = max(int(self%xmin), floor(x1))
    i2 = min(int(self%xmax), ceiling(x2))

    ierr = nf90_inq_dimid ( self%ncid, "south_north", south_north_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'south_north'" )
    ierr = nf90_inquire_dimension ( self%ncid, south_north_dimid, len=self%nj )
    call error_handler ( ierr, "Problem inquiring dimension south_north." )

    ierr = nf90_inq_dimid ( self%ncid, "south_north_stag", south_north_stag_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'south_north_stag'" )

    self%ymin = 1
    self%ymax = self%nj
    j1 = max(int(self%ymin), floor(y1))
    j2 = min(int(self%ymax), ceiling(y2))

    ierr = nf90_inq_dimid ( self%ncid, "bottom_top", bottom_top_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'bottom_top_dimid'" )
    ierr = nf90_inquire_dimension ( self%ncid, bottom_top_dimid, len=self%nk )
    call error_handler ( ierr, "Problem inquiring dimension bottom_top_dimid." )

    ierr = nf90_inq_dimid ( self%ncid, "bottom_top_stag", bottom_top_stag_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'bottom_top_stag_dimid'" )

    self%zmin = 1
    self%zmax = self%nk

    ierr = nf90_inq_dimid ( self%ncid, "Time", time_dimid )
    call error_handler ( ierr, "Problem inq_dimid for 'Time'" )
    ierr = nf90_inquire_dimension ( self%ncid, time_dimid, len=self%ntime )
    call error_handler ( ierr, "Problem inquiring dimension time." )

    !
    ! Build geopotential height field on full levels.
    !

    ! allocate ( self%zf ( self%ni, self%nj, self%nk+1 ) )
    allocate ( self%zf ( i1:i2, j1:j2, self%nk+1 ) )
    ierr = nf90_inq_varid ( self%ncid, "PHB", varid )
    call error_handler ( ierr, "Problem inquiring varid PHB" )
    ierr = nf90_get_var ( self%ncid, varid, self%zf, start=(/i1,j1,1,time_frame/), count=(/i2+1-i1,j2+1-j1,self%nk+1,1/)  )
    call error_handler ( ierr, "Problem get var PHB" )

    allocate ( ph ( i1:i2, j1:j2, self%nk+1 ) )
    ierr = nf90_inq_varid ( self%ncid, "PH", varid )
    call error_handler ( ierr, "Problem inquiring varid PH" )
    ierr = nf90_get_var ( self%ncid, varid, ph, start=(/i1,j1,1,time_frame/), count=(/i2+1-i1,j2+1-j1,self%nk+1,1/)    )
    call error_handler ( ierr, "Problem get var PH" )

    self%zf = ( self%zf + ph ) / 9.81
    deallocate(ph)

    !
    !  Interpolate to half levels.
    !

    ! allocate(self%zh(self%ni,self%nj,self%nk))
    allocate(self%zh(i1:i2,j1:j2,self%nk))
    do k = 1, self%nk
        self%zh(:,:,k) = 0.5*(self%zf(:,:,k) + self%zf(:,:,k+1))
    enddo

    if (flight_level_coordinate == "P") then
        !
        !  Build pressure field on half levels.
        !

        allocate ( self%p ( i1:i2, j1:j2, self%nk ) )
        ierr = nf90_inq_varid ( self%ncid, "P", varid )
        call error_handler ( ierr, "Problem inquiring varid P" )
        ierr = nf90_get_var ( self%ncid, varid, self%p, start=(/i1,j1,1,time_frame/), count=(/i2+1-i1,j2+1-j1,self%nk,1/)    )
        call error_handler ( ierr, "Problem get var P" )

        allocate ( pb ( i1:i2, j1:j2, self%nk ) )
        ierr = nf90_inq_varid ( self%ncid, "PB", varid )
        call error_handler ( ierr, "Problem inquiring varid PB" )
        ierr = nf90_get_var ( self%ncid, varid, pb, start=(/i1,j1,1,time_frame/), count=(/i2+1-i1,j2+1-j1,self%nk,1/)    )
        call error_handler ( ierr, "Problem get var PB" )

        self%p = ( self%p + pb )
        deallocate(pb)
    endif
    
    !
    !  Get list of 3-d fields ( defined with 4d arrays (including time) in wrfout netcdf)
    !

    ierr = nf90_inquire ( self%ncid, nVariables=nVariables )
    call error_handler ( ierr, "Problem inquiring nVariables" )

    self%nfields = 0
    do varid = 1, nVariables
        ierr = nf90_inquire_variable ( self%ncid, varid, xtype=xType, name=varname, ndims=nDims, dimids=dimid )
        call error_handler ( ierr, "Problem inquire variable" )
        if ( xType /= NF90_FLOAT ) cycle  ! Only process variables of type NF90_FLOAT
        if ( ndims /= 4          ) cycle  ! Only process 3d variables (plus time dimension)
        if ( ( dimid(1) /= west_east_dimid   ) .and. ( dimid(1) /= west_east_stag_dimid   ) ) cycle ! Make sure x dimension is as expected
        if ( ( dimid(2) /= south_north_dimid ) .and. ( dimid(2) /= south_north_stag_dimid ) ) cycle ! Make sure y dimension is as expected
        if ( ( dimid(3) /= bottom_top_dimid  ) .and. ( dimid(3) /= bottom_top_stag_dimid  ) ) cycle ! Make sure z dimension is as expected
        if (   dimid(4) /= time_dimid  ) cycle ! Make sure time dimension is as expected
        self%nfields = self%nfields + 1
        self%list_of_fields(self%nfields) = trim(varname)
    enddo

    ierr = nf90_inq_varid ( self%ncid, "U", self%varid_u )
    call error_handler ( ierr, "Problem inquiring varid u" )

    ierr = nf90_inq_varid ( self%ncid, "V", self%varid_v )
    call error_handler ( ierr, "Problem inquiring varid v" )

    ierr = nf90_inq_varid ( self%ncid, "W", self%varid_w )
    call error_handler ( ierr, "Problem inquiring varid w" )

    ! allocate ( self%u ( self%ni+1, self%nj, self%nk ) )
    allocate ( self%u ( i1:i2+1, j1:j2, self%nk ) )
    ierr = nf90_get_var ( self%ncid, self%varid_u, self%u, start=(/i1,j1,1,time_frame/), count=(/i2+2-i1,j2+1-j1,self%nk,1/)  )
    call error_handler ( ierr, "Problem get_var u" )

    ! allocate ( self%v ( self%ni, self%nj+1, self%nk ) )
    allocate ( self%v ( i1:i2, j1:j2+1, self%nk ) )
    ierr = nf90_get_var ( self%ncid, self%varid_v, self%v, start=(/i1,j1,1,time_frame/), count=(/i2+1-i1,j2+2-j1,self%nk,1/) )
    call error_handler ( ierr, "Problem get_var v" )

    ! allocate ( self%w ( self%ni, self%nj, self%nk+1 ) )
    allocate ( self%w ( i1:i2, j1:j2, self%nk+1 ) )
    ierr = nf90_get_var ( self%ncid, self%varid_w, self%w, start=(/i1,j1,1,time_frame/), count=(/i2+1-i1, j2+1-j1, self%nk+1,1/) )
    call error_handler ( ierr, "Problem get_var w" )

  end subroutine wrf_open_precompute
  
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_print(self)
    implicit none
    class(wrf_metadata_type), intent(in) :: self
    write(*,'("ncid  = ", I10)') self%ncid
    write(*,'("flnm  = ", A)') trim(self%flnm)
    write(*,'("ni    = ", I10)') self%ni
    write(*,'("nj    = ", I10)') self%nj
    write(*,'("nk    = ", I10)') self%nk
    write(*,'("ntime = ", I10)') self%ntime
  end subroutine wrf_print

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_close(self, silent)
    use netcdf, only : nf90_close
    implicit none
    class(wrf_metadata_type), intent(inout) :: self
    logical, optional, intent(in) :: silent
    integer :: ierr

    if (self%ncid > 0) then
        if (present(silent)) then
            if (.not. silent) then
                write(*,'("Closing WRF file ",A, " and time frame ", I3)') trim(self%flnm), self%time_frame
            endif
        else
            write(*,'("Closing WRF file ",A, " and time frame ", I3)') trim(self%flnm), self%time_frame
        endif
        ierr = nf90_close(self%ncid)
        call error_handler ( ierr , "Problem closing file." )
    endif

    self%ncid = -1
    self%flnm = ""
    self%time_frame = -1

  end subroutine wrf_close

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_debug_on(self)
    implicit none
    class(wrf_metadata_type), intent(inout) :: self
    self%debug_flag = .true.
  end subroutine wrf_debug_on

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_debug_off(self)
    implicit none
    class(wrf_metadata_type), intent(inout) :: self
    self%debug_flag = .false.
  end subroutine wrf_debug_off

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  real(kind=RKIND) function wrf_find_z_index(self, x_index, y_index, z) result ( z_index )
    implicit none
    class(wrf_metadata_type), intent(in) :: self
    real(kind=RKIND), intent(in) :: x_index, y_index, z
    integer :: k,im,ip,jm,jp
    real(kind=RKIND) :: xa, xb, ya, yb, zm, zp

    z_index = -1.E36
    if ( x_index < 1 ) return
    if ( y_index < 1 ) return
    if ( x_index > self%ni ) return
    if ( y_index > self%nj ) return

    ! Interpolate z's to a column at x_index, y_index.

    im = int(x_index)
    ip = im+1

    xb = ( x_index - im )
    xa = 1.0-xb

    jm = int(y_index)
    jp = jm+1

    yb = ( y_index - jm )
    ya = 1.0-yb

    ! Find bracketing indices of z and interpolate.
    zm = xa*ya*self%zh(im,jm,1) + xa*yb*self%zh(im,jp,1) + xb*ya*self%zh(ip,jm,1) + xb*yb*self%zh(ip,jp,1)
    if ( zm > z ) return
    KSEEK : do k = 2, self%nk
        zp = xa*ya*self%zh(im,jm,k) + xa*yb*self%zh(im,jp,k) + xb*ya*self%zh(ip,jm,k) + xb*yb*self%zh(ip,jp,k)
        if ( zp > z ) then
            z_index = (k-1) + ( z - zm ) / ( zp - zm )
            return
        endif
        zm = zp
    enddo KSEEK

  end function wrf_find_z_index

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  real(kind=RKIND) function wrf_get_z_at_p(self, x_index, y_index, p, z_index) result ( z )
    implicit none
    class(wrf_metadata_type), intent(in) :: self
    real(kind=RKIND),           intent(in)  :: x_index, y_index, p
    real(kind=RKIND), optional, intent(out) :: z_index
    integer :: k,im,ip,jm,jp
    real(kind=RKIND) :: xa, xb, ya, yb, zm, zp, pm, pp

    z = -1.E36
    if ( present ( z_index ) ) z_index = -1.E36
    if ( x_index < 1 ) return
    if ( y_index < 1 ) return
    if ( x_index > self%ni - 1 ) return
    if ( y_index > self%nj - 1 ) return

    ! Interpolate z's in a column at x_index, y_index.

    im = int(x_index)
    ip = im+1

    xb = ( x_index - im )
    xa = 1.0-xb

    jm = int(y_index)
    jp = jm+1

    yb = ( y_index - jm )
    ya = 1.0-yb

    !
    !  Find bracketing indices of p and interpolate.
    !

    !  Start with finding pressure at the lowest model level.
    pm = xa*ya*self%p (im,jm,1) + xa*yb*self%p (im,jp,1) + xb*ya*self%p (ip,jm,1) + xb*yb*self%p (ip,jp,1)
    if ( p > pm ) return ! desired level is below the lowest model level (i.e., desired P is higher than the lowest model level P)
    KSEEK : do k = 2, self%nk
        pp = xa*ya*self%p(im,jm,k) + xa*yb*self%p(im,jp,k) + xb*ya*self%p(ip,jm,k) + xb*yb*self%p(ip,jp,k)
        if ( pp < p ) then
            zp = xa*ya*self%zh(im,jm,k  ) + xa*yb*self%zh(im,jp,k  ) + xb*ya*self%zh(ip,jm,k  ) + xb*yb*self%zh(ip,jp,k  )
            zm = xa*ya*self%zh(im,jm,k-1) + xa*yb*self%zh(im,jp,k-1) + xb*ya*self%zh(ip,jm,k-1) + xb*yb*self%zh(ip,jp,k-1)
            z = zm + ( (zp-zm) * (log(p/pm))/(log(pp/pm)) )
            if ( present ( z_index ) ) z_index = (k-1) + ( z - zm ) / ( zp - zm )
            ! write(*,'(6F15.5)') pm,p,pp,zm,z,zp
            return
        endif
        pm = pp
    enddo KSEEK

  end function wrf_get_z_at_p


!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_interp ( self, x , y , z , varname, out_array )
    use netcdf, only : nf90_inq_varid
    use netcdf, only : nf90_get_var
    use netcdf, only : nf90_get_att
    use netcdf, only : nf90_inquire_variable
    use netcdf, only : nf90_inquire_dimension
    implicit none
    class(wrf_metadata_type), intent(in) :: self
    character(len=*), intent(in) :: varname
    real(kind=RKIND), dimension(:,:), intent(in) :: x, y, z
    real(kind=RKIND), dimension(:,:), intent(out) :: out_array
    real(kind=RKIND), allocatable, dimension(:,:,:) :: source_array, p, pb, t, qv

    integer :: iray, igate, im, ip, jm, jp, km, kp, kidx
    real(kind=RKIND)    :: xa, xb, ya, yb, za, zb
    integer :: ierr, varid, t_varid, p_varid, pb_varid, q_varid, numdims, zdimid, nzd, nlevs, xsz, ysz
    integer, dimension(4) :: dims
    real(kind=RKIND) :: xoffs, yoffs, zoffs
    character(len=8) :: stagger
    real(kind=RKIND), parameter :: rgas = 287.04
    real(kind=RKIND), parameter :: cp = 1004.
    real(kind=RKIND), parameter :: gamma = rgas/cp

    real(kind=RKIND), parameter :: Rd = 287.058 ! Value used in crsim
    real(kind=RKIND), parameter :: eps = 0.622

    real(kind=RKIND) :: sttime, entime, actime, thr1_stime, thr1_etime
    real(kind=RKIND) :: thr2_time, fstart_omp, fend_omp
    real(kind=RKIND) :: fst_netcdf, fen_netcdf, fst_ncread, fen_ncread

    select case ( varname ) 
    case default

        call cpu_time(fst_netcdf)
        ierr = nf90_inq_varid ( self%ncid, trim(varname), varid )
        call error_handler(ierr, "Problem inquire varid for " // trim(varname) )

        ierr = nf90_inquire_variable ( self%ncid, varid, ndims = numdims, dimids = dims )
        call error_handler(ierr, "Problem inquire variable " // trim(varname) )

        ierr = nf90_inquire_dimension ( self%ncid, dims(3), len = nzd)
        call error_handler(ierr, "Problem inquire bottom-top dimension " // trim(varname) )

        ierr = nf90_get_att ( self%ncid, varid, "stagger", stagger )
        call error_handler(ierr, "Problem get stagger attribute " // trim(varname) )
        call cpu_time(fen_netcdf)
 
        !
        ! Special treatment for the usual WRF grid staggering if we find the "stagger" attribute.
        !

        call cpu_time(sttime)
        select case ( stagger )
        case default
            allocate ( source_array(self%ni,self%nj,self%nk) )
            xoffs = 0.0
            yoffs = 0.0
            zoffs = 0.0
        case ("X")
            allocate ( source_array(self%ni+1,self%nj,self%nk) )
            xoffs = 0.5
            yoffs = 0.0
            zoffs = 0.0
        case ("Y")
            allocate ( source_array(self%ni,self%nj+1,self%nk) )
            xoffs = 0.0
            yoffs = 0.5
            zoffs = 0.0
        case ("Z")
            allocate ( source_array(self%ni,self%nj,self%nk+1) )
            xoffs = 0.0
            yoffs = 0.0
            zoffs = 0.5
        end select
        call cpu_time(entime)
        actime = entime-sttime
        
        call cpu_time(fst_ncread)

        ierr = nf90_get_var ( self%ncid, varid, source_array, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get " // trim(varname) )
        !!$OMP END PARALLEL

        call cpu_time(fen_ncread)

    case ( "TEMPERATURE" )
        stop "TEMPERATURE"
        ! Special treatment if "TEMPERATURE" is called for
        stagger = ""
        xoffs = 0.0
        yoffs = 0.0
        zoffs = 0.0

        ierr = nf90_inq_varid ( self%ncid, "T", t_varid )
        call error_handler(ierr, "Problem inquire T" )

        ierr = nf90_inq_varid ( self%ncid, "P", p_varid )
        call error_handler(ierr, "Problem inquire P" )

        ierr = nf90_inq_varid ( self%ncid, "PB", pb_varid )
        call error_handler(ierr, "Problem inquire PB" )

        allocate ( source_array ( self%ni , self%nj , self%nk ) )
        allocate ( p            ( self%ni , self%nj , self%nk ) )
        allocate ( pb           ( self%ni , self%nj , self%nk ) )

        ierr = nf90_get_var ( self%ncid, t_varid, source_array, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get T")

        ierr = nf90_get_var ( self%ncid, p_varid, p, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get P" )

        ierr = nf90_get_var ( self%ncid, pb_varid, pb, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PB" )

        source_array = (source_array + 300.) * ( (p+pb) / 1.E5 ) ** gamma

        deallocate(p,pb)

    case ( "RHO_D")
        stop "RHO_D"
        stagger = ""
        xoffs = 0.0
        yoffs = 0.0
        zoffs = 0.0

        ierr = nf90_inq_varid ( self%ncid, "QVAPOR", q_varid )
        call error_handler(ierr, "Problem inquire QVAPOR" )

        ierr = nf90_inq_varid ( self%ncid, "P", p_varid )
        call error_handler(ierr, "Problem inquire P" )

        ierr = nf90_inq_varid ( self%ncid, "PB", pb_varid )
        call error_handler(ierr, "Problem inquire PB" )

        ierr = nf90_inq_varid ( self%ncid, "T", t_varid )
        call error_handler(ierr, "Problem inquire T" )

        allocate ( source_array ( self%ni , self%nj , self%nk ) )
        allocate ( p            ( self%ni , self%nj , self%nk ) )
        allocate ( pb           ( self%ni , self%nj , self%nk ) )
        allocate ( t            ( self%ni , self%nj , self%nk ) )
        allocate ( qv           ( self%ni , self%nj , self%nk ) )

        ierr = nf90_get_var ( self%ncid, t_varid, t, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get T")

        ierr = nf90_get_var ( self%ncid, p_varid, p, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get P" )

        ierr = nf90_get_var ( self%ncid, pb_varid, pb, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PB" )

        p = p+pb
        deallocate(pb)

        t = (t + 300.) * (p/1.E5) ** gamma
        source_array=(p*qv)/(eps+qv)      !  e water vapor pressure in Pa 
        source_array=(p-source_array)/(Rd*t)  !  kg/m^3    

        deallocate(p)
        deallocate(t)
        deallocate(qv)

    case ( "PRESSURE" )
        stop "PRESSURE"
        ! Special treatment if "PRESSURE" is called for
        stagger = ""
        xoffs = 0.0
        yoffs = 0.0
        zoffs = 0.0

        ierr = nf90_inq_varid ( self%ncid, "P", p_varid )
        call error_handler(ierr, "Problem inquire P" )

        ierr = nf90_inq_varid ( self%ncid, "PB", pb_varid )
        call error_handler(ierr, "Problem inquire PB" )

        allocate ( source_array ( self%ni , self%nj , self%nk ) )
        allocate ( pb           ( self%ni , self%nj , self%nk ) )

        ierr = nf90_get_var ( self%ncid, p_varid, source_array, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get P")

        ierr = nf90_get_var ( self%ncid, pb_varid, pb, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PB" )

        source_array = source_array + pb

        deallocate(pb)

    case ( "HT" )
        stop "HT"
        ! Special treatment if "HT", the geopotential height, is called for
        stagger = "Z"
        xoffs = 0.0
        yoffs = 0.0
        zoffs = 0.5

        ierr = nf90_inq_varid ( self%ncid, "PH", p_varid )
        call error_handler(ierr, "Problem inquire PH" )

        ierr = nf90_inq_varid ( self%ncid, "PHB", pb_varid )
        call error_handler(ierr, "Problem inquire PHB" )

        allocate ( source_array ( self%ni , self%nj , self%nk+1 ) )
        allocate ( pb           ( self%ni , self%nj , self%nk+1 ) )

        ierr = nf90_get_var ( self%ncid, p_varid, source_array, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PH")

        ierr = nf90_get_var ( self%ncid, pb_varid, pb, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PHB" )

        source_array = ( source_array + pb ) / 9.81

        deallocate(pb)
    ! Adding this for the Mibrandt and Yau MP scheme (BWK, 3/14/2022)
    case ( "RHO_DS")
        stop "RHO_DS"
        stagger = ""
        xoffs = 0.0
        yoffs = 0.0
        zoffs = 0.0

        ierr = nf90_inq_varid ( self%ncid, "QVAPOR", q_varid )
        call error_handler(ierr, "Problem inquire QVAPOR" )

        ierr = nf90_inq_varid ( self%ncid, "P", p_varid )
        call error_handler(ierr, "Problem inquire P" )

        ierr = nf90_inq_varid ( self%ncid, "PB", pb_varid )
        call error_handler(ierr, "Problem inquire PB" )

        ierr = nf90_inq_varid ( self%ncid, "T", t_varid )
        call error_handler(ierr, "Problem inquire T" )

        allocate ( source_array ( self%ni , self%nj , self%nk ) )
        allocate ( p            ( self%ni , self%nj , self%nk ) )
        allocate ( pb           ( self%ni , self%nj , self%nk ) )
        allocate ( t            ( self%ni , self%nj , self%nk ) )
        allocate ( qv           ( self%ni , self%nj , self%nk ) )

        ierr = nf90_get_var ( self%ncid, t_varid, t, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get T")

        ierr = nf90_get_var ( self%ncid, p_varid, p, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get P" )

        ierr = nf90_get_var ( self%ncid, pb_varid, pb, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PB" )

        p = p+pb
        deallocate(pb)

        do kidx = 1, self%nk
           if ( kidx == 1 ) then

              t(:,:,kidx) = (t(:,:,kidx) + 300.) * (p(:,:,kidx)/1.E5) ** gamma
              source_array(:,:,kidx)=(p(:,:,kidx)*qv(:,:,kidx))/(eps+qv(:,:,kidx))      !  e water vapor pressure in Pa
              source_array(:,:,kidx)=(p(:,:,kidx)-source_array(:,:,kidx))/(Rd*t(:,:,kidx))  !  kg/m^3
           else
              source_array(:,:,kidx) = source_array(:,:,1)
           end if
        end do

        deallocate(p)
        deallocate(t)
        deallocate(qv)
    end select

    out_array = -9999.

    do iray = 1, size(x,2)
        do igate = 1, size(x,1)
            if ( z(igate,iray) < -9998 ) cycle

            im = int ( x(igate,iray) + xoffs )
            ip = im+1
            if ( ip > self%ni ) cycle

            jm = int ( y(igate,iray) + yoffs )
            jp = jm+1

            km = int ( z(igate,iray) + zoffs )
            kp = km+1

            xb = ( x(igate,iray) + xoffs - im )
            xa = 1.0-xb

            yb = ( y(igate,iray) + yoffs - jm )
            ya = 1.0-yb

            zb = ( z(igate,iray) + zoffs - km )
            za = 1.0-zb

            out_array(igate,iray) = xa*ya*za*source_array(im,jm,km) + xa*ya*zb*source_array(im,jm,kp) + &
                 &                  xa*yb*za*source_array(im,jp,km) + xa*yb*zb*source_array(im,jp,kp) + &
                 &                  xb*ya*za*source_array(ip,jm,km) + xb*ya*zb*source_array(ip,jm,kp) + &
                 &                  xb*yb*za*source_array(ip,jp,km) + xb*yb*zb*source_array(ip,jp,kp)

        enddo
    enddo

  end subroutine wrf_interp

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_Ku ( self, x , y , z , out_Ku, out_u )
    implicit none
    class(wrf_metadata_type), intent(in) :: self
    real(kind=RKIND), dimension(:,:), intent(in) :: x, y, z
    real(kind=RKIND), dimension(:,:), intent(out) :: out_Ku
    real(kind=RKIND), dimension(:,:), intent(out) :: out_u

    real(kind=RKIND) :: xoffs
    integer :: iray, igate, im, ip, jm, jp, km, kp
    real(kind=RKIND)    :: xa, xb, ya, yb, za, zb
    integer :: ierr, varid

    real(kind=RKIND) :: dummm, dummp, dumpm, dumpp, dupmm, dupmp, duppm, duppp
    integer :: i2

    ! write(*,'(" WRF   build Ku, U")')

    out_Ku = -9999.
    out_u  = -9999.
    xoffs = 0.5

    do iray = 1, size(x,2)
        do igate = 1, size(x,1)
            if ( z(igate,iray) < -9998 ) cycle

            im = int ( x(igate,iray) )
            ip = im+1
            if ( ip > self%ni ) cycle
            i2 = ip+1

            jm = int ( y(igate,iray) )
            jp = jm+1

            km = int ( z(igate,iray) )
            kp = km+1

            xb = ( x(igate,iray) - im )
            xa = 1.0-xb

            yb = ( y(igate,iray) - jm )
            ya = 1.0-yb

            zb = ( z(igate,iray) - km )
            za = 1.0-zb
            
            dummm = self%u(ip,jm,km)-self%u(im,jm,km)
            dumpm = self%u(ip,jp,km)-self%u(im,jp,km)
            dummp = self%u(ip,jm,kp)-self%u(im,jm,kp)
            dumpp = self%u(ip,jp,kp)-self%u(im,jp,kp)

            dupmm = self%u(i2,jm,km)-self%u(ip,jm,km)
            duppm = self%u(i2,jp,km)-self%u(ip,jp,km)
            dupmp = self%u(i2,jm,kp)-self%u(ip,jm,kp)
            duppp = self%u(i2,jp,kp)-self%u(ip,jp,kp)
            
            out_Ku(igate,iray) = xa*ya*za*dummm + xa*ya*zb*dummp + &
                 &               xa*yb*za*dumpm + xa*yb*zb*dumpp + &
                 &               xb*ya*za*dupmm + xb*ya*zb*dupmp + &
                 &               xb*yb*za*duppm + xb*yb*zb*duppp

            out_Ku(igate,iray) = out_Ku(igate,iray) * self%rdx

            ! Get new im, ip, xa, xb for the u field itself.  We can use the jm, jp, km, kp, etc. computed above
            im = int ( x(igate,iray) + xoffs )
            ip = im+1
            xb = ( x(igate,iray) + xoffs - im )
            xa = 1.0-xb

            out_u(igate,iray) = xa*ya*za*self%u(im,jm,km) + xa*ya*zb*self%u(im,jm,kp) + &
                 &              xa*yb*za*self%u(im,jp,km) + xa*yb*zb*self%u(im,jp,kp) + &
                 &              xb*ya*za*self%u(ip,jm,km) + xb*ya*zb*self%u(ip,jm,kp) + &
                 &              xb*yb*za*self%u(ip,jp,km) + xb*yb*zb*self%u(ip,jp,kp)

        enddo
    enddo
  end subroutine wrf_Ku

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_Kv ( self, x , y , z , out_Kv, out_v )

    implicit none
    class(wrf_metadata_type), intent(in) :: self
    real(kind=RKIND), dimension(:,:), intent(in) :: x, y, z
    real(kind=RKIND), dimension(:,:), intent(out) :: out_Kv
    real(kind=RKIND), dimension(:,:), intent(out) :: out_v

    integer :: iray, igate, im, ip, jm, jp, km, kp
    real(kind=RKIND)    :: xa, xb, ya, yb, za, zb
    integer :: ierr, varid
    real(kind=RKIND) :: yoffs

    real(kind=RKIND) :: dvmmm, dvmmp, dvmpm, dvmpp, dvpmm, dvpmp, dvppm, dvppp
    integer :: j2

    ! write(*,'(" WRF   build Kv, V")')

    out_Kv = -9999.
    out_v  = -9999.
    yoffs = 0.5

    do iray = 1, size(x,2)
        do igate = 1, size(x,1)
            if ( z(igate,iray) < -9998 ) cycle

            im = int ( x(igate,iray) )
            ip = im+1

            jm = int ( y(igate,iray) )
            jp = jm+1
            if ( jp > self%nj ) cycle
            j2 = jp+1

            km = int ( z(igate,iray) )
            kp = km+1

            xb = ( x(igate,iray) - im )
            xa = 1.0-xb

            yb = ( y(igate,iray) - jm )
            ya = 1.0-yb

            zb = ( z(igate,iray) - km )
            za = 1.0-zb
            
            dvmmm = self%v(im,jp,km)-self%v(im,jm,km)
            dvmpm = self%v(im,j2,km)-self%v(im,jp,km)
            dvmmp = self%v(im,jp,kp)-self%v(im,jm,kp)
            dvmpp = self%v(im,j2,kp)-self%v(im,jp,kp)

            dvpmm = self%v(ip,jp,km)-self%v(ip,jm,km)
            dvppm = self%v(ip,j2,km)-self%v(ip,jp,km)
            dvpmp = self%v(ip,jp,kp)-self%v(ip,jm,kp)
            dvppp = self%v(ip,j2,kp)-self%v(ip,jp,kp)
            
            out_Kv(igate,iray) = xa*ya*za*dvmmm + xa*ya*zb*dvmmp + &
                 &               xa*yb*za*dvmpm + xa*yb*zb*dvmpp + &
                 &               xb*ya*za*dvpmm + xb*ya*zb*dvpmp + &
                 &               xb*yb*za*dvppm + xb*yb*zb*dvppp

            out_Kv(igate,iray) = out_Kv(igate,iray) * self%rdy
            
            ! Get new jm, jp, ya, yb for the u field itself.  We can use the jm, jp, km, kp, etc. computed above
            jm = int ( y(igate,iray) + yoffs )
            jp = jm+1
            yb = ( y(igate,iray) + yoffs - jm )
            ya = 1.0-yb

            out_v(igate,iray) = xa*ya*za*self%v(im,jm,km) + xa*ya*zb*self%v(im,jm,kp) + &
                 &              xa*yb*za*self%v(im,jp,km) + xa*yb*zb*self%v(im,jp,kp) + &
                 &              xb*ya*za*self%v(ip,jm,km) + xb*ya*zb*self%v(ip,jm,kp) + &
                 &              xb*yb*za*self%v(ip,jp,km) + xb*yb*zb*self%v(ip,jp,kp)

        enddo
    enddo
  end subroutine wrf_Kv

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_Kw ( self, x , y , z , out_Kw, out_w, out_HT )

    implicit none
    class(wrf_metadata_type), intent(in) :: self
    real(kind=RKIND), dimension(:,:), intent(in) :: x, y, z
    real(kind=RKIND), dimension(:,:), intent(out) :: out_Kw
    real(kind=RKIND), dimension(:,:), intent(out) :: out_W
    real(kind=RKIND), dimension(:,:), intent(out) :: out_HT
    ! real(kind=RKIND), allocatable, dimension(:,:,:) :: wgrid, geopt, phb

    integer :: iray, igate, im, ip, jm, jp, km, kp
    real(kind=RKIND)    :: xa, xb, ya, yb, za, zb
    integer :: ierr, varid
    real(kind=RKIND) :: zoffs

    real(kind=RKIND) :: dwmmm, dwmmp, dwmpm, dwmpp, dwpmm, dwpmp, dwppm, dwppp
    integer :: k2

    ! write(*,'(" WRF   build Kw, W, HT")')

    out_Kw = -9999.
    out_w  = -9999.
    out_HT = -9999.

    zoffs = 0.5
    do iray = 1, size(x,2)
        do igate = 1, size(x,1)
            if ( z(igate,iray) < -9998 ) cycle

            im = int ( x(igate,iray) )
            ip = im+1

            jm = int ( y(igate,iray) )
            jp = jm+1

            km = int ( z(igate,iray) )
            kp = km+1
            if ( kp > self%nk ) cycle
            k2 = kp+1

            xb = ( x(igate,iray) - im )
            xa = 1.0-xb

            yb = ( y(igate,iray) - jm )
            ya = 1.0-yb

            zb = ( z(igate,iray) - km )
            za = 1.0-zb
            
            dwmmm = (self%w(im,jm,kp)-self%w(im,jm,km))/(self%zf(im,jm,kp)-self%zf(im,jm,km))
            dwmpm = (self%w(im,jp,kp)-self%w(im,jp,km))/(self%zf(im,jp,kp)-self%zf(im,jp,km))
            dwmmp = (self%w(im,jm,k2)-self%w(im,jm,kp))/(self%zf(im,jm,k2)-self%zf(im,jm,kp))
            dwmpp = (self%w(im,jp,k2)-self%w(im,jp,kp))/(self%zf(im,jp,k2)-self%zf(im,jp,kp))

            dwpmm = (self%w(ip,jm,kp)-self%w(ip,jm,km))/(self%zf(ip,jm,kp)-self%zf(ip,jm,km))
            dwppm = (self%w(ip,jp,kp)-self%w(ip,jp,km))/(self%zf(ip,jp,kp)-self%zf(ip,jp,km))
            dwpmp = (self%w(ip,jm,k2)-self%w(ip,jm,kp))/(self%zf(ip,jm,k2)-self%zf(ip,jm,kp))
            dwppp = (self%w(ip,jp,k2)-self%w(ip,jp,kp))/(self%zf(ip,jp,k2)-self%zf(ip,jp,kp))
            
            out_Kw(igate,iray) = xa*ya*za*dwmmm + xa*ya*zb*dwmmp + &
                 &               xa*yb*za*dwmpm + xa*yb*zb*dwmpp + &
                 &               xb*ya*za*dwpmm + xb*ya*zb*dwpmp + &
                 &               xb*yb*za*dwppm + xb*yb*zb*dwppp

            km = int ( z(igate,iray) + zoffs )
            kp = km+1
            zb = ( z(igate,iray) + zoffs - km )
            za = 1.0-zb

            out_W(igate,iray) = xa*ya*za*self%w(im,jm,km) + xa*ya*zb*self%w(im,jm,kp) + &
                 &              xa*yb*za*self%w(im,jp,km) + xa*yb*zb*self%w(im,jp,kp) + &
                 &              xb*ya*za*self%w(ip,jm,km) + xb*ya*zb*self%w(ip,jm,kp) + &
                 &              xb*yb*za*self%w(ip,jp,km) + xb*yb*zb*self%w(ip,jp,kp)

            out_HT(igate,iray)= xa*ya*za*self%zf(im,jm,km) + xa*ya*zb*self%zf(im,jm,kp) + &
                 &              xa*yb*za*self%zf(im,jp,km) + xa*yb*zb*self%zf(im,jp,kp) + &
                 &              xb*ya*za*self%zf(ip,jm,km) + xb*ya*zb*self%zf(ip,jm,kp) + &
                 &              xb*yb*za*self%zf(ip,jp,km) + xb*yb*zb*self%zf(ip,jp,kp)

        enddo
    enddo
  end subroutine wrf_Kw

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_RHO_d ( self, x , y , z , out_RHO_d, out_T, out_RHO_ds )
    !
    !  Given x,y,z as arrays in [gate,ray] space, 
    !     read needed fields from an open NetCDF file (self%ncid),
    !     and compute RHO_D and TEMPERATURE fields in [gate,ray] space.
    !
    use netcdf, only : nf90_inq_varid
    use netcdf, only : nf90_get_var
    use netcdf, only : nf90_get_att
    implicit none
    class(wrf_metadata_type), intent(in) :: self
    real(kind=RKIND), dimension(:,:), intent(in) :: x, y, z
    real(kind=RKIND), dimension(:,:), intent(out) :: out_RHO_d, out_RHO_ds ! Added RHO_ds (BWK, 3/14/2022)
    real(kind=RKIND), dimension(:,:), intent(out) :: out_T
    real(kind=RKIND), allocatable, dimension(:,:,:) :: RHO_D_grid, tgrid, qv, RHO_DS_grid ! Added RHO_DS (BWK, 3/14/2022)
    ! real(kind=RKIND), allocatable, dimension(:,:,:) :: pgrid, pb

    integer :: iray, igate, im, ip, jm, jp, km, kp, kidx ! Added kidx as part of the RHO_DS work (BWK, 3/14/2022)
    real(kind=RKIND)    :: xa, xb, ya, yb, za, zb
    integer :: ierr, t_varid, q_varid
    ! integer :: p_varid, pb_varid
    real(kind=RKIND) :: xoffs, yoffs, zoffs
    character(len=8) :: stagger
    real(kind=RKIND), parameter :: rgas = 287.04
    real(kind=RKIND), parameter :: cp = 1004.
    real(kind=RKIND), parameter :: gamma = rgas/cp

    real(kind=RKIND), parameter :: Rd = 287.058 ! Value used in crsim
    real(kind=RKIND), parameter :: eps = 0.622

    ! write(*,'(" WRF   build RHO_D, TEMPERATURE")')

    xoffs = 0.0
    yoffs = 0.0
    zoffs = 0.0

    ierr = nf90_inq_varid ( self%ncid, "QVAPOR", q_varid )
    call error_handler(ierr, "Problem inquire QVAPOR" )

    ierr = nf90_inq_varid ( self%ncid, "T", t_varid )
    call error_handler(ierr, "Problem inquire T" )

    allocate ( RHO_D_grid ( self%ni , self%nj , self%nk ) )
    allocate ( tgrid      ( self%ni , self%nj , self%nk ) )
    allocate ( qv         ( self%ni , self%nj , self%nk ) )
    allocate ( RHO_DS_grid ( self%ni , self%nj , self%nk ) )

    ierr = nf90_get_var ( self%ncid, t_varid, tgrid, start=(/1,1,1,self%time_frame/) )
    call error_handler(ierr, "Problem get T")

    ierr = nf90_get_var ( self%ncid, q_varid, qv, start=(/1,1,1,self%time_frame/) )
    call error_handler(ierr, "Problem get QVAPOR")

    tgrid = (tgrid + 300.) * (self%p/1.E5) ** gamma
    RHO_D_grid=(self%p*qv)/(eps+qv)      !  e water vapor pressure in Pa
    RHO_D_grid=(self%p-RHO_D_grid)/(Rd*tgrid)  !  kg/m^3

    ! Added by BWK, 3/14/2022
    do kidx = 1, self%nk
       if ( kidx == 1 ) then
         RHO_DS_grid(:,:,kidx)=(self%p(:,:,kidx)*qv(:,:,kidx))/(eps+qv(:,:,kidx))      !  e water vapor pressure in Pa
         RHO_DS_grid(:,:,kidx)=(self%p(:,:,kidx)-RHO_DS_grid(:,:,kidx))/(Rd*tgrid(:,:,kidx))  !  kg/m^3
       else
         RHO_DS_grid(:,:,kidx) = RHO_DS_grid(:,:,1)
       endif
    enddo

    out_RHO_d  = -9999.
    out_T      = -9999.
    out_RHO_ds = -9999.

    do iray = 1, size(x,2)
        do igate = 1, size(x,1)
            if ( z(igate,iray) < -9998 ) cycle

            im = int ( x(igate,iray) + xoffs )
            ip = im+1
            if ( ip > self%ni ) cycle

            jm = int ( y(igate,iray) + yoffs )
            jp = jm+1

            km = int ( z(igate,iray) + zoffs )
            kp = km+1

            xb = ( x(igate,iray) + xoffs - im )
            xa = 1.0-xb

            yb = ( y(igate,iray) + yoffs - jm )
            ya = 1.0-yb

            zb = ( z(igate,iray) + zoffs - km )
            za = 1.0-zb


            out_RHO_d(igate,iray) = xa*ya*za*RHO_D_grid(im,jm,km) + xa*ya*zb*RHO_D_grid(im,jm,kp) + &
                 &                  xa*yb*za*RHO_D_grid(im,jp,km) + xa*yb*zb*RHO_D_grid(im,jp,kp) + &
                 &                  xb*ya*za*RHO_D_grid(ip,jm,km) + xb*ya*zb*RHO_D_grid(ip,jm,kp) + &
                 &                  xb*yb*za*RHO_D_grid(ip,jp,km) + xb*yb*zb*RHO_D_grid(ip,jp,kp)


            out_t(igate,iray) = xa*ya*za*tgrid(im,jm,km) + xa*ya*zb*tgrid(im,jm,kp) + &
                 &              xa*yb*za*tgrid(im,jp,km) + xa*yb*zb*tgrid(im,jp,kp) + &
                 &              xb*ya*za*tgrid(ip,jm,km) + xb*ya*zb*tgrid(ip,jm,kp) + &
                 &              xb*yb*za*tgrid(ip,jp,km) + xb*yb*zb*tgrid(ip,jp,kp)

            ! Added by BWK, 3/14/2022
            out_RHO_ds(igate,iray) = xa*ya*za*RHO_DS_grid(im,jm,km) + xa*ya*zb*RHO_DS_grid(im,jm,kp) + &
                 &                  xa*yb*za*RHO_DS_grid(im,jp,km) + xa*yb*zb*RHO_DS_grid(im,jp,kp) + &
                 &                  xb*ya*za*RHO_DS_grid(ip,jm,km) + xb*ya*zb*RHO_DS_grid(ip,jm,kp) + &
                 &                  xb*yb*za*RHO_DS_grid(ip,jp,km) + xb*yb*zb*RHO_DS_grid(ip,jp,kp)

        enddo
    enddo
  end subroutine wrf_RHO_d

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! THESE ARE THE VARIABLE BEAMWIDTH VERSIONS OF THE INTERPOLATIONS ABOVE !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrf_interp_varbw ( self, which, bwtype, volume, scan, grid_dx, varname, out_array )
    use netcdf, only : nf90_inq_varid
    use netcdf, only : nf90_get_var
    use netcdf, only : nf90_get_att
    use netcdf, only : nf90_inquire_variable
    use module_scanning,   only : scan_type
    use module_cfradial_output, only : volume_type
    !use module_llxy, only : PI
    implicit none
    class(wrf_metadata_type), intent(in) :: self
    character(len=1),             intent(in)    :: which
    integer,                      intent(in)    :: bwtype
    type(volume_type),            target        :: volume
    type(scan_type),              intent(in)    :: scan
    character(len=*), intent(in) :: varname
    real(kind=RKIND), dimension(:,:), intent(out) :: out_array
    real(kind=RKIND), intent(in) :: grid_dx   ! model horizontal spacing
    !real(kind=RKIND), intent(in) :: ac_xgrid, ac_ygrid, ac_z, vol_range, vol_azm, vol_elev, vol_bw_h, vol_bw_v, scan_dr  ! information from the volume structure
    real(kind=RKIND), dimension(:,:), pointer :: x, y, z    ! This is the index of the range gate
    real(kind=RKIND), dimension(:,:), pointer :: x_ll, x_lr, x_ur, x_ul    ! New x-index variables associated with the beamwidth
    real(kind=RKIND), dimension(:,:), pointer :: y_ll, y_lr, y_ur, y_ul    ! New y-index variables associated with the beamwidth
    real(kind=RKIND), dimension(:,:), pointer :: z_ll, z_lr, z_ur, z_ul    ! New z-index variables associated with the beamwidth
    real(kind=RKIND), allocatable, dimension(:,:,:) :: source_array, p, pb, t, qv

    integer :: iray, igate, im, ip, jm, jp, km, kp, kidx ! Added kidx to go with RHO_DS (BWK, 3/14/2022)
    integer :: dx, dy, dz, dxl, dyl, dzl
    real(kind=RKIND)    :: xa, xb, ya, yb, za, zb, tmp_gate_array, tmp_gate_array_cr, wgt_use1, wgt_tot1, wgt_use_cr, wgt_tot_cr
    integer :: ierr, varid, t_varid, p_varid, pb_varid, q_varid, numdims
    integer, dimension(5) :: dims
    real(kind=RKIND)      :: maxx, maxy, maxz, minx, miny, minz
    real(kind=RKIND) :: xoffs, yoffs, zoffs
    character(len=8) :: stagger
    real(kind=RKIND), parameter :: rgas = 287.04
    real(kind=RKIND), parameter :: cp = 1004.
    real(kind=RKIND), parameter :: gamma = rgas/cp

    real(kind=RKIND), parameter :: Rd = 287.058 ! Value used in crsim
    real(kind=RKIND), parameter :: eps = 0.622

    ! Variables from the CR-SIM post-processing - Added by B. Klotz (12/14/2023)
    real(kind=RKIND)       :: tmp_dx_meters, tmp_dy_meters, tmp_dz_meters
    real(kind=RKIND)       :: tmp_rad_out, tmp_azm_out, tmp_elev_out ! model point in radar coord.
    real(kind=RKIND)       :: wfac,wfacr
!   integer                :: nmeth
!
    real(kind=RKIND)       :: dr, d_r,d_az,d_el
    real(kind=RKIND)       :: d_r2,d_az2,d_el2
    real(kind=RKIND)       :: dr2,daz2,del2
    real(kind=RKIND)       :: Wr,Wa,We
    real(kind=RKIND)       :: fac

    x => volume%point_to_data("VX")
    y => volume%point_to_data("VY")
    x_ll => volume%point_to_data("VX_LL")
    y_ll => volume%point_to_data("VY_LL")
    x_lr => volume%point_to_data("VX_LR")
    y_lr => volume%point_to_data("VY_LR")
    x_ur => volume%point_to_data("VX_UR")
    y_ur => volume%point_to_data("VY_UR")
    x_ul => volume%point_to_data("VX_UL")
    y_ul => volume%point_to_data("VY_UL")
    select case ( which )
    case default
        stop
    case ("A")
        z => volume%point_to_data("ZINDXA")
        z_ll => volume%point_to_data("ZINDXA_LL")
        z_lr => volume%point_to_data("ZINDXA_LR")
        z_ur => volume%point_to_data("ZINDXA_UR")
        z_ul => volume%point_to_data("ZINDXA_UL")
    case ("B")
        z => volume%point_to_data("ZINDXB")
        z_ll => volume%point_to_data("ZINDXB_LL")
        z_lr => volume%point_to_data("ZINDXB_LR")
        z_ur => volume%point_to_data("ZINDXB_UR")
        z_ul => volume%point_to_data("ZINDXB_UL")
    end select


    dr = scan%meters_between_gates

    select case ( varname )
    case default

        ierr = nf90_inq_varid ( self%ncid, trim(varname), varid )
        call error_handler(ierr, "Problem inquire " // trim(varname) )

        ierr = nf90_get_att ( self%ncid, varid, "stagger", stagger )
        call error_handler(ierr, "Problem get stagger attribute " // trim(varname) )

        ierr = nf90_inquire_variable ( self%ncid, varid, ndims = numdims, dimids = dims )
        !print *, "Num Dims, Dim IDs = ", numdims, dims

        !
        ! Special treatment for the usual WRF grid staggering if we find the "stagger" attribute.
        !

        select case ( stagger )
        case default
            allocate ( source_array(self%ni,self%nj,self%nk) )
            xoffs = 0.0
            yoffs = 0.0
            zoffs = 0.0
        case ("X")
            allocate ( source_array(self%ni+1,self%nj,self%nk) )
            xoffs = 0.5
            yoffs = 0.0
            zoffs = 0.0
        case ("Y")
            allocate ( source_array(self%ni,self%nj+1,self%nk) )
            xoffs = 0.0
            yoffs = 0.5
            zoffs = 0.0
        case ("Z")
            allocate ( source_array(self%ni,self%nj,self%nk+1) )
            xoffs = 0.0
            yoffs = 0.0
            zoffs = 0.5
        end select

        ierr = nf90_get_var ( self%ncid, varid, source_array, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get " // trim(varname) )

    case ( "TEMPERATURE" )
        stop "varbw_TEMPERATURE"
        ! Special treatment if "TEMPERATURE" is called for
        stagger = ""
        xoffs = 0.0
        yoffs = 0.0
        zoffs = 0.0

        ierr = nf90_inq_varid ( self%ncid, "T", t_varid )
        call error_handler(ierr, "Problem inquire T" )

        ierr = nf90_inq_varid ( self%ncid, "P", p_varid )
        call error_handler(ierr, "Problem inquire P" )

        ierr = nf90_inq_varid ( self%ncid, "PB", pb_varid )
        call error_handler(ierr, "Problem inquire PB" )

        allocate ( source_array ( self%ni , self%nj , self%nk ) )
        allocate ( p            ( self%ni , self%nj , self%nk ) )
        allocate ( pb           ( self%ni , self%nj , self%nk ) )

        ierr = nf90_get_var ( self%ncid, t_varid, source_array, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get T")

        ierr = nf90_get_var ( self%ncid, p_varid, p, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get P" )

        ierr = nf90_get_var ( self%ncid, pb_varid, pb, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PB" )

        source_array = (source_array + 300.) * ( (p+pb) / 1.E5 ) ** gamma

        deallocate(p,pb)

    case ( "RHO_D")
        stop "varbw_RHO_D"
        stagger = ""
        xoffs = 0.0
        yoffs = 0.0
        zoffs = 0.0

        ierr = nf90_inq_varid ( self%ncid, "QVAPOR", q_varid )
        call error_handler(ierr, "Problem inquire QVAPOR" )

        ierr = nf90_inq_varid ( self%ncid, "P", p_varid )
        call error_handler(ierr, "Problem inquire P" )

        ierr = nf90_inq_varid ( self%ncid, "PB", pb_varid )
        call error_handler(ierr, "Problem inquire PB" )

        ierr = nf90_inq_varid ( self%ncid, "T", t_varid )
        call error_handler(ierr, "Problem inquire T" )

        allocate ( source_array ( self%ni , self%nj , self%nk ) )
        allocate ( p            ( self%ni , self%nj , self%nk ) )
        allocate ( pb           ( self%ni , self%nj , self%nk ) )
        allocate ( t            ( self%ni , self%nj , self%nk ) )
        allocate ( qv           ( self%ni , self%nj , self%nk ) )

        ierr = nf90_get_var ( self%ncid, t_varid, t, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get T")

        ierr = nf90_get_var ( self%ncid, p_varid, p, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get P" )

        ierr = nf90_get_var ( self%ncid, pb_varid, pb, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PB" )

        p = p+pb
        deallocate(pb)

        t = (t + 300.) * (p/1.E5) ** gamma
        source_array=(p*qv)/(eps+qv)      !  e water vapor pressure in Pa
        source_array=(p-source_array)/(Rd*t)  !  kg/m^3

        deallocate(p)
        deallocate(t)
        deallocate(qv)

    case ( "PRESSURE" )
        stop "varbw_PRESSURE"
        ! Special treatment if "PRESSURE" is called for
        stagger = ""
        xoffs = 0.0
        yoffs = 0.0
        zoffs = 0.0

        ierr = nf90_inq_varid ( self%ncid, "P", p_varid )
        call error_handler(ierr, "Problem inquire P" )

        ierr = nf90_inq_varid ( self%ncid, "PB", pb_varid )
        call error_handler(ierr, "Problem inquire PB" )

        allocate ( source_array ( self%ni , self%nj , self%nk ) )
        allocate ( pb           ( self%ni , self%nj , self%nk ) )

        ierr = nf90_get_var ( self%ncid, p_varid, source_array, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get P")

        ierr = nf90_get_var ( self%ncid, pb_varid, pb, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PB" )

        source_array = source_array + pb

        deallocate(pb)

    case ( "HT" )
        stop "varbw_HT"
        ! Special treatment if "HT", the geopotential height, is called for
        stagger = "Z"
        xoffs = 0.0
        yoffs = 0.0
        zoffs = 0.5

        ierr = nf90_inq_varid ( self%ncid, "PH", p_varid )
        call error_handler(ierr, "Problem inquire PH" )

        ierr = nf90_inq_varid ( self%ncid, "PHB", pb_varid )
        call error_handler(ierr, "Problem inquire PHB" )

        allocate ( source_array ( self%ni , self%nj , self%nk+1 ) )
        allocate ( pb           ( self%ni , self%nj , self%nk+1 ) )

        ierr = nf90_get_var ( self%ncid, p_varid, source_array, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PH")

        ierr = nf90_get_var ( self%ncid, pb_varid, pb, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PHB" )

        source_array = ( source_array + pb ) / 9.81

        deallocate(pb)
    case ( "RHO_DS")
        stop "varbw_DS"
        ! For the Mibrandt and Yau MP scheme
        stagger = ""
        xoffs = 0.0
        yoffs = 0.0
        zoffs = 0.0

        ierr = nf90_inq_varid ( self%ncid, "QVAPOR", q_varid )
        call error_handler(ierr, "Problem inquire QVAPOR" )

        ierr = nf90_inq_varid ( self%ncid, "P", p_varid )
        call error_handler(ierr, "Problem inquire P" )

        ierr = nf90_inq_varid ( self%ncid, "PB", pb_varid )
        call error_handler(ierr, "Problem inquire PB" )

        ierr = nf90_inq_varid ( self%ncid, "T", t_varid )
        call error_handler(ierr, "Problem inquire T" )

        allocate ( source_array ( self%ni , self%nj , self%nk ) )
        allocate ( p            ( self%ni , self%nj , self%nk ) )
        allocate ( pb           ( self%ni , self%nj , self%nk ) )
        allocate ( t            ( self%ni , self%nj , self%nk ) )
        allocate ( qv           ( self%ni , self%nj , self%nk ) )

        ierr = nf90_get_var ( self%ncid, t_varid, t, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get T")

        ierr = nf90_get_var ( self%ncid, p_varid, p, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get P" )

        ierr = nf90_get_var ( self%ncid, pb_varid, pb, start=(/1,1,1,self%time_frame/) )
        call error_handler(ierr, "Problem get PB" )

        p = p+pb
        deallocate(pb)

        do kidx = 1, self%nk
            if ( kidx == 1 ) then

                t(:,:,kidx) = (t(:,:,kidx) + 300.) * (p(:,:,kidx)/1.E5) ** gamma
                source_array(:,:,kidx)=(p(:,:,kidx)*qv(:,:,kidx))/(eps+qv(:,:,kidx))      !  e water vapor pressure in Pa
                source_array(:,:,kidx)=(p(:,:,kidx)-source_array(:,:,kidx))/(Rd*t(:,:,kidx))  !  kg/m^3
            else
                source_array(:,:,kidx) = source_array(:,:,1)
            endif
        enddo

        deallocate(p)
        deallocate(t)
        deallocate(qv)
    end select

    out_array = -9999.

    ! From CR-SIM PostProcessing - Added by B. Klotz (12/14/2023)

    wfac=-2.e0*log(2.e0)
    wfacr=pi2/wfac

    ! radar range resolution - from CR-SIM PostProcessing (B. Klotz, 12/14/2023)
    dr2  = dr*dr

    do iray = 1, size(x,2)

        ! angular resolution in az, el (radians) - Added by B. Klotz from CR-SIM PostProcessing (12/14/2023)
        daz2 = (volume%beamwidth_h(iray))*(volume%beamwidth_h(iray))
        del2 = (volume%beamwidth_v(iray))*(volume%beamwidth_v(iray))

        do igate = 1, size(x,1)
            if ( z(igate,iray) < -9998 ) cycle

            ! Creating new filter to account for the beamwidth (B. Klotz, 12/04/2023)

            tmp_gate_array    = 0.0
            tmp_gate_array_cr = 0.0
         
            ! This is the standard version as previously done - really only used if no points found with
            ! the alternative method described in the CR-SIM postprocessing tool
            wgt_use1 = 0.0
            wgt_tot1 = 0.0

            wgt_use_cr = 0.0
            wgt_tot_cr = 0.0


            ! Get the minimum and maximum indices in each direction
            maxx = MAX(x_ul(igate,iray),x_ur(igate,iray),x_lr(igate,iray),x_ll(igate,iray))
            maxy = MAX(y_ul(igate,iray),y_ur(igate,iray),y_lr(igate,iray),y_ll(igate,iray))
            maxz = MAX(z_ul(igate,iray),z_ur(igate,iray),z_lr(igate,iray),z_ll(igate,iray))

            minx = MIN(x_ul(igate,iray),x_ur(igate,iray),x_lr(igate,iray),x_ll(igate,iray))
            miny = MIN(y_ul(igate,iray),y_ur(igate,iray),y_lr(igate,iray),y_ll(igate,iray))
            minz = MIN(z_ul(igate,iray),z_ur(igate,iray),z_lr(igate,iray),z_ll(igate,iray))

            if (minz < 0) minz = 1
            if (maxz < 0) maxz = 1

            ! Distance in grid points between max and min positions
            dx = (int ( maxx + xoffs ) + 1) - (int ( minx + xoffs ) )
            dy = (int ( maxy + yoffs ) + 1) - (int ( miny + yoffs ) )
            dz = (int ( maxz + zoffs ) + 1) - (int ( minz + zoffs ) )

            !dx = 1
            !dy = 1
            !dz = 1
            ! Now let's determine the weight for each associated grid point
            ! Loops over the number of total indices (dx * dy * dz)
            !total_inds = dx * dy * dz

            im = int ( minx + xoffs )
            jm = int ( miny + yoffs )
            km = int ( minz + zoffs )
 
            !im = int (x(igate,iray) + xoffs)
            !jm = int (y(igate,iray) + yoffs)
            !km = int (z(igate,iray) + zoffs)

            !print *, "Igate, iray", igate, iray
            !print *, "Min x, min y, min z", minx, miny, minz
            !print *, "Max x, max y, max z", maxx, maxy, maxz
            !print *, "Dx, dy, dz", dx, dy, dz
            !print *, "X, Y, Z",(x(igate,iray) + xoffs),(y(igate,iray) + yoffs),(z(igate,iray) + zoffs)
            !print *, "im, jm, km", im, jm, km

            DX_LOOP : do dxl = 1, dx+1

                if (dxl .eq. 1) then
                    !im = int (x(igate,iray) + xoffs)
                    im = int (minx + xoffs)
                endif
 
                if (im < 1 .or. im > self%ni) then
                    im=im+1
                    cycle
                endif

                tmp_dx_meters  = (im - volume%aircraft_xgrid(iray)) * grid_dx

                xb = ABS( (x(igate,iray) + xoffs) - im) / (dx)
                xa = 1.0 - xb

                DY_LOOP : do dyl = 1, dy+1

                    if (dyl .eq. 1) then
                        !jm = int (y(igate,iray) + yoffs)
                        jm = int (miny + yoffs)
                    endif
       
                    if (jm < 1 .or. jm > self%nj) then
                        jm=jm+1
                        cycle
                    endif

                    tmp_dy_meters  = (jm - volume%aircraft_ygrid(iray)) * grid_dx

                    yb = ABS( (y(igate,iray) + yoffs) - jm) / (dy)
                    ya = 1.0 - yb

                    DZ_LOOP : do dzl = 1, dz+1

                        if (dzl .eq. 1) then
                            km = int (z(igate,iray) + zoffs)
                        endif
                   

                        if (km < 1 .or. km > self%nk) then
                            km=km+1
                            cycle
                        endif
        
                        tmp_dz_meters = (self%zf(im,jm,km)) - volume%aircraft_z(iray) ! Added by B. Klotz (12/14/2023)
                        ! print *, "im, jm, km, dx_m, dy_m, dz_m",im,jm,km,tmp_dx_meters,tmp_dy_meters,tmp_dz_meters

                        call WeightFuncCRS(tmp_dx_meters,tmp_dy_meters,tmp_dz_meters,volume%range(igate),volume%azimuth(iray), &
                             & volume%elevation(iray),volume%beamwidth_h(iray),volume%beamwidth_v(iray),dr,wfacr,wfac,wgt_use_cr)
                        
                   ! Getting the location of the grid point relative to the aircraft in spherical coordiantes
                   ! Added by B. Klotz (12/14/2023)
                   !!!! This was moved to the above function - can probably be removed
                   !call cart2sph(tmp_dx_meters, tmp_dy_meters, tmp_dz_meters, tmp_rad_out, tmp_azm_out, tmp_elev_out)
                   !tmp_azm_out = 90-tmp_azm_out

                   ! Added from the CR-SIM post processing code (B. Klotz, 12/14/2023)
                   ! distances between centres of radar volume and model grid
                   !d_r=tmp_rad_out-volume%range(igate)
                   !d_az=tmp_azm_out-(volume%azimuth(iray))
                   !d_el=tmp_elev_out-(volume%elevation(iray))


                   !if (abs(d_az) > 180) d_az = 360-abs(d_az)
                   !print *, "vrad,vaz,vel,rad_out,azm_out,elev_out",volume%range(igate),volume%azimuth(iray),volume%elevation(iray),tmp_rad_out,tmp_azm_out,tmp_elev_out
                   !print *, "d_r, d_az, d_el",d_r,d_az,d_el
!
                   ! squared distances between centres of radar volume and model grid
                   !d_r2=d_r*d_r
                   !d_az2=d_az*d_az
                   !d_el2=d_el*d_el

                   ! radar range, azimuth and elevation weighting function
                   !if ((abs(d_r)<=(dr+(0.25*dr))) .and. (abs(d_az)<=(volume%beamwidth_h(iray)/2.0)+(0.25*volume%beamwidth_h(iray)/2.0)) .and. (abs(d_el)<=(volume%beamwidth_v(iray)/2.0)+(0.25*volume%beamwidth_v(iray)/2.0))) then
                   !   fac=d_r2/dr2
                   !   Wr = exp(wfacr*fac)
!
                   !   fac=d_az2/daz2
                   !   Wa = exp(wfac*fac)
                   !   fac=d_el2/del2 !by oue?
                   !   We = exp(wfac*fac)

                   !   zb = ABS( (z(igate,iray) + zoffs) - km) / (dz)
                   !   za = 1.0 - zb
                  !else
                  !    Wr = 0.e0
                  !    Wa = 0.e0
                  !    We = 0.e0

                   !   za = 0.0;
                  !endif
!
                  !wgt_use_cr = Wr*Wa*We
                        if (wgt_use_cr<=1.e-5) wgt_use_cr=0.e0

                        zb = ABS( (z(igate,iray) + zoffs) - km) / (dz)
                        za = 1.0 - zb

                        ! Get the aircraft relative x,y,z position

                        wgt_use1 = xa*ya*za
                        wgt_tot1 = wgt_tot1+wgt_use1

                        wgt_tot_cr = wgt_tot_cr+wgt_use_cr

                        tmp_gate_array = tmp_gate_array + (wgt_use1*source_array(im,jm,km))
                        tmp_gate_array_cr = tmp_gate_array_cr + (wgt_use_cr*source_array(im,jm,km))
                        !print *, "im, jm, km, xa, ya, za", im, jm, km, xa, ya, za, wgt_tot
                        !print *, "Wr, We, Wa, Wtot, Wfacr", Wr, We, Wa, wgt_tot, wfacr
                  

                        km = km+1
                    enddo DZ_LOOP

                    jm = jm+1
                enddo DY_LOOP

                im = im+1
            enddo DX_LOOP
    
            ! This is a similar procedure to the CR-SIM post-processing code
            if (wgt_tot_cr > 0) then
                out_array(igate,iray) = tmp_gate_array_cr/wgt_tot_cr
            elseif (wgt_tot_cr <=0 .and. wgt_tot1 > 0) then
                out_array(igate,iray) = tmp_gate_array/wgt_tot1
            else
                out_array(igate,iray) = -9999.
            endif

        enddo
    enddo

  end subroutine wrf_interp_varbw

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_Ku_varbw ( self, which, bwtype, volume, scan, grid_dx, out_Ku, out_u )
    use module_scanning,   only : scan_type
    use module_cfradial_output, only : volume_type
    implicit none
    class(wrf_metadata_type),     intent(in)    :: self
    character(len=1),             intent(in)    :: which ! "A" or "B"
    integer,                      intent(in)    :: bwtype
    type(volume_type), target                        :: volume
    type(scan_type),              intent(in)    :: scan
    real(kind=RKIND), intent(in) :: grid_dx   ! model horizontal spacing

    real(kind=RKIND), pointer, dimension(:,:) :: x, y, z    ! This is the index of the range gate
    real(kind=RKIND), pointer, dimension(:,:) :: x_ll, x_lr, x_ur, x_ul    ! New x-index variables associated with the beamwidth
    real(kind=RKIND), pointer, dimension(:,:) :: y_ll, y_lr, y_ur, y_ul    ! New y-index variables associated with the beamwidth
    real(kind=RKIND), pointer, dimension(:,:) :: z_ll, z_lr, z_ur, z_ul    ! New z-index variables associated with the beamwidth

    real(kind=RKIND), dimension(:,:), intent(out) :: out_Ku
    real(kind=RKIND), dimension(:,:), intent(out) :: out_u

    !real(kind=RKIND) :: xoffs
    !integer :: iray, igate, im, ip, jm, jp, km, kp
    !real(kind=RKIND)    :: xa, xb, ya, yb, za, zb
    !integer :: ierr, varid

    real(kind=RKIND) :: dummm, dummp, dumpm, dumpp, dupmm, dupmp, duppm, duppp
    integer :: i2

    integer :: iray, igate, im, ip, jm, jp, km, kp, kidx ! Added kidx to go with RHO_DS (BWK, 3/14/2022)
    integer :: dx, dy, dz, dxl, dyl, dzl
    real(kind=RKIND)    :: xa, xb, ya, yb, za, zb, tmp_gate_array_u, tmp_gate_array_Ku,tmp_gate_array_cr_u, tmp_gate_array_cr_Ku, wgt_use1, wgt_tot1, wgt_use_cr, wgt_tot_cr
    integer :: ierr, varid, t_varid, p_varid, pb_varid, q_varid, numdims
    integer, dimension(5) :: dims
    real(kind=RKIND)      :: maxx, maxy, maxz, minx, miny, minz
    real(kind=RKIND)      :: xoffs, yoffs, zoffs
    character(len=8)      :: stagger


    ! Variables from the CR-SIM post-processing
    real(kind=RKIND)       :: tmp_dx_meters, tmp_dy_meters, tmp_dz_meters
    real(kind=RKIND)       :: tmp_rad_out, tmp_azm_out, tmp_elev_out ! model point in radar coord.
    real(kind=RKIND)       :: wfac,wfacr

    real(kind=RKIND)       :: dr, d_r,d_az,d_el
    real(kind=RKIND)       :: d_r2,d_az2,d_el2
    real(kind=RKIND)       :: dr2,daz2,del2
    real(kind=RKIND)       :: Wr,Wa,We
    real(kind=RKIND)       :: fac

    x => volume%point_to_data("VX")
    y => volume%point_to_data("VY")
    if (bwtype > 0) then
        x_ll => volume%point_to_data("VX_LL")
        y_ll => volume%point_to_data("VY_LL")
        x_lr => volume%point_to_data("VX_LR")
        y_lr => volume%point_to_data("VY_LR")
        x_ur => volume%point_to_data("VX_UR")
        y_ur => volume%point_to_data("VY_UR")
        x_ul => volume%point_to_data("VX_UL")
        y_ul => volume%point_to_data("VY_UL")
    endif
    select case ( which )
    case default
        stop
    case ("A")
        z => volume%point_to_data("ZINDXA")
        if (bwtype > 0) then
            z_ll => volume%point_to_data("ZINDXA_LL")
            z_lr => volume%point_to_data("ZINDXA_LR")
            z_ur => volume%point_to_data("ZINDXA_UR")
            z_ul => volume%point_to_data("ZINDXA_UL")
        endif
    case ("B")
        z => volume%point_to_data("ZINDXB")
        if (bwtype > 0) then
            z_ll => volume%point_to_data("ZINDXB_LL")
            z_lr => volume%point_to_data("ZINDXB_LR")
            z_ur => volume%point_to_data("ZINDXB_UR")
            z_ul => volume%point_to_data("ZINDXB_UL")
        endif
    end select

    dr = scan%meters_between_gates

    ! write(*,'(" WRF   build Ku, U")')

    out_Ku = -9999.
    out_u  = -9999.
    xoffs = 0.5
    yoffs = 0.0
    zoffs = 0.0

! From CR-SIM PostProcessing - Added by B. Klotz (12/14/2023)

    wfac=-2.e0*log(2.e0)
    wfacr=pi2/wfac

! radar range resolution - from CR-SIM PostProcessing (B. Klotz, 12/14/2023)
    dr2  = dr*dr

    do iray = 1, size(x,2)

! angular resolution in az, el (radians) - Added by B. Klotz from CR-SIM PostProcessing (12/14/2023)
        daz2 = (volume%beamwidth_h(iray))*(volume%beamwidth_h(iray))
        del2 = (volume%beamwidth_v(iray))*(volume%beamwidth_v(iray))

        do igate = 1, size(x,1)

! Creating new filter to account for the beamwidth (B. Klotz, 12/04/2023)

            tmp_gate_array_u     = 0.0
            tmp_gate_array_Ku    = 0.0
            tmp_gate_array_cr_u  = 0.0
            tmp_gate_array_cr_Ku = 0.0

! This is the standard version as previously done - really only used if no points found with
! the alternative method described in the CR-SIM postprocessing tool
            wgt_use1 = 0.0
            wgt_tot1 = 0.0

            wgt_use_cr = 0.0
            wgt_tot_cr = 0.0


            ! Get the minimum and maximum indices in each direction
            maxx = MAX(x_ul(igate,iray),x_ur(igate,iray),x_lr(igate,iray),x_ll(igate,iray))
            maxy = MAX(y_ul(igate,iray),y_ur(igate,iray),y_lr(igate,iray),y_ll(igate,iray))
            maxz = MAX(z_ul(igate,iray),z_ur(igate,iray),z_lr(igate,iray),z_ll(igate,iray))

            minx = MIN(x_ul(igate,iray),x_ur(igate,iray),x_lr(igate,iray),x_ll(igate,iray))
            miny = MIN(y_ul(igate,iray),y_ur(igate,iray),y_lr(igate,iray),y_ll(igate,iray))
            minz = MIN(z_ul(igate,iray),z_ur(igate,iray),z_lr(igate,iray),z_ll(igate,iray))

            if (minz < 0) minz = 1
            if (maxz < 0) maxz = 1

            ! Distance in grid points between max and min positions
            dx = (int ( maxx + xoffs ) + 1) - (int ( minx + xoffs ) )
            dy = (int ( maxy + yoffs ) + 1) - (int ( miny + yoffs ) )
            dz = (int ( maxz + zoffs ) + 1) - (int ( minz + zoffs ) )

            im = int ( minx + xoffs )
            jm = int ( miny + yoffs )
            km = int ( minz + zoffs )

            DX_LOOP : do dxl = 1, dx+1

                if (dxl .eq. 1) then
                    !im = int (x(igate,iray) + xoffs)
                    im = int (minx + xoffs)
                endif

                if (im < 1 .or. im > self%ni) then
                    im=im+1
                    cycle
                endif

                tmp_dx_meters  = (im - volume%aircraft_xgrid(iray)) * grid_dx ! Added by B. Klotz (12/14/2023)

                xb = ABS( (x(igate,iray) + xoffs) - im) / (dx)
                xa = 1.0 - xb

                DY_LOOP : do dyl = 1, dy+1

                    if (dyl .eq. 1) then
                        !jm = int (y(igate,iray) + yoffs)
                        jm = int (miny + yoffs)
                    endif

                    if (jm < 1 .or. jm > self%nj) then
                        jm=jm+1
                        cycle
                    endif

                    tmp_dy_meters  = (jm - volume%aircraft_ygrid(iray)) * grid_dx ! Added by B. Klotz (12/14/2023)

                    yb = ABS( (y(igate,iray) + yoffs) - jm) / (dy)
                    ya = 1.0 - yb

                    DZ_LOOP : do dzl = 1, dz+1

                        if (dzl .eq. 1) then
                            km = int (z(igate,iray) + zoffs)
                        endif
         

                        if (km < 1 .or. km > self%nk) then
                            km=km+1
                            cycle
                        endif

                        tmp_dz_meters = (self%zf(im,jm,km)) - volume%aircraft_z(iray) ! Added by B. Klotz (12/14/2023)
                        ! print *, "im, jm, km, dx_m, dy_m, dz_m",im,jm,km,tmp_dx_meters,tmp_dy_meters,tmp_dz_meters

                        call WeightFuncCRS(tmp_dx_meters,tmp_dy_meters,tmp_dz_meters,volume%range(igate),volume%azimuth(iray), &
                             & volume%elevation(iray),volume%beamwidth_h(iray),volume%beamwidth_v(iray),dr,wfacr,wfac,wgt_use_cr)
!
                        !wgt_use_cr = Wr*Wa*We
                        if (wgt_use_cr<=1.e-5) wgt_use_cr=0.e0

                        zb = ABS( (z(igate,iray) + zoffs) - km) / (dz)
                        za = 1.0 - zb

                        ! Get the aircraft relative x,y,z position

                        wgt_use1 = xa*ya*za
                        wgt_tot1 = wgt_tot1+wgt_use1

                        wgt_tot_cr = wgt_tot_cr+wgt_use_cr

                        dummm = self%u(im,jm,km)-self%u(int(x(igate,iray)),int(y(igate,iray)),int(z(igate,iray)))

                        tmp_gate_array_Ku = tmp_gate_array_Ku + (wgt_use1*dummm)
                        tmp_gate_array_cr_Ku = tmp_gate_array_cr_Ku + (wgt_use_cr*dummm)

                        tmp_gate_array_u = tmp_gate_array_u + (wgt_use1*self%u(im,jm,km))
                        tmp_gate_array_cr_u = tmp_gate_array_cr_u + (wgt_use_cr*self%u(im,jm,km))

                        km = km+1
                    enddo DZ_LOOP

                    jm = jm+1
                enddo DY_LOOP

                im = im+1
            enddo DX_LOOP

            ! This is a similar procedure to the CR-SIM post-processing code
            if (wgt_tot_cr > 0) then
                out_Ku(igate,iray) = (tmp_gate_array_cr_Ku * self%rdx)/wgt_tot_cr
                out_u(igate,iray)  = (tmp_gate_array_cr_u)/wgt_tot_cr
            elseif (wgt_tot_cr <=0 .and. wgt_tot1 > 0) then
                out_Ku(igate,iray) = (tmp_gate_array_Ku * self%rdx)/wgt_tot1
                out_u(igate,iray)  = (tmp_gate_array_u)/wgt_tot1
            else
                out_Ku(igate,iray) = -9999.
                out_u(igate,iray)  = -9999.
            endif

        enddo
    enddo
  end subroutine wrf_Ku_varbw

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_Kv_varbw ( self, which, bwtype, volume, scan, grid_dx, out_Kv, out_v )
    use module_scanning,   only : scan_type
    use module_cfradial_output, only : volume_type
    implicit none
    class(wrf_metadata_type), intent(in) :: self
    character(len=1),             intent(in)    :: which
    integer,                      intent(in)    :: bwtype
    type(volume_type),            target        :: volume
    type(scan_type),              intent(in)    :: scan
    real(kind=RKIND), intent(in) :: grid_dx   ! model horizontal spacing
    real(kind=RKIND), dimension(:,:), pointer :: x, y, z    ! This is the index of the range gate
    real(kind=RKIND), dimension(:,:), pointer :: x_ll, x_lr, x_ur, x_ul    ! New x-index variables associated with the beamwidth
    real(kind=RKIND), dimension(:,:), pointer :: y_ll, y_lr, y_ur, y_ul    ! New y-index variables associated with the beamwidth
    real(kind=RKIND), dimension(:,:), pointer :: z_ll, z_lr, z_ur, z_ul    ! New z-index variables associated with the beamwidth

    real(kind=RKIND), dimension(:,:), intent(out) :: out_Kv
    real(kind=RKIND), dimension(:,:), intent(out) :: out_v

    real(kind=RKIND) :: dummm, dummp, dumpm, dumpp, dupmm, dupmp, duppm, duppp
    integer :: j2

    integer :: iray, igate, im, ip, jm, jp, km, kp, kidx ! Added kidx to go with RHO_DS (BWK, 3/14/2022)
    integer :: dx, dy, dz, dxl, dyl, dzl  !(Added by B. Klotz, 12/07/2023)
    real(kind=RKIND)    :: xa, xb, ya, yb, za, zb, tmp_gate_array_v, tmp_gate_array_Kv, tmp_gate_array_cr_v, tmp_gate_array_cr_Kv, wgt_use1, wgt_tot1, wgt_use_cr, wgt_tot_cr
    integer :: ierr, varid, t_varid, p_varid, pb_varid, q_varid, numdims
    integer, dimension(5) :: dims
    real(kind=RKIND)      :: maxx, maxy, maxz, minx, miny, minz
    real(kind=RKIND)      :: xoffs, yoffs, zoffs
    character(len=8)      :: stagger


    ! Variables from the CR-SIM post-processing - Added by B. Klotz (12/14/2023)
    real(kind=RKIND)       :: tmp_dx_meters, tmp_dy_meters, tmp_dz_meters
    real(kind=RKIND)       :: tmp_rad_out, tmp_azm_out, tmp_elev_out ! model point in radar coord.
    real(kind=RKIND)       :: wfac,wfacr

    real(kind=RKIND)       :: dr, d_r,d_az,d_el
    real(kind=RKIND)       :: d_r2,d_az2,d_el2
    real(kind=RKIND)       :: dr2,daz2,del2
    real(kind=RKIND)       :: Wr,Wa,We
    real(kind=RKIND)       :: fac

    x => volume%point_to_data("VX")
    y => volume%point_to_data("VY")
    x_ll => volume%point_to_data("VX_LL")
    y_ll => volume%point_to_data("VY_LL")
    x_lr => volume%point_to_data("VX_LR")
    y_lr => volume%point_to_data("VY_LR")
    x_ur => volume%point_to_data("VX_UR")
    y_ur => volume%point_to_data("VY_UR")
    x_ul => volume%point_to_data("VX_UL")
    y_ul => volume%point_to_data("VY_UL")
    select case ( which )
    case default
        stop
    case ("A")
        z => volume%point_to_data("ZINDXA")
        z_ll => volume%point_to_data("ZINDXA_LL")
        z_lr => volume%point_to_data("ZINDXA_LR")
        z_ur => volume%point_to_data("ZINDXA_UR")
        z_ul => volume%point_to_data("ZINDXA_UL")
    case ("B")
        z => volume%point_to_data("ZINDXB")
        z_ll => volume%point_to_data("ZINDXB_LL")
        z_lr => volume%point_to_data("ZINDXB_LR")
        z_ur => volume%point_to_data("ZINDXB_UR")
        z_ul => volume%point_to_data("ZINDXB_UL")
    end select


    dr = scan%meters_between_gates

    out_Kv = -9999.
    out_v  = -9999.
    xoffs = 0.0
    yoffs = 0.5
    zoffs = 0.0

    ! From CR-SIM PostProcessing

    wfac=-2.e0*log(2.e0)
    wfacr=pi2/wfac

    ! radar range resolution - from CR-SIM PostProcessing
    dr2  = dr*dr

    do iray = 1, size(x,2)

        ! angular resolution in az, el (radians) - from CR-SIM PostProcessing
        daz2 = (volume%beamwidth_h(iray))*(volume%beamwidth_h(iray))
        del2 = (volume%beamwidth_v(iray))*(volume%beamwidth_v(iray))

        do igate = 1, size(x,1)

            ! Creating new filter to account for the beamwidth

            tmp_gate_array_v     = 0.0
            tmp_gate_array_Kv    = 0.0
            tmp_gate_array_cr_v  = 0.0
            tmp_gate_array_cr_Kv = 0.0

            ! This is the standard version as previously done - really only used if no points found with
            ! the alternative method described in the CR-SIM postprocessing tool
            wgt_use1 = 0.0
            wgt_tot1 = 0.0

            wgt_use_cr = 0.0
            wgt_tot_cr = 0.0

            ! Get the minimum and maximum indices in each direction
            maxx = MAX(x_ul(igate,iray),x_ur(igate,iray),x_lr(igate,iray),x_ll(igate,iray))
            maxy = MAX(y_ul(igate,iray),y_ur(igate,iray),y_lr(igate,iray),y_ll(igate,iray))
            maxz = MAX(z_ul(igate,iray),z_ur(igate,iray),z_lr(igate,iray),z_ll(igate,iray))

            minx = MIN(x_ul(igate,iray),x_ur(igate,iray),x_lr(igate,iray),x_ll(igate,iray))
            miny = MIN(y_ul(igate,iray),y_ur(igate,iray),y_lr(igate,iray),y_ll(igate,iray))
            minz = MIN(z_ul(igate,iray),z_ur(igate,iray),z_lr(igate,iray),z_ll(igate,iray))

            if (minz < 0) minz = 1
            if (maxz < 0) maxz = 1

            ! Distance in grid points between max and min positions
            dx = (int ( maxx + xoffs ) + 1) - (int ( minx + xoffs ) )
            dy = (int ( maxy + yoffs ) + 1) - (int ( miny + yoffs ) )
            dz = (int ( maxz + zoffs ) + 1) - (int ( minz + zoffs ) )

            im = int ( minx + xoffs )
            jm = int ( miny + yoffs )
            km = int ( minz + zoffs )

            DX_LOOP : do dxl = 1, dx+1

                if (dxl .eq. 1) then
                    !im = int (x(igate,iray) + xoffs)
                    im = int (minx + xoffs)
                endif

                if (im < 1 .or. im > self%ni) then
                    im=im+1
                    cycle
                endif

                tmp_dx_meters  = (im - volume%aircraft_xgrid(iray)) * grid_dx ! Added by B. Klotz (12/14/2023)

                xb = ABS( (x(igate,iray) + xoffs) - im) / (dx)
                xa = 1.0 - xb

                DY_LOOP : do dyl = 1, dy+1

                    if (dyl .eq. 1) then
                        !jm = int (y(igate,iray) + yoffs)
                        jm = int (miny + yoffs)
                    endif

                    if (jm < 1 .or. jm > self%nj) then
                        jm=jm+1
                        cycle
                    endif

                    tmp_dy_meters  = (jm - volume%aircraft_ygrid(iray)) * grid_dx ! Added by B. Klotz (12/14/2023)

                    yb = ABS( (y(igate,iray) + yoffs) - jm) / (dy)
                    ya = 1.0 - yb

                    DZ_LOOP : do dzl = 1, dz+1

                        if (dzl .eq. 1) then
                            km = int (z(igate,iray) + zoffs)
                        endif

                        if (km < 1 .or. km > self%nk) then
                            km=km+1
                            cycle
                        endif

                        tmp_dz_meters = (self%zf(im,jm,km)) - volume%aircraft_z(iray) ! Added by B. Klotz (12/14/2023)

                        call WeightFuncCRS(tmp_dx_meters,tmp_dy_meters,tmp_dz_meters,volume%range(igate),volume%azimuth(iray), &
                             & volume%elevation(iray),volume%beamwidth_h(iray),volume%beamwidth_v(iray),dr,wfacr,wfac,wgt_use_cr)

                        !wgt_use_cr = Wr*Wa*We
                        if (wgt_use_cr<=1.e-5) wgt_use_cr=0.e0

                        zb = ABS( (z(igate,iray) + zoffs) - km) / (dz)
                        za = 1.0 - zb

                        ! Get the aircraft relative x,y,z position

                        wgt_use1 = xa*ya*za
                        wgt_tot1 = wgt_tot1+wgt_use1

                        wgt_tot_cr = wgt_tot_cr+wgt_use_cr

                        dummm = self%v(im,jm,km)-self%v(int(x(igate,iray)),int(y(igate,iray)),int(z(igate,iray)))

                        tmp_gate_array_Kv = tmp_gate_array_Kv + (wgt_use1*dummm)
                        tmp_gate_array_cr_Kv = tmp_gate_array_cr_Kv + (wgt_use_cr*dummm)

                        tmp_gate_array_v = tmp_gate_array_v + (wgt_use1*self%v(im,jm,km))
                        tmp_gate_array_cr_v = tmp_gate_array_cr_v + (wgt_use_cr*self%v(im,jm,km))

                        km = km+1
                    enddo DZ_LOOP

                    jm = jm+1
                enddo DY_LOOP

                im = im+1
            enddo DX_LOOP

            ! This is a similar procedure to the CR-SIM post-processing code
            if (wgt_tot_cr > 0) then
                out_Kv(igate,iray) = (tmp_gate_array_cr_Kv * self%rdy)/wgt_tot_cr
                out_v(igate,iray)  = (tmp_gate_array_cr_v)/wgt_tot_cr
            elseif (wgt_tot_cr <=0 .and. wgt_tot1 > 0) then
                out_Kv(igate,iray) = (tmp_gate_array_Kv * self%rdy)/wgt_tot1
                out_v(igate,iray)  = (tmp_gate_array_v)/wgt_tot1
            else
                out_Kv(igate,iray) = -9999.
                out_v(igate,iray)  = -9999.
            endif

        enddo
    enddo
  end subroutine wrf_Kv_varbw

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_Kw_varbw ( self, which, bwtype, volume, scan, grid_dx, out_Kw, out_w, out_HT )
    use module_scanning,   only : scan_type
    use module_cfradial_output, only : volume_type
    implicit none
    class(wrf_metadata_type), intent(in) :: self
    character(len=1),             intent(in)    :: which
    integer,                      intent(in)    :: bwtype
    type(volume_type),            target        :: volume
    type(scan_type),              intent(in)    :: scan
    real(kind=RKIND), intent(in) :: grid_dx   ! model horizontal spacing
    real(kind=RKIND), dimension(:,:), pointer :: x, y, z    ! This is the index of the range gate
    real(kind=RKIND), dimension(:,:), pointer :: x_ll, x_lr, x_ur, x_ul    ! New x-index variables associated with the beamwidth
    real(kind=RKIND), dimension(:,:), pointer :: y_ll, y_lr, y_ur, y_ul    ! New y-index variables associated with the beamwidth
    real(kind=RKIND), dimension(:,:), pointer :: z_ll, z_lr, z_ur, z_ul    ! New z-index variables associated with the beamwidth

    real(kind=RKIND), dimension(:,:), intent(out) :: out_Kw
    real(kind=RKIND), dimension(:,:), intent(out) :: out_w
    real(kind=RKIND), dimension(:,:), intent(out) :: out_HT

    real(kind=RKIND) :: dummm, dummp, dumpm, dumpp, dupmm, dupmp, duppm, duppp, zfdif
    integer :: k2

    integer :: iray, igate, im, ip, jm, jp, km, kp, kidx ! Added kidx to go with RHO_DS (BWK, 3/14/2022)
    integer :: dx, dy, dz, dxl, dyl, dzl  !(Added by B. Klotz, 12/07/2023)
    real(kind=RKIND)    :: xa, xb, ya, yb, za, zb, tmp_gate_array_w, tmp_gate_array_Kw, tmp_gate_array_cr_w, tmp_gate_array_cr_Kw, tmp_gate_array_HT, tmp_gate_array_cr_HT, wgt_use1, wgt_tot1, wgt_use_cr, wgt_tot_cr
    integer :: ierr, varid, t_varid, p_varid, pb_varid, q_varid, numdims
    integer, dimension(5) :: dims
    real(kind=RKIND)      :: maxx, maxy, maxz, minx, miny, minz
    real(kind=RKIND)      :: xoffs, yoffs, zoffs
    character(len=8)      :: stagger


    ! Variables from the CR-SIM post-processing - Added by B. Klotz (12/14/2023)
    real(kind=RKIND)       :: tmp_dx_meters, tmp_dy_meters, tmp_dz_meters
    real(kind=RKIND)       :: tmp_rad_out, tmp_azm_out, tmp_elev_out ! model point in radar coord.
    real(kind=RKIND)       :: wfac,wfacr

    real(kind=RKIND)       :: dr, d_r,d_az,d_el
    real(kind=RKIND)       :: d_r2,d_az2,d_el2
    real(kind=RKIND)       :: dr2,daz2,del2
    real(kind=RKIND)       :: Wr,Wa,We
    real(kind=RKIND)       :: fac

    x => volume%point_to_data("VX")
    y => volume%point_to_data("VY")
    x_ll => volume%point_to_data("VX_LL")
    y_ll => volume%point_to_data("VY_LL")
    x_lr => volume%point_to_data("VX_LR")
    y_lr => volume%point_to_data("VY_LR")
    x_ur => volume%point_to_data("VX_UR")
    y_ur => volume%point_to_data("VY_UR")
    x_ul => volume%point_to_data("VX_UL")
    y_ul => volume%point_to_data("VY_UL")
    select case ( which )
    case default
        stop
    case ("A")
        z => volume%point_to_data("ZINDXA")
        z_ll => volume%point_to_data("ZINDXA_LL")
        z_lr => volume%point_to_data("ZINDXA_LR")
        z_ur => volume%point_to_data("ZINDXA_UR")
        z_ul => volume%point_to_data("ZINDXA_UL")
    case ("B")
        z => volume%point_to_data("ZINDXB")
        z_ll => volume%point_to_data("ZINDXB_LL")
        z_lr => volume%point_to_data("ZINDXB_LR")
        z_ur => volume%point_to_data("ZINDXB_UR")
        z_ul => volume%point_to_data("ZINDXB_UL")
    end select

    dr = scan%meters_between_gates

    out_Kw = -9999.
    out_w  = -9999.
    out_HT = -9999.

    zoffs = 0.5
    xoffs = 0.0
    yoffs = 0.0

    ! From CR-SIM PostProcessing - Added by B. Klotz (12/14/2023)

    wfac=-2.e0*log(2.e0)
    wfacr=pi2/wfac

    ! radar range resolution - from CR-SIM PostProcessing (B. Klotz, 12/14/2023)
    dr2  = dr*dr

    do iray = 1, size(x,2)

        ! angular resolution in az, el (radians) - Added by B. Klotz from CR-SIM PostProcessing (12/14/2023)
        daz2 = (volume%beamwidth_h(iray))*(volume%beamwidth_h(iray))
        del2 = (volume%beamwidth_v(iray))*(volume%beamwidth_v(iray))

        do igate = 1, size(x,1)

            ! Creating new filter to account for the beamwidth (B. Klotz, 12/04/2023)

            tmp_gate_array_w     = 0.0
            tmp_gate_array_Kw    = 0.0
            tmp_gate_array_cr_w  = 0.0
            tmp_gate_array_cr_Kw = 0.0

            ! This is the standard version as previously done - really only used if no points found with
            ! the alternative method described in the CR-SIM postprocessing tool
            wgt_use1 = 0.0
            wgt_tot1 = 0.0

            wgt_use_cr = 0.0
            wgt_tot_cr = 0.0

            ! Get the minimum and maximum indices in each direction
            maxx = MAX(x_ul(igate,iray),x_ur(igate,iray),x_lr(igate,iray),x_ll(igate,iray))
            maxy = MAX(y_ul(igate,iray),y_ur(igate,iray),y_lr(igate,iray),y_ll(igate,iray))
            maxz = MAX(z_ul(igate,iray),z_ur(igate,iray),z_lr(igate,iray),z_ll(igate,iray))

            minx = MIN(x_ul(igate,iray),x_ur(igate,iray),x_lr(igate,iray),x_ll(igate,iray))
            miny = MIN(y_ul(igate,iray),y_ur(igate,iray),y_lr(igate,iray),y_ll(igate,iray))
            minz = MIN(z_ul(igate,iray),z_ur(igate,iray),z_lr(igate,iray),z_ll(igate,iray))

            if (minz < 0) minz = 1
            if (maxz < 0) maxz = 1

            ! Distance in grid points between max and min positions
            dx = (int ( maxx + xoffs ) + 1) - (int ( minx + xoffs ) )
            dy = (int ( maxy + yoffs ) + 1) - (int ( miny + yoffs ) )
            dz = (int ( maxz + zoffs ) + 1) - (int ( minz + zoffs ) )

            im = int ( minx + xoffs )
            jm = int ( miny + yoffs )
            km = int ( minz + zoffs )

            DX_LOOP : do dxl = 1, dx+1

                if (dxl .eq. 1) then
                    !im = int (x(igate,iray) + xoffs)
                    im = int (minx + xoffs)
                endif

                if (im < 1 .or. im > self%ni) then
                    im=im+1
                    cycle
                endif

                tmp_dx_meters  = (im - volume%aircraft_xgrid(iray)) * grid_dx ! Added by B. Klotz (12/14/2023)

                xb = ABS( (x(igate,iray) + xoffs) - im) / (dx)
                xa = 1.0 - xb

                DY_LOOP : do dyl = 1, dy+1

                    if (dyl .eq. 1) then
                        !jm = int (y(igate,iray) + yoffs)
                        jm = int (miny + yoffs)
                    endif

                    if (jm < 1 .or. jm > self%nj) then
                        jm=jm+1
                        cycle
                    endif

                    tmp_dy_meters  = (jm - volume%aircraft_ygrid(iray)) * grid_dx ! Added by B. Klotz (12/14/2023)

                    yb = ABS( (y(igate,iray) + yoffs) - jm) / (dy)
                    ya = 1.0 - yb

                    DZ_LOOP : do dzl = 1, dz+1

                        if (dzl .eq. 1) then
                            km = int (z(igate,iray) + zoffs)
                        endif


                        if (km < 1 .or. km > self%nk) then
                            km=km+1
                            cycle
                        endif

                        tmp_dz_meters = (self%zf(im,jm,km)) - volume%aircraft_z(iray) ! Added by B. Klotz (12/14/2023)

                        call WeightFuncCRS(tmp_dx_meters,tmp_dy_meters,tmp_dz_meters,volume%range(igate),volume%azimuth(iray), &
                             & volume%elevation(iray),volume%beamwidth_h(iray),volume%beamwidth_v(iray),dr,wfacr,wfac,wgt_use_cr)

                        !wgt_use_cr = Wr*Wa*We
                        if (wgt_use_cr<=1.e-5) wgt_use_cr=0.e0

                        zb = ABS( (z(igate,iray) + zoffs) - km) / (dz)
                        za = 1.0 - zb

                        ! Get the aircraft relative x,y,z position

                        wgt_use1 = xa*ya*za
                        wgt_tot1 = wgt_tot1+wgt_use1

                        wgt_tot_cr = wgt_tot_cr+wgt_use_cr

                        zfdif = self%zf(im,jm,km)-self%zf(int(x(igate,iray)),int(y(igate,iray)),int(z(igate,iray)))
                        !print *,'zfdif:  ',zfdif
                        if (zfdif <= 0.0) zfdif = 1.0
                        !dummm = (self%w(im,jm,km)-self%w(x(igate,iray),y(igate,iray),z(igate,iray))) / (self%zf(im,jm,km)-self%zf(x(igate,iray),y(igate,iray),z(igate,iray)))
                        tmp_gate_array_Kw = tmp_gate_array_Kw + (wgt_use1*dummm)
                        tmp_gate_array_cr_Kw = tmp_gate_array_cr_Kw + (wgt_use_cr*dummm)

                        tmp_gate_array_w = tmp_gate_array_w + (wgt_use1*self%w(im,jm,km))
                        tmp_gate_array_cr_w = tmp_gate_array_cr_w + (wgt_use_cr*self%w(im,jm,km))

                        tmp_gate_array_HT = tmp_gate_array_HT + (wgt_use1*self%zf(im,jm,km))
                        tmp_gate_array_cr_HT = tmp_gate_array_cr_HT + (wgt_use_cr*self%zf(im,jm,km))
                        km = km+1
                    enddo DZ_LOOP

                    jm = jm+1
                enddo DY_LOOP

                im = im+1
            enddo DX_LOOP

            ! This is a similar procedure to the CR-SIM post-processing code
            if (wgt_tot_cr > 0) then
                out_Kw(igate,iray) = (tmp_gate_array_cr_Kw)/wgt_tot_cr
                out_w(igate,iray)  = (tmp_gate_array_cr_w)/wgt_tot_cr
                out_HT(igate,iray)  = (tmp_gate_array_cr_HT)/wgt_tot_cr
            elseif (wgt_tot_cr <=0 .and. wgt_tot1 > 0) then
                out_Kw(igate,iray) = (tmp_gate_array_Kw)/wgt_tot1
                out_w(igate,iray)  = (tmp_gate_array_w)/wgt_tot1
                out_HT(igate,iray)  = (tmp_gate_array_HT)/wgt_tot1
            else
                out_Kw(igate,iray) = -9999.
                out_w(igate,iray)  = -9999.
                out_HT(igate,iray)  = -9999.
            endif

  
        enddo
    enddo
  end subroutine wrf_Kw_varbw

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_RHO_d_varbw ( self, which, bwtype, volume, scan, grid_dx, out_RHO_d, out_T, out_RHO_ds)
    use netcdf, only : nf90_inq_varid
    use netcdf, only : nf90_get_var
    use netcdf, only : nf90_get_att
    use netcdf, only : nf90_inquire_variable
    use module_scanning,   only : scan_type
    use module_cfradial_output, only : volume_type
    implicit none
    class(wrf_metadata_type), intent(in) :: self
    character(len=1),             intent(in)    :: which
    integer,                      intent(in)    :: bwtype
    type(volume_type),            target        :: volume
    type(scan_type),              intent(in)    :: scan
    real(kind=RKIND), intent(in) :: grid_dx   ! model horizontal spacing
    real(kind=RKIND), dimension(:,:), pointer :: x, y, z    ! This is the index of the range gate
    real(kind=RKIND), dimension(:,:), pointer :: x_ll, x_lr, x_ur, x_ul    ! New x-index variables associated with the beamwidth
    real(kind=RKIND), dimension(:,:), pointer :: y_ll, y_lr, y_ur, y_ul    ! New y-index variables associated with the beamwidth
    real(kind=RKIND), dimension(:,:), pointer :: z_ll, z_lr, z_ur, z_ul    ! New z-index variables associated with the beamwidth
    real(kind=RKIND), dimension(:,:), intent(out) :: out_RHO_d, out_RHO_ds, out_T
    real(kind=RKIND), allocatable, dimension(:,:,:) :: RHO_D_grid, tgrid, qv, RHO_DS_grid ! Added RHO_DS (BWK, 3/14/2022)

    integer :: iray, igate, im, ip, jm, jp, km, kp, kidx ! Added kidx to go with RHO_DS (BWK, 3/14/2022)
    integer :: dx, dy, dz, dxl, dyl, dzl  !(Added by B. Klotz, 12/07/2023)
    real(kind=RKIND)    :: xa, xb, ya, yb, za, zb, tmp_gate_array_rho_d, tmp_gate_array_rho_ds, tmp_gate_array_t
    real(kind=RKIND)    :: tmp_gate_array_cr_rho_d, tmp_gate_array_cr_rho_ds, tmp_gate_array_cr_t
    real(kind=RKIND)    :: wgt_use1, wgt_tot1, wgt_use_cr, wgt_tot_cr
    integer :: ierr, varid, t_varid, p_varid, pb_varid, q_varid, numdims
    integer, dimension(5) :: dims
    real(kind=RKIND)      :: maxx, maxy, maxz, minx, miny, minz
    real(kind=RKIND) :: xoffs, yoffs, zoffs
    character(len=8) :: stagger
    real(kind=RKIND), parameter :: rgas = 287.04
    real(kind=RKIND), parameter :: cp = 1004.
    real(kind=RKIND), parameter :: gamma = rgas/cp

    real(kind=RKIND), parameter :: Rd = 287.058 ! Value used in crsim
    real(kind=RKIND), parameter :: eps = 0.622

    ! Variables from the CR-SIM post-processing - Added by B. Klotz (12/14/2023)
    real(kind=RKIND)       :: tmp_dx_meters, tmp_dy_meters, tmp_dz_meters
    real(kind=RKIND)       :: tmp_rad_out, tmp_azm_out, tmp_elev_out ! model point in radar coord.
    real(kind=RKIND)       :: wfac,wfacr

    real(kind=RKIND)       :: dr, d_r,d_az,d_el
    real(kind=RKIND)       :: d_r2,d_az2,d_el2
    real(kind=RKIND)       :: dr2,daz2,del2
    real(kind=RKIND)       :: Wr,Wa,We
    real(kind=RKIND)       :: fac

    x => volume%point_to_data("VX")
    y => volume%point_to_data("VY")
    x_ll => volume%point_to_data("VX_LL")
    y_ll => volume%point_to_data("VY_LL")
    x_lr => volume%point_to_data("VX_LR")
    y_lr => volume%point_to_data("VY_LR")
    x_ur => volume%point_to_data("VX_UR")
    y_ur => volume%point_to_data("VY_UR")
    x_ul => volume%point_to_data("VX_UL")
    y_ul => volume%point_to_data("VY_UL")
    select case ( which )
    case default
        stop
    case ("A")
        z => volume%point_to_data("ZINDXA")
        z_ll => volume%point_to_data("ZINDXA_LL")
        z_lr => volume%point_to_data("ZINDXA_LR")
        z_ur => volume%point_to_data("ZINDXA_UR")
        z_ul => volume%point_to_data("ZINDXA_UL")
    case ("B")
        z => volume%point_to_data("ZINDXB")
        z_ll => volume%point_to_data("ZINDXB_LL")
        z_lr => volume%point_to_data("ZINDXB_LR")
        z_ur => volume%point_to_data("ZINDXB_UR")
        z_ul => volume%point_to_data("ZINDXB_UL")
    end select

    !
    !  Given x,y,z as arrays in [gate,ray] space,
    !     read needed fields from an open NetCDF file (self%ncid),
    !     and compute RHO_D and TEMPERATURE fields in [gate,ray] space.
    !

    ! write(*,'(" WRF   build RHO_D, TEMPERATURE")')

    xoffs = 0.0
    yoffs = 0.0
    zoffs = 0.0

    ierr = nf90_inq_varid ( self%ncid, "QVAPOR", q_varid )
    call error_handler(ierr, "Problem inquire QVAPOR" )

    ierr = nf90_inq_varid ( self%ncid, "T", t_varid )
    call error_handler(ierr, "Problem inquire T" )

    allocate ( RHO_D_grid ( self%ni , self%nj , self%nk ) )
    allocate ( tgrid      ( self%ni , self%nj , self%nk ) )
    allocate ( qv         ( self%ni , self%nj , self%nk ) )
    allocate ( RHO_DS_grid ( self%ni , self%nj , self%nk ) )

    ierr = nf90_get_var ( self%ncid, t_varid, tgrid, start=(/1,1,1,self%time_frame/) )
    call error_handler(ierr, "Problem get T")

    ierr = nf90_get_var ( self%ncid, q_varid, qv, start=(/1,1,1,self%time_frame/) )
    call error_handler(ierr, "Problem get QVAPOR")

    tgrid = (tgrid + 300.) * (self%p/1.E5) ** gamma
    RHO_D_grid=(self%p*qv)/(eps+qv)      !  e water vapor pressure in Pa
    RHO_D_grid=(self%p-RHO_D_grid)/(Rd*tgrid)  !  kg/m^3

    ! Added by BWK, 3/14/2022
    do kidx = 1, self%nk
       if ( kidx == 1 ) then
         RHO_DS_grid(:,:,kidx)=(self%p(:,:,kidx)*qv(:,:,kidx))/(eps+qv(:,:,kidx))      !  e water vapor pressure in Pa
         RHO_DS_grid(:,:,kidx)=(self%p(:,:,kidx)-RHO_DS_grid(:,:,kidx))/(Rd*tgrid(:,:,kidx))  !  kg/m^3
       else
         RHO_DS_grid(:,:,kidx) = RHO_DS_grid(:,:,1)
       endif
    enddo

    out_RHO_d  = -9999.
    out_T      = -9999.
    out_RHO_ds = -9999.

    ! From CR-SIM PostProcessing - Added by B. Klotz (12/14/2023)

    dr = scan%meters_between_gates
    wfac=-2.e0*log(2.e0)
    wfacr=pi2/wfac

    do iray = 1, size(x,2)
        do igate = 1, size(x,1)
            if ( z(igate,iray) < -9998 ) cycle

            tmp_gate_array_rho_d = 0.0
            tmp_gate_array_t = 0.0
            tmp_gate_array_rho_ds = 0.0
            wgt_use1 = 0.0
            wgt_tot1 = 0.0

            tmp_gate_array_cr_rho_d = 0.0
            tmp_gate_array_cr_t = 0.0
            tmp_gate_array_cr_rho_ds = 0.0
            wgt_use_cr = 0.0
            wgt_tot_cr = 0.0

            ! Get the minimum and maximum indices in each direction
            maxx = MAX(x_ul(igate,iray),x_ur(igate,iray),x_lr(igate,iray),x_ll(igate,iray))
            maxy = MAX(y_ul(igate,iray),y_ur(igate,iray),y_lr(igate,iray),y_ll(igate,iray))
            maxz = MAX(z_ul(igate,iray),z_ur(igate,iray),z_lr(igate,iray),z_ll(igate,iray))

            minx = MIN(x_ul(igate,iray),x_ur(igate,iray),x_lr(igate,iray),x_ll(igate,iray))
            miny = MIN(y_ul(igate,iray),y_ur(igate,iray),y_lr(igate,iray),y_ll(igate,iray))
            minz = MIN(z_ul(igate,iray),z_ur(igate,iray),z_lr(igate,iray),z_ll(igate,iray))

            if (minz < 0) minz = 1
            if (maxz < 0) maxz = 1

            ! Distance in grid points between max and min positions
            dx = (int ( maxx + xoffs ) + 1) - (int ( minx + xoffs ) ) + 3
            dy = (int ( maxy + yoffs ) + 1) - (int ( miny + yoffs ) ) + 3
            dz = (int ( maxz + zoffs ) + 1) - (int ( minz + zoffs ) ) + 3

            ! Now let's determine the weight for each associated grid point
            ! Loops over the number of total indices (dx * dy * dz)
            !total_inds = dx * dy * dz

            im = int ( minx + xoffs )
            jm = int ( miny + yoffs )
            km = int ( minz + zoffs )

            DX_LOOP : do dxl = 1, dx+1

                if (dxl .eq. 1) then
                    !im = int (x(igate,iray) + xoffs)
                    im = int (minx + xoffs)
                endif

                if (im < 1 .or. im > self%ni) then
                    im=im+1
                    cycle
                endif

                xb = ABS( (x(igate,iray) + xoffs) - im) / (dx)
                xa = 1.0 - xb

                tmp_dx_meters  = (im - volume%aircraft_xgrid(iray)) * grid_dx ! Added by B. Klotz (12/14/2023)

                DY_LOOP : do dyl = 1, dy+1

                    if (dyl .eq. 1) then
                        !jm = int (y(igate,iray) + yoffs)
                        jm = int (miny + yoffs)
                    endif

                    if (jm < 1 .or. jm > self%nj) then
                        jm=jm+1
                        cycle
                    endif

                    yb = ABS( (y(igate,iray) + yoffs) - jm) / (dy)
                    ya = 1.0 - yb

                    tmp_dy_meters  = (jm - volume%aircraft_ygrid(iray)) * grid_dx ! Added by B. Klotz (12/14/2023)

                    DZ_LOOP : do dzl = 1, dz+1

                        if (dzl .eq. 1) then
                            km = int (z(igate,iray) + zoffs)
                        endif
         

                        if (km < 1 .or. km > self%nk) then
                            km=km+1
                            cycle
                        endif

                        tmp_dz_meters = (self%zf(im,jm,km)) - volume%aircraft_z(iray) ! Added by B. Klotz (12/14/2023)
                        ! print *, "im, jm, km, dx_m, dy_m, dz_m",im,jm,km,tmp_dx_meters,tmp_dy_meters,tmp_dz_meters

                        call WeightFuncCRS(tmp_dx_meters,tmp_dy_meters,tmp_dz_meters,volume%range(igate),volume%azimuth(iray), &
                             & volume%elevation(iray),volume%beamwidth_h(iray),volume%beamwidth_v(iray),dr,wfacr,wfac,wgt_use_cr)
        
                        !wgt_use_cr = Wr*Wa*We
                        if (wgt_use_cr<=1.e-5) wgt_use_cr=0.e0

                        zb = ABS( (z(igate,iray) + zoffs) - km) / (dz)
                        za = 1.0 - zb

                        ! Get the aircraft relative x,y,z position

                        wgt_use1 = xa*ya*za
                        wgt_tot1 = wgt_tot1+wgt_use1

                        wgt_tot_cr = wgt_tot_cr+wgt_use_cr

                        tmp_gate_array_rho_d = tmp_gate_array_rho_d + (wgt_use1*RHO_D_grid(im,jm,km))
                        tmp_gate_array_cr_rho_d = tmp_gate_array_cr_rho_d + (wgt_use_cr*RHO_D_grid(im,jm,km))
                        tmp_gate_array_rho_ds = tmp_gate_array_rho_ds + (wgt_use1*RHO_DS_grid(im,jm,km))
                        tmp_gate_array_cr_rho_ds = tmp_gate_array_cr_rho_ds + (wgt_use_cr*RHO_DS_grid(im,jm,km))
                        tmp_gate_array_t = tmp_gate_array_t + (wgt_use1*tgrid(im,jm,km))
                        tmp_gate_array_cr_t = tmp_gate_array_cr_t + (wgt_use_cr*tgrid(im,jm,km))
                        !print *, "im, jm, km, xa, ya, za", im, jm, km, xa, ya, za, wgt_tot
                        !print *, "Wr, We, Wa, Wtot, Wfacr", Wr, We, Wa, wgt_tot, wfacr

                        km = km+1
                    enddo DZ_LOOP

                    jm = jm+1
                enddo DY_LOOP

                im = im+1
            enddo DX_LOOP

            ! This is a similar procedure to the CR-SIM post-processing code
            if (wgt_tot_cr > 0) then
                out_RHO_d(igate,iray)  = tmp_gate_array_cr_rho_d/wgt_tot_cr
                out_RHO_ds(igate,iray) = tmp_gate_array_cr_rho_ds/wgt_tot_cr
                out_T(igate,iray)      = tmp_gate_array_cr_t/wgt_tot_cr
            elseif (wgt_tot_cr <=0 .and. wgt_tot1 > 0) then
                out_RHO_d(igate,iray)  = tmp_gate_array_rho_d/wgt_tot1
                out_RHO_ds(igate,iray) = tmp_gate_array_rho_ds/wgt_tot1
                out_T(igate,iray)      = tmp_gate_array_t/wgt_tot1
            else
                out_RHO_d(igate,iray)  = -9999.
                out_RHO_ds(igate,iray) = -9999.
                out_T(igate,iray)      = -9999.
            endif

            ! End addition of loop filter

        enddo
    enddo
  end subroutine wrf_RHO_d_varbw

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
! End variable beamwidth versions
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine wrf_value ( self , x , y , z , u , v , w )
    implicit none
    class(wrf_metadata_type), intent(in) :: self
    real(kind=RKIND), intent(in) :: x, y, z
    real(kind=RKIND), intent(out) :: u, v, w
    
    integer          :: i,j,k
    integer          :: im, jm, km
    integer          :: ip, jp, kp
    real(kind=RKIND) :: xa, ya, za, xb, yb, zb
    real(kind=RKIND) :: xu, yv

    u = -9999.
    v = -9999.
    w = -9999.
    if ( z < 1 ) return

    km = int(z)
    kp = km+1
    zb = ( z - int(z) )
    za = 1.0-zb


    !  Get bracketing i indices for the u field (staggered u in x gives us the 0.5 offset).
    xu = x+0.5
    im = int(xu)
    ip = im+1
    if (ip > self%ni+1) return
    xb = ( xu - int(xu) )
    xa = 1.0-xb

    !  Get bracketing j indices for the u field (no stagger in y for u)
    jm = int(y)
    jp = jm+1
    if (jp > self%nj) return
    yb = ( y - int(y) )
    ya = 1.0-yb

    u =  xa*ya*za*self%u(im,jm,km) + xa*ya*zb*self%u(im,jm,kp) + &
         xa*yb*za*self%u(im,jp,km) + xa*yb*zb*self%u(im,jp,kp) + &
         xb*ya*za*self%u(ip,jm,km) + xb*ya*zb*self%u(ip,jm,kp) + &
         xb*yb*za*self%u(ip,jp,km) + xb*yb*zb*self%u(ip,jp,kp)

    !  Get bracketing i indices for the v field (no stagger in x for v)
    im = int(x)
    ip = im+1
    if (ip > self%ni) then
       u = -9999.
       return
    endif
    xb = ( x - int(x) )
    xa = 1.0-xb

    !  Get bracketing j indices for the v field (staggered v in y gives us the 0.5 offset).
    yv = y+0.5
    jm = int(yv)
    jp = jm+1
    if (jp > self%nj+1) then
       u = -9999.
       return
    endif
    yb = ( yv - int(yv) )
    ya = 1.0-yb

    v =  xa*ya*za*self%v(im,jm,km) + xa*ya*zb*self%v(im,jm,kp) + &
         xa*yb*za*self%v(im,jp,km) + xa*yb*zb*self%v(im,jp,kp) + &
         xb*ya*za*self%v(ip,jm,km) + xb*ya*zb*self%v(ip,jm,kp) + &
         xb*yb*za*self%v(ip,jp,km) + xb*yb*zb*self%v(ip,jp,kp)

    !  No stagger in x or y for w
    im = int(x)
    ip = im+1
    xb = ( x - int(x) )
    xa = 1.0-xb

    jm = int(y)
    jp = jm+1
    yb = ( y - int(y) )
    ya = 1.0-yb

    w =  xa*ya*za*self%w(im,jm,km) + xa*ya*zb*self%w(im,jm,kp) + &
         xa*yb*za*self%w(im,jp,km) + xa*yb*zb*self%w(im,jp,kp) + &
         xb*ya*za*self%w(ip,jm,km) + xb*ya*zb*self%w(ip,jm,kp) + &
         xb*yb*za*self%w(ip,jp,km) + xb*yb*zb*self%w(ip,jp,kp)

  end subroutine wrf_value

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine error_handler ( ierr , message )
#ifdef _PARALLEL_
    use mpi, only: MPI_ABORT
    use mpi, only: MPI_COMM_WORLD
#endif
    use netcdf, only : NF90_NOERR
    use netcdf, only : nf90_strerror
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: message
    integer :: xerr
    integer :: ierrmpi
    if ( ierr == NF90_NOERR ) return
    write(*,'("NETCDF >> ", A)') nf90_strerror(ierr)
    write(*,'("MODULE_ACCESS_WRF: ",A)') message
#ifdef _PARALLEL_
    call MPI_ABORT(MPI_COMM_WORLD, xerr, ierrmpi)
#endif
    stop
  end subroutine error_handler

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

end module module_access_wrf
