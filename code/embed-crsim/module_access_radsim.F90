module module_access_radsim
  implicit none
  private

  integer, parameter :: MAX_RADAR_FILES = 1

  type, public :: radar_simulator_type
      integer :: ncid = -999999
      integer :: nx = -999999
      integer :: ny = -999999
      integer :: nz = -999999
      integer :: time_index = -999999
      real :: x_index = -999999
      real :: y_index = -999999
      real :: z_index = -999999
      real :: height = -999999.
      real :: beamwidth = -999999.
      real :: range_resolution = -999999.
      integer :: nfields = -999999
      character(len=256), dimension(1000) :: list_of_fields
    contains
      procedure :: open => radar_open
      procedure :: find_info => radar_find_info
      procedure :: close => radar_close
      procedure :: interp => radar_interp
  end type radar_simulator_type

contains

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine radar_open ( self , flnm )
    use netcdf, only : NF90_NOWRITE
    use netcdf, only : NF90_MAX_NAME
    use netcdf, only : NF90_MAX_VAR_DIMS
    use netcdf, only : NF90_FLOAT
    use netcdf, only : nf90_open
    use netcdf, only : nf90_inq_dimid
    use netcdf, only : nf90_inquire
    use netcdf, only : nf90_inquire_dimension
    use netcdf, only : nf90_inquire_variable
    implicit none
    class(radar_simulator_type), intent(inout) :: self
    character(len=*), intent(in) :: flnm
    integer :: ierr
    integer :: dimid_nx, dimid_ny, dimid_nz
    integer :: nVariables, nDims, xType
    integer :: varid, vindx
    character(len=NF90_MAX_NAME) :: varname
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimid

    !
    !  Open file
    !
    ierr = nf90_open(trim(flnm), NF90_NOWRITE, self%ncid)
    call error_handler ( ierr , "Problem opening file "//trim(flnm) )

    !
    !  Get dimensions.
    !

    ierr = nf90_inq_dimid ( self%ncid, "nx", dimid_nx )
    call error_handler ( ierr, "Problem inq_dimid for 'nx'" )
    ierr = nf90_inquire_dimension ( self%ncid, dimid_nx, len=self%nx )
    call error_handler ( ierr, "Problem inquiring dimension nx." )

    ierr = nf90_inq_dimid ( self%ncid, "ny", dimid_ny )
    call error_handler ( ierr, "Problem inq_dimid for 'ny'" )
    ierr = nf90_inquire_dimension ( self%ncid, dimid_ny, len=self%ny )
    call error_handler ( ierr, "Problem inquiring dimension ny." )

    ierr = nf90_inq_dimid ( self%ncid, "nz", dimid_nz )
    call error_handler ( ierr, "Problem inq_dimid for 'nz'" )
    ierr = nf90_inquire_dimension ( self%ncid, dimid_nz, len=self%nz )
    call error_handler ( ierr, "Problem inquiring dimension nz." )

    !
    !  Get list of 3-d fields
    !

    ierr = nf90_inquire ( self%ncid, nVariables=nVariables )
    call error_handler ( ierr, "Problem inquiring nVariables" )

    self%nfields = 0
    do varid = 1, nVariables
        ierr = nf90_inquire_variable ( self%ncid, varid, xtype=xType, name=varname, ndims=nDims, dimids=dimid )
        call error_handler ( ierr, "Problem inquire variable" )
        if ( xType /= NF90_FLOAT ) cycle  ! Only process variables of type NF90_FLOAT
        if ( ndims /= 3          ) cycle  ! Only process 3d variables
        ! write(*,'(A20, 10I6)') trim(varname), dimid(1:nDims)
        self%nfields = self%nfields + 1
        self%list_of_fields(self%nfields) = trim(varname)
    enddo

  end subroutine radar_open
    
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine radar_find_info(self, radar_glob_pattern)
    use netcdf, only : NF90_NOWRITE
    use netcdf, only : NF90_GLOBAL
    use netcdf, only : nf90_open
    use netcdf, only : nf90_inq_varid
    use netcdf, only : nf90_get_var
    use netcdf, only : nf90_close
    use netcdf, only : nf90_get_att
    implicit none
    class(radar_simulator_type), intent(inout) :: self
    character(len=*), intent(in) :: radar_glob_pattern
    character(len=1024) :: radar_file, string, str
    integer :: i
    integer :: ncid
    integer :: ierr
    integer :: varid

    ! print*, 'radar_glob_pattern = ', trim(radar_glob_pattern)
    call execute_command_line("ls "//trim(radar_glob_pattern)//" > .rthing")

    open(10,file=".rthing", status='old', form='formatted', action='read')

    do i = 1, MAX_RADAR_FILES
        read(10,'(A)',end=999) string
        print*, trim(string)
        radar_file = trim(string)

        ierr = nf90_open(trim(string), NF90_NOWRITE, ncid)
        call error_handler ( ierr , "Problem opening file "//trim(string) )

        ! ierr = nf90_inq_varid ( ncid, "XTIME", varid )
        ! call error_handler ( ierr, "Problem inquiring varid XTIME" )

        ierr = nf90_get_att ( ncid, NF90_GLOBAL, "scene_extracted_at_time_step", str)
        call error_handler ( ierr, "Problem get att scene_extracted_at_time_step" )
        read(str,*) self%time_index
        write(*,'("radar time index = ", I10)') self%time_index

        ierr = nf90_inq_varid ( ncid, "rad_ixc", varid )
        call error_handler ( ierr, "Problem inquiring varid rad_ixc" )

        ierr = nf90_get_var ( ncid, varid, self%x_index )
        call error_handler ( ierr, "Problem get var rad_ixc" )
        write(*,'("radar position x = ", F10.4)') self%x_index

        ierr = nf90_inq_varid ( ncid, "rad_iyc", varid )
        call error_handler ( ierr, "Problem inquiring varid rad_iyc" )

        ierr = nf90_get_var ( ncid, varid, self%y_index )
        call error_handler ( ierr, "Problem get var rad_iyc" )
        write(*,'("radar position y = ", F10.4)') self%y_index

        ierr = nf90_inq_varid ( ncid, "rad_zc", varid )
        call error_handler ( ierr, "Problem inquiring varid rad_zc" )

        ierr = nf90_get_var ( ncid, varid, self%height )
        call error_handler ( ierr, "Problem get var rad_zc" )
        write(*,'("radar height = ", F10.4)') self%height

        ierr = nf90_inq_varid ( ncid, "rad_beamwidth", varid )
        call error_handler ( ierr, "Problem inquiring varid rad_beamwidth" )

        ierr = nf90_get_var ( ncid, varid, self%beamwidth )
        call error_handler ( ierr, "Problem get var rad_beamwidth" )
        write(*,'("radar beamwidth = ", F10.4)') self%beamwidth

        ierr = nf90_inq_varid ( ncid, "rad_range_resolution", varid )
        call error_handler ( ierr, "Problem inquiring varid rad_range_resolution" )

        ierr = nf90_get_var ( ncid, varid, self%range_resolution )
        call error_handler ( ierr, "Problem get var rad_range_resolution" )
        write(*,'("radar range resolution = ", F10.4)') self%range_resolution

        ierr = nf90_close(ncid)
        call error_handler ( ierr, "Problem close")

    enddo
999 continue
    close(10, status='delete')

  end subroutine radar_find_info

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine radar_interp ( self, x , y , z , varname, out_array )
    use netcdf, only : nf90_inq_varid
    use netcdf, only : nf90_get_var
    use netcdf, only : nf90_get_att
    implicit none
    class(radar_simulator_type), intent(in) :: self
    character(len=*), intent(in) :: varname
    real, dimension(:,:) :: x, y, z
    real, dimension(:,:) :: out_array
    real, allocatable, dimension(:,:,:) :: source_array

    integer :: iray, igate, im, ip, jm, jp, km, kp
    real    :: xa, xb, ya, yb, za, zb
    integer :: ierr, varid
    real :: xoffs, yoffs, zoffs
    character(len=8) :: stagger

    ierr = nf90_inq_varid ( self%ncid, trim(varname), varid )
    call error_handler(ierr, "Problem inquire " // trim(varname) )
    xoffs = 0.0
    yoffs = 0.0
    zoffs = 0.0

    allocate ( source_array ( self%nx, self%ny, self%nz ) )

    ierr = nf90_get_var ( self%ncid, varid, source_array )
    call error_handler(ierr, "Problem get " // trim(varname) )

    out_array = -9999.

    do iray = 1, size(x,2)
        do igate = 1, size(x,1)
            if ( z(igate,iray) < -9998 ) cycle

            im = int ( x(igate,iray) + xoffs )
            ip = im+1
            if ( ip > size(x,1) ) cycle

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
            
            !
            !  If any one of the eight source values looks like it's a bad value, 
            !  don't interpolate.
            !
            
            if ( minval(source_array(im:ip,jm:jp,km:kp)) < -998 ) cycle

            !
            !  All seems well.... Interpolate.
            !
            out_array(igate,iray) = xa*ya*za*source_array(im,jm,km) + xa*ya*zb*source_array(im,jm,kp) + &
                 &                  xa*yb*za*source_array(im,jp,km) + xa*yb*zb*source_array(im,jp,kp) + &
                 &                  xb*ya*za*source_array(ip,jm,km) + xb*ya*zb*source_array(ip,jm,kp) + &
                 &                  xb*yb*za*source_array(ip,jp,km) + xb*yb*zb*source_array(ip,jp,kp)

        enddo
    enddo
  end subroutine radar_interp

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine radar_close(self)
    use netcdf, only : nf90_close
    implicit none
    class(radar_simulator_type), intent(inout) :: self
    integer :: ierr

    ierr = nf90_close(self%ncid)
    call error_handler ( ierr , "Problem closing file." )

    self%ncid = -999999

  end subroutine radar_close

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine error_handler ( ierr , message )
    use netcdf, only : NF90_NOERR
    use netcdf, only : nf90_strerror
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: message
    if ( ierr == NF90_NOERR ) return
    write(*,'("NETCDF >> ", A)') nf90_strerror(ierr)
    write(*,'("MODULE_ACCESS_RADSIM: ",A)') message
    stop
  end subroutine error_handler

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
end module module_access_radsim
