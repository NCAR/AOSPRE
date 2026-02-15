module module_scanning
  use module_configuration, only : RKIND
  implicit none

  private

  type, public :: scan_type

      character(len=24) :: primary_axis
      character(len=24), dimension(2) :: table_provides
      character(len=24) :: sweep_mode
      real(kind=RKIND) :: meters_between_gates
      real(kind=RKIND) :: meters_to_center_of_first_gate
      real(kind=RKIND) :: max_range_in_meters
      real(kind=RKIND) :: snr_mask_threshold

      real(kind=RKIND) :: seconds_for_scan_cycle
      real(kind=RKIND) :: seconds_plus_skip
      real(kind=RKIND) :: pulse_repetition_frequency
      integer :: pulses_per_pulse_set
      integer :: revisits_per_acquisition_time
      integer :: beams_per_acquisition_time

      real(kind=RKIND) :: pulse_repetition_time
      real(kind=RKIND) :: revisit_time
      real(kind=RKIND) :: data_acquisition_time
      integer          :: pulses_per_revisit_time
      integer          :: pulses_per_acquisition_time

      real(kind=RKIND) :: lambda
      real(kind=RKIND) :: Nyquist_Velocity

      integer          :: scan_cycle
      real(kind=RKIND) :: skip_seconds_between_scans

      real(kind=RKIND) :: fold_limit_lower
      real(kind=RKIND) :: fold_limit_upper
      integer :: sweep_count
      integer :: beam_count
      real(kind=RKIND), allocatable, dimension(:) :: beam_time_from_start_of_scan_cycle
      integer, allocatable, dimension(:) :: sweep_start_index
      integer, allocatable, dimension(:) :: sweep_end_index
      real(kind=RKIND), allocatable, dimension(:) :: rotation
      real(kind=RKIND), allocatable, dimension(:) :: tilt
    contains
      procedure :: read_scanning_table_file
  end type scan_type

contains

  subroutine read_scanning_table_file ( self, filename, herky_jerky )
    implicit none
    class ( scan_type )             , intent(inout) :: self
    character(len=*)                , intent(in)    :: filename
    logical                         , intent(in)    :: herky_jerky
    integer :: ierr
    character(len=1024) :: string
    character(len=24) :: lhs,rhs
    integer :: k
    integer :: indx
    real(kind=RKIND), allocatable, dimension(:) :: temporary_rotation
    real(kind=RKIND), allocatable, dimension(:) :: temporary_tilt
    real(kind=RKIND), allocatable, dimension(:) :: temporary_start
    real(kind=RKIND), allocatable, dimension(:) :: temporary_end
    real(kind=RKIND) :: rdrota
    real(kind=RKIND) :: rdtilt
    integer :: line_number
    real(kind=RKIND) :: offset1
    real(kind=RKIND) :: offset2
    integer :: ibi

    allocate(temporary_rotation(1000000))
    allocate(temporary_tilt    (1000000))
    allocate(temporary_start   (10000))
    allocate(temporary_end     (10000))
    temporary_rotation = -1.E36
    temporary_tilt  = -1.E36
    temporary_start = -1.E36
    temporary_end   = -1.E36

    open(12,file=trim(filename), status='old', form='formatted', action='read', iostat=ierr)
    if ( ierr /= 0 ) then
        write(*,'("Problem opening scanning table file ''", A, "''")') trim(filename)
        stop "MODULE_SCANNING : READ_SCANNING_TABLE_FILE"
    endif

    self%beam_count = 0
    self%sweep_count = 0
    self%primary_axis = "NULL"
    self%table_provides = "NULL"
    self%sweep_mode = "NULL"

    line_number = 0
    READLOOP_K : do k = 1, 1000000

        read(12,'(A)',end=200,err=300) string
        line_number = line_number + 1

        ! Remove any leading white space
        string = adjustl(string)

        ! Skip blank lines
        if ( string == " " ) cycle READLOOP_K

        ! Skip lines that begin with "#"
        if ( string(1:1) == "#") cycle READLOOP_K

        ! Remove anything trailing a "#"
        indx = index(string,"#")
        if (indx > 0) then
            string = trim(string(1:indx-1))
            if ( string == " " ) cycle READLOOP_K
        endif

        ! Split lines that have "=".  We recognize "PRIMARY_AXIS", "PARAMETERS", and "SWEEP_MODE"
        indx = index(string,"=")
        if (indx>0) then
            lhs = trim(adjustl(string(1:indx-1)))
            rhs = trim(adjustl(string(indx+1:)))
            select case (lhs)
            case default
                write(*,'("string = ''",A,"''")') trim(string)
                write(*,'("line number = ",I10)') line_number
                stop "MODULE_SCANNING:READ_SCANNING_TABLE_FILE:  Failure to interpret string."
            case ( "PRIMARY_AXIS", "Primary_Axis", "Primary_axis", "primary_axis", &
                 & "PRIMARY-AXIS", "Primary-Axis", "Primary-axis", "primary-axis", &
                 & "PRIMARY AXIS", "Primary Axis", "Primary axis", "primary axis" )
                self%primary_axis = rhs
                cycle READLOOP_K
            case ("SWEEP_MODE", "Sweep_Mode","Sweep_mode","sweep_mode")
                select case (rhs)
                case default
                    write(*,'("sweep_mode = ",A)') trim(rhs)
                    stop "Unrecognized sweep_mode.  Check scanning table setup, (or maybe a coding TODO)."
                case ("rhi","RHI")
                    self%sweep_mode = "rhi"
                case ("sector","Sector")
                    self%sweep_mode = "sector"
                case ("ppi","PPI")
                    self%sweep_mode = "ppi"
                case ("azimuth_surveillance")
                    self%sweep_mode = "azimuth_surveillance"
                case ("elevation_surveillance")
                    self%sweep_mode = "elevation_surveillance"
                end select
                cycle READLOOP_K
            case ("PARAMETERS", "Parameters","parameters")
                ! Recognize "ROT", "TILT"
                indx = index(rhs,",")
                self%table_provides(1) = trim(adjustl(rhs(1:indx-1)))
                self%table_provides(2) = trim(adjustl(rhs(indx+1:)))

                select case (self%table_provides(1))
                case default
                    stop "Not ROT"
                case("ROT","Rot","rot","ROTATION","Rotation","rotation")
                    continue
                end select

                select case (self%table_provides(2))
                case default
                    stop "Not TILT"
                case("TILT","Tilt","tilt")
                    continue
                end select

                cycle READLOOP_K
            end select
        endif

        ! Everything else should be the set of angles that build up a sweep.
        if (string == "<sweep>") then

            self%sweep_count = self%sweep_count + 1

            if ( self%sweep_count > size(temporary_start) ) then
                write(*,'("Too many sweeps.  Increase size of work arrays ''temporary_start'' and ''temporary_end''")')
                stop "MODULE_SCANNING : READ_SCANNING_TABLE_FILE"
            endif

            temporary_start(self%sweep_count) = self%beam_count
            cycle READLOOP_K

        endif

        if (string == "</sweep>") then
            temporary_end(self%sweep_count) = self%beam_count-1
            cycle READLOOP_K
        endif



        read(string,*,iostat=ierr) rdrota, rdtilt
        if ( ierr /= 0 ) then
            write(*,'("Problem reading scanning table file ''", A, "'', line ", I8)') trim(filename), k
            stop "MODULE_SCANNING : READ_SCANNING_TABLE_FILE"
        endif
        line_number = line_number + 1

        self%beam_count = self%beam_count + 1
        if ( self%beam_count > size(temporary_rotation) ) then
            write(*,'("Too many entries.  Increase size of work arrays ''temporary_rotation'' and ''temporary_tilt''")')
            stop "MODULE_SCANNING : READ_SCANNING_TABLE_FILE"
        endif
        temporary_rotation(self%beam_count) = rdrota
        temporary_tilt    (self%beam_count) = rdtilt

    enddo READLOOP_K

200 continue

    close(12)

    write(*,'("module_scanning: Sensor PRIMARY_AXIS: ''",A,"''")') trim(self%primary_axis)
    write(*,'("module_scanning: sweep_mode: ''",A,"''")') trim(self%sweep_mode)
    write(*,'("module_scanning: Table Provides: ''",A,"'' , ''",A,"''")') trim(self%table_provides(1)), trim(self%table_provides(2))
    
    if (self%sweep_mode == "NULL") stop "sweep_mode not found.  Check scanning table setup."
    if (self%primary_axis == "NULL") stop "sensor PRIMARY_AXIS not found.  Check scanning table setup."
    if (self%table_provides(1) == "NULL") stop "sweep PARAMETERS(1) not found.  Check scanning table setup."
    if (self%table_provides(2) == "NULL") stop "sweep PARAMETERS(2) not found.  Check scanning table setup."

    allocate(self%rotation(self%beam_count))
    allocate(self%tilt(self%beam_count))
    allocate(self%sweep_start_index(self%sweep_count))
    allocate(self%sweep_end_index  (self%sweep_count))
    self%rotation = temporary_rotation(1:self%beam_count)
    self%tilt     = temporary_tilt    (1:self%beam_count)
    self%sweep_start_index = temporary_start(1:self%sweep_count)
    self%sweep_end_index = temporary_end(1:self%sweep_count)

    !  Now that we've read the scanning table, we can compute more timing details
    !  of individual beams.

    allocate(self%beam_time_from_start_of_scan_cycle(self%beam_count))

    if (herky_jerky) then

        ! For intermittent scans
        ! Assume instantaneous scans, so the only time increment is going to be the time to skip between scans.
        self%seconds_for_scan_cycle = 0
        self%seconds_plus_skip = self%skip_seconds_between_scans
        self%beam_time_from_start_of_scan_cycle = 0

    else

        !  Calculate seconds_for_scan_cycle, padding to a whole number of dwells
        if ( mod(self%beam_count, self%beams_per_acquisition_time)>0) then
            ! Integer division; add 1 for a partial set of beams in acquisition time.
            self%seconds_for_scan_cycle = (1 + (self%beam_count / self%beams_per_acquisition_time)) * self%data_acquisition_time
        else
            ! Integer division
            self%seconds_for_scan_cycle = (self%beam_count / self%beams_per_acquisition_time) * self%data_acquisition_time
        endif

        do k = 1, self%beam_count

            !  Offset1 is the time of the start of a data acquisition cycle.
            !  Note the integer division to put each beam into the right data acquisition cycle.

            offset1 = ((k-1)/self%beams_per_acquisition_time) * self%data_acquisition_time

            !  offset2 is the time of the beam from the start of a revisit cycle.
            !  ibi is the index of the beam (1,2,3,etc) in a revisit cycle (or data acquisition cycle, for that matter)
            ibi = 1+mod(k-1,self%beams_per_acquisition_time)
            offset2 = ( ( (self%revisits_per_acquisition_time-1) * self%pulses_per_revisit_time ) + &
                 &      ( ibi * self%pulses_per_pulse_set ) ) * self%pulse_repetition_time

            !  Time of the beam from the start of a scan cycle.
            self%beam_time_from_start_of_scan_cycle(k) = offset1  + offset2

        enddo

        self%seconds_plus_skip = self%seconds_for_scan_cycle + self%skip_seconds_between_scans

    endif

    write(*,'(" beam_count                               = ", I12)') self%beam_count
    write(*,'(" sweep_count                              = ", I12)') self%sweep_count
    write(*,'(" seconds_for_scan_cycle                   = ", F12.6)') self%seconds_for_scan_cycle
    write(*,'(" skip_seconds_between_scan                = ", F12.6)') self%skip_seconds_between_scans
    write(*,'(" seconds_plus_skip                        = ", F12.6)') self%seconds_plus_skip

    do k = 1, self%sweep_count
        write(*,'(6x,"Sweep ",I6, " :: beams ", I6, " through ", I6, " :: timing ", F12.6, " through ", F12.6, " seconds")') &
             k, self%sweep_start_index(k), self%sweep_end_index(k), &
             self%beam_time_from_start_of_scan_cycle(self%sweep_start_index(k)+1), &
             self%beam_time_from_start_of_scan_cycle(self%sweep_end_index(k)+1)
    enddo
    write(*,*)

    return

300 continue
    write(*,'("Problem reading scanning table file ''", A, "''")') trim(filename)
    stop "MODULE_SCANNING : READ_SCANNING_TABLE_FILE"

  end subroutine read_scanning_table_file

end module module_scanning
