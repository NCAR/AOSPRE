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
      real(kind=RKIND) :: MDS_1km

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

    type, public :: scan_type_idt

      character(len=24) :: primary_axis
      character(len=24), dimension(2) :: table_provides
      character(len=24) :: sweep_mode
      real(kind=RKIND) :: meters_between_gates
      real(kind=RKIND) :: meters_to_center_of_first_gate
      real(kind=RKIND) :: max_range_in_meters
      real(kind=RKIND) :: snr_mask_threshold
      real(kind=RKIND) :: MDS_1km

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
      integer          :: nlobes              ! This is the number of lobes = main + up to 4 positional sidelobes (BK, 5/6/26)
      real(kind=RKIND) :: skip_seconds_between_scans

      real(kind=RKIND) :: fold_limit_lower
      real(kind=RKIND) :: fold_limit_upper
      integer :: sweep_count
      integer :: beam_count
      real(kind=RKIND), allocatable, dimension(:) :: beam_time_from_start_of_scan_cycle
      integer, allocatable, dimension(:) :: sweep_start_index
      integer, allocatable, dimension(:) :: sweep_end_index
      real(kind=RKIND), allocatable, dimension(:,:) :: rotation  ! rot_main, rot_sl1, rot_sl2, rot_sl3, rot_sl4
      real(kind=RKIND), allocatable, dimension(:,:) :: tilt      ! tilt_main, tilt_sl1, tilt_sl2, tilt_sl3, tilt_sl4
      real(kind=RKIND), allocatable, dimension(:,:) :: bwth      ! bw_main, bw_sl1, bw_sl2, bw_sl3, bw_sl4
      real(kind=RKIND), allocatable, dimension(:,:) :: direct    ! dir_main, dir_sl1, dir_sl2, dir_sl3, dir_sl4
      real(kind=RKIND), allocatable, dimension(:,:) :: azim      ! will be updated in the 'beam_geometry' subroutine; This provides the main+SL structure
      real(kind=RKIND), allocatable, dimension(:,:) :: elev      ! will be updated in the 'beam_geometry' subroutine
      real(kind=RKIND), allocatable, dimension(:,:) :: beamwidth_h      ! will be updated in the 'beam_geometry' subroutine; This provides the main+SL structure
      real(kind=RKIND), allocatable, dimension(:,:) :: beamwidth_v      ! will be updated in the 'beam_geometry' subroutine
      real(kind=RKIND), allocatable, dimension(:,:) :: sll       ! will be updated in the 'beam_geometry' subroutine; This is the side lobe level, relative to the main lobe directivity
      real(kind=RKIND), allocatable, dimension(:) :: tilt_err, rot_err, mds
    contains
      procedure :: read_scanning_table_file_idt
  end type scan_type_idt

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
                self%table_provides(2) = trim(adjustl(rhs(indx+1:indx+8)))

                select case (self%table_provides(1))
                case default
                    stop "Not ROT"
                case("ROT","Rot","rot","ROTATION","Rotation","rotation")
                    continue
                end select

                select case (self%table_provides(2))
                case default
                    stop "Not TILT"
                case("TILT","Tilt","tilt","TIL","Til","til")
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

  subroutine read_scanning_table_file_idt ( self, filename, herky_jerky )
    implicit none
    class ( scan_type_idt )             , intent(inout) :: self
    character(len=*)                , intent(in)    :: filename
    logical                         , intent(in)    :: herky_jerky
    integer :: ierr
    character(len=1024) :: string
    character(len=24) :: lhs,rhs, tmprhs
    integer :: k
    integer :: indx, indx2, bmnum, tmp_nlobe
    real(kind=RKIND), allocatable, dimension(:,:) :: temporary_rotation
    real(kind=RKIND), allocatable, dimension(:,:) :: temporary_tilt
    real(kind=RKIND), allocatable, dimension(:) :: temporary_roterr
    real(kind=RKIND), allocatable, dimension(:) :: temporary_tilterr
    real(kind=RKIND), allocatable, dimension(:,:) :: temporary_dir
    real(kind=RKIND), allocatable, dimension(:) :: temporary_mds
    real(kind=RKIND), allocatable, dimension(:,:) :: temporary_bw
    
    ! Testing to see if this makes sense (BWK, 04/21/26)
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_rotation_sl1
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_tilt_sl1
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_dir_sl1
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_bw_sl1
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_rotation_sl2
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_tilt_sl2
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_dir_sl2
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_bw_sl2
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_rotation_sl3
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_tilt_sl3
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_dir_sl3
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_bw_sl3
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_rotation_sl4
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_tilt_sl4
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_dir_sl4
    !real(kind=RKIND), allocatable, dimension(:) :: temporary_bw_sl4
    
    real(kind=RKIND), allocatable, dimension(:) :: temporary_start
    real(kind=RKIND), allocatable, dimension(:) :: temporary_end
    real(kind=RKIND), dimension(5) :: rdrot      ! rdrot_main, rdrot_sl1, rdrot_sl2, rdrot_sl3, rdrot_sl4
    real(kind=RKIND), dimension(5) :: rdtilt     ! rdtilt_main, rdtilt_sl1, rdtilt_sl2, rdtilt_sl3, rdtilt_sl4
    real(kind=RKIND), dimension(5) :: rddirect        ! dir_main, dir_sl1, dir_sl2, dir_sl3, dir_sl4
    real(kind=RKIND), dimension(5) :: rdbwth         ! bw_main, bw_sl1, bw_sl2, bw_sl3, bw_sl4
    real(kind=RKIND) :: rdrot_err, rdtilt_err, mds
    integer :: line_number, ioffset, offcnt, iloop
    integer, allocatable, dimension(:) :: ipos
    real(kind=RKIND) :: offset1
    real(kind=RKIND) :: offset2
    integer :: ibi

    allocate(ipos(1:3))
    allocate(temporary_rotation(1000000,5))
    allocate(temporary_tilt   (1000000,5))
    allocate(temporary_roterr(1000000))
    allocate(temporary_tilterr    (1000000))
    allocate(temporary_dir(1000000,5))
    allocate(temporary_mds(1000000))
    allocate(temporary_bw    (1000000,5))
    
    !allocate(temporary_rotation_sl1(1000000))
    !allocate(temporary_tilt_sl1    (1000000))
    !allocate(temporary_dir_sl1(1000000))
    !allocate(temporary_bw_sl1    (1000000))
    !allocate(temporary_rotation_sl2(1000000))
    !allocate(temporary_tilt_sl2    (1000000))
    !allocate(temporary_dir_sl2(1000000))
    !allocate(temporary_bw_sl2    (1000000))
    !allocate(temporary_rotation_sl3(1000000))
    !allocate(temporary_tilt_sl3    (1000000))
    !allocate(temporary_dir_sl3(1000000))
    !allocate(temporary_bw_sl3    (1000000))
    !allocate(temporary_rotation_sl4(1000000))
    !allocate(temporary_tilt_sl4    (1000000))
    !allocate(temporary_dir_sl4(1000000))
    !allocate(temporary_bw_sl4    (1000000))
    
    allocate(temporary_start   (10000))
    allocate(temporary_end     (10000))

    temporary_rotation = -1.E36
    temporary_tilt  = -1.E36
    temporary_roterr = -1.E36
    temporary_tilterr  = -1.E36
    temporary_dir = -1.E36
    temporary_mds  = -1.E36
    
    !temporary_rotation_sl1 = -1.E36
    !temporary_tilt_sl1  = -1.E36
    !temporary_dir_sl1 = -1.E36
    !temporary_bw_sl1  = -1.E36
    !temporary_rotation_sl2 = -1.E36
    !temporary_tilt_sl2  = -1.E36
    !temporary_dir_sl2 = -1.E36
    !temporary_bw_sl2  = -1.E36
    !temporary_rotation_sl3 = -1.E36
    !temporary_tilt_sl3  = -1.E36
    !temporary_dir_sl3 = -1.E36
    !temporary_bw_sl3  = -1.E36
    !temporary_rotation_sl4 = -1.E36
    !temporary_tilt_sl4  = -1.E36
    !temporary_dir_sl4 = -1.E36
    !temporary_bw_sl4  = -1.E36
    
    temporary_start = -1.E36
    temporary_end   = -1.E36
    tmp_nlobe = 0

    open(12,file=trim(filename), status='old', form='formatted', action='read', iostat=ierr)
    if ( ierr /= 0 ) then
        write(*,'("Problem opening scanning table file ''", A, "''")') trim(filename)
        stop "MODULE_SCANNING : READ_SCANNING_TABLE_FILE"
    endif

    self%beam_count = 0
    self%sweep_count = 0
    self%nlobes = 0
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
            print *, "RHS = ", rhs
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

                ! Need to loop over the string
                ioffset = 0
                offcnt = 0
                !print *, LBOUND(ipos,1)
                do iloop = 1, 2
                   !print *,"Loop count = ",iloop
                   !print *, rhs(ioffset+1:)
                   tmprhs = trim(adjustl(rhs(ioffset+1:)))
                   indx2 = index(rhs(ioffset+1:),",")
                   !indx2 = index(tmprhs,",")
                   !print *, ioffset, indx2
                   !print *, rhs(ioffset+1:ioffset+indx2-1)
                   
                   if (indx2 > 0) then
                      offcnt = offcnt+1
                   
                      ipos(iloop) = indx2
                      ioffset = indx2
                      !print *, "IPos, Ioff = ",ipos(iloop),ioffset
                   else
                      if (offcnt == 1) then
                         ipos(2) = len_trim(rhs(ioffset+1:))
                      endif
                   endif
                enddo

                print*, "IPOS 1, IPOS 2 = ", ipos(1), ipos(2)
                !print *, rhs(ipos(1)+1:ipos(2)-1)
                self%table_provides(1) = trim(adjustl(rhs(1:ipos(1)-1)))   
                self%table_provides(2) = trim(adjustl(rhs(ipos(1)+1:ipos(1)+ipos(2)-1)))

                print *, self%table_provides(1), self%table_provides(2)

                select case (self%table_provides(1))
                case default
                    stop "Not ROT"
                case("ROT","Rot","rot","ROTATION","Rotation","rotation")
                    continue
                end select

                select case (self%table_provides(2))
                case default
                    stop "Not TILT"
                case("TILT","Tilt","tilt","TIL","Til","til")
                    continue
                end select

                cycle READLOOP_K
            end select
        endif

        ! Everything else should be the set of angles that build up a sweep.
        if (string == "<sweep>") then

            self%sweep_count = self%sweep_count + 1
            print *, 'STRING = ', string

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

        print *, 'STRING = ', len(string)

        !read(string,*,iostat=ierr) rdrot, rdtilt, rdrot_err, rdtilt_err, dir, mds, bw, rdrot_sl1, rdtilt_sl1, &
        !    dir_sl1, bw_sl1, rdrot_sl2, rdtilt_sl2, dir_sl2, bw_sl2, rdrot_sl3, rdtilt_sl3, dir_sl3, bw_sl3, rdrot_sl4, rdtilt_sl4, &
        !    dir_sl4, bw_sl4
        read(string,*,iostat=ierr) rdrot(1), rdtilt(1), rdrot_err, rdtilt_err, rddirect(1), rdbwth(1), mds, rdrot(2), rdtilt(2), &
            rddirect(2), rdbwth(2), rdrot(3), rdtilt(3), rddirect(3), rdbwth(3), rdrot(4), rdtilt(4), rddirect(4), rdbwth(4), &
            rdrot(5), rdtilt(5), rddirect(5), rdbwth(5)    
        if ( ierr /= 0 ) then
            write(*,'("Problem reading scanning table file ''", A, "'', line ", I8)') trim(filename), k
            stop "MODULE_SCANNING : READ_SCANNING_TABLE_FILE"
        endif
        line_number = line_number + 1

        ! Determine the number of lobes in the scan file. Minimum number should be 1 to account for the main lobe (BWK, 5/6/26)
        if (self%beam_count == 1) then
             tmp_nlobe = count(rdrot(:) > -999.)
        endif

        self%beam_count = self%beam_count + 1
        if ( self%beam_count > size(temporary_rotation) ) then
            write(*,'("Too many entries.  Increase size of work arrays ''temporary_rotation'' and ''temporary_tilt''")')
            stop "MODULE_SCANNING : READ_SCANNING_TABLE_FILE"
        endif
        temporary_rotation(self%beam_count,1) = rdrot(1)
        temporary_tilt    (self%beam_count,1) = rdtilt(1)
        temporary_roterr(self%beam_count) = rdrot_err
        temporary_tilterr    (self%beam_count) = rdtilt_err
        temporary_dir(self%beam_count,1) = rddirect(1)
        temporary_mds(self%beam_count) = mds
        temporary_bw    (self%beam_count,1) = rdbwth(1)
        temporary_rotation(self%beam_count,2) = rdrot(2)
        temporary_tilt    (self%beam_count,2) = rdtilt(2)
        temporary_dir(self%beam_count,2) = rddirect(2)
        temporary_bw    (self%beam_count,2) = rdbwth(2)
        temporary_rotation(self%beam_count,3) = rdrot(3)
        temporary_tilt    (self%beam_count,3) = rdtilt(3)
        temporary_dir(self%beam_count,3) = rddirect(3)
        temporary_bw    (self%beam_count,3) = rdbwth(3)
        temporary_rotation(self%beam_count,4) = rdrot(4)
        temporary_tilt    (self%beam_count,4) = rdtilt(4)
        temporary_dir(self%beam_count,4) = rddirect(4)
        temporary_bw    (self%beam_count,4) = rdbwth(4)
        temporary_rotation(self%beam_count,5) = rdrot(5)
        temporary_tilt    (self%beam_count,5) = rdtilt(5)
        temporary_dir(self%beam_count,5) = rddirect(5)
        temporary_bw    (self%beam_count,5) = rdbwth(5)

        ! Part of the old setup - delete after testing new setup
        !temporary_rotation_sl1(self%beam_count) = rdrot_sl1
        !temporary_tilt_sl1    (self%beam_count) = rdtilt_sl1
        !temporary_dir_sl1(self%beam_count) = dir_sl1
        !temporary_bw_sl1    (self%beam_count) = bw_sl1
        !temporary_rotation_sl2(self%beam_count) = rdrot_sl2
        !temporary_tilt_sl2    (self%beam_count) = rdtilt_sl2
        !temporary_dir_sl2(self%beam_count) = dir_sl2
        !temporary_bw_sl2    (self%beam_count) = bw_sl2
        !temporary_rotation_sl3(self%beam_count) = rdrot_sl3
        !temporary_tilt_sl3    (self%beam_count) = rdtilt_sl3
        !temporary_dir_sl3(self%beam_count) = dir_sl3
        !temporary_bw_sl3    (self%beam_count) = bw_sl3
        !temporary_rotation_sl4(self%beam_count) = rdrot_sl4
        !temporary_tilt_sl4    (self%beam_count) = rdtilt_sl4
        !temporary_dir_sl4(self%beam_count) = dir_sl4
        !temporary_bw_sl4    (self%beam_count) = bw_sl4

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

    allocate(self%rotation(self%beam_count,5))
    allocate(self%tilt(self%beam_count,5))
    allocate(self%direct(self%beam_count,5))
    allocate(self%bwth(self%beam_count,5))
    allocate(self%rot_err(self%beam_count))
    allocate(self%tilt_err(self%beam_count))
    allocate(self%mds(self%beam_count))
    allocate(self%azim(self%beam_count,5))
    allocate(self%elev(self%beam_count,5))
    allocate(self%beamwidth_h(self%beam_count,5))
    allocate(self%beamwidth_v(self%beam_count,5))
    allocate(self%sll(self%beam_count,5))
    allocate(self%sweep_start_index(self%sweep_count))
    allocate(self%sweep_end_index  (self%sweep_count))

    self%nlobes = tmp_nlobe
    self%rotation = temporary_rotation(1:self%beam_count,1:5)
    self%tilt     = temporary_tilt    (1:self%beam_count,1:5)
    self%rot_err = temporary_roterr(1:self%beam_count)
    self%tilt_err     = temporary_tilterr    (1:self%beam_count)
    self%direct = temporary_dir(1:self%beam_count,1:5)
    self%mds = temporary_mds(1:self%beam_count)
    self%bwth     = temporary_bw    (1:self%beam_count,1:5)
    self%azim(:,:)     = -9999.
    self%elev(:,:)     = -9999.
    self%beamwidth_h     = -9999.
    self%beamwidth_v     = -9999.
    self%sll             = 0.
    self%sweep_start_index = temporary_start(1:self%sweep_count)
    self%sweep_end_index = temporary_end(1:self%sweep_count)
    
    ! Part of the old setup - to be deleted once tested
    !self%rot_sl1 = temporary_rotation_sl1(1:self%beam_count)
    !self%tilt_sl1     = temporary_tilt_sl1    (1:self%beam_count)
    !self%dir_sl1 = temporary_dir_sl1(1:self%beam_count)
    !self%bw_sl1     = temporary_bw_sl1    (1:self%beam_count)
    !self%rot_sl2 = temporary_rotation_sl2(1:self%beam_count)
    !self%tilt_sl2     = temporary_tilt_sl2    (1:self%beam_count)
    !self%dir_sl2 = temporary_dir_sl2(1:self%beam_count)
    !self%bw_sl2     = temporary_bw_sl2    (1:self%beam_count)
    !self%rot_sl3 = temporary_rotation_sl3(1:self%beam_count)
    !self%tilt_sl3     = temporary_tilt_sl3    (1:self%beam_count)
    !self%dir_sl3 = temporary_dir_sl3(1:self%beam_count)
    !self%bw_sl3     = temporary_bw_sl3    (1:self%beam_count)
    !self%rot_sl4 = temporary_rotation_sl4(1:self%beam_count)
    !self%tilt_sl4     = temporary_tilt_sl4    (1:self%beam_count)
    !self%dir_sl4 = temporary_dir_sl4(1:self%beam_count)
    !self%bw_sl4     = temporary_bw_sl4    (1:self%beam_count)
   

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

  end subroutine read_scanning_table_file_idt

end module module_scanning
