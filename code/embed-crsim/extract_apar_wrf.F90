program extract_apar
  use module_configuration, only : RKIND
  use module_access_wrf, only : MAX_WRF_FILES
  use module_access_wrf, only : MAX_WRF_TIMES
  use module_access_wrf, only : find_all_wrf_times
  use module_access_wrf, only : wrf_metadata_type
  use module_access_wrf, only : options_type
  use module_crsim_wrapper, only : crsim_wrapper
  use module_flightpath, only : aircraft_metadata_type
  use module_flightpath, only : flightpath_class
  use module_scanning,   only : scan_type
  use module_cfradial_output, only : cfradial_type
  use module_cfradial_output, only : volume_type
  use module_cfradial_output, only : volume_field_type
  use module_llxy, only : PROJ_MERC
  use module_llxy, only : PROJ_LC
  use module_llxy, only : proj_info
  use module_llxy, only : map_init
  use module_llxy, only : map_set
  use kwm_date_utilities, only : geth_newdate
  use module_external_attitude, only : use_external_attitudes
  use module_external_attitude, only : read_attitude_file
  use module_external_attitude, only : attitude_details_at_time
  use module_timing_utilities, only : timer_pause
  use module_timing_utilities, only : timer_resume
  use module_timing_utilities, only : timer_print
  use module_timing_utilities, only : timer_value

  use crsim_mod, only : conf_var
#ifdef _PARALLEL_
  ! use omp_lib
  use mpi
#endif

  implicit none

  ! Get the OMP variables
  integer             :: THREADS, THRD, RANK, ik, nperiods, k2
  integer             :: CPU_COUNT

#ifdef _PARALLEL_
  ! MPI initializations
  integer :: ierrmpi
  integer :: ierr
  integer :: mpi_real_selection
#endif

  integer :: scount
  integer :: icount
  real(kind=RKIND), dimension(100000) :: scantime_vector = 1.E36
  integer,          dimension(100000) :: ipanel_vector = -999999
  real(kind=RKIND)        :: scantime_offset_from_initial
  real(kind=RKIND)        :: total_cpu_time
  character(len=1024)     :: flnm, outflnm
  type(options_type)      :: options
  type(proj_info)         :: proj
  integer                 :: error_flag, iflag
  integer                 :: seedlen

  ! integer :: max_procs

  type(conf_var) :: conf ! CRSim config options

  type(aircraft_metadata_type)     :: aircraft
  type(volume_field_type), pointer :: ptr

  type (wrf_metadata_type) :: metaA! (40)  ! Assigned this value to account for maximum number of processors allowed
  type (wrf_metadata_type) :: metaB! (40)

  type panel_type
      character(len=1024)     :: output_filename_format_string
      type(scan_type)         :: scan
      type(volume_type)       :: volume
      integer                 :: sweep_number = 0
      integer                 :: scan_cycles_per_volume
      logical                 :: update_flag
      integer                 :: bwtype
      integer, dimension(256) :: seed
  end type panel_type
  type(panel_type), dimension(4), target :: panel
  type(volume_type), pointer :: volume => NULL()
  type(scan_type),   pointer :: scan => NULL()
  integer,           pointer :: scan_cycles_per_volume => NULL()

  character(len=1024) :: namelist_file = ""
  character(len=1024) :: flightpath_name = ""
  character(len=1024) :: fp_format_string = ""

  integer :: igate
  integer :: iray

  type(cfradial_type) :: cf

  character(len=1024), dimension(MAX_WRF_FILES) :: wrf_files
  real(kind=RKIND),    dimension(MAX_WRF_FILES*MAX_WRF_TIMES) :: wrf_xtimes
  character(len=19),   dimension(MAX_WRF_FILES*MAX_WRF_TIMES) :: wrf_times
  integer,             dimension(MAX_WRF_FILES*MAX_WRF_TIMES) :: wrf_file_index
  integer,             dimension(MAX_WRF_FILES*MAX_WRF_TIMES) :: wrf_time_index

  real(kind=RKIND) :: leg_current_time
  integer          :: i
  integer          :: j
  integer          :: WRFindexA
  integer          :: WRFindexB
  real(kind=RKIND) :: time_interpolation_factor
  real(kind=RKIND) :: dx
  integer          :: idim
  integer          :: jdim

  integer          :: isweep
  integer          :: ibeam
  integer          :: nwindow

  character(len=19) :: wrf_reference_time
  real(kind=RKIND), pointer, dimension(:,:) :: vx, vy, vz, zindxA, zindxB, Kwork, Kwork2, Kwork3

  !  Beamwidth vertex points
  real(kind=RKIND), pointer, dimension(:,:) :: vx_ul, vy_ul, vz_ul, zindxA_ul, zindxB_ul
  real(kind=RKIND), pointer, dimension(:,:) :: vx_ur, vy_ur, vz_ur, zindxA_ur, zindxB_ur
  real(kind=RKIND), pointer, dimension(:,:) :: vx_lr, vy_lr, vz_lr, zindxA_lr, zindxB_lr
  real(kind=RKIND), pointer, dimension(:,:) :: vx_ll, vy_ll, vz_ll, zindxA_ll, zindxB_ll

  type(flightpath_class) :: flightpath
  real(kind=RKIND)       :: beamtime

  character(len=256) :: l1, l2
  character(len=19)  :: ldate1 ! YYYY-mm-dd_HH:MM:SS
  character(len=19)  :: ldate2 ! YYYY-mm-dd_HH:MM:SS

  integer :: ipanel
  integer :: num_panels
  integer :: bwtype

  num_panels = command_argument_count()
  if (num_panels == 0) then
     write(*,'(80("*"))')
     write(*, '("*****",70X,"*****")')
     write(*, '("*****    No namelist found.  One or more namelist files must be specified  *****")')
     write(*, '("*****    as command-line arguments                                         *****")')
     write(*, '("*****",70X,"*****")')
     write(*,'(80("*"))')
     stop "ERROR EXIT: No namelist found.  One or more namelist files must be provided as command line arguments"
  endif
  
#ifdef _PARALLEL_
  call MPI_INIT(ierrmpi)
  call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, ierrmpi)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, cpu_count, ierrmpi)
  print *,'Processor id: ',RANK,' of ',cpu_count,' processors.'
  ! Choose an appropriate MPI data type based on the configurred Fortran REAL (kind=4) or (kind=8)
  if (RKIND == 8) then
      mpi_real_selection= MPI_REAL8
  else if (RKIND == 4) then
      mpi_real_selection= MPI_REAL
  endif
#else
  RANK = 0
  CPU_COUNT = 1
#endif

  call timer_resume("TOTAL")

  do ipanel = 1, num_panels
      call get_command_argument(ipanel, namelist_file)
      call namelist_options(namelist_file, seedlen, options, panel(ipanel)%scan, conf)
      panel(ipanel)%output_filename_format_string = options%output_filename_format_string
      panel(ipanel)%bwtype = options%bwtype
      call random_seed(get=panel(ipanel)%seed(1:seedlen))
  enddo

  call timer_resume("WRFTIMES")
  ! Added the conv_minute term to the following call (B. Klotz, 4/5/2021)  
  if (RANK == 0) then
      call find_all_wrf_times ( trim(options%wrf_glob_pattern), options%conv_minute, wrf_reference_time, &
           &                        wrf_files, wrf_times, wrf_xtimes, wrf_file_index, wrf_time_index, dx, idim, jdim )
  endif

#ifdef _PARALLEL_
  call MPI_BARRIER(MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(wrf_reference_time, 19, MPI_INTEGER1, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(wrf_files, 1024*MAX_WRF_FILES, MPI_INTEGER1, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(wrf_xtimes, MAX_WRF_FILES*MAX_WRF_TIMES, mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(wrf_times, 19*MAX_WRF_FILES*MAX_WRF_TIMES, MPI_INTEGER1, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(wrf_file_index, MAX_WRF_FILES*MAX_WRF_TIMES, MPI_INTEGER, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(wrf_time_index, MAX_WRF_FILES*MAX_WRF_TIMES, MPI_INTEGER, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(dx, 1, mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(idim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(jdim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BARRIER(MPI_COMM_WORLD, ierrmpi)
#endif

  !
  !  <wrf_reference_time> will be a 19-character string of the form "YYYY-MM-DD_hh:mm:ss".
  ! 

  !
  !    Idealized simulations, so set up an arbitrary map projection
  !

  call timer_pause("WRFTIMES")
  call timer_resume("MAP")

  call map_init(proj)
  call map_set ( PROJ_LC, proj, truelat1=33.0_RKIND, truelat2=37.0_RKIND, lat1=35.0_RKIND, lon1=-100.0_RKIND, &
       &         knowni=real(idim,kind=RKIND)/2.0, knownj=real(jdim,kind=RKIND)/2.0, stdlon=-100.0_RKIND, dx=dx )

  call timer_pause("MAP")

  do ipanel = 1, num_panels
      panel(ipanel)%scan_cycles_per_volume = 1
      panel(ipanel)%scan%scan_cycle = 0
  enddo

  if ( use_external_attitudes ) then
      ! Read an external netcdf file of aircraft attitudes

      ! call read_attitude_file ( "attitude.nc" , leg_initial_time, start_at=93.0, end_at=346.0 )
      call read_attitude_file ( options%leg_initial_time, start_at=1.0_RKIND )
  endif

  flightpath = flightpath_class ( options%waypoint, options%leg_initial_time, proj%dx )
  aircraft%track = flightpath%track
  aircraft%air_speed = options%air_speed
  aircraft%drift = 0.0
  aircraft%heading = aircraft%track - aircraft%drift

  if ( use_external_attitudes ) then
      call attitude_details_at_time ( options%leg_initial_time, aircraft%roll, aircraft%pitch, aircraft%drift, aircraft%heading, &
           &                                                    aircraft%dxdt, aircraft%dydt, aircraft%dzdt, error_flag )
      if ( error_flag  > 0 ) stop "Attitude details?"

      call attitude_details_at_time ( options%leg_initial_time, aircraft%roll, aircraft%pitch, aircraft%drift, aircraft%heading, &
           &                                                    aircraft%dxdt, aircraft%dydt, aircraft%dzdt, error_flag )
      if ( error_flag  > 0 ) stop "Attitude details?"
  endif

  !
  !  Loop over panels and scans for each panel, to created a list of times we'll make a scan, ordered by time, and the 
  !  panel associated with each scan.
  !

  if (RANK==0) then
      scount = 0
      do ipanel = 1, num_panels
          SCANTIMES: do icount = 0, 9999999
              scantime_offset_from_initial = icount * panel(ipanel)%scan%seconds_plus_skip
              ! Test on the end time of the scan
              ! If the end time is beyond leg_time_seconds, do ont include that scan.
              if (scantime_offset_from_initial + panel(ipanel)%scan%seconds_for_scan_cycle - panel(ipanel)%scan%skip_seconds_between_scans > options%leg_time_seconds) exit SCANTIMES
              leg_current_time = options%leg_initial_time + scantime_offset_from_initial
              ! Insert, probably pretty slow, but our arrays are small enough.
              JLOOP : do j = 1, scount+1
                  if (leg_current_time < scantime_vector(j)) then
                      scantime_vector(j+1:scount+1) = scantime_vector(j:scount)
                      ipanel_vector(j+1:scount+1) = ipanel_vector(j:scount)
                      scantime_vector(j) = leg_current_time
                      ipanel_vector(j) = ipanel
                      exit JLOOP
                  endif
              enddo JLOOP
              scount = scount + 1
          enddo SCANTIMES
      enddo
  endif

  call timer_resume("FLIGHTPATH")
  RANK_IF: if (RANK==0) then
      call flightpath%precompute(options, proj%dx, wrf_xtimes, wrf_times, wrf_files, wrf_file_index, wrf_time_index)
      call timer_print("OpenA")
      call timer_print("OpenB")
      flightpath_name = "flightpath.ascii"
      call flightpath%dump(trim(flightpath_name), rank)
  endif  RANK_IF
  call timer_pause("FLIGHTPATH")

#ifdef _PARALLEL_

  call timer_resume("BARRIER")
  call MPI_BARRIER(MPI_COMM_WORLD, ierrmpi)
  call timer_pause("BARRIER")
  call MPI_BCAST(scount, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(scantime_vector, scount, mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(ipanel_vector, scount, MPI_INTEGER, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(flightpath%pathcount, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(flightpath%time_history, (flightpath%pathcount+1), mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(flightpath%track_history, (flightpath%pathcount+1), mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(flightpath%tmp_track_history, (flightpath%pathcount+1), mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(flightpath%location_history, (3*(flightpath%pathcount+1)), mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(flightpath%ground_velocity_history, (3*(flightpath%pathcount+1)), mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(flightpath%wind_history, (3*(flightpath%pathcount+1)), mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(flightpath%time, 1, mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(flightpath%location, 3, mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(flightpath%wind, 3, mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call MPI_BCAST(flightpath%track, 1, mpi_real_selection, 0, MPI_COMM_WORLD, ierrmpi)
  call timer_resume("BARRIER")
  call MPI_BARRIER(MPI_COMM_WORLD, ierrmpi)
  call timer_pause("BARRIER")
#endif

  ! Traverse the flightpath, using the info in flightpath as a lookup table for the location of the aircraft.

  call timer_resume("LEG_LOOP")
  LEG: do icount = 1, scount
      
      if (rank /= (icount-1)*cpu_count/scount) then
          panel(ipanel)%scan%scan_cycle = panel(ipanel)%scan%scan_cycle + 1
          cycle LEG
      endif
      ipanel = ipanel_vector(icount)
      scan => panel(ipanel)%scan
      bwtype = panel(ipanel)%bwtype
      print*, "seed for panel = ", panel(ipanel)%seed(1:seedlen)
      call random_seed(put=panel(ipanel)%seed(1:seedlen))
      ! <leg_current_time> is seconds after <wrf_reference_time>
      leg_current_time = scantime_vector(icount)
      print*, rank, "leg_current_time = ", leg_current_time

      if (RANK+1 > CPU_COUNT) stop "RANK,CPU_COUNT"

      call timer_resume("WRFSEARCH")
      WRFindexA = -1
      WRFindexB = -1

      time_interpolation_factor = -999999  ! Default dummy value to be replaced in WRF_TIME_SEARCH loop
      
      WRF_TIME_SEARCH : do i = 1, MAX_WRF_FILES*MAX_WRF_TIMES
          if (wrf_file_index(i) < 0) exit WRF_TIME_SEARCH
          if (abs(wrf_xtimes(i) - leg_current_time) < -1.E3) then
              ! Time match
              time_interpolation_factor = 1.0
              WRFindexA = i
              WRFindexB = -1
              exit WRF_TIME_SEARCH
          else if (wrf_xtimes(i) > leg_current_time) then
              WRFindexB = i
              time_interpolation_factor = ( leg_current_time - wrf_xtimes(WRFindexA) ) / &
                   ( wrf_xtimes(WRFindexB) - wrf_xtimes(WRFindexA) )
              time_interpolation_factor = 1.0-time_interpolation_factor
              exit WRF_TIME_SEARCH
          endif
          WRFindexA = i
      enddo WRF_TIME_SEARCH
      if (time_interpolation_factor == -999999) then
         write(*,'("***** ERROR *****")')
         write(*,'("***** time_interpolation factor not set")')
         write(*,'("***** Bracketing WRF output times for time ", F20.10, " not found")') leg_current_time
         write(*,'("***** First model output time = ", F20.10, ":", 3X, A)') wrf_xtimes(1), wrf_times(1)
         do i = 1, MAX_WRF_FILES*MAX_WRF_TIMES
            if (wrf_xtimes(i) < -1.E25) then
               write(*,'("***** Final model output time = ", F20.10, ":", 3X, A)') wrf_xtimes(i-1), wrf_times(i-1)
               exit
            endif
         enddo
#ifdef _PARALLEL_
         call MPI_ABORT(MPI_COMM_WORLD, ierr, ierrmpi)
#endif
         exit LEG
      endif
      call timer_pause("WRFSEARCH")

      call timer_resume("WRFOPEN")
      if ((metaA%flnm /= wrf_files(wrf_file_index(WRFindexA))).or.(metaA%time_frame /= wrf_time_index(WRFindexA))) then
          call metaA%close()
          call metaA%open(wrf_files(wrf_file_index(WRFindexA)), wrf_time_index(WRFindexA))
      endif
      if (WRFindexB > 0) then
          if ((metaB%flnm /= wrf_files(wrf_file_index(WRFindexB))).or.(metaB%time_frame /= wrf_time_index(WRFindexB))) then
              call metaB%close()
              call metaB%open(wrf_files(wrf_file_index(WRFindexB)), wrf_time_index(WRFindexB))
          endif
      endif
      call timer_pause("WRFOPEN")

      call timer_resume("PRESWEEP")

      !
      !  So now we know the aircraft location.  Define all the rays of one scan cycle.
      !

      volume => panel(ipanel)%volume
      scan_cycles_per_volume => panel(ipanel)%scan_cycles_per_volume

      volume = volume_type ( volume%number+1, wrf_reference_time , int(leg_current_time) , &
           &                 scan%max_range_in_meters ,  scan%meters_between_gates ,       &
           &                 scan%fold_limit_lower, scan%fold_limit_upper,                 &
           &                 metaA%ncid, metaA%list_of_fields(1:metaA%nfields) )

      if (scan%sweep_count > size(volume%sweep_number)) then
          write(*,'("*** Problem:  sweep_count > max_sweeps:  ", I6, " > ",I6)') scan%sweep_count, size(volume%sweep_number)
          stop "*** increase MAX_SWEEPS in module_cfradial_output"
      endif

      volume%unambiguous_range(:) = 0.5 * ( 2.9979E8 * scan%pulse_repetition_time )

      volume%nrays = scan%beam_count
      if (volume%nrays > size(volume%time)) then
          write(*, '(/,"   *****")')
          write(*, '("   *****   Number of rays in scanning table exceeds size allocated in ")')
          write(*, '("   *****   module_cfradial_output:MAX_RAYS.")')
          write(*, '("   *****        NRAYS    = ", I10)') volume%nrays
          write(*, '("   *****        MAX_RAYS = ", I10)') size(volume%time)
          write(*, '("   *****   Increase parameter MAX_RAYS in module_cfradial_output, ")')
          write(*, '("   *****   or decrease number of rays in scanning table")')
          write(*, '("   *****",/)')
          stop "MAX_RAYS EXCEEDED"
      endif

      zindxA => volume%new_array("ZINDXA")
      zindxB => volume%new_array("ZINDXB")
      vx   => volume%new_array("VX")
      vy   => volume%new_array("VY")
      vz   => volume%new_array("VZ")

      if (bwtype > 0) then
          !   BW x, y, z indices, UL, UR, LR, LL
          zindxA_ul  => volume%new_array("ZINDXA_UL")
          zindxB_ul  => volume%new_array("ZINDXB_UL")
          vx_ul => volume%new_array("VX_UL")
          vy_ul => volume%new_array("VY_UL")
          vz_ul => volume%new_array("VZ_UL")

          zindxA_ur  => volume%new_array("ZINDXA_UR")
          zindxB_ur  => volume%new_array("ZINDXB_UR")
          vx_ur => volume%new_array("VX_UR")
          vy_ur => volume%new_array("VY_UR")
          vz_ur => volume%new_array("VZ_UR")

          zindxA_lr  => volume%new_array("ZINDXA_LR")
          zindxB_lr  => volume%new_array("ZINDXB_LR")
          vx_lr => volume%new_array("VX_LR")
          vy_lr => volume%new_array("VY_LR")
          vz_lr => volume%new_array("VZ_LR")

          zindxA_ll  => volume%new_array("ZINDXA_LL")
          zindxB_ll  => volume%new_array("ZINDXB_LL")
          vx_ll => volume%new_array("VX_LL")
          vy_ll => volume%new_array("VY_LL")
          vz_ll => volume%new_array("VZ_LL")
      endif

      call timer_pause("PRESWEEP")
      call timer_resume("SWEEP_LOOP")

      iray = 0
      beamtime = leg_current_time
      SWEEP_LOOP : do isweep = 1, scan%sweep_count
          volume%sweep = volume%sweep + 1
          volume%sweep_number(volume%sweep) = volume%sweep-1
          volume%sweep_mode(volume%sweep) = scan%sweep_mode
          volume%sweep_start_ray_index(volume%sweep) = scan%sweep_start_index(isweep)
          volume%sweep_end_ray_index(volume%sweep) = scan%sweep_end_index(isweep)

          call timer_resume("BEAM_LOOP")
          BEAM_LOOP : do ibeam = scan%sweep_start_index(isweep), scan%sweep_end_index(isweep)
              iray = iray + 1

              if (options%herky_jerky) then
                  !  Time stands still for beams and sweeps in a particular scan cycle
                  volume%time(iray) = int(leg_current_time) - volume%start_time_integer_seconds
                  ! print *, leg_current_time, volume%time_coverage_start, volume%start_time_integer_seconds, volume%time(1)
                  ! volume%time(iray) = leg_current_time
              else
                  !  Each beam gets an updated time

                  beamtime = leg_current_time + scan%beam_time_from_start_of_scan_cycle(iray)

                  !
                  !  volume%time for CFRadial output is described in the CFRadial documentation as 
                  !  the time of the center of the beam, but for now I'll keep it as the time of 
                  !  the last pulse in the beam.
                  !  It is in (fractional) seconds after <volume%time_interval_start>
                  !

                  volume%time(iray) = beamtime - volume%start_time_integer_seconds

                  ! write(*,'("ibeam, iray, beamtime = ", I5, I5, F12.6, F12.6)') ibeam, iray, beamtime, volume%time(iray)

              endif

              ! Extract instantaneous aircraft details from the flightpath information.
              call flightpath%getloc(beamtime, proj, aircraft, error_flag)
              if ( error_flag > 0 ) exit LEG

              ! Fill in beam-specific details of the volume structure.
              call beam_geometry ( aircraft , scan , iray , volume , proj%dx, conf%Theta1, bwtype, options%ref_angle )

              call timer_resume("GATE_LOOP")
              GATE_LOOP : do igate = 1, volume%ngates

                  zindxA(igate,iray) = metaA%find_z_index(vx(igate,iray), vy(igate,iray), vz(igate,iray))
                  if (bwtype > 0) then
                      zindxA_ul(igate,iray) = metaA%find_z_index(vx_ul(igate,iray), vy_ul(igate,iray), vz_ul(igate,iray))
                      zindxA_ur(igate,iray) = metaA%find_z_index(vx_ur(igate,iray), vy_ur(igate,iray), vz_ur(igate,iray))
                      zindxA_lr(igate,iray) = metaA%find_z_index(vx_lr(igate,iray), vy_lr(igate,iray), vz_lr(igate,iray))
                      zindxA_ll(igate,iray) = metaA%find_z_index(vx_ll(igate,iray), vy_ll(igate,iray), vz_ll(igate,iray))
                  endif
                  if ( time_interpolation_factor < 1.0 ) then
                      ZindxB(igate,iray) = metaB%find_z_index(vx(igate,iray), vy(igate,iray), vz(igate,iray))
                      if (bwtype > 0) then
                          zindxB_ul(igate,iray) = metaB%find_z_index(vx_ul(igate,iray), vy_ul(igate,iray), vz_ul(igate,iray))
                          zindxB_ur(igate,iray) = metaB%find_z_index(vx_ur(igate,iray), vy_ur(igate,iray), vz_ur(igate,iray))
                          zindxB_lr(igate,iray) = metaB%find_z_index(vx_lr(igate,iray), vy_lr(igate,iray), vz_lr(igate,iray))
                          zindxB_ll(igate,iray) = metaB%find_z_index(vx_ll(igate,iray), vy_ll(igate,iray), vz_ll(igate,iray))
                      endif
                  endif

              enddo GATE_LOOP
              call timer_pause("GATE_LOOP")

          enddo BEAM_LOOP
          call timer_pause("BEAM_LOOP")

          !
          ! fixed_angle depends on the sweep_mode, and on the coordinate system in use (Z, Y-Prime, etc.)
          !

          select case ( scan%sweep_mode )
          case ( "rhi" )
              select case (scan%primary_axis)
              case ("Z", "axis_z")
                  volume%fixed_angle(volume%sweep) = volume%azimuth( scan%sweep_start_index(isweep) + 1 )
              case ("Y-Prime", "axis_y_prime")
                   volume%fixed_angle(volume%sweep) = volume%elevation( scan%sweep_start_index(isweep) + 1 )
                  !volume%fixed_angle(volume%sweep) = volume%tilt( scan%sweep_start_index(isweep) + 1 )
              case ("Y", "axis_y")
                  volume%fixed_angle(volume%sweep) = volume%azimuth( scan%sweep_start_index(isweep) + 1 ) 
              case ("X","axis_x")    
                  volume%fixed_angle(volume%sweep) = volume%elevation( scan%sweep_start_index(isweep) + 1 ) 
              case default
                  write(*,'("Unrecognized primary axis (",A,") for sweep_moode = "A)') trim(scan%primary_axis), trim(scan%sweep_mode)
                  stop
              end select
          case ( "sector" )
              select case (scan%primary_axis)
              case ("Z", "axis_z")
                  volume%fixed_angle(volume%sweep) = scan%tilt( scan%sweep_start_index(isweep) + 1 )
              case default
                  write(*,'("Unrecognized primary axis (",A,") for sweep_moode = "A)') trim(scan%primary_axis), trim(scan%sweep_mode)
                  stop
              end select
          case ( "elevation_surveillance" )
                select case (scan%primary_axis)
                case ("Z", "axis_z")
                   volume%fixed_angle(volume%sweep) = scan%tilt( scan%sweep_start_index(isweep) + 1 )
                case ("Y-Prime", "axis_y_prime")
                   volume%fixed_angle(volume%sweep) = volume%tilt( scan%sweep_start_index(isweep) + 1 )
                case default
                   write(*,'("Unrecognized primary axis (",A,") for sweep_moode = "A)') trim(scan%primary_axis), trim(scan%sweep_mode)
                   stop
                end select
          case default
              write(*,'("Unrecognized sweep_mode:  ",A)') trim(scan%sweep_mode)
              stop
          end select
          ! volume%fixed_angle(volume%sweep) = scan%tilt( scan%sweep_start_index(isweep) + 1 )

      enddo SWEEP_LOOP

      nullify(vx_ul, vy_ul, vz_ul, zindxA_ul, zindxB_ul)
      nullify(vx_ur, vy_ur, vz_ur, zindxA_ur, zindxB_ur)
      nullify(vx_lr, vy_lr, vz_lr, zindxA_lr, zindxB_lr)
      nullify(vx_ll, vy_ll, vz_ll, zindxA_ll, zindxB_ll)

      call timer_pause("SWEEP_LOOP")

      !
      ! Having found the X/Y/Z/Time coordinates of the entire scan, interpolate other fields to those locations.
      !

      call timer_resume("POSTSWEEP")

      call timer_resume("PTR_INTERP")
      ptr => volume%next_point ( NULL() )
      do while ( associated ( ptr ) )
          select case ( ptr%source )
          case ( "WRF" )
              call timer_resume("PTRALLOC")
              if ( allocated ( ptr%data ) ) then
                  write(*, '("ptr%field_name ", A, " already allocated.")') trim(ptr%field_name)
              else
                  allocate ( ptr%data ( volume%ngates, volume%nrays ) )
              endif
              call timer_pause("PTRALLOC")
              if (ptr%field_name == "U") then
                  continue
              elseif (ptr%field_name == "V") then
                  continue
              elseif (ptr%field_name == "W") then
                  continue
              elseif (ptr%field_name == "HT") then
                  continue
              elseif (ptr%field_name == "TEMPERATURE") then
                  continue
              elseif (ptr%field_name == "RHO_D") then
                  continue
              elseif (ptr%field_name == "RHO_DS") then    ! Added by BWK, 3/14/2022
                  continue
              else
                  call timer_resume("VARINTERP")
                  if (bwtype == 0) then
                      call metaA%interp ( vx , vy , zindxA , ptr%field_name , ptr%data )
                      if (time_interpolation_factor < 1.0 ) then
                          allocate(Kwork(volume%ngates, volume%nrays))
                          call metaB%interp ( vx , vy , zindxB , ptr%field_name , Kwork )
                          ptr%data = ptr%data*time_interpolation_factor + Kwork*(1.0-time_interpolation_factor)
                          deallocate(Kwork)
                      endif
                  else
                      call metaA%interp_varbw ("A", bwtype, volume, scan, proj%dx, ptr%field_name, ptr%data )
                      if (time_interpolation_factor < 1.0 ) then
                          allocate(Kwork(volume%ngates, volume%nrays))
                          call metaB%interp_varbw ("B", bwtype, volume, scan, proj%dx, ptr%field_name, Kwork )
                          ptr%data = ptr%data*time_interpolation_factor + Kwork*(1.0-time_interpolation_factor)
                          deallocate(Kwork)
                      endif
                  endif
                  call timer_pause("VARINTERP")
              endif
          end select
          ptr => volume%next_point(ptr)
      enddo
      call timer_pause("PTR_INTERP")

      !  Treat time evolution of derived fields between WRF history snapshots (metaA and metaB).
      !  Put this section into a subroutine (in a module or class somewhere, or maybe standalone?)
      call timer_resume("KU_INTERP")

      ! Compute Ku at time A
      Kwork => volume%new_array("Ku")
      Kwork2 => volume%point_to_data("U")
      if (bwtype == 0) then
          call metaA%Ku(vx, vy, zindxA, Kwork, Kwork2)
      else
          call metaA%Ku_varbw("A", bwtype, volume, scan, proj%dx, Kwork, Kwork2)
      endif

      nullify(Kwork)
      nullify(Kwork2)

      ! Interpolate U and Ku to appropriate time
      if (time_interpolation_factor < 1.0 ) then

          ! Compute Ku at time B
          Kwork => volume%new_array("Ku_timeB")
          Kwork2 => volume%new_array("U_timeB")
          if (bwtype == 0) then
              call metaB%Ku(vx, vy, zindxB, Kwork, Kwork2)
          else
              call metaB%Ku_varbw("B", bwtype, volume, scan, proj%dx, Kwork, Kwork2)
          endif
          nullify(Kwork)
          nullify(Kwork2)

          Kwork => volume%point_to_data("U")
          Kwork2 => volume%point_to_data("U_timeB")
          where (Kwork > -9998 .and. Kwork2 > -9998)
              Kwork = Kwork*time_interpolation_factor + Kwork2*(1.0-time_interpolation_factor)
          elsewhere
              Kwork = -9999.0
          end where
          nullify(Kwork)
          nullify(Kwork2)

          Kwork => volume%point_to_data("Ku")
          Kwork2 => volume%point_to_data("Ku_timeB")
          where (Kwork > -9998 .and. Kwork2 > -9998)
              Kwork = Kwork*time_interpolation_factor + Kwork2*(1.0-time_interpolation_factor)
          elsewhere
              Kwork = -9999.0
          end where
          nullify(Kwork)
          nullify(Kwork2)
      endif
      call timer_pause("KU_INTERP")

      !  Treat time evolution between WRF history snapshots for other derived fields.
      call timer_resume("KV_INTERP")
      Kwork => volume%new_array("Kv")
      Kwork2 => volume%point_to_data("V")
      if (bwtype == 0) then
         call metaA%Kv(vx, vy, zindxA, Kwork, Kwork2)
      else
         call metaA%Kv_varbw("A", bwtype, volume, scan, proj%dx, Kwork, Kwork2)   ! Added by B. Klotz (12/18/2023), part of the variable beamwidth
      endif
      nullify(Kwork)
      nullify(Kwork2)

      ! Interpolate V and Kv to appropriate time
      if (time_interpolation_factor < 1.0 ) then

          ! Compute Kv at time B
          Kwork => volume%new_array("Kv_timeB")
          Kwork2 => volume%new_array("V_timeB")
          if (bwtype == 0) then
             call metaB%Kv(vx, vy, zindxB, Kwork, Kwork2)
          else
             call metaB%Kv_varbw("B", bwtype, volume, scan, proj%dx, Kwork, Kwork2)
          endif
          nullify(Kwork)
          nullify(Kwork2)

          Kwork => volume%point_to_data("V")
          Kwork2 => volume%point_to_data("V_timeB")
          where (Kwork > -9998 .and. Kwork2 > -9998)
              Kwork = Kwork*time_interpolation_factor + Kwork2*(1.0-time_interpolation_factor)
          elsewhere
              Kwork = -9999.0
          end where
          nullify(Kwork)
          nullify(Kwork2)

          Kwork => volume%point_to_data("Kv")
          Kwork2 => volume%point_to_data("Kv_timeB")
          where (Kwork > -9998 .and. Kwork2 > -9998)
              Kwork = Kwork*time_interpolation_factor + Kwork2*(1.0-time_interpolation_factor)
          elsewhere
              Kwork = -9999.0
          end where
          nullify(Kwork)
          nullify(Kwork2)

      endif
      call timer_pause("KV_INTERP")

      ! treat time evolution with multiple WRF datasets (metaA and metaB).
      call timer_resume("KW_INTERP")
      Kwork => volume%new_array("Kw")
      Kwork2 => volume%point_to_data("W")
      Kwork3 => volume%point_to_data("HT")
      if (bwtype == 0) then
          call metaA%Kw(vx, vy, zindxA, Kwork, Kwork2, Kwork3)
      else
          call metaA%Kw_varbw("A", bwtype, volume, scan, proj%dx, Kwork, Kwork2, Kwork3)
      endif
      nullify(Kwork)
      nullify(Kwork2)
      nullify(Kwork3)

      ! Interpolate W and Kw to appropriate time
      if (time_interpolation_factor < 1.0 ) then
          ! Compute Kw at time B
          Kwork => volume%new_array("Kw_timeB")
          Kwork2 => volume%new_array("W_timeB")
          Kwork3 => volume%new_array("HT_timeB")
          if (bwtype == 0) then
              call metaB%Kw(vx, vy, zindxB, Kwork, Kwork2, Kwork3)
          else
              call metaB%Kw_varbw("B", bwtype, volume, scan, proj%dx, Kwork, Kwork2, Kwork3)
          endif
          nullify(Kwork)
          nullify(Kwork2)
          nullify(Kwork3)

          Kwork => volume%point_to_data("W")
          Kwork2 => volume%point_to_data("W_timeB")
          where (Kwork > -9998 .and. Kwork2 > -9998)
              Kwork = Kwork*time_interpolation_factor + Kwork2*(1.0-time_interpolation_factor)
          elsewhere
              Kwork = -9999.0
          end where
          nullify(Kwork)
          nullify(Kwork2)

          Kwork => volume%point_to_data("Kw")
          Kwork2 => volume%point_to_data("Kw_timeB")
          where (Kwork > -9998 .and. Kwork2 > -9998)
              Kwork = Kwork*time_interpolation_factor + Kwork2*(1.0-time_interpolation_factor)
          elsewhere
              Kwork = -9999.0
          end where
          nullify(Kwork)
          nullify(Kwork2)

          Kwork => volume%point_to_data("HT")
          Kwork2 => volume%point_to_data("HT_timeB")
          where (Kwork > -9998 .and. Kwork2 > -9998)
              Kwork = Kwork*time_interpolation_factor + Kwork2*(1.0-time_interpolation_factor)
          elsewhere
              Kwork = -9999.0
          end where
          nullify(Kwork)
          nullify(Kwork2)
      endif
      call timer_pause("KW_INTERP")

      ! treat time evolution with multiple WRF datasets (metaA and metaB).
      ! All associations of KWork3 are added by BWK, 3/14/2022
      call timer_resume("OTHER_INTERP")
      Kwork  => volume%point_to_data("RHO_D")
      Kwork2 => volume%point_to_data("TEMPERATURE")
      Kwork3 => volume%point_to_data("RHO_DS")
      if (bwtype == 0) then
          call metaA%RHO_D(vx, vy, zindxA, Kwork, Kwork2, Kwork3) ! Added Kwork3 (BWK, 3/14/2022)
      else
          call metaA%RHO_D_varbw("A", bwtype, volume, scan, proj%dx, Kwork, Kwork2, Kwork3)
      endif
      nullify(Kwork)
      nullify(Kwork2)

      if (time_interpolation_factor < 1.0) then
          ! Compute RHO_D and TEMPERATURE at time B
          Kwork  => volume%new_array("RHO_D_timeB")
          Kwork2 => volume%new_array("TEMPERATURE_timeB")
          Kwork3 => volume%new_array("RHO_DS_timeB")

          if (bwtype == 0) then
              call metaB%RHO_D(vx, vy, zindxB, Kwork, Kwork2, Kwork3)
          else
              call metaB%RHO_D_varbw("B", bwtype, volume, scan, proj%dx, Kwork, Kwork2, Kwork3)
          endif
          nullify(Kwork)
          nullify(Kwork2)
          nullify(Kwork3)

          Kwork => volume%point_to_data("RHO_D")
          Kwork2 => volume%point_to_data("RHO_D_timeB")
          where (Kwork > -9998 .and. Kwork2 > -9998)
              Kwork = Kwork*time_interpolation_factor + Kwork2*(1.0-time_interpolation_factor)
          elsewhere
              Kwork = -9999.0
          end where
          nullify(Kwork)
          nullify(Kwork2)

          Kwork => volume%point_to_data("TEMPERATURE")
          Kwork2 => volume%point_to_data("TEMPERATURE_timeB")
          where (Kwork > -9998 .and. Kwork2 > -9998)
              Kwork = Kwork*time_interpolation_factor + Kwork2*(1.0-time_interpolation_factor)
          elsewhere
              Kwork = -9999.0
          end where
          nullify(Kwork)
          nullify(Kwork2)

          Kwork => volume%point_to_data("RHO_DS")
          Kwork2 => volume%point_to_data("RHO_DS_timeB")
          where (Kwork > -9998 .and. Kwork2 > -9998)
              Kwork = Kwork*time_interpolation_factor + Kwork2*(1.0-time_interpolation_factor)
          elsewhere
              Kwork = -9999.0
          end where
          nullify(Kwork)
          nullify(Kwork2)

      endif
      call timer_pause("OTHER_INTERP")
      nullify(vx   , vy   , vz   , zindxA   , zindxB   )
      nullify(vx_ul, vy_ul, vz_ul, zindxA_ul, zindxB_ul)
      nullify(vx_ur, vy_ur, vz_ur, zindxA_ur, zindxB_ur)
      nullify(vx_lr, vy_lr, vz_lr, zindxA_lr, zindxB_lr)
      nullify(vx_ll, vy_ll, vz_ll, zindxA_ll, zindxB_ll)

      call timer_pause("POSTSWEEP")

      !
      !  Compute radial velocity.
      !

      call timer_resume("RADVEL")
      call volume%radial_velocity ( proj%dx )
      call timer_pause("RADVEL")

      !
      !  Call CRSIM routines to compute radar moments etc. on our defined rays.
      !

      call timer_resume("CRSIM")
      call crsim_wrapper(conf, volume)
      call timer_pause("CRSIM")

      !
      ! Compute new velocity with fall velocites included
      !

      call volume%radial_velocity_fvel ( proj%dx )

      !
      !  Consider attenuation.
      !

      call timer_resume("ATTENUATE")
      call volume%attenuate(scan%meters_between_gates)
      call timer_pause("ATTENUATE")

      !
      !  Compute Signal-to-Noise Ratio (SNR)
      !

      call timer_resume("SNR")
      call volume%compute_snr()
      call timer_pause("SNR")

      !
      !  Compute Independent Pulse Sampling relations.
      !

      call timer_resume("IPS")
      call volume%compute_ips(scan%revisits_per_acquisition_time, scan%pulse_repetition_time, conf%freq)
      call timer_pause("IPS")

      !
      !  Add noise to Zhh_attenuated
      !

      call timer_resume("NOISE")
      call volume%add_zhh_noise(scan%snr_mask_threshold)
      call timer_pause("NOISE")

      !
      ! Compute PHIdp from Kdp and Differential Back Phase
      !

      call timer_resume("PHIDP")
      call volume%compute_phidp( scan%snr_mask_threshold )
      call timer_pause("PHIDP")

      !
      !  Apply radial velocity folding, add noise to velocity, etc.
      !

      call timer_resume("FOLD")
      call volume%fold_radial_velocity(scan%snr_mask_threshold)
      call timer_pause("FOLD")

      !
      !  Get beamwidth from CRSIM configuration.  Assuming horizontal and vertical beamwidths are the same.
      !

      volume%beamwidth_h = conf%Theta1
      volume%beamwidth_v = conf%Theta1

      !
      !  Write some number of scan cycles as a single volume.
      !

      if ( mod ( scan%scan_cycle + 1, scan_cycles_per_volume ) == 0 ) then
          call timer_resume("CFOUTPUT")
          write(outflnm, trim(panel(ipanel)%output_filename_format_string) ) volume%number

          call geth_newdate(ldate1, volume%time_coverage_start, int(volume%time(1)))
          ! print*, "ldate1 = ", volume%time_coverage_start, " + ", (volume%time(1)), " = ", ldate1
          call geth_newdate(ldate2, volume%time_coverage_start, int(volume%time(volume%nrays)))
          ! print*, "ldate2 = ", volume%time_coverage_start, " + ", (volume%time(volume%nrays)), " = ", ldate2

          ! print*, "Fold upper limit = ", volume%fold_limit_upper
          ! print*, "Fold lower limit = ", volume%fold_limit_lower
          ! call geth_newdate(ldate1, volume%reference_time, int(volume%time(1)))
          ! print*, "ldate1 = ", volume%reference_time, " + ", (volume%time(1)), " = ", ldate1
          ! call geth_newdate(ldate2, volume%reference_time, int(volume%time(volume%nrays)))
          ! print*, "ldate2 = ", volume%reference_time, " + ", (volume%time(volume%nrays)), " = ", ldate2

          ! call geth_newdate(ldate1, volume%time_coverage_start, int(volume%time(1)-options%leg_initial_time))
          ! print*, "ldate1 = ", volume%time_coverage_start, " + ", volume%time(1), "-",options%leg_initial_time, " = ", ldate1
          ! call geth_newdate(ldate2, volume%time_coverage_start, int(volume%time(volume%nrays)-options%leg_initial_time))
          ! print*, "ldate2 = ", volume%time_coverage_start, " + ", volume%time(volume%nrays), "-",options%leg_initial_time, " = ", ldate2

          write(l1,'(A,"_",A,F4.3)') ldate1(1:4)//ldate1(6:7)//ldate1(9:10),ldate1(12:13)//ldate1(15:16)//ldate1(18:19), volume%time(1)-floor(volume%time(1))
          write(l2,'(A,"_",A,F4.3)') ldate2(1:4)//ldate2(6:7)//ldate2(9:10),ldate2(12:13)//ldate2(15:16)//ldate2(18:19), volume%time(volume%nrays)-floor(volume%time(volume%nrays))
          write(outflnm, trim(panel(ipanel)%output_filename_format_string) ) trim(adjustl(l1)), trim(adjustl(l2))
          call cf%open ( trim(outflnm), trim(namelist_file), scan%meters_between_gates, scan%meters_to_center_of_first_gate, &
               &         volume%nrays, volume%ngates, volume%sweep, volume%time_coverage_start, &
               &         volume%fold_limit_lower, volume%fold_limit_upper )
          call cf%prepare_metadata(volume)
          call cf%write_volume(volume, scan%primary_axis)
          call cf%close()
          call timer_pause("CFOUTPUT")
      endif

      scan%scan_cycle = scan%scan_cycle + 1
      call panel(ipanel)%volume%final()
      call random_seed(get=panel(ipanel)%seed(1:seedlen))
      nullify(scan, volume, scan_cycles_per_volume)
  enddo LEG

  call timer_pause("LEG_LOOP")
  call timer_pause("TOTAL")
  total_cpu_time = timer_value("TOTAL")

  do i = 0, cpu_count-1
#ifdef _PARALLEL_
      call MPI_BARRIER(MPI_COMM_WORLD, ierrmpi)
#endif
      if (i /= RANK) cycle
      print*, "RANK = ", RANK, "----------------------------"
      call timer_print("TOTAL"       , label="total time (s):")
#ifdef _PARALLEL_
      call timer_print("BARRIER"     , label="total wait time (s):"              , percent_of=total_cpu_time)
#endif
      call timer_print("FLIGHTPATH"  , label="Flightpath time (s):"              , percent_of=total_cpu_time)
      call timer_print("WRFTIMES"    , label="WRFtimes time (s):"                , percent_of=total_cpu_time)
      call timer_print("MAP"         , label="Map time (s):"                     , percent_of=total_cpu_time)
      call timer_print("LEG_LOOP"    , label="Leg loop time (s):"                , percent_of=total_cpu_time)
      call timer_print("WRFSEARCH"   , label="WRF File search total time (s):"    , percent_of=total_cpu_time)
      call timer_print("WRFOPEN"     , label="WRF File open total time (s):"      , percent_of=total_cpu_time)
      call timer_print("PRESWEEP"    , label="Pre sweep loop total time (s):"     , percent_of=total_cpu_time)
      call timer_print("GATE_LOOP"   , label="Total gate loop time (s):"          , percent_of=total_cpu_time)
      call timer_print("BEAM_LOOP"   , label="Total beam loop time (s):"          , percent_of=total_cpu_time)
      call timer_print("SWEEP_LOOP"  , label="Total sweep loop time (s):"         , percent_of=total_cpu_time)
      call timer_print("PTRALLOC"    , label="Total Pointer allocation time (s):" , percent_of=total_cpu_time)
      call timer_print("VARINTERP"   , label="Total var interp time (s):"         , percent_of=total_cpu_time)
      call timer_print("PTR_INTERP"  , label="Total data pointer interp time (s):", percent_of=total_cpu_time)
      call timer_print("KU_INTERP"   , label="Total Ku interp time (s):"          , percent_of=total_cpu_time)
      call timer_print("KV_INTERP"   , label="Total Kv interp time (s):"          , percent_of=total_cpu_time)
      call timer_print("KW_INTERP"   , label="Total Kw interp time (s):"          , percent_of=total_cpu_time)
      call timer_print("OTHER_INTERP", label="Total Other interp time (s):"       , percent_of=total_cpu_time)
      call timer_print("POSTSWEEP"   , label="Post sweep loop total time (s):"    , percent_of=total_cpu_time)
      call timer_print("RADVEL"      , label="Total radial velocity time (s):"    , percent_of=total_cpu_time)
      call timer_print("CRSIM"       , label='Total crsim time (s):'              , percent_of=total_cpu_time)
      call timer_print("ATTENUATE"   , label="Total attenuate time (s):"          , percent_of=total_cpu_time)
      call timer_print("SNR"         , label="Total snr time (s):"                , percent_of=total_cpu_time)
      call timer_print("IPS"         , label="Total ips time (s):"                , percent_of=total_cpu_time)
      call timer_print("NOISE"       , label="Total noise time (s):"              , percent_of=total_cpu_time)
      call timer_print("PHIDP"       , label="Total phidp time (s):"              , percent_of=total_cpu_time)
      call timer_print("FOLD"        , label="Total fold time (s):"               , percent_of=total_cpu_time)
      call timer_print("CFOUTPUT"    , label="Total CfRadial output time (s):"    , percent_of=total_cpu_time)
  enddo

#ifdef _PARALLEL_
  call MPI_Finalize(ierrmpi)
#endif

end program extract_apar
   
!
!------------------------------------------------------------------------------------------------------------------
!

subroutine beam_geometry(aircraft, scan, ibeam, volume, grid_dx, theta1, bw_type, pref_ang)
  use module_configuration, only : RKIND
  use module_geometry,   only : coordinate_transformation
  use module_flightpath, only : aircraft_metadata_type
  use module_scanning,   only : scan_type
  use module_cfradial_output, only : volume_type
  implicit none
  type(aircraft_metadata_type), intent(in)    :: aircraft
  type(scan_type),              intent(in)    :: scan
  integer,                      intent(in)    :: ibeam
  type(volume_type),            intent(inout) :: volume
  real(kind=RKIND),             intent(in)    :: grid_dx
  real(kind=RKIND),             intent(in)    :: theta1
  integer,                      intent(in)    :: bw_type
  real(kind=RKIND),             intent(in)    :: pref_ang

  integer :: igate
  integer :: verts

  real(kind=RKIND), pointer, dimension(:,:) :: vx
  real(kind=RKIND), pointer, dimension(:,:) :: vy
  real(kind=RKIND), pointer, dimension(:,:) :: vz

  real(kind=RKIND), pointer, dimension(:,:) :: vx_ul, vy_ul, vz_ul   ! Upper left corner of beam position
  real(kind=RKIND), pointer, dimension(:,:) :: vx_ur, vy_ur, vz_ur   ! Upper right corner of beam position
  real(kind=RKIND), pointer, dimension(:,:) :: vx_lr, vy_lr, vz_lr   ! Lower right corner of beam position
  real(kind=RKIND), pointer, dimension(:,:) :: vx_ll, vy_ll, vz_ll   ! Lower left corner of beam position

  real(kind=RKIND)    :: azimuth
  real(kind=RKIND)    :: elevation
  real(kind=RKIND)    :: x_unit_dist
  real(kind=RKIND)    :: y_unit_dist
  real(kind=RKIND)    :: z_unit_dist
  real(kind=RKIND)    :: bwidthuse_h
  real(kind=RKIND)    :: bwidthuse_v
  real(kind=RKIND)    :: gate_dist
  real(kind=RKIND)    :: gate_x_index
  real(kind=RKIND)    :: gate_y_index
  real(kind=RKIND)    :: gate_z_meters

  !  For incorporation of variable beamwidth
  real(kind=RKIND), dimension(4) :: azm_bw_tmp, elev_bw_tmp, x_bw, y_bw, z_bw
  real(kind=RKIND), dimension(4) :: gate_xind_bw, gate_yind_bw, gate_zmeters_bw

  volume%latitude  ( ibeam ) = aircraft%lat
  volume%longitude ( ibeam ) = aircraft%lon
  volume%altitude  ( ibeam ) = aircraft%z
  volume%heading   ( ibeam ) = aircraft%heading
  volume%roll      ( ibeam ) = aircraft%roll
  volume%pitch     ( ibeam ) = aircraft%pitch
  volume%drift     ( ibeam ) = aircraft%drift
  volume%rotation  ( ibeam ) = scan%rotation(ibeam)
  volume%tilt      ( ibeam ) = scan%tilt(ibeam)
  volume%eastward_velocity  ( ibeam ) = aircraft%dxdt
  volume%northward_velocity ( ibeam ) = aircraft%dydt
  volume%vertical_velocity  ( ibeam ) = aircraft%dzdt
  volume%eastward_wind      ( ibeam ) = aircraft%u
  volume%northward_wind     ( ibeam ) = aircraft%v
  volume%vertical_wind      ( ibeam ) = aircraft%w
  volume%heading_rate       ( ibeam ) = 0.0
  volume%roll_rate          ( ibeam ) = 0.0
  volume%pitch_rate         ( ibeam ) = 0.0

  volume%aircraft_xgrid     ( ibeam ) = aircraft%xgrid
  volume%aircraft_ygrid     ( ibeam ) = aircraft%ygrid
  volume%aircraft_z         ( ibeam ) = aircraft%z

  ! Convert radar rotation and tilt (which are relative to the airframe) and aircraft roll and pitch (zero for now),
  ! to azimuth and elevation (which are relative to true north and true horizontal)

  ! call coordinate_transformation ( trim(scan%primary_axis), volume%tilt(ibeam), volume%rotation(ibeam), &
  !      &                           volume%heading(ibeam), volume%roll(ibeam), volume%pitch(ibeam), &
  !      &                           azimuth, elevation, x_unit_dist, y_unit_dist, z_unit_dist )

  call coordinate_transformation ( trim(scan%primary_axis), volume%tilt(ibeam), volume%rotation(ibeam), &
       &                           volume%heading(ibeam), volume%roll(ibeam), volume%pitch(ibeam), theta1, pref_ang, &
       &                           azimuth, elevation, x_unit_dist, y_unit_dist, z_unit_dist, bwidthuse_h, bwidthuse_v )

  volume%azimuth         ( ibeam ) = azimuth
  volume%elevation       ( ibeam ) = elevation
  volume%georefs_applied ( ibeam ) = 1

  !
  ! Assign the beamwidth to the beam based on the sine of the elevation
  !

  if (bw_type < 2) then
      volume%beamwidth_h( ibeam ) = theta1 !+ sin(abs(azimuth * RAD_PER_DEG))   !!! This needs to be edited to account for the other panels, right now it assumes the side panels only
      volume%beamwidth_v( ibeam ) = theta1 !+ sin(abs(elevation * RAD_PER_DEG))
  elseif (bw_type == 2) then
      !volume%beamwidth_h( ibeam ) = theta1 + sin(abs(azimuth * RAD_PER_DEG))   !!! This needs to be edited to account for the other panels, right now it assumes the side panels only
      !volume%beamwidth_v( ibeam ) = theta1 + sin(abs(elevation * RAD_PER_DEG))
      volume%beamwidth_h( ibeam ) = bwidthuse_h
      volume%beamwidth_v( ibeam ) = bwidthuse_v
  endif

  vx   => volume%point_to_data("VX")
  vy   => volume%point_to_data("VY")
  vz   => volume%point_to_data("VZ")

  if (bw_type > 0) then
     ! Upper left indices of beam
     vx_ul => volume%point_to_data("VX_UL")
     vy_ul => volume%point_to_data("VY_UL")
     vz_ul => volume%point_to_data("VZ_UL")

     ! Upper right indices of beam
     vx_ur => volume%point_to_data("VX_UR")
     vy_ur => volume%point_to_data("VY_UR")
     vz_ur => volume%point_to_data("VZ_UR")

     ! Lower right indices of beam
     vx_lr => volume%point_to_data("VX_LR")
     vy_lr => volume%point_to_data("VY_LR")
     vz_lr => volume%point_to_data("VZ_LR")

     ! Lower left indices of beam
     vx_ll => volume%point_to_data("VX_LL")
     vy_ll => volume%point_to_data("VY_LL")
     vz_ll => volume%point_to_data("VZ_LL")
  endif


  ! Fill in x,y,z locations of each gate.
  GATE_LOOP : do igate = 1, volume%ngates

      ! 3d distance of gate from aircraft
      gate_dist = scan%meters_to_center_of_first_gate + (igate-1) * scan%meters_between_gates

      ! Put gate_x, gate_y, and gate_z all in grid coordinates
      gate_x_index  = aircraft%xgrid + ( x_unit_dist * gate_dist / grid_dx )
      gate_y_index  = aircraft%ygrid + ( y_unit_dist * gate_dist / grid_dx )
      gate_z_meters = aircraft%z     + ( z_unit_dist * gate_dist )

      volume%range(igate)    = gate_dist
      vx(igate,ibeam) = gate_x_index
      vy(igate,ibeam) = gate_y_index
      vz(igate,ibeam) = gate_z_meters

      if (bw_type > 0) then

         !
         ! Determine the min/max x,y,z vertices based on the beamwidth
         !

         ! First determine the width in meters based on the radius - arc length
         !   bwidth_h_radians = volume%beamwidth_h(ibeam) * RAD_PER_DEG
         !   bwidth_v_radians = volume%beamwidth_v(ibeam) * RAD_PER_DEG
         !
         ! Now convert the spherical coordinates to a cartesian coordinate.
         ! This assumes that elevation rotates upward from the horizontal (x,y) plane and
         ! azimuth rotates positively CCW from the positive x vector. This requires recomputing
         ! azimuth = 90-azimuth to get the correct position
         !
         ! These are temporary variables that are only applicable to the range gate in question
         ! This is creating the rectangular vertices of the beam, starting lower left (beam-relative)
         ! and then clockwise around - so its LL, UL, UR, LR
         ! Note: we may want this to be circular in the future, but for now, the beam is rectangular
         !

         azm_bw_tmp(1) = (90-volume%azimuth(ibeam))-(volume%beamwidth_h(ibeam) / 2.0)
         azm_bw_tmp(2) = azm_bw_tmp(1)
         azm_bw_tmp(3) = (90-volume%azimuth(ibeam))+(volume%beamwidth_h(ibeam) / 2.0)
         azm_bw_tmp(4) = azm_bw_tmp(3)

         elev_bw_tmp(1) = (volume%elevation(ibeam))-(volume%beamwidth_v(ibeam) / 2.0)
         elev_bw_tmp(2) = (volume%elevation(ibeam))+(volume%beamwidth_v(ibeam) / 2.0)
         elev_bw_tmp(3) = elev_bw_tmp(2)
         elev_bw_tmp(4) = elev_bw_tmp(1)


         ! Let's loop over the bounds to get the x, y, z position in meters - this is relative to the radar
         BW_LOOP: do verts = 1, 4
            ! Here's the combined horizontal and vertical beamwidth
            call sph2cart(azm_bw_tmp(verts),elev_bw_tmp(verts),gate_dist,x_bw(verts),y_bw(verts),z_bw(verts))

            ! Now compute the x, y grid index as above
            gate_xind_bw(verts) = aircraft%xgrid + ( x_bw(verts) / grid_dx )
            gate_yind_bw(verts) = aircraft%ygrid + ( y_bw(verts) / grid_dx )
            gate_zmeters_bw(verts) = aircraft%z + z_bw(verts)


         enddo BW_LOOP

         ! Define the points of the vertices that make up the square around the gate

         ! Upper left
         vx_ul(igate,ibeam) = gate_xind_bw(2)
         vy_ul(igate,ibeam) = gate_yind_bw(2)
         vz_ul(igate,ibeam) = gate_zmeters_bw(2)

         ! Upper right
         vx_ur(igate,ibeam) = gate_xind_bw(3)
         vy_ur(igate,ibeam) = gate_yind_bw(3)
         vz_ur(igate,ibeam) = gate_zmeters_bw(3)

         ! Lower right
         vx_lr(igate,ibeam) = gate_xind_bw(4)
         vy_lr(igate,ibeam) = gate_yind_bw(4)
         vz_lr(igate,ibeam) = gate_zmeters_bw(4)

         ! Lower left
         vx_ll(igate,ibeam) = gate_xind_bw(1)
         vy_ll(igate,ibeam) = gate_yind_bw(1)
         vz_ll(igate,ibeam) = gate_zmeters_bw(1)
         ! End of additions by B. Klotz (11/30/2023)
      endif

  enddo GATE_LOOP

contains

  subroutine sph2cart(azm_in, elev_in, rad_in, x_out, y_out, z_out)

    !
    !  Convert from spherical to Cartesian coordinates.
    !  Assumes azimuth is 0 degrees in the positive x axis and 90 degrees in the positive y-axis
    !

    use module_configuration, only : RKIND
    use module_llxy, only : RAD_PER_DEG
    use module_llxy, only : DEG_PER_RAD
    implicit none

    real(kind=RKIND), intent(in)  :: azm_in, elev_in, rad_in
    real(kind=RKIND), intent(out) :: x_out, y_out, z_out

    x_out = rad_in * cos(elev_in * RAD_PER_DEG) * cos(azm_in * RAD_PER_DEG)
    y_out = rad_in * cos(elev_in * RAD_PER_DEG) * sin(azm_in * RAD_PER_DEG)
    z_out = rad_in * sin(elev_in * RAD_PER_DEG)

  end subroutine sph2cart

end subroutine beam_geometry

subroutine WeightFuncCRS(dxm,dym,dzm,rng,az,el,bwh,bwv,dru,wfr,wf,wgtu)
   use module_configuration, only : RKIND
   implicit none

   real(kind=RKIND), intent(in)   :: dxm, dym, dzm, rng, az, el, bwh, bwv, dru, wfr, wf
   real(kind=RKIND), intent(out)  :: wgtu
   real(kind=RKIND)               :: d_r, d_az, d_el, d_r2, d_az2, d_el2, fac, Wr, We, Wa
   real(kind=RKIND)               :: rad_out, azm_out, elev_out

   !
   ! Get the location of the grid point relative to the aircraft in spherical coordiantes.
   !

   call cart2sph(dxm, dym, dzm, rad_out, azm_out, elev_out)
   azm_out = 90-azm_out

   !
   !  Added from the CR-SIM post processing code,
   !  distances between centres of radar volume and model grid
   !

   d_r=rad_out-rng
   d_az=azm_out-az
   d_el=elev_out-el

  if (abs(d_az) > 180) d_az = 360-abs(d_az)
  !print *, "vrad,vaz,vel,rad_out,azm_out,elev_out",volume%range(igate),volume%azimuth(iray),volume%elevation(iray),tmp_rad_out,tmp_azm_out,tmp_elev_out
  !print *, "d_r, d_az, d_el",d_r,d_az,d_el

  !
  !  Squared distances between centres of radar volume and model grid
  !

  d_r2=d_r*d_r
  d_az2=d_az*d_az
  d_el2=d_el*d_el

  !
  !  Radar range, azimuth, and elevation weighting function
  !

  if ((abs(d_r)<=(dru+(0.25*dru))) .and. (abs(d_az)<=(bwh/2.0)+(0.25*bwh/2.0)) .and. (abs(d_el)<=(bwv/2.0)+(0.25*bwv/2.0))) then
      fac=d_r2/(dru*dru)
      Wr = exp(wfr*fac)
      fac=d_az2/(bwh*bwh)
      Wa = exp(wf*fac)
      fac=d_el2/(bwv*bwv) !by oue?
      We = exp(wf*fac)
  else
      Wr = 0.e0
      Wa = 0.e0
      We = 0.e0
  endif

  wgtu = Wr*Wa*We

contains
  subroutine cart2sph(x_in, y_in, z_in, rad_out, azm_out, elev_out)
    use module_configuration, only : RKIND
    use module_llxy, only : RAD_PER_DEG
    use module_llxy, only : DEG_PER_RAD
    implicit none

    real(kind=RKIND), intent(in)  :: x_in, y_in, z_in
    real(kind=RKIND), intent(out) :: rad_out, azm_out, elev_out

    rad_out  = SQRT(x_in*x_in + y_in*y_in + z_in*z_in)
    azm_out  = ATAN2(y_in,x_in) * DEG_PER_RAD
    elev_out = ATAN2(z_in,SQRT(x_in*x_in + y_in*y_in)) * DEG_PER_RAD

  end subroutine cart2sph
end subroutine WeightFuncCRS

!
!------------------------------------------------------------------------------------------------------------------
!

subroutine namelist_options ( namelist_file , seedlen , opts , scan , conf )
  use module_configuration, only : RKIND
  use module_access_wrf, only : options_type
  use module_access_wrf, only : waypoint_type
  use module_scanning,   only : scan_type
  use crsim_mod, only : conf_var
  use module_external_attitude, only : use_external_attitudes
  use module_external_attitude, only : attitude_file
  use module_external_attitude, only : attitude_orientation_rotate_degrees
#ifdef _PARALLEL_
  use mpi, only : MPI_ABORT
  use mpi, only : MPI_COMM_WORLD
#endif
  implicit none
  character(len=*),   intent(in)  :: namelist_file
  integer,            intent(out) :: seedlen
  type(options_type), intent(out) :: opts
  type(scan_type),    intent(out) :: scan
  type(conf_var),     intent(out) :: conf ! CRSim config options



#ifdef _NOPE_SURVEILLANCE_
  type(scan_type)  :: surveillance_scan
#endif

  logical             :: lexist
  character(len=1024) :: wrf_glob_pattern
  character(len=1024) :: output_filename_format_string
  character(len=1)    :: flight_level_coordinate

  real(kind=RKIND), dimension(0:999) :: flight_waypoints_x
  real(kind=RKIND), dimension(0:999) :: flight_waypoints_y
  real(kind=RKIND), dimension(0:999) :: flight_waypoints_vert

  ! heading:  The aircraft heading for the route ( degrees N, clockwise )
  real(kind=RKIND) :: heading

  ! air_speed:  The aircraft air speed for the route ( m/s )
  real(kind=RKIND) :: air_speed

  ! When to start the leg, in seconds.
  real(kind=RKIND) :: leg_initial_time

  ! How long the leg last, in seconds.
  real(kind=RKIND) :: leg_time_seconds

  ! Whether to consider time_evolution from multiple files.
  logical :: time_evolution

  ! Whether to convert minutes to seconds or not for XTIME
  integer :: conv_minute

  ! Individual panel broadside frame of reference; aircraft relative.  Valid values:
  !         0: Top
  !        90: Starboard (right) side
  !       180: Bottom
  !       270: Port (left) side
  real(kind=RKIND) :: ref_angle   

  ! Whether to to use variable or constant beamwidth - broadside beamwidth is set in CONFIG file.  Valid values:
  !    0: Infinitesimal beamwidth, i.e, no beamwidth considerations applied.
  !    1: Constant beamwidth
  !    2: Variable beamwidth
  integer :: bwtype

  ! Helicopter mode:  Aircraft is stationary with respect to the ground.
  logical :: helicopter

  !  Fixed aircraft location for a full sweep (herky_jerky=.TRUE. )
  !           or each beam gets a new aircraft location (herky_jerky=.FALSE.) (default)
  logical :: herky_jerky

  integer, dimension(256) :: seed

  namelist/options/ wrf_glob_pattern, time_evolution, output_filename_format_string,  &
       &            flight_level_coordinate,                                              &
       &            flight_waypoints_x, flight_waypoints_y, flight_waypoints_vert, heading, air_speed, &
       &            leg_initial_time, leg_time_seconds, conv_minute, &
       &            bwtype, ref_angle, &
       &            herky_jerky, helicopter, seed

  real(kind=RKIND) :: meters_between_gates
  real(kind=RKIND) :: meters_to_center_of_first_gate
  real(kind=RKIND) :: max_range_in_meters
  real(kind=RKIND) :: SNR_mask_threshold

  real(kind=RKIND) :: pulse_repetition_frequency
  integer :: pulses_per_pulse_set
  integer :: revisits_per_acquisition_time
  integer :: beams_per_acquisition_time
  real(kind=RKIND) :: pulse_repetition_time
  real(kind=RKIND) :: data_acquisition_time
  integer :: pulses_per_revisit_time
  integer :: pulses_per_acquisition_time
  real(kind=RKIND) :: skip_seconds_between_scans

  real(kind=RKIND) :: fold_limit_lower
  real(kind=RKIND) :: fold_limit_upper
  character(len=1024) :: scanning_table
  real(kind=RKIND) :: offset1
  real(kind=RKIND) :: offset2
  integer :: ibi

  character(len=1024) :: CRSIM_Config

#ifdef _PARALLEL_
  integer :: ierrmpi
#endif
  

  namelist/scanning/ meters_between_gates, meters_to_center_of_first_gate,   &
       &             max_range_in_meters, pulse_repetition_frequency,        &
       &             pulses_per_pulse_set, revisits_per_acquisition_time,    &
       &             beams_per_acquisition_time,                             &
       &             skip_seconds_between_scans,                             &
       &             SNR_mask_threshold,                                     &
       &             CRSIM_Config,                                           &
       &             scanning_table

#ifdef _NOPE_SURVEILLANCE_
  namelist/surveillance_scanning/ meters_between_gates, meters_to_center_of_first_gate,   &
       &             max_range_in_meters, seconds_for_scan_cycle,            &
       &             scan_interval_seconds,                                  &
       &             scanning_table
#endif

  namelist/attitude_external_source/ use_external_attitudes, attitude_file, attitude_orientation_rotate_degrees

  integer :: namelist_unit
  integer :: k
  integer :: waypoint_count
  integer :: ierr

  output_filename_format_string = '("test_",I6.6,".nc")'
  ! Defaul values for &options namelist
  leg_time_seconds = 600
  leg_initial_time = -1.E36
  flight_level_coordinate = "Z"
  time_evolution = .FALSE.
  conv_minute = -9999
  bwtype = 0
  ref_angle = -1.E36
  herky_jerky = .FALSE.
  helicopter = .FALSE.
  flight_waypoints_x    = -1.E36
  flight_waypoints_y    = -1.E36
  flight_waypoints_vert = -1.E36
  seed = -999999
  
  if (len_trim(namelist_file) == 0) then
      write(*,'(/," ***** ERROR *****")')
      write(*,'(" ***** User-provided namelist filename must be included as the sole command-line argument.")')
      write(*,'(" ***** ERROR EXIT",/)')
      stop
  endif

  open(newunit=namelist_unit, form='formatted', action='read', file=trim(namelist_file), iostat=ierr )
  if (ierr/=0) then
      write(*,'(80("*"))')
      write(*,'("*****   ERROR ",61x,"*****")')
      write(*,'("*****   Problem opening user-specified namelist file                       *****")')
      write(*,'("*****   File: ",A)') trim(namelist_file)
      write(*,'("*****   ERROR EXIT",57x,"*****")')
      write(*,'(80("*"))')
      stop "Problem opening user-specified namelist file"
  endif

  
  read(namelist_unit, options)
      
  ! write(*,options)

  call random_seed(size=seedlen)
  if (any(seed/=-999999)) then
      write(*,'(/,"Using provided seed array for random number generator.")')
      write(*,'("seed array may have up to", I4, " integers",/)') seedlen
      call random_seed(put=seed(1:seedlen))
  else
      write(*,'(/,"Creating seed array on the fly for random number generator.",/)')
      call init_random_seed()
  endif
  
  call random_seed(get=seed)
  write(*,'("seed array = ")', advance="no")
  write(*,'(" ", 5I15)') seed(1:seedlen)
  write(*,*)

  COUNTLOOP : do k = 0, size(flight_waypoints_x)-1
      if (flight_waypoints_x(k) < -1.E25) then
          waypoint_count = k
          exit COUNTLOOP
      endif
  enddo COUNTLOOP

  allocate(opts%waypoint(0:waypoint_count-1))
  do k = 0, waypoint_count-1
      select case (flight_level_coordinate)
      case ("Z")
          ! Z in meters
          opts%waypoint(k) = waypoint_type( flight_waypoints_x(k), flight_waypoints_y(k), flight_waypoints_vert(k), -1.E36 )
      case ("P")
          ! P in hPa (convert to Pa when we store it as a waypoint type).
          opts%waypoint(k) = waypoint_type( flight_waypoints_x(k), flight_waypoints_y(k), -1.E36, 1.E2*flight_waypoints_vert(k) )
      end select
  enddo

  opts%wrf_glob_pattern = wrf_glob_pattern
  opts%output_filename_format_string = output_filename_format_string
  opts%flight_level_coordinate = flight_level_coordinate
  opts%air_speed = air_speed
  opts%leg_initial_time = leg_initial_time
  opts%leg_time_seconds = leg_time_seconds
  opts%time_evolution = time_evolution
  opts%conv_minute = conv_minute
!KWM  opts%fpar = fpar
  opts%bwtype = bwtype
  opts%ref_angle  = ref_angle
  opts%herky_jerky = herky_jerky
  opts%helicopter = helicopter

  !
  !  Sanity checks on /options/
  !

  if (opts%leg_initial_time < -1.E25) then
      write(*,'(/,"Option ''leg_initial_time'' not set in namelist &options.",/)')
      stop
  end if

  select case (opts%flight_level_coordinate)
  case default
      write(*,'("Invalid namelist value flight_level_coordinate = ''",A,"''")') opts%flight_level_coordinate
      stop
  case ("Z","z")
      opts%flight_level_coordinate = 'Z'
  case ("P","p")
      opts%flight_level_coordinate = 'P'
  end select

  if ((opts%ref_angle == 0.0) .or. (opts%ref_angle == 90.0) .or. (opts%ref_angle == 180.0) .or. (opts%ref_angle == 270.0)) then
      ! Ok
  else if (opts%ref_angle < -1.E35) then
      write(*,'(" ***** ERROR ***** namelist value refs_angle not set.")')
      stop
  else
      write(*,'(" ***** ERROR ***** Invalid namelist value refs_angle = ''",F15.5,"''")') opts%ref_angle
      stop
  endif

  select case (opts%bwtype)
  case default
      write(*,'("Invalid namelist value bwtype = ''",I10,"''")') opts%bwtype
      stop
  case ( 0, 1, 2 )
      ! Ok
  end select

  !
  !  Default values for namelist /attitude_external_source/
  !

  use_external_attitudes = .FALSE.
  attitude_file = "nonexistant_attitude_file"
  read(namelist_unit, attitude_external_source)

  !
  ! Default values for namelist/scanning/
  !

  meters_between_gates = -1.E36
  meters_to_center_of_first_gate = -1.E36
  max_range_in_meters = -1.E36
  pulse_repetition_frequency = -1.E36
  pulses_per_pulse_set = -999999
  revisits_per_acquisition_time = -999999
  beams_per_acquisition_time = -999999
  skip_seconds_between_scans = 0.0
  scanning_table = ""
  CRSIM_Config = ""
  snr_mask_threshold = 0.0

  read(namelist_unit, scanning)

  close(namelist_unit)

  ! Check that variables are set
  if ( meters_between_gates < -1.E25           ) stop "namelist/scanning/: set METERS_BETWEEN_GATES"
  if ( meters_to_center_of_first_gate < -1.E25 ) stop "namelist/scanning/: set METERS_TO_CENTER_OF_FIRST_GATE"
  if ( max_range_in_meters < -1.E25            ) stop "namelist/scanning/: set MAX_RANGE_IN_METERS"
  if ( pulse_repetition_frequency < -1.E25     ) stop "namelist/scanning/:  set PULSE_REPETITION_FREQUENCY"
  if ( pulses_per_pulse_set < 0                ) stop "namelist/scanning/:  set PULSES_PER_PULSE_SET"
  if ( revisits_per_acquisition_time < 0       ) stop "namelist/scanning/:  set REVISITS_PER_ACQUISITION_TIME"
  if ( beams_per_acquisition_time < 0          ) stop "namelist/scanning/:  set BEAMS_PER_ACQUISITION_TIME"
  if ( scanning_table == " "                   ) stop "namelist/scanning/: set SCANNING_TABLE"
  if ( CRSIM_Config == " "                     ) stop "namelist/scanning/: set CRSIM_Config"

  !
  ! Read a CRSim config file:
  !

  write(*,'("CRSIM Configuration options from file ",A)') trim(CRSIM_Config)
  inquire(file=trim(CRSIM_Config), exist=lexist)
  if (.not. lexist) then
      write(*,'(80("*"))')
      write(*,'("*****   ERROR ",61x,"*****")')
      write(*,'("*****   Problem opening CRSIM_Config file as specified in the namelist     *****")')
      write(*,'("*****   File: ",A)') trim(CRSIM_Config)
      write(*,'("*****   ERROR EXIT",57x,"*****")')
      write(*,'(80("*"))')
      stop "Problem opening CRSIM_Config file as specified in the namelist"
  endif
  call ReadConfParameters(trim(CRSIM_Config),conf)

  scan%snr_mask_threshold             = snr_mask_threshold
  scan%meters_between_gates           = meters_between_gates
  scan%meters_to_center_of_first_gate = meters_to_center_of_first_gate
  scan%max_range_in_meters            = max_range_in_meters
  scan%pulse_repetition_frequency     = pulse_repetition_frequency
  scan%pulses_per_pulse_set           = pulses_per_pulse_set
  scan%revisits_per_acquisition_time  = revisits_per_acquisition_time
  scan%beams_per_acquisition_time     = beams_per_acquisition_time

  scan%seconds_for_scan_cycle         = -1.E36 !  Computed after we've read the Scanning Table and 
  !  know how many beams we're going to have in our sweeps.
  scan%skip_seconds_between_scans     = skip_seconds_between_scans

  write(*,'(/,"SPECIFIED:  pulse_repetition_frequency    = ", F12.6)') scan%pulse_repetition_frequency
  write(*,'("SPECIFIED:  pulses_per_pulse_set          = ", I12)'  ) scan%pulses_per_pulse_set
  write(*,'("SPECIFIED:  revisits_per_acquisition_time = ", I12)'  ) scan%revisits_per_acquisition_time
  write(*,'("SPECIFIED:  beams_per_acquisition_time    = ", I12)'  ) scan%beams_per_acquisition_time
  write(*,*)

  !
  !  And compute a few more time-related parameters
  !

  scan%pulses_per_revisit_time = scan%pulses_per_pulse_set * scan%beams_per_acquisition_time
  scan%pulses_per_acquisition_time = scan%pulses_per_revisit_time * scan%revisits_per_acquisition_time
  scan%pulse_repetition_time = 1.0 / scan%pulse_repetition_frequency
  scan%revisit_time = scan%pulses_per_revisit_time * scan%pulse_repetition_time
  scan%data_acquisition_time = scan%pulses_per_acquisition_time  * scan%pulse_repetition_time
  write(*,'("COMPUTED:   pulses_per_revisit_time       = ", I12)') scan%pulses_per_revisit_time
  write(*,'("COMPUTED:   pulses_per_acquisition_time   = ", I12)') scan%pulses_per_acquisition_time
  write(*,'("COMPUTED:   pulse_repetition_time         = ", F12.6)') scan%pulse_repetition_time
  write(*,'("COMPUTED:   revisit_time                  = ", F12.6)') scan%revisit_time
  write(*,'("COMPUTED:   data_acquisition_time         = ", F12.6)') scan%data_acquisition_time

  !
  ! Having read both the namelist and the CRSim configuration, we can compute a couple of other things.
  !
  
  scan%lambda = 2.9979E8 / (conf%freq * 1.E9)
  scan%Nyquist_Velocity = scan%lambda / (4.0*scan%pulse_repetition_time)
  scan%fold_limit_lower = -scan%Nyquist_Velocity
  scan%fold_limit_upper =  scan%Nyquist_Velocity

  write(*,'("COMPUTED:   Wavelength(lambda)            = ", F12.6)') scan%lambda
  write(*,'("COMPUTED:   Nyquist Velocity              = ", F12.6)') scan%Nyquist_Velocity

  write(*,*)

  

  call scan%read_scanning_table_file(trim(scanning_table))

  !
  !  Now that we've read the scanning table, we can compute more timing details
  !  of individual beams.
  !

  if ( mod(scan%beam_count, scan%beams_per_acquisition_time)>0) then
      ! Integer division; add 1 for a partial set of beams in acquisition time.
      scan%seconds_for_scan_cycle = (1 + (scan%beam_count / scan%beams_per_acquisition_time)) * scan%data_acquisition_time
  else
      ! Integer division
      scan%seconds_for_scan_cycle = (scan%beam_count / scan%beams_per_acquisition_time) * scan%data_acquisition_time
  endif

  allocate(scan%beam_time_from_start_of_scan_cycle(scan%beam_count))
  do k = 1, scan%beam_count

      !  Offset1 is the time of the start of a data acquisition cycle.
      !  Note the integer division to put each beam into the right data acquisition cycle.

      offset1 = ((k-1)/scan%beams_per_acquisition_time) * scan%data_acquisition_time

      !  offset2 is the time of the beam from the start of a revisit cycle.
      !  ibi is the index of the beam (1,2,3,etc) in a revisit cycle (or data acquisition cycle, for that matter)
      ibi = 1+mod(k-1,scan%beams_per_acquisition_time) 
      offset2 = ( ( (scan%revisits_per_acquisition_time-1) * scan%pulses_per_revisit_time ) + &
           &      ( ibi * scan%pulses_per_pulse_set ) ) * scan%pulse_repetition_time

      !  Time of the beam from the start of a scan cycle.
      scan%beam_time_from_start_of_scan_cycle(k) = offset1  + offset2

  enddo

  if (opts%herky_jerky) then
      ! For intermittent scans
      ! Assume instantaneous scans, so the only time increment is going to be the time to skip between scans.
      scan%seconds_plus_skip = scan%skip_seconds_between_scans
  else
      scan%seconds_plus_skip = scan%seconds_for_scan_cycle + scan%skip_seconds_between_scans
  endif


  write(*,*)
  write(*,'(" sweep_count = ", I6)') scan%sweep_count
  write(*,'(" beam_count  = ", I6)') scan%beam_count
  write(*,'(" seconds_for_scan_cycle  = ", F12.8)') scan%seconds_for_scan_cycle
  write(*,'(" skip_seconds_between_scan  = ", F12.8)') scan%skip_seconds_between_scans
  write(*,'(" seconds_plus_skip  = ", F12.8)') scan%seconds_plus_skip
  do k = 1, scan%sweep_count
      write(*,'(6x,"Sweep ",I6, " :: beams ", I6, " through ", I6, " :: timing ", F12.6, " through ", F12.6, " seconds")') &
           k, scan%sweep_start_index(k), scan%sweep_end_index(k), &
           scan%beam_time_from_start_of_scan_cycle(scan%sweep_start_index(k)+1), &
           scan%beam_time_from_start_of_scan_cycle(scan%sweep_end_index(k)+1)
  enddo
  write(*,*)

end subroutine namelist_options

subroutine init_random_seed()
  !
  ! From gcc.gnu.org example of using random_seed.
  !
#ifdef __INTEL_COMPILER
  use ifport, only : getpid
#endif
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t
#ifdef __NVCOMPILER_LLVM__
  integer, external :: getpid
#endif

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
      read(un) seed
      close(un)
  else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.
      call system_clock(t)
      if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
      end if
      pid = getpid()
      print*, 'pid = ', pid
      t = ieor(t, int(pid, kind(t)))
      do i = 1, n
          seed(i) = lcg(t)
      end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
        s = 104729
    else
        s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed
