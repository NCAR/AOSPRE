&options
 seed = 875270836, 1299229797
 wrf_glob_pattern = "./wrfout_test_????.nc"

 output_filename_format_string = '("LHS_rhi_",A,"_to_",A,".nc")'

 leg_initial_time = 4651
 leg_time_seconds = 8
 time_evolution = .TRUE.

  flight_waypoints_x    =       1,     10
  flight_waypoints_y    =     105,    105
  flight_waypoints_vert =    2000,   2000

 flight_level_coordinate   = "Z"

 air_speed = 120. ! meters per second
 herky_jerky = .FALSE.
 namelist_iq = "../iq/namelist.iq"

 bwtype = 0,
 ref_angle = 270,

/

&attitude_external_source
  use_external_attitudes = .False.,
  attitude_file = "attitude.nc",
  attitude_orientation_rotate_degrees = 182.0,
/

&scanning

 CRSIM_Config = "CONFIG_crsim"

 scanning_table = "scanning_lhs.txt"
 scanning_table_iq = "../iq/Scanning_IQ_LHS.txt"
 pulse_repetition_frequency     = 2500
 pulses_per_pulse_set           = 2
 revisits_per_acquisition_time  = 20
 beams_per_acquisition_time     = 6

skip_seconds_between_scans = 0.0

 meters_between_gates           = 150.
 meters_to_center_of_first_gate = 150.
 max_range_in_meters            = 75000.
 
 ! meters_between_gates           = 75.
 !meters_to_center_of_first_gate = 1.
 !max_range_in_meters            = 75000.
 
/

&config_output
/
