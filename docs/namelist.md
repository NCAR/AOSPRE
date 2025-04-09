#**Namelist file**

Fortran namelist to control aspects of the AOSPRE simulation.  Organized in sections [&options](#options), [&scanning](#scanning), and [&config_output](#config_output), each representing a Fortran *namelist record*.

## &options

### **seed**
*integer array*

* integer seeds to initialize the Fortran random number generator.  The number of integers needed to initialize the random number generator is compiler-dependent.
* optional, defaults to system-dependent random initialization.
```
seed = 198276, 98120419, 166, -16189763, 66389692, 1470800236
```

---------------------------------------------------------------------------------------------------------------------------------------------------

### **wrf_glob_pattern**
*string*

* pattern to find model output files, according to conventional POSIX globbing (wildcard) rules.
```
 wrf_glob_pattern   = "/data/wrf/supercell/100m/wrfout_d01_????-??-??_??:??:??"
```

---------------------------------------------------------------------------------------------------------------------------------------------------

### **output_filename_format_string**
*string*

* Fortran format string to build output file names.
```
output_filename_format_string = '("supercell_",A,"_",A,".nc")'
```

---------------------------------------------------------------------------------------------------------------------------------------------------

### **leg_initial_time**
*integer*

* time (seconds past model initialization time) to start the flight.
* The model initialization time is the time of the model start, and not necessarily the time of the first model output. In WRF, model start time can be found using `ncdump -h /path/to/wrfout.nc | grep :SIMULATION_START_DATE`.
```
 leg_initial_time = 4500
```
---------------------------------------------------------------------------------------------------------------------------------------------------

### **leg_time_seconds**
*integer*

* duration (in seconds) of the entire flight.
```
 leg_time_seconds = 300
```
---------------------------------------------------------------------------------------------------------------------------------------------------

### **flight_waypoints_x**
*float array*

* x-coordinates (model grid units) of the waypoints of the flight.
```
flight_waypoints_x    =    150.5,  350.0,  350.0,
```
---------------------------------------------------------------------------------------------------------------------------------------------------

### **flight_waypoints_y**
*float array*

* y-coordinates (model grid units) of the waypoints of the flight.
```
flight_waypoints_y    =    150.75,  150.25,  350.0,
```
---------------------------------------------------------------------------------------------------------------------------------------------------

### **flight_waypoints_vert**
*float array*

* vertical coordinates (in units of the flight_level_coordinate) of the waypoints of the flight.
```
 flight_waypoints_vert =  1500.0,  1500.0,  1500.0,
```

---------------------------------------------------------------------------------------------------------------------------------------------------

### **flight_level_coordinate**
*single character*

* Vertical coordinate used for flight waypoints.  Valid values are <ul><ul> - 'Z' (height in meters)</ul> <ul> - 'P' (Pressure in hPa).</ul></ul>
```
flight_level_coordinate   = "Z"
```
---------------------------------------------------------------------------------------------------------------------------------------------------

### **air_speed**
*float*

* constant air speed (m/s) for the flight.
```
 air_speed = 120.0
```

---------------------------------------------------------------------------------------------------------------------------------------------------

### **herky_jerky**
*logical*

* explain this.
```
 herky_jerky = .FALSE.
```
---------------------------------------------------------------------------------------------------------------------------------------------------
<!--
## &attitude_external_source
*Needs work; probably not functioning well*
---------------------------------------------------------------------------------------------------------------------------------------------------
-->

## &scanning

### **crsim_config**
*string*

* file name for CR-SIM configuration
```
crsim_config = "CONFIG"
```
---------------------------------------------------------------------------------------------------------------------------------------------------
### **scanning_table**
*string* 

* file name for the table that defines the scanning strategy
```
scanning_table = "Scanning_Table_LHS"
```
---------------------------------------------------------------------------------------------------------------------------------------------------
### **pulse_repetition_frequency**
*integer*

* describe this
```
pulse_repetition_frequency = 2000
```
---------------------------------------------------------------------------------------------------------------------------------------------------
### **pulses_per_pulse_set**
*integer*

* describe this
```
pulses_per_pulse_set = 2
```
---------------------------------------------------------------------------------------------------------------------------------------------------
### **revisits_per_acquisition_time**
*integer*

* describe this
```
revisits_per_acquisition_time = 7
```
---------------------------------------------------------------------------------------------------------------------------------------------------
### **beams_per_acquisition_time**
*integer*

* describe this
```
beams_per_acquisition_time     = 3
```
---------------------------------------------------------------------------------------------------------------------------------------------------
### **seconds_for_scan_cycle**
*float*

* describe this
```
 ! seconds_for_scan_cycle         = 2.0
```

---------------------------------------------------------------------------------------------------------------------------------------------------

### **meters_between_gates**
*float*

* describe this
```
 meters_between_gates           = 150.
```
---------------------------------------------------------------------------------------------------------------------------------------------------
### **meters_to_center_of_first_gate**
*float*

* describe this
```
 meters_to_center_of_first_gate = 150.
```
---------------------------------------------------------------------------------------------------------------------------------------------------
### **max_range_in_meters**
*float*

* describe this
```
 max_range_in_meters            = 75000.
```
---------------------------------------------------------------------------------------------------------------------------------------------------
### **fold_limit_lower**
*float*

* describe this
* optional (*is this used anymore or can we always determine the limits from the pulse configuration?*)
```
fold_limit_lower = -16.905
```
---------------------------------------------------------------------------------------------------------------------------------------------------
### **fold_limit_upper**
*float*

* describe this
* optional (*is this used anymore or can we always determine the limits from the pulse configuration?*)
```
fold_limit_upper =  16.905
```
---------------------------------------------------------------------------------------------------------------------------------------------------
 


## &config_output
