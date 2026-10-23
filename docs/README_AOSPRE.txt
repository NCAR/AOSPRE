%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
% README for the NCAR EOL APAR Observation Simulation, Processing,    %
% and Research Environment (AOSPRE)                                   %
%                                                                     %
% The purpose of this document is to describe the process of running  %
% AOSPRE system. This process includes installation of all            %
% appropriate libraries needed to run the code, designing             %
% flight tracks, ensuring you have the correct model output, and      %
% visualizing the AOSPRE output.                                      %
%                                                                     %
% AOSPRE creators: Kevin Manning and Scott Ellis                      %
% AOSPRE contributors: Brad Klotz and George Bryan                    %
% Scientific Lead: Wen-Chau Lee                                       %
%                                                                     %
% Document version: 1.0.1                                             %
% Date:             04/02/2021                                        %
% Author:           Brad Klotz                                        %
% Contact:          bradklotz@ucar.edu                                %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Version History %%%%

04/02/2021 - 1.0: Document initiated
09/25/2024 - 1.0.1: Updated in preparation for repository release (B. Klotz)



%%%% README Contents %%%%

1) AOSPRE Purpose
2) Pre-Install checklist
3) Installation instructions
4) Running the AOSPRE
5) Visualizing AOSPRE output                             


%%%% 1) AOSPRE Purpose %%%%

The intent of the AOSPRE system is to provide a glimpse of the APAR
observations in the context of certain weather phenomena. It is 
promoted as a risk reduction tool and for investigative purposes to
determine limits of the APAR and appropriate scan strategies prior
to installation and operation on the NSF/NCAR C-130. It is up to
the user how to operate the aircraft in this simulated environment
while providing some restraint for known operational guidelines.


%%%% 2) Pre-install Checklist %%%%

The main AOSPRE code is setup in a Github repository with specific instructions
for downloading, building, testing, and usage. The build and operation of the AOSPRE
set to run on Linux/Unix systems. The current version has been successfully tested 
on is CentOS Stream release 8, 4.18.0-277.el8.x86_64, AlmaLinux 8.10 (Cerulean Leopard), 
and MacOS Ventura, Version 13.6.9.

The code requires NetCDF libraries - make sure you are correctly linked
to those libraries. I make sure this is always correct by adding the
following line to my .bashrc file:

export LD_LIBRARY_PATH={Enter path to NetCDF lib/ directory}

The AOSPRE runs in a bash shell using Fortran as its main programming
language. Ensure that you have the GNU Fortran compiler installed.
For Linux, this should be installed automatically, but its always good
to check. If you have admin privileges, you can check by typing:

yum install gcc-gfortran

The third requirement is a dependency on a third-party software called
the Cloud Resolving Model Radar Simulator (CR-SIM). There are look-up
tables (LUTS) in the software that are required for AOSPRE. To download, go
to: https://you.stonybrook.edu/radar/research/radar-simulators/. Once you
have installed CR-SIM, you can link the "aux" directory to your AOSPRE code
directory. This will be discussed in a little more detail later in this
document.

The last pre-requisite to run the AOSPRE system is WRF or CM1 model output.
The functions are expecting certain variables, so a list of these are
provided in a separate document called: WRF_VARIABLES.txt. It is 
effectively an "ncdump" of the variables in one example output file
used in the AOSPRE testing. If your model output has at least all of 
these variables, it should be able to run in the AOSPRE system. 

An optional component of this software is to run with parallel
computing resources. It is currently set up to use OpenMPI. If
the user chooses to install with this capability, the correct
libraries and software should be properly located. This will be 
discussed in the installation instructions.

To recap:
- Set your .bashrc file to explicitly define the NetCDF lib path
- Set your .bashrc to set unlimited for memory constraints
- Install GNU Fortran compiler and NetCDF, if not already done
- Install CR-SIM software
- Have WRF or CM1 model output ready to use
- (Optional) Have OpenMPI installed with proper Fortran compilers


%%%% 3) Installation instructions (with Makefile) %%%%

There are several ways to install or build the code. The main way to do this
would be to clone the repository by entering the command:

git clone https://github.com/NCAR/APAR-Observing-Simulator

A second option would be dowloading the tar'ed and zipped file to
the user's system. Unzip and untar this file in the desired
operating location.

The structure of the repository is such that all of the source Fortran
code is stored in the directory is "code/embed-crsim". It contains all the functions
to run the main portion of the code and generate output files. The first expectation
is that the user links the lookup tables from CR-SIM for access by the 
AOSPRE main code. To do this from the main AOSPRE directory, locate the "aux"
directory from the CR-SIM installation: 

crsim/src/crsim-{version}/share/crsim/aux

Then link it to the current AOSPRE main directory:

ln -s crsim/src/crsim-{version}/share/crsim/aux

Now the user should go to the "code/embed-crsim" directory. The current 
expectation for installation is to use the "Makefile" within this directory.
The user will have to determine their system compilers and uncomment the 
selected lines for the appropriate situation. It is also important to
declare the location of the NetCDF libraries. The variable "NETCDF" should
be explicitly declared with that location so the software has access
to its "/lib" directory. After ensuring the setup in the Makefile is
correct type:

make clean
make all

Barring any compilation errors, this will clean any old versions and 
reinstall all the scripts with the correct library locations. The executable
is called aospre. You are now able to run the AOSPRE code!


%%%% 4) Running the AOSPRE %%%%

Before you can run the AOSPRE, there are four things you must do. These
include verifying the CONFIG file, determining your flight path,
updating your namelist file, and defining a scan file. 

a) The CONFIG file comes from CR-SIM, and the only variable that might
need changing is the section that says "Specify radar frequency".
For APAR purposes, this is set to 5.5d0 for C-band, but the other
frequencies listed are also viable.

b) The next item to set up or check is your scan file. These files will
be listed as "Scanning_Table_..." and they define your scan strategy.
There are examples of these scan files in the "Scanning_Table_Library" directory
found in the main directory. The scan tables expect columns of rotation and tilt 
angles and the specified axis of rotation (primary axis).

c) The third item is the namelist file. Currently, there are two namelist
files used - one for PPI and one for RHI scans. Several lines will need
to be adjusted to make it work for your particular flight.

   i)    'wrf_glob_pattern' - this should reflect the location and naming
			style of your model output files
   ii)   'output_filename_...' - the start of the name can be changed, but
                          the A,"_to_",A,".nc" portion should stay the same
   iii)   'leg_initial_time' - this is the start time of your flight in 
                          seconds
   iv)  'leg_time_seconds' - this is the total length of the flight in 
                          seconds
   v)   'time_evolution' - this should always be kept at .TRUE.

   vi)  'flight_waypoints...' - these describe the points of the flight 
                          path, which is discussed below.
   vii) 'air_speed'      - this is the speed of the aircraft and may
			vary depending on your specific plane
   viii) 'herky_jerky'    - should be kept at .FALSE.
   ix)   ' bwtype'    - this is defining whether to run without consideration
			of beamwidth (0), constant beamwidth (1), and variable
                        beamwidth (2)
   x)    ' ref_angle' - defines the reference angle for the beamwidth, where
                        270 = LHS, 90 = RHS, 0 = TOP, and 180 = BOT
   xi)   'use_external_attitudes' - if you have actual attitude information
                           then you can set this to true and provide the
                           file name.
   xii)    'CRSIM_Config'   - list the name of the CONFIG file you used
   xiii)   'scanning_table' - list the name of the scanning table you used
   
   The remaining variables are specific to APAR and may need to be adjusted
   to your specific radar.

d) The final item to take care of is determining the waypoints that get
included in the namelist file. A bash script is provided to determine
those values. It is called "flight_planner.sh" and has a configuration
file called "flight_config". Instructions are provided in that script for
each of the entries. The main thing is to make sure you determine a good
initial starting location and time. The waypoints are defined in terms
of model grid index (x,y). To get the waypoints, type ./flight_planner.sh,
and copy the output to your namelist file.

e) Once you have these files set up, you can run the AOSPRE. The typical
procedure is create a directory for running tests associated with a 
particular weather simulation in the main directory, such as "run_tests/storm1".
Make sure all your necessaryfiles are in that directory, including CONFIG, 
namelist, and ScanningTable. All output will go in that directory initially
but can be moved to a more appropriate location. The expectation is that the user
could use this storm directory as a based for running different flights
for the same simulation.

AOSPRE output are in NetCDF format. A file called flightpath.ascii is also
output as the aircraft flight tracker during the simulated flight.

To invoke the code, you should be in your output directory, and then
type: ./../code/embed-crsim/aospre <namelist>

If running with parallel processing, type:
mpirun -n N ../code/embed-crsim/aospre <namelist>


%%%% Visualizing the AOSPRE output %%%%

There are options to visualize the data with 3rd party tools,  including
ncview, LROSE, or other tools. An existing set of Matlab-based scripts
can be used to create plots, but they are not currently part of this code
release. If the user would like access to these scripts, a request can
be made through the Github repository.



 
