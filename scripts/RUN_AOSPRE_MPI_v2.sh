#!/bin/bash

# RUN_AOSPRE_MPI.sh
#
# Brad Klotz, NCAR Project Scientist
# Created: 06/16/2021
# Modified: 08/07/2025 - for testing with the parallel computing mode

# This program  will initiate the APAR Observing Simulation (AOS)
# system code. You will need to specify a few items in order to complete the
# code requirements. You will need to configure your namelist
# file prior to running this script for now. That functionality
# may change in the future.

# The script will read information that you sepcify at the top of
# this script, which include: a base directory, operating
# directory, scan type, namelist file, experiment_directory,
# test directory, and output directory.

# Recommended script execution - creates a log file and runs in the background
# ./RUN_AOS.sh fpar ftime > out_aos.log &
#
# Even though the <fpar> and <ftime> variables are listed elsewhere, they need to be included
# in the command as well. So <fpar> indicates whether you want to run parallel or not and <ftime>
# is the flight time in seconds.
#
# fpar: Enter 1 for parallel or 0 for serial
# ftime: Enter the length of the flight in seconds
#
# Other entries of interest
#
# Base Directory: This should be your main AOS directory you've created
#                 example: /export/wind1/bradklotz/AOSPRE
#
# Operating Directory: This is the directory containing your namelist and config files
#                 example: ${basedir}/code/WRF_supercell_100m
#
# Scan Mode: This is simply the type of scan - surveillance or rhi
#
# Namelist file: This is the file that controls the functionality of the AOS code
#                 example: namelist.surveillance
#
# Experiment Directory: If you are performing several runs as part of a larger
#                       experiment, enter the experiment directory here. Otherwise leave
#                       as an empty string.
#                       example: experiments/test_01
#
# Test Directory: This is the name of the test-specific directory you
#                 want to create. It will typically have information about the
#		  particular AOS run you are executing.
#		  example: supercell_100m_10s_Cband_1kmAGL
#
# Output Directory: This is the combined strings of the Base Directory, the
#                   Experiment Directory, and the Test Directory.
#No user input is needed to create this directory.
#

##### USER-DEFINED VARIABLES #####

# Base directory
base_dir="/export/wind1/bradklotz/AOSPRE_projection"

# Operating directory
operating_dir=${base_dir}/WRF_AR_20211024

# Scan mode: rhi or surveillance
scan_mode="rhi"
#scan_mode="surveillance"

# Namelist file
#namelist_file=namelist_0.5_AK
#namelist_file=full_ppi_namelist_reduced_-126to-122_37to42.txt
namelist_file=rhs_rhi_namelist_reduced_-126to-122_37to42.txt
#namelist_file=namelist_rhi_supercell_100m_0.1kmAGL_WestEast_x215_y655_singleleg_10minute

# Experiment directory - if none, leave as empty string
experiment_dir="individual_tests"
#experiment_dir="experiments/reduced_panel_tests"

# Experiment resolution (100m or 500m)
exp_res=600

# Test directory
# This will typically provide some specific information
# about your run of the AOS - such as the storm type, resolution
# aircraft initial height, microwave band, etc.
test_dir="atmospheric_river_600m_10s_2.0kmAGL_Cband_2.0BW_x380_y220_43211-43871s_FVEL_beamwidth2.0_const_ver2.0"

#### No user input needed beyond this point ####

# Output directory
# First check that the experiment directory exists (if it is not an
# empty string)

if [ ! -z "$experiment_dir" ]; then
   
   if [ ! -d ${base_dir}/${experiment_dir} ]; then
      mkdir -p ${base_dir}/${experiment_dir}
   fi

   output_dir=${base_dir}/${experiment_dir}/${test_dir}
else
   output_dir=${base_dir}/${test_dir}
fi

# Get your current working directory
current_dir=$(pwd)

echo "Running the AOS code..."
echo "Your current working directory is: " ${current_dir}
echo "Your operating directory is: " ${operating_dir}
echo "Your scan mode is: " ${scan_mode}
echo "Your namelist file is: " ${namelist_file}
echo "Your output directory is: " ${output_dir}

##### END USER VARIABLES #####

##### EXECUTION OF AOS CODE #####

# Verifying the base and operating directories  exist
if [[ ! -d ${base_dir} || -z "$base_dir" ]]; then
	echo "Your base directory is not correct...please fix and re-run."
	echo "Exiting..."
	exit
fi

if [[ ! -d ${operating_dir} || -z "$operating_dir" ]]; then
	echo "Your operating directory is not correct...please fix and re-run."
	echo "Exiting..."
        exit
fi

# Create the output directory if it doesn't exist
echo ""
echo "All specified directories are verified...Proceeding with AOS execution..."

if [ ! -d ${output_dir} ]; then
	   mkdir ${output_dir}
fi

# Go to the operating directory
echo ""
echo "The current working directory is now the operating directory: " ${operating_dir}
cd ${operating_dir}

echo ""
echo "Command line entries: (par1, flighttime)   $1 $2"

# Get the starting date and time of the AOS run
start_time=`date +%s.%N`

# Pull some information out of the namelist file?
# We need to know whether we will need to operate in parallel and the time designated
# for the flight

#fp_line=`grep -nw fpar ${namelist_file}`
#echo ${fp_line}
numprocs=`nproc`
nphalf=$(( $numprocs / 2 ))

if [ ${exp_res} -eq 500 ]; then
    max_procs=19;
else
    max_procs=7;
fi

#echo "Do you need to run the AOS in parallel? (1=yes,0=no): "
par1=$1

if [ $par1 -eq 0 ]; then
   npuse=1;
   runver=1
   echo ""
   echo "You chose to run in serial mode, so only 1 processor will be used..."
   echo ""

else
   echo ""
   echo "You chose to run in parallel mode...determining the suggested number of processors..."
   #echo "Enter the length of your flight in seconds: "
   flighttime=$2
   runver=0

   echo ""
   echo "The current breakdown time window is set to 30 seconds for each processor..."
   echo "You have ${numprocs} CPUs on this machine...but the recommendation is to use no more than half of them..."
   echo "This means you should use no more than ${nphalf} processors..."
   echo "The maximum allowed for this resolution is ${max_procs}...this will be the number of processors allowed..."

   # Compute the expected number of processors needed based on time
   nwin=30   # seconds for each flight window
   npuse=$(( $flighttime / $nwin ))
   ii=1   # number of 30 second time windows


   while [ $npuse -gt ${max_procs} ];
   do
   # if [ $npuse -gt $nphalf ]; then
      ii=$(( ${ii}+1 ))
      nwin=$(( $nwin * $ii ))
      echo "The maximum number of processors was exceeded...using the ${nwin} second time window..."
      #nwin=$(( $nwin * 2 ))
      #npuse=$(( $flighttime / $nwin ))
      
      #if [ $nphalf -gt ${max_procs} ]; then
      #   nwin=$(( $flighttime / ${max_procs} ))
      #   echo "The maximum number of processors was exceeded...using ${nwin} time window..."
         npuse=${max_procs}
      #fi
         
   #fi
   done
   
   echo "Given your flight time specified, you should use ${npuse} processors with $nwin second windows"
   echo ""
   echo "Moving forward with ${npuse} processors and $nwin second windows"
   
fi
echo ""
#echo "Sleeping..."
#sleep 20


# Check to make sure mpirun exists, if it doesn't, then you have to run using the traditional commands
# Assign the check variable to 0, if the command isn't found, set it to 1
#runver=0
if ! command -v mpirun &> /dev/null
then
   echo "mpirun does not exist on this computer...running AOS in the traditional manner..."
   echo ""
   runver=1
fi 

# Execute the code
echo ""
echo "Initiating AOS code..."
echo ""
echo "And off we go..."
echo ""
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "%%        AOS SPECIFIC OUTPUT BEGINS        %%"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

if [ ${runver} -eq 1 ]; then
#   ./../embed-crsim/a.out ${namelist_file}
   ./../code/embed-crsim/aospre ./namelist_files/${namelist_file}
else
#   mpirun --mca bil openid,self,vader -np ${npuse} -report-child-jobs-separately -tag-output ../embed-crsim_openmpi/atest.out ${namelist_file}
   mpirun  -np 1  ./../code/embed-crsim/aospre ./namelist_files/${namelist_file} 
fi

echo ""
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "%%         AOS SPECIFIC OUTPUT ENDS         %%"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

# Get the end date and time of the AOS run
end_time=`date +%s.%N`

runtime=$( echo "$end_time - $start_time" | bc -l )

echo ""
echo "AOSPRE processing time = ${runtime} seconds"

# Create a new flightpath file from the multiple thread file version
#numfp=`ls -1 flightpath_* | wc -l`
#numfp2=$(expr ${numfp} - 1)
#echo $numfp2

#for (( i=0 ; i<=$numfp2 ; i++ ));
#do
#    echo $i
#    if [ $i -eq 0 ]; then
#       if [ -f "flightpath.ascii" ]; then
#          rm flightpath.ascii
#       fi
#       cat flightpath_${i}.ascii > flightpath.ascii
#    else
#       cat flightpath_${i}.ascii >> flightpath.ascii
#    fi
    
#    rm flightpath_${i}.ascii
#done

# Once complete, move created files to output directory
echo ""
echo "Processing complete...moving files to output directory..."
if [ ${scan_mode} == 'rhi' ]; then

   # Creating a directory within your output for the particular scan mode
   if [ ! -d ${output_dir}/${scan_mode} ]; then
       mkdir ${output_dir}/${scan_mode}
   fi

   mv RHI.*.nc ${output_dir}/${scan_mode}
   cp 'flightpath.ascii' ${output_dir}/${scan_mode}
   cp ./namelist_files/${namelist_file} ${output_dir}/${scan_mode}
else

   # Creating a directory within your output for the particular scan mode
   if [ ! -d ${output_dir}/${scan_mode} ]; then
       mkdir ${output_dir}/${scan_mode}
   fi

   mv surveillance_*.nc ${output_dir}/${scan_mode}
   cp 'flightpath.ascii' ${output_dir}/${scan_mode}
   cp ./namelist_files/${namelist_file} ${output_dir}/${scan_mode}
fi

# Return to the original directory when you started
echo ""
echo "Returning to your original working directory: " ${current_dir}
cd ${current_dir}

echo ""
echo "AOSPRE code execution is complete...Done!"


