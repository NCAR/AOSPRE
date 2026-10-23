project(aospre)
cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)

option(PARALLEL "Enable parallel AOSPRE" OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}") 

add_executable(${PROJECT_NAME}
    crsim_luts_mod.F90
    crsim_mod.F90
    extract_apar_wrf.F90
    kwm_date_utilities.F90
    module_access_wrf.F90
    module_aircraft.F90
    module_cfradial_output.F90
    module_configuration.F90
    module_crsim_wrapper.F90
    module_external_attitude.F90
    module_geometry.F90
    module_llxy.F90
    module_quicksort.F90
    module_scanning.F90
    module_timing_utilities.F90
    phys_param_mod.F90
    ReadConfParameters.F90
    wrf_var_mod.F90
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    Fortran_FORMAT FREE
)

target_compile_options(${PROJECT_NAME} PRIVATE $<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none>)

# MacOS specific flags to prevent an "Illegal instruction" error from occurring.
# TODO: Investigate why this is necessary.
set_property(SOURCE "extract_apar_wrf.F90" PROPERTY COMPILE_FLAGS $<$<PLATFORM_ID:Darwin>:"-Og">)
set_property(SOURCE "extract_apar_wrf.F90" PROPERTY COMPILE_FLAGS $<$<PLATFORM_ID:Darwin>:"-ftree-sra">)

find_package(NetCDF REQUIRED COMPONENTS Fortran)
target_link_libraries(${PROJECT_NAME} PRIVATE NetCDF::NetCDF_Fortran NetCDF::NetCDF_C)

message(STATUS "Parallel AOSPRE: ${PARALLEL}")
if(PARALLEL)
    if(NOT MPI_HOME AND NOT ENV{MPI_HOME})
        find_path(MPI_HOME "bin/mpiexec" PATH_SUFFIXES mpich openmpi PATHS /usr/lib /usr/lib64)
    endif()
    find_package(MPI
        REQUIRED
        COMPONENTS
            Fortran
    )

    # gfortran 10+ is more strict on argument mismatches which causes MPI versions earlier than 4 to throw an error
    # an argument mismatch. Add the '-fallow-argument-mismatch' flag to turn that into a warning.
    set(is_gnu "$<Fortran_COMPILER_ID:GNU>")
    set(fortran_v10_or_later "$<VERSION_GREATER_EQUAL:$<Fortran_COMPILER_VERSION>,10>")
    set(mpi_earlier_than_v4 "$<VERSION_LESS:${MPI_Fortran_VERSION_MAJOR},4>")
    set(meets_requirements "$<AND:${is_gnu},${fortran_v10_or_later},${mpi_earlier_than_v4}>")
    set_property(SOURCE "extract_apar_wrf.F90" PROPERTY COMPILE_FLAGS $<${meets_requirements}:-fallow-argument-mismatch>)

    target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_Fortran)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _PARALLEL_)
endif()