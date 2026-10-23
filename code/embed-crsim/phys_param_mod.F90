  !! ----------------------------------------------------------------------------
  !! ----------------------------------------------------------------------------
  !!  *MODULE* phys_param_mod
  !!  @version crsim 3.x
  !!
  !!  *SRC_FILE*
  !!  crsim/src/phys_param_mod.f90
  !!
  !!
  !!  *LAST CHANGES*
  !!
  !! Aug 20, 2019       AT      Conversion real to double precission for two constants 
  !!                            Addition of integer value of zero
  !!
  !!  *DESCRIPTION* 
  !!
  !!  This module contains the physical parameters and constants
  !!  used throughout the CR-SIM. 
  !!
  !!  ==========================
  !!  CRSIM for APAR BSD LICENSE
  !!  ==========================
  !!  
  !!  Copyright (c) 2024, Stony Brook University - Brookhaven National 
  !!  Laboratory - McGill University Radar Science Group 
  !!  
  !!  Redistribution and use in source and binary forms, with or without
  !!  modification, are permitted provided that the following conditions are
  !!  met:
  !!  
  !!  1) If the software is modified to produce derivative works, such
  !!  modified software should be clearly marked, so as not to confuse it
  !!  with the version available from the Radar Science Group.
  !!  
  !!  2) Redistributions of source code must retain the above copyright
  !!  notice, this list of conditions and the following disclaimer.
  !!  
  !!  3) Redistributions in binary form must reproduce the above copyright
  !!  notice, this list of conditions and the following disclaimer in the
  !!  documentation and/or other materials provided with the distribution.
  !!  
  !!  4) Neither the name of the Stony Brook University - Brookhaven National
  !!  Laboratory - McGill University Radar Science Group nor the names of its 
  !!  contributors, if any, may be used to endorse or promote products derived 
  !!  from this software without specific prior written permission.
  !!  
  !!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  !!  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  !!  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  !!  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  !!  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  !!  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  !!  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  !!  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  !!  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  !!  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  !!  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  !!
  !!
  !!----------------------------------------------------------------------------------------------------------------- 
  !!-----------------------------------------------------------------------------------------------------------------
  !!
  module phys_param_mod
  Implicit None


  Integer, parameter :: maxNumOfHydroTypes = 8  ! the maximal number of Hydrometeor Types
                                              !  Allowed (currrently 8 for RAMS, MP=40)


  real*8,Parameter   :: pi=dacos(-1.d0) ! value of pi
  real*8,Parameter   :: r2d=180.d0/pi  ! radians to degrees conversion factor
  real*8,Parameter   :: d2r=pi/180.d0  ! degree to radians conversion factor
  
  real*8,Parameter   :: K2=0.92d0 ! reflective index factor (for water)
  
  Real*8,Parameter   :: Rearth = 6371.d0  !km
  Real*8,Parameter   :: c=2.9979d0*1.d+8 !!! Speed of light m/s
  
  
  
  
  
  Real*8,Parameter   :: Rd=287.058d0   !  J kg-1 K-1  the individual gas constant for dry air   
  Real*8,Parameter   :: eps=0.622d0 ! eps=Rd/Rv=mv/mv           
  Real*8,Parameter   :: p0=1000.d0*1.d+2  ! Pa  pressure level of 1000 mb
  Real*8,Parameter   :: p0_850=850.d0*1.d+2  ! Pa  pressure level of 850 mb
  Real*8,Parameter   :: p0_600=600.d0*1.d+2  ! Pa  pressure level of 600 mb
  Real*8,Parameter   :: cp=1004.d0  ! J kg-1 K-1
  Real*8,Parameter   :: grav=9.80665d0 ! m s-2conventional standard value forgravitational acceleration 
  
  
  
  real*8,parameter                 :: rhow=1000.d0 ! kg m-3 density of water
  real*8,parameter                 :: rhow_25C=997.d0 ! kg m-3 density of water at 25 deg C
  real*8,parameter                 :: T0K =273.15d0 ! K Temperature of 0 deg in K
  
  
  
  real*8,parameter                 :: oneOverThree=1.d0/3.d0
  real*8,parameter                 :: twoOverThree=2.d0/3.d0
  
  
  ! OTHER PARAMETERS
  Integer,Parameter     :: N_layers=10 ! maximum number of cloud layers used for ARSCL
  
  
  
  ! SOME CONSTANTS USED IN THE CODE
  Real*8,Parameter                   :: m999=-999.d0
  Integer,Parameter                  :: izero = 0
  Real*8,Parameter                   :: zero=0.d0
  Real*8,Parameter                   :: Zthr=1.d-10 ! mm^6/m^3  or  -100 dBZ ! threshold value of reflectivity;
                                                    ! computed reflectivity set to zero if larger than Zthr
  
  
  !  Default setting close to ARM radars for Doppler spectra simulation 
  
  real*8, parameter  :: TimeSampling     = 2.d0                 ! integration time in secs
  real*8, parameter  :: PRF_S            = 600.d0               ! Pulse Repetition Frequency Hz for S band
  real*8, parameter  :: PRF_C            = 1240.d0              ! Pulse Repetition Frequency Hz for C band
  real*8, parameter  :: PRF_X            = 2000.d0              ! Pulse Repetition Frequency Hz for X band
  real*8, parameter  :: PRF_Ka           = 2700.d0              ! Pulse Repetition Frequency Hz for Ka band
  real*8, parameter  :: PRF_Ku           = 2000.d0              ! Pulse Repetition Frequency Hz for Ku, but use it for X temporarily
  real*8, parameter  :: PRF_W            = 7500.d0              ! Pulse Repetition Frequency Hz for W band
  real*8, parameter  :: NOISE_1km        = 1.d-4                ! Linear units of noise power in dBZ at 1-km
  integer, parameter :: NFFT             = 256                  ! number of FFT points
  
  
  end module phys_param_mod
