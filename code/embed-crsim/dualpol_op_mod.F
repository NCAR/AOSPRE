module dualpol_op_mod
!--------------------------
! Parameterized polarimetric radar operator following
! Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators 
!   for simulation and assimilation of polarimetric radar data with 
!   numerical weather predictions. Adv. Atmos. Sci., 38(5), 737−754.
!
! Author: Zhiquan (Jake) Liu, NCAR/MMM
! Initial version: July 2022
!--------------------------------------
  implicit none


  ! Setting up some various contants

   real (kind=8), parameter :: density_rain    = 1.0, & ! g/cm^3=1000kg/m^3
                              density_snow    = 0.1, & ! dry density, g/cm^3
                              density_graupel = 0.5, & ! dry density, g/cm^3
                              density_hail    = 0.917  ! dry density, g/cm^3

   real (kind=8), parameter :: density_air     = 1.225  ! kg/m^3=0.001 g/cm^3
   real (kind=8), parameter :: dpi             = 3.1415926

   real (kind=8), dimension(0:12,0:3) :: sband_snow_coefs, sband_graupel_coefs, sband_hail_coefs
   real (kind=8), dimension(0:12,0:3) :: cband_snow_coefs, cband_graupel_coefs, cband_hail_coefs
   real (kind=8), dimension(1:4,0:4)  :: sband_rain_coefs, cband_rain_coefs

!  private

 contains

  ! ------------------------------------------------------------------------------
  subroutine read_coefs_icephase(filename, coefs)
  !------------------------------------------------------------------------
  ! Read coef C as defined in Table 1 (snow), 2 (hail) and 3 (graupel) 
  ! of Zhang et al., 2021, 13 x 4 array
  !------------------------------------------------------------------------
    implicit none
    character (len=*), intent(in)      :: filename
    real (kind=8),     intent(inout)   :: coefs(0:12,0:3)   ! coefs

    integer :: i

    open(unit=15,file=trim(filename),status='old',form='formatted',action='read')
    do i=0,12
       read(15,*) coefs(i,:)
       !write(*,'(4f12.8)') coefs(i,:)
    enddo
    close(15)

  end subroutine read_coefs_icephase

  ! ------------------------------------------------------------------------------
  subroutine read_coefs_rain(filename, coefs)
  !------------------------------------------------------------------------
  ! Read coef  in Eq. (13) - (16) of Zhang et al., 2021
  !  4 x 5 array
  !------------------------------------------------------------------------
    implicit none
    character (len=*), intent(in)      :: filename
    real (kind=8),     intent(inout)   :: coefs(1:4,0:4)   ! coefs

    integer :: i

    open(unit=16,file=trim(filename),status='old',form='formatted',action='read')
    do i=1,4
       read(16,*) coefs(i,:)
       !write(*,'(4f12.8)') coefs(i,:)
    enddo
    close(16)

  end subroutine read_coefs_rain

  ! ------------------------------------------------------------------------------
  subroutine init_coefs
  !------------------------------------------------------------------------
  ! Read coef of rain/snow/graupel/hail of Zhang et al., 2021
  !------------------------------------------------------------------------
    implicit none

    !character (len=*) :: filedir_sbr !,filedir_sbs,filedir_sbg,filedir_sbh
    !character (len=*) :: filedir_cbr,filedir_cbs,filedir_cbg,filedir_cbh
    
    !filedir_sbr='/export/wind1/bradklotz/AOS/code/embed-crsim_jakeliu_dualpol/sband_rain_coefs.txt'

    !real (kind=8), dimension(0:12,0:3) :: sband_snow_coefs, sband_graupel_coefs, sband_hail_coefs
    !real (kind=8), dimension(0:12,0:3) :: cband_snow_coefs, cband_graupel_coefs, cband_hail_coefs
    !real (kind=8), dimension(1:4,0:4)  :: sband_rain_coefs, cband_rain_coefs
    

    call read_coefs_rain('../code/embed-crsim/coefs/sband_rain_coefs.txt',sband_rain_coefs)
        !print *,"Sband rain coefs: ",sband_rain_coefs
        !call sleep(90)
    !call read_coefs_rain(filedir_sbr,sband_rain_coefs)
    call read_coefs_icephase('../code/embed-crsim/coefs/sband_snow_coefs.txt',sband_snow_coefs)
    call read_coefs_icephase('../code/embed-crsim/coefs/sband_graupel_coefs.txt',sband_graupel_coefs)
    call read_coefs_icephase('../code/embed-crsim/coefs/sband_hail_coefs.txt',sband_hail_coefs)

    call read_coefs_rain('../code/embed-crsim/coefs/cband_rain_coefs.txt',cband_rain_coefs)
    call read_coefs_icephase('../code/embed-crsim/coefs/cband_snow_coefs.txt',cband_snow_coefs)
    call read_coefs_icephase('../code/embed-crsim/coefs/cband_graupel_coefs.txt',cband_graupel_coefs)
    call read_coefs_icephase('../code/embed-crsim/coefs/cband_hail_coefs.txt',cband_hail_coefs)

  end subroutine init_coefs

  ! ------------------------------------------------------------------------------
  subroutine dualpol_op_rain(w, dm, a, zh, zdr, kdp, phv)
  !--------------------------------------------------------
  ! Eqs. (13) - (16) of Zhang et al., 2021
  !------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: w   ! density_air * mixing ratio , g/m^3
    real(kind=8), intent(in)   :: dm  ! mass/volume-weighted mean diameter, mm
    real(kind=8), intent(in)   :: a(1:4,0:4) ! coefs for rain
    real(kind=8), intent(out)  :: zh  ! horizontal reflectivity, linear units
    real(kind=8), intent(out)  :: zdr ! differential reflectivity, linear units
    real(kind=8), intent(out)  :: kdp ! specific differential phase, degree/km^-1
    real(kind=8), intent(out)  :: phv ! co-polar correlation coefficient, 0-1

    zh  = w*(a(1,0) + a(1,1)*dm + a(1,2)*dm**2 + a(1,3)*dm**3 + a(1,4)*dm**4)**2
    !print *,'a(1,0)= ',a(1,0),' a(1,1)= ',a(1,1),' a(1,2)= ',a(1,2)
    !print *,'zhr_tmp= ',zh
    !print *,'w= ',w
    !print *,'dm= ',dm

    zdr =    a(2,0) + a(2,1)*dm + a(2,2)*dm**2 + a(2,3)*dm**3 + a(2,4)*dm**4

    kdp = w*(a(3,0) + a(3,1)*dm + a(3,2)*dm**2 + a(3,3)*dm**3 + a(3,4)*dm**4)

    phv =    a(4,0) + a(4,1)*dm + a(4,2)*dm**2 + a(4,3)*dm**3 + a(4,4)*dm**4
    !print *,'Zdr= ',zdr,'Kdp= ',kdp,'Phv= ',phv

  end subroutine dualpol_op_rain

  ! ------------------------------------------------------------------------------
  subroutine dualpol_op_icephase(z, w, rho, dm, a, zh, zdr, kdp, phv)
  !-----------------------------------------------------------------
  ! Eqs. (17) - (20) of Zhang et al., 2021 for snow/graupel/hail
  !---------------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: z   ! 6th moment of PSD, Eq. (6)
    real(kind=8), intent(in)   :: w   ! density_air * mixing ratio , g/m^3
    real(kind=8), intent(in)   :: rho ! density of melting snow/graupel/hail, g/cm^3
    real(kind=8), intent(in)   :: dm  ! mass/volume-weighted mean diameter, mm
    real(kind=8), intent(in)   :: a(0:12) ! coef

    real(kind=8), intent(out)  :: zh  ! horizontal reflectivity, linear units
    real(kind=8), intent(out)  :: zdr ! differential reflectivity, linear units
    real(kind=8), intent(out)  :: kdp ! specific differential phase, degree/km^-1
    real(kind=8), intent(out)  :: phv ! co-polar correlation coefficient, 0-1

    zh  = z*(a(0) + a(1)*dm + a(2)*dm**2 + a(3)*dm**3)**2

    zdr =    a(4) + a(5)*dm + a(6)*dm**2

    kdp = w*(a(7) + a(8)*dm + a(9)*dm**2)/rho

    phv =   a(10) + a(11)*dm + a(12)*dm**2

  end subroutine dualpol_op_icephase

  ! --------------------------------------------------------------
  subroutine dualpol_op_total ( zh_rain,  zh_snow,  zh_graupel, &
                               zdr_rain, zdr_snow, zdr_graupel, &
                               kdp_rain, kdp_snow, kdp_graupel, &
                               phv_rain, phv_snow, phv_graupel, &
                               zh, zdr, kdp, phv)
  !---------------------------------------------------------------
  ! compute total zh, zdr, kdp, phv
  ! Eqs. 22 - 25 of Zhang et al., 2021
  !---------------------------------------------------------------

    implicit none

    real(kind=8), intent(in)  :: zh_rain, zh_snow, zh_graupel    ! horizontal reflectivity, dBZ
    real(kind=8), intent(in)  :: zdr_rain, zdr_snow, zdr_graupel ! differential reflectivity, dB
    real(kind=8), intent(in)  :: kdp_rain, kdp_snow, kdp_graupel ! specific differential phase, degree/km^-1
    real(kind=8), intent(in)  :: phv_rain, phv_snow, phv_graupel ! co-polar correlation coefficient, 0-1

    real(kind=8), intent(out)  :: zh  ! horizontal reflectivity, linear units
    real(kind=8), intent(out)  :: zdr ! differential reflectivity, linear units
    real(kind=8), intent(out)  :: kdp ! specific differential phase, degree/km^-1
    real(kind=8), intent(out)  :: phv ! co-polar correlation coefficient, 0-1

    real(kind=8) :: zhzdr_rain, zhzdr_snow, zhzdr_graupel
    real(kind=8) :: dzhrainzdrrain, dzhsnowzdrsnow, dzhgraupelzdrgraupel

    zh  = zh_rain + zh_snow + zh_graupel ! Eq. 22

    if (zdr_rain > 0.) then
       dzhrainzdrrain = zh_rain/zdr_rain
       zhzdr_rain    = zh_rain    / sqrt(zdr_rain)
    else
       dzhrainzdrrain =0.
       zhzdr_rain =0.
    end if
    if (zdr_snow > 0.) then
       dzhsnowzdrsnow = zh_snow/zdr_snow
       zhzdr_snow = zh_snow    / sqrt(zdr_snow)
    else
       dzhsnowzdrsnow = 0.
       zhzdr_snow = 0.
    end if
    if (zdr_graupel > 0.) then
       dzhgraupelzdrgraupel = zh_graupel/zdr_graupel
       zhzdr_graupel = zh_graupel / sqrt(zdr_graupel)
    else
       dzhgraupelzdrgraupel = 0.
       zhzdr_graupel =0.
    end if

    if (dzhrainzdrrain+dzhsnowzdrsnow+dzhgraupelzdrgraupel > 0.) then
       zdr = zh/(zh_rain/zdr_rain + zh_snow/zdr_snow + zh_graupel/zdr_graupel) ! Eq. 23
    else
       zdr =0.
    end if

    kdp = kdp_rain + kdp_snow + kdp_graupel ! Eq. 24

    if (zhzdr_rain + zhzdr_snow + zhzdr_graupel /= 0.) then
       phv = (zhzdr_rain*phv_rain + zhzdr_snow*phv_snow + zhzdr_graupel*phv_graupel) &
           /(zhzdr_rain + zhzdr_snow + zhzdr_graupel) ! Eq. 25
    else
       phv =0.
    end if
  end subroutine dualpol_op_total

  ! ------------------------------------------------------------------------------
  subroutine coef_a(c, gamma, a)
  !---------------------------------------
  ! Eq. (21) of Zhang et al., 2021
  !---------------------------------------
    implicit none
    real(kind=8), intent(in)   :: c(0:3)  ! coef c in Table
    real(kind=8), intent(in)   :: gamma ! coef for Zdr
    real(kind=8), intent(out)  :: a  ! horizontal reflectivity, dBZ

    integer  :: i

    a = c(0)
    do i = 1, 3
      a = a + c(i) * gamma**i
    end do

  end subroutine coef_a

  !-------------------------------------------------------------------
  subroutine rainfrac(qr, qx, ratio)
  !------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: qr    ! mixing ratio of rain, kg/kg
    real(kind=8), intent(in)   :: qx    ! mixing ratio of hydrometeor, kg/kg
    real(kind=8), intent(out)  :: ratio ! qr/(qr+qx)

    ratio = qr/(qr+qx)

  end subroutine rainfrac

  !-------------------------------------------------------------------
  subroutine weticephase_density(precip_type, ratio, rho)
  !------------------------------------------------------
  ! Eq. (7) of Zhang et al., 2021
  !--------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type
    real(kind=8), intent(in)   :: ratio  ! qr/(qr+qx)
    real(kind=8), intent(out)  :: rho    ! density of melting precip, g/cm^3

    select case (precip_type)
       case ('snow')
          rho =  density_snow * (1.0 - ratio**2) + ratio**2
       case ('graupel')
          rho =  density_graupel * (1.0 - ratio**2) + ratio**2
       case ('hail')
          rho =  density_hail * (1.0 - ratio**2) + ratio**2
       case default
          write(*,*) 'Unknow precipitation type: ', precip_type
          stop
    end select

  end subroutine weticephase_density

  !-------------------------------------------------------------------
  subroutine watercontent(air_density, qx, w)
  !------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: qx            ! mixing ratio of precip, g/kg
    real(kind=8), intent(out)  :: w             ! water content, g/m^3

    w = air_density * qx

  end subroutine watercontent

  !------------------------------------------------------------------------------
  subroutine dm_z_2moment(air_density, qx, rhox, ntx, dm, z)
  !-----------------------------------------------------------
  ! Eqs. (5) and (6) of Zhang et al., 2021
  ! Assume 2-moment MP scheme with ntx and qx as prognostic variables
  ! Exponential PSD: N(D) = N0 * exp(-lambda D)
  !-------------------------------------------------------------------
    implicit none
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: rhox          ! density of wet hydrometeor,g/cm^3 
                                                ! rain density = 1 g/cm^3
    real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
    real(kind=8), intent(in)   :: ntx           ! number concentration, ?/m^3
    real(kind=8), intent(out)  :: dm            ! mean diameter, mm
    real(kind=8), intent(out)  :: z             ! 6th moment of PSD

    if (qx > 0.0) then
       dm = 4.0 * ( (air_density*qx)/(dpi*rhox*ntx) )**(1.0/3.0)
        z = 11250.0 * ( (air_density*qx)/(dpi*rhox) )* dm**3
        print *,' Z_lin= ',z
    else
       dm = 0.0
        z = 0.0
    endif
    
  end subroutine dm_z_2moment

  !------------------------------------------------------------------------------
  subroutine dm_z_wsm6(precip_type, air_density, temp, qx, rhox, dm, z)
  !------------------------------------------------------------------------
  ! Single moment scheme assuming Exponential PSD: N(D) = N0 * exp(-lambda D)
  ! N0 is fixed except for snow, then calcuate lambda from N0
  !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type  ! rain/snow/graupel
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: rhox          ! density of wet hydrometeor,g/cm^3 
                                                ! rain density = 1 g/cm^3
    real(kind=8), intent(in)   :: temp          ! temperature, Kelvin
    real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
    real(kind=8), intent(out)  :: dm            ! mean diameter, mm
    real(kind=8), intent(out)  :: z             ! 6th moment of PSD

    real(kind=8)  :: n0  ! intercept parameter of exponential PSD
    real(kind=8)  :: lambda ! slope parameter of exponential PSD

    select case (precip_type)
       case ('rain')
           n0 = 8.0e6 ! m^(-4)
       case ('snow')
           n0 = 2.0e6*exp(0.12*(273.15-temp))
       case ('graupel')
           n0 = 4.0e6 ! m^(-4)
       case default
          write(*,*) 'Unknow precipitation type: ', precip_type
          stop
    end select

      !rhox = rhox * 1000.0 ! g/cm^3 = 10^-3 kg/ (10-2 m)^3 = 1000 kg/m^3
      !  qx = qx * 0.001   ! g/kg -> kg/kg
    lambda = ((dpi*1000.0*rhox*n0)/(air_density*qx*0.001))**0.25 ! m^-1, Year2-report
        dm = 4000.0 / lambda  ! Eq. (5), should be in mm, *1000?
         z = 11250.0 * ( (air_density*qx)/(dpi*rhox) )* dm**3 ! Eq. (6)

  end subroutine dm_z_wsm6

  !------------------------------------------------------------------------------
  subroutine dm_z_gceop(precip_type, air_density, qx, rhox, dm, z)
  !------------------------------------------------------------------------
  ! Single moment scheme assuming Exponential PSD: N(D) = N0 * exp(-lambda D)
  ! N0 is fixed, calcuate lambda from N0
  ! Goddard Cumulus Ensemble scheme used in CWB WRF operation
  ! only difference from wsm6: n0 for snow is a fixed value
  !-------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type  ! rain/snow/graupel
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: rhox          ! density of wet hydrometeor,g/cm^3 
                                                ! rain density = 1 g/cm^3
    real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
    real(kind=8), intent(out)  :: dm            ! mean diameter, mm
    real(kind=8), intent(out)  :: z             ! 6th moment of PSD

    real(kind=8)  :: n0  ! intercept parameter of exponential PSD
    real(kind=8)  :: lambda ! slope parameter of exponential PSD

    select case (precip_type)
       case ('rain')
           n0 = 8.0e6 ! m^(-4)
       case ('snow')
           n0 = 1.6e7
       case ('graupel')
           n0 = 4.0e6 ! m^(-4)
       case default
          write(*,*) 'Unknow precipitation type: ', precip_type
          stop
    end select

      !rhox = rhox * 1000.0 ! 1 g/cm^3 = 10^-3 kg/ (10-2 m)^3 = 1000 kg/m^3
      !  qx = qx * 0.001    ! 1 g/kg -> 0.001 * kg/kg
    lambda = ((dpi*1000.0*rhox*n0)/(air_density*qx*0.001))**0.25 ! m^-1, Year2-report
        dm = 4000.0 / lambda  ! Eq. (5), should be in mm, *1000?
         z = 11250.0 * ( (air_density*qx)/(dpi*rhox) )* dm**3 ! Eq. (6)

  end subroutine dm_z_gceop

  !------------------------------------------------------------------------------
  subroutine dm_z_gcecwb(precip_type, air_density, temp, qx, frac, dm, z, rhox)
  !-------------------------------------------------------------------------------
  ! Single moment scheme assuming Gamma PSD: N(D) = N0 * D^alpha * exp(-lambda D)
  ! Goddard Cumulus Ensemble scheme modified by CWB WRF
  ! hydrometeor particle density and mean diameter are diagnosed.
  !
  ! T.-C. Tsai, S.-Y. Jiang, J.-P. Chen, H.-L. Huang, Yi.-C. Lo, J.-S. Hong, 2022:
  !  A Gamma-type Single Moment Bulk Microphysics Scheme for Operational Forecasting.
  !----------------------------------------------------------------------------------
    implicit none
    character (len=*), intent(in) :: precip_type  ! rain/snow/graupel
    real(kind=8), intent(in)   :: air_density   ! Air Density, 1.225 kg/m^3 
    real(kind=8), intent(in)   :: temp          ! temperature, Kelvin
    real(kind=8), intent(in)   :: qx            ! mixing ratio of hydrometeor, g/kg
    real(kind=8), intent(in)   :: frac          ! melting fraction=qr/(qr+qx)
    real(kind=8), intent(out)  :: dm            ! mean diameter, mm
    real(kind=8), intent(out)  :: z             ! 6th moment of PSD
    real(kind=8), intent(out)  :: rhox          ! density of melting snow/graupel

    real(kind=8)  :: logw  ! ln(water content), kg/m^3
    real(kind=8)  :: logt  ! ln(temperature)

    
    logw = log(air_density * (0.001*qx)) ! kg/m^3
    logt = log(temp)  ! temp in Kelvin

  ! Diagnose density
    select case (precip_type)
       case ('rain')
          rhox = 1000.0 ! kg/m^3
            dm = exp(-20.189916 + 6.888453*logt + 0.25*logw  &
                         - 0.3132*logt**2 + 1.43e-18*logw**2  &
                         - 5.27932e-16*logt*logw) ! in micron meter
       case ('snow')
          rhox = exp(65.84696 - 17.124*logt + 0.84975*logw  &
                         + 1.069132*logt**2 - 0.00614*logw**2 &
                         - 0.208331*logt*logw) ! kg/m^3
            dm = exp(-87.705062 + 26.74956*logt - 1.327406*logw  &
                         - 1.6701*logt**2 + 0.009598*logw**2  &
                         + 0.32543632*logt*logw) ! in micron meter
       case ('graupel')
          rhox = 1000.0*frac + &  ! kg/m^3
                 (1-frac)*max(250.0,exp(0.056433*logw+6.608048))
            dm = exp(-74.804847 + 16.08519*logt - 2.530305*logw  &
                         - 0.19094*logt**2 + 0.005846*logw**2  &
                         + 0.51917156*logt*logw) ! in micron meter
       case default
          write(*,*) 'Unknow precipitation type: ', precip_type
          stop
    end select

    ! air_density(kg/m^3), qx(g/kg), rhox(g/cm^3), dm(mm)
    rhox = rhox * 0.001 ! kg/m^3 -> g/cm^3
      dm = dm * 0.001   ! um -> mm
       z = 11250.0 * ( (air_density*qx)/(dpi*rhox) ) * dm**3 ! Eq. (6)

  end subroutine dm_z_gcecwb

  ! ------------------------------------------------------------------------------

end module dualpol_op_mod
