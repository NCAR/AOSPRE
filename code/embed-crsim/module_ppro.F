module module_ppro
  use module_configuration, only : RKIND
  use module_cfradial_output, only : volume_type
  use phys_param_mod, only: d2r
  use dualpol_op_mod
  
  implicit none

contains

  subroutine ppro_wrapper(volume)
   use dualpol_op_mod

  implicit none

  class (volume_type) :: volume
  integer :: igate
  integer :: iray

  integer       :: i
  real (kind=RKIND) :: air_density ! kg/m^3
  real (kind=RKIND) :: temp ! temperature, Kelvin
  real (kind=RKIND) :: qrain, qsnow, qgraupel ! mixing ratio, kg/kg
  real (kind=RKIND) :: rho_rain, rho_snow, rho_graupel ! g/cm^3
  real (kind=RKIND) :: dm ! mean diameter, mm
  real (kind=RKIND) :: w   ! g/m^3
  real (kind=RKIND) :: ratio   ! water fraction in melting snow/graupel, 0-1
  real (kind=RKIND) :: rho ! density, g/cm^3
  real (kind=RKIND) :: z, a(0:12)
  real (kind=RKIND) :: zh_rain, zdr_rain, kdp_rain, phv_rain
  real (kind=RKIND) :: zh_snow, zdr_snow, kdp_snow, phv_snow
  real (kind=RKIND) :: zh_graupel, zdr_graupel, kdp_graupel, phv_graupel
  real (kind=RKIND) :: zhzdr_rain, zhzdr_snow, zhzdr_graupel
  real (kind=RKIND) :: zh, zdr, kdp, phv

  ! Added by B. Klotz (10/08/2025)
  ! Variables in CR-SIM but are not computed in P-PRO
  real (kind=RKIND) :: fvel, Ah, Av, Adp, Zvv, Zvh, LDRh 
  real (kind=RKIND) :: Dopp, DPhase, SWtot, SWh, SWs, SWt, SWv
  ! End addition (10/08/2025)

  !real(kind=8), dimension(1:4,0:4)  :: sband_rain_coefs, cband_rain_coefs
  !real(kind=8), dimension(0:12,0:3) :: sband_snow_coefs, cband_snow_coefs, sband_hail_coefs, cband_hail_coefs
  !real(kind=8), dimension(0:12,0:3) :: sband_graupel_coefs, cband_graupel_coefs

  real(kind=RKIND), pointer, dimension(:,:) :: rhodptr, tptr
  real(kind=RKIND), pointer, dimension(:,:) :: qrptr, qsptr, qgptr
  real(kind=RKIND), pointer, dimension(:,:) :: qnrptr, qnsptr, qngptr

  real(kind=RKIND), pointer, dimension(:,:) :: zh_ptr
  real(kind=RKIND), pointer, dimension(:,:) :: phv_ptr
  real(kind=RKIND), pointer, dimension(:,:) :: zdr_ptr
  real(kind=RKIND), pointer, dimension(:,:) :: kdp_ptr

  real(kind=RKIND), pointer, dimension(:,:) :: zhrain_ptr
  real(kind=RKIND), pointer, dimension(:,:) :: zhsnow_ptr
  real(kind=RKIND), pointer, dimension(:,:) :: zhgraup_ptr

  ! Added by B. Klotz (10/08/2025)
  ! Variables included in CR-SIM but not computed in P-PRO
  real(kind=RKIND), pointer, dimension(:,:) :: fvelptr
  real(kind=RKIND), pointer, dimension(:,:) :: ahptr
  real(kind=RKIND), pointer, dimension(:,:) :: avptr
  real(kind=RKIND), pointer, dimension(:,:) :: adpptr
  real(kind=RKIND), pointer, dimension(:,:) :: zvvptr
  real(kind=RKIND), pointer, dimension(:,:) :: zvhptr
  real(kind=RKIND), pointer, dimension(:,:) :: ldrhptr
  real(kind=RKIND), pointer, dimension(:,:) :: Doppptr
  real(kind=RKIND), pointer, dimension(:,:) :: DPhaseptr
  !real(kind=RKIND), pointer, dimension(:,:) :: PHIdp_DPh_ptr
  real(kind=RKIND), pointer, dimension(:,:) :: swtotptr
  real(kind=RKIND), pointer, dimension(:,:) :: SWsptr
  real(kind=RKIND), pointer, dimension(:,:) :: swhptr
  real(kind=RKIND), pointer, dimension(:,:) :: swtptr
  real(kind=RKIND), pointer, dimension(:,:) :: swvptr
  real(kind=RKIND), pointer, dimension(:,:) :: zhcloud_ptr
  real(kind=RKIND), pointer, dimension(:,:) :: zhice_ptr

  ! End Addition of CR-SIM variables

  real(kind=RKIND) :: rho_d, qsum, qruse, qsuse, qguse, qnruse, qnsuse, qnguse
  real(kind=RKIND) :: Zthr = 1.d-10     ! in dBZ not linear units

  air_density = density_air ! use 1.225 kg/m^3 fixed value for this example
  temp = 260.0 ! Kelvin

  call init_coefs
  print *,"Cband rain coefs: ",cband_rain_coefs
  !call sleep(5)

    tptr     => volume%point_to_data("TEMPERATURE")
    rhodptr  => volume%point_to_data("RHO_D")
    qrptr    => volume%point_to_data("QRAIN")
    qsptr    => volume%point_to_data("QSNOW")
    qgptr    => volume%point_to_data("QGRAUP")
    qnrptr   => volume%point_to_data("QNRAIN")
    qnsptr   => volume%point_to_data("QNSNOW")
    qngptr   => volume%point_to_data("QNGRAUPEL")

    zh_ptr   => volume%new_array("Zhh")
    zhrain_ptr   => volume%new_array("Zhh_rain")
    zhsnow_ptr   => volume%new_array("Zhh_snow")
    zhgraup_ptr  => volume%new_array("Zhh_graupel")

    phv_ptr  => volume%new_array("RHOhv")
    zdr_ptr  => volume%new_array("Zdr")
    kdp_ptr  => volume%new_array("Kdp")

    ! Added by B. Klotz (10/08/2025)
    fvelptr => volume%new_array("fvel")
    ahptr   => volume%new_array("Ah")
    avptr   => volume%new_array("Av")
    adpptr  => volume%new_array("Adp")
    zvvptr  => volume%new_array("Zvv")
    zvhptr  => volume%new_array("Zvh")
    ldrhptr => volume%new_array("LDRh")
    Doppptr => volume%new_array("Dopp")
    DPhaseptr => volume%new_array("DPhase")
    !PHIdp_DPh_ptr   => volume%new_array("PHIdp_DPh")
    swtotptr => volume%new_array("SWtot")
    SWsptr   => volume%new_array("SWs")
    swhptr   => volume%new_array("SWh")
    swtptr   => volume%new_array("SWt")
    swvptr   => volume%new_array("SWv")
    zhcloud_ptr   => volume%new_array("Zhh_cloud")
    zhice_ptr   => volume%new_array("Zhh_ice")

    ! End addition of CR-SIM variables

    zh_ptr        = -9999.
    zhrain_ptr    = -9999.
    zhsnow_ptr    = -9999.
    zhgraup_ptr   = -9999.

    phv_ptr   = -9999.
    zdr_ptr   = -9999.
    kdp_ptr   = -9999.

    ! Added by B. Klotz (10/08/2025)
    fvelptr  = -9999.
    ahptr    = -9999.
    avptr   = -9999.
    adpptr  = -9999.
    zvvptr  = -9999.
    zvhptr  = -9999.
    ldrhptr = -9999.
    Doppptr = -9999.
    DPhaseptr = -9999.
    !PHIdp_DPh_ptr   = -9999.
    swtotptr = -9999.
    SWsptr   = -9999.
    swhptr   = -9999.
    swtptr   = -9999.
    swvptr   = -9999.
    zhcloud_ptr   = -9999.
    zhice_ptr   = -9999.

    ! End CR-SIM variables

    qruse = 0.0
    qsuse = 0.0
    qguse = 0.0
    qnruse = 0.0
    qnsuse = 0.0
    qnguse = 0.0
    
    print *, "NRAYS = ", volume%nrays
    do iray = 1, volume%nrays
        IGATE_LOOP : do igate = 1, volume%ngates
            !print *, "IRAY = ",iray," IGATE = ",igate

            qsum = 0.0
            if ( qrptr(igate,iray) > 1.d-10 ) then
               qsum = qsum + qrptr(igate,iray)
               qruse = qrptr(igate,iray)
               qnruse = qnrptr(igate,iray)
            else
               qruse = 0.0
               qnruse = 0.0
            endif
            if ( qsptr(igate,iray) > 1.d-15 ) then
               qsum = qsum + qsptr(igate,iray)
               qsuse = qsptr(igate,iray)
               qnsuse = qnsptr(igate,iray)
            else
               qsuse = 0.0
               qnsuse = 0.0
            endif
            if ( qrptr(igate,iray) > 1.d-15 ) then
               qsum = qsum + qgptr(igate,iray)
               qguse = qgptr(igate,iray)
               qnguse = qngptr(igate,iray)
            else
               qguse = 0.0
               qnguse = 0.0
            endif

            if (qsum <= 0.0) then
               cycle IGATE_LOOP
            else
               rho_d = rhodptr(igate,iray)

            ! Convert to the units expected in scheme_2moment
               print *, 'Gate= ',igate,' Ray= ',iray
               print *, 'QR= ',qruse,' QS= ',qsuse,' QG= ',qguse
               print *, 'QNR= ',qnruse,' QNS= ',qnsuse,' QNG= ',qnguse

              qruse = 1000.0 * qruse  ! to convert to g/kg
              qsuse = 1000.0 * qsuse  ! to convert to g/kg
              qguse = 1000.0 * qguse  ! to convert to g/kg

              qnruse = qnruse / 1000.0 ! convert #/kg to #/m^3; 1000 kg/1 m^3
              qnsuse = qnsuse / 1000.0
              qnguse = qnguse / 1000.0

              call scheme_2moment(2, air_density, rho_d, qruse, qsuse, qguse, &
                   qnruse, qnsuse, qnguse, &
                   zh, zdr, kdp, phv, &
                   zh_rain, zdr_rain, kdp_rain, phv_rain, &
                   zh_snow, zdr_snow, kdp_snow, phv_snow, &
                   zh_graupel, zdr_graupel, kdp_graupel, phv_graupel )

           ! sband_rain_coefs(1:4,0:4), sband_snow_coefs(0:12,0:3), sband_graupel_coefs(0:12,0:3), &
           ! cband_rain_coefs(1:4,0:4), cband_snow_coefs(0:12,0:3), cband_graupel_coefs(0:12,0:3), &

           ! print *,' Zh_AM = ',zh
            ! print *, 'Zh_AM_rain = ',zh_rain, ' Zh_AM_snow= ', zh_snow, ' Zh_graup= ',zh_graupel

             
             if (zh > Zthr) then
                zh_ptr(igate,iray)    = zh
                phv_ptr(igate,iray)   = phv
                kdp_ptr(igate,iray)   = kdp
                zdr_ptr(igate,iray)   = zdr
                
                !if (zdr_ptr(igate,iray) < 0) then
                !   zdr_ptr(igate,iray)   = -9999.9
                !endif
                
             else
                zh_ptr(igate,iray)    = -9999.
                phv_ptr(igate,iray)   = -9999.
                zdr_ptr(igate,iray)   = -9999.
                kdp_ptr(igate,iray)   = -9999.
             endif

            if (zh_rain > Zthr) then
                zhrain_ptr(igate,iray)    = zh_rain
            else
                zhrain_ptr(igate,iray)    = -9999.
            endif

            if (zh_snow > Zthr) then
               zhsnow_ptr(igate,iray)    = zh_snow
            else
               zhsnow_ptr(igate,iray)    = -9999.
            endif

            if (zh_graupel > Zthr) then
               zhgraup_ptr(igate,iray)    = zh_graupel
            else
               zhgraup_ptr(igate,iray)    = -9999.
            endif
          endif

          ! Added by B. Klotz (10/08/2025)
          ! CR-SIM variables not computed in P-PRO
          fvelptr(igate,iray) = -9999.
          ahptr(igate,iray) = -9999.
          avptr(igate,iray) = -9999.
          adpptr(igate,iray) = -9999.
          zvvptr(igate,iray) = -9999.
          zvhptr(igate,iray) = -9999.
          ldrhptr(igate,iray) = -9999.
          Doppptr(igate,iray) = -9999.
          DPhaseptr(igate,iray) = -9999.
          !PHIdp_DPh_ptr(igate,iray) = -9999.
          swtotptr(igate,iray) = -9999.
          SWsptr(igate,iray) = -9999.
          swhptr(igate,iray) = -9999.
          swtptr(igate,iray) = -9999.
          swvptr(igate,iray) = -9999.
          zhcloud_ptr(igate,iray) = -9999.
          zhice_ptr(igate,iray) = -9999.

          ! End addition of CR-SIM variables
        enddo IGATE_LOOP
        !call sleep(1)
    enddo

    where (zh_ptr >= Zthr)
       zh_ptr = 10*log10(zh_ptr) ! zh is provided in linear units; same for hydrometeor classes
       zdr_ptr = 10*log10(zdr_ptr)  ! zdr is provided in linear units
       phv_ptr = phv_ptr
       kdp_ptr = kdp_ptr
    elsewhere
       zh_ptr = -9999.
       zdr_ptr = -9999.
       phv_ptr = -9999.
       kdp_ptr = -9999.
    endwhere

    where (zhrain_ptr >= Zthr)
       zhrain_ptr = 10*log10(zhrain_ptr)
    elsewhere
       zhrain_ptr = -9999.
    endwhere

    where (zhsnow_ptr >= Zthr)
       zhsnow_ptr = 10*log10(zhsnow_ptr)
    elsewhere
       zhsnow_ptr = -9999.
    endwhere

    where (zhgraup_ptr >= Zthr)
       zhgraup_ptr = 10*log10(zhgraup_ptr)
    elsewhere
       zhgraup_ptr = -9999.
    endwhere

  end subroutine ppro_wrapper


!---------------------------------------------------
! example code to compute polarimetric radar variables
! input from Morrison2Moment scheme: air_density,qrain,qsnow,qgraupel,
!                                    ntx_rain,ntx_snow,ntx_graupel
!---------------------------------------------------
  subroutine scheme_2moment(iband, air_density, rhod, qrain, qsnow, qgraupel, &
                     ntx_rain, ntx_snow, ntx_graupel, &
                     zh, zdr, kdp, phv, &
                     zh_rain, zdr_rain, kdp_rain, phv_rain, &
                     zh_snow, zdr_snow, kdp_snow, phv_snow, &
                     zh_graupel, zdr_graupel, kdp_graupel, phv_graupel )

    implicit none
    integer       :: i,iband
    real (kind=8), intent(in) :: air_density, rhod ! kg/m^3
    ! real (kind=8), intent(in) :: density_rain, density_snow, density_graupel ! g/cm^3 = 1000 km/m^3
    real (kind=8), intent(in) :: qrain, qsnow, qgraupel ! mixing ratio, kg/kg
    real (kind=8), intent(in) :: ntx_rain, ntx_snow, ntx_graupel  ! number concentration, ?/m^3
    !real (kind=8), intent(in) :: sband_rain_coefs, sband_snow_coefs, sband_graupel_coefs ! coefficients read by calling init_coefs
    !real (kind=8), intent(in) :: cband_rain_coefs, cband_snow_coefs, cband_graupel_coefs ! coefficients read by calling init_coefs
    real (kind=8), intent(out):: zh_rain, zdr_rain, kdp_rain, phv_rain
    real (kind=8), intent(out) :: zh_snow, zdr_snow, kdp_snow, phv_snow
    real (kind=8), intent(out) :: zh_graupel, zdr_graupel, kdp_graupel,phv_graupel
    real (kind=8), intent(out) :: zh, zdr, kdp, phv
    real (kind=8) :: rho_rain, rho_snow, rho_graupel ! g/cm^3
    real (kind=8) :: dm ! mean diameter, mm
    real (kind=8) :: w   ! g/m^3
    real (kind=8) :: ratio   ! water fraction in melting snow/graupel, 0-1
    real (kind=8) :: rho ! density, g/cm^3
    real (kind=8) :: z, a(0:12)
    real (kind=8) :: zhzdr_rain, zhzdr_snow, zhzdr_graupel
    real(kind=8), dimension(1:4,0:4)  :: rain_coefs !, sband_rain_coefs, cband_rain_coefs
    real(kind=8), dimension(0:12,0:3) :: snow_coefs !, sband_snow_coefs, cband_snow_coefs
    real(kind=8), dimension(0:12,0:3) :: graupel_coefs !, sband_graupel_coefs, cband_graupel_coefs

    !call init_coefs
    !print *,"Cband rain coefs = ",cband_rain_coefs

! 1.0 RAIN
!-------------
    if (iband == 1) then
       rain_coefs = sband_rain_coefs
       snow_coefs = sband_snow_coefs
       graupel_coefs = sband_graupel_coefs
    elseif (iband == 2) then
       rain_coefs = cband_rain_coefs
       snow_coefs = cband_snow_coefs
       graupel_coefs = cband_graupel_coefs

       !print *,'rain_coefs= ',rain_coefs
    endif
    w = air_density * qrain ! kg/m^3 * g/kg = g/m^3
    !w = rhod * qrain
    !rho_rain = 1.0 ! 1 g/cm^3
    rho_rain = density_rain ! 1 g/cm^3
    if (qrain > 10e-10) then
       call dm_z_2moment(air_density, qrain, rho_rain, ntx_rain, dm, z)
       call dualpol_op_rain(w, dm, rain_coefs, zh_rain, zdr_rain, kdp_rain, phv_rain)

     !write(*,'(a,3f10.5)') 'Rain: q(g/kg) w(g/m^3) dm(mm) = ',qrain,w,dm
     !write(*,'(a,4f12.5)') 'ZH ZDR KDP PHV = ', 10.0*log10(zh_rain), &
     !                           10.0*log10(zdr_rain), kdp_rain, phv_rain
    else
      zh_rain = 0.
    endif

! 2.0 SNOW
!-----------------------
    w = air_density * qsnow ! kg/m^3 * g/kg = g/m^3
    if ((qrain+qsnow) > 0.0 .and. qsnow > 0.0) then
      ratio = qrain/(qrain+qsnow) ! water fraction
      call weticephase_density('snow', ratio, rho_snow)
      call dm_z_2moment(air_density, qsnow, rho_snow, ntx_snow, dm, z)

    ! calcuate a coefficients
      do i = 0,12
        call coef_a(snow_coefs(i,0:3), ratio, a(i))
      !print *, coef_a(i)
      enddo

!   calcuate dual-pol variables
      call dualpol_op_icephase(z, w, rho_snow, dm, a(0:12), &
                        zh_snow, zdr_snow, kdp_snow, phv_snow)

    !write(*,'(a,4f10.5)') 'Snow q(g/kg) w(g/m^3) dm(mm) ratio = ',qsnow,w,dm,ratio
    !write(*,'(a,4f12.5)') 'ZH ZDR KDP PHV = ',10.0*log10(zh_snow), &
    !                        10.0*log10(zdr_snow), kdp_snow, phv_snow

    else
      zh_snow = 0.
    endif

! 3.0 Graupel
!-----------------------
    w = air_density * qgraupel ! kg/m^3 * g/kg = g/m^3
    if ((qrain+qgraupel) > 0.0 .and. qgraupel > 0.0) then
      ratio = qrain/(qrain+qgraupel) ! water fraction
      call weticephase_density('graupel', ratio, rho_graupel) ! g/cm^3
      call dm_z_2moment(air_density, qgraupel, rho_graupel, ntx_graupel, dm, z)

    ! calcuate a coefficients
      do i = 0,12
        call coef_a(graupel_coefs(i,0:3), ratio, a(i))
      !print *, coef_a(i)
      enddo

!   calcuate dual-pol variables
      call dualpol_op_icephase(z, w, rho_graupel, dm, a(0:12), &
                        zh_graupel, zdr_graupel, kdp_graupel, phv_graupel)

    !write(*,'(a,4f10.5)') 'Graupel q(g/kg) w(g/m^3) dm(mm) ratio = ',qgraupel,w,dm,ratio
    !write(*,'(a,4f12.5)') 'ZH ZDR KDP PHV = ',10.0*log10(zh_graupel), &
    !                        10.0*log10(zdr_graupel), kdp_graupel, phv_graupel
    else
      zh_graupel = 0.
    endif
! 4.0 total zh, zdr, kdp, phv
! Eqs. 22, 23, 24, 25 of Zhang et al., 2021

    call dualpol_op_total( &
         zh_rain,  zh_snow,  zh_graupel, &
        zdr_rain, zdr_snow, zdr_graupel, &
        kdp_rain, kdp_snow, kdp_graupel, &
        phv_rain, phv_snow, phv_graupel, &
         zh, zdr, kdp, phv)

  end subroutine scheme_2moment

end module module_ppro

  
