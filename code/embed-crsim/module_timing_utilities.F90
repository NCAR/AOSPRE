module module_timing_utilities
  use module_configuration, only : RKIND
  implicit none

  private

  integer, parameter :: MAX_TIMERS = 100
  type :: timer_list
      character(len=1024) :: name = ""
      real(kind=RKIND)    :: start
      real(kind=RKIND)    :: total
  end type timer_list
  type(timer_list), dimension(MAX_TIMERS) :: timers

  public :: timer_start
  public :: timer_pause
  public :: timer_print
  public :: timer_resume
  public :: timer_value

contains

  !
  !-------------------------------------------------------------------------------
  !

  subroutine timer_start(name)
    implicit none
    character(len=*),   intent(in)    :: name
    integer :: k
    do k = 1, MAX_TIMERS
        if (timers(k)%name == "") then
            timers(k)%name = trim(name)
        endif
        if (timers(k)%name == name) then
            timers(k)%total = 0.0
            call cpu_time(timers(k)%start)
            write(*,'("Start timer ", A, " at cpu_time ", F12.4)') name, timers(k)%start
            return
        endif
    enddo
    write(*,'( "***** ERROR in TIMER_INITIALIZE ***** Timer name ", A, " not created.")') trim(name)
    stop "TIMER INITIALIZE"
  end subroutine timer_start

  !
  !-------------------------------------------------------------------------------
  !

  subroutine timer_pause(name)
    implicit none
    character(len=*), intent(in) :: name
    real :: instant
    integer :: k

    do k = 1, MAX_TIMERS
        if (timers(k)%name == name) then
            call cpu_time(instant)
            timers(k)%total = timers(k)%total + (instant - timers(k)%start)
            timers(k)%start = -1.E36
            return
        endif
    enddo
    write(*,'( "***** ERROR in TIMER_PAUSE ***** Timer name ", A, " not found.")') trim(name)
    stop "TIMER_PAUSE"
  end subroutine timer_pause

  !
  !-------------------------------------------------------------------------------
  !

  subroutine timer_resume(name)
    implicit none
    character(len=*), intent(in) :: name
    integer :: k
    do k = 1, MAX_TIMERS

        if (timers(k)%name == "") then
            ! Timer of name <name> is a new timer, so initialize the total to zero.
            timers(k)%name = trim(name)
            timers(k)%total = 0.0
            call cpu_time(timers(k)%start)
            write(*,'("New timer start ", A, " at cpu_time ", F12.4)') name, timers(k)%start
            return
        endif

        if (timers(k)%name == name) then
            ! Timer of name <name> is an existing timer, so do not reinitialize total.
            call cpu_time(timers(k)%start)
            ! write(*,'("Resume timer ", A, "from ", F12.4, " at cpu_time ", F12.4)') name, timers(k)%total, timers(k)%start
            return
        endif
    enddo
    write(*,'( "***** ERROR in TIMER_RESUME ***** Timer name ", A, " not found.")') trim(name)
    stop "TIMER_RESUME"
  end subroutine timer_resume

  !
  !-------------------------------------------------------------------------------
  !

  subroutine timer_print(name, label, percent_of)
    implicit none
    character(len=*), intent(in) :: name
    character(len=*), optional, intent(in) :: label
    real(kind=RKIND), optional, intent(in) :: percent_of
    real                            :: instant
    integer :: k
    do k = 1, MAX_TIMERS
        if (timers(k)%name == name) then
            call cpu_time(instant)
            if ( present(label) ) then
                if (present(percent_of)) then
                    write(*,'(A35, F12.5, " (",F8.4,"%)")') label, timers(k)%total, 1.E2*timers(k)%total/percent_of
                else
                    write(*,'(A35, F12.5)') label, timers(k)%total
                endif
            else
                if (present(percent_of)) then
                    write(*,'("TIMER: ", A, " = ", F12.5, " (", F8.4, "%)")') trim(name), timers(k)%total, 1.E2*timers(k)%total/percent_of
                else
                    write(*,'("TIMER: ", A, " = ", F12.5)') trim(name), timers(k)%total
                endif
            endif
            return
        endif
    enddo
    write(*,'( "***** ERROR in TIMER_PRINT ***** Timer name ", A, " not found.")') trim(name)
    stop "TIMER_PRINT"

  end subroutine timer_print

  !
  !-------------------------------------------------------------------------------
  !

  real function timer_value(name)
    implicit none
    character(len=*), intent(in) :: name
    real                            :: instant
    integer :: k
    do k = 1, MAX_TIMERS
        if (timers(k)%name == name) then
            timer_value = timers(k)%total
            return
        endif
    enddo
    write(*,'( "***** ERROR in TIMER_VALUE ***** Timer name ", A, " not found.")') trim(name)
    stop "TIMER_PRINT"

  end function timer_value

  !
  !-------------------------------------------------------------------------------
  !

end module module_timing_utilities
