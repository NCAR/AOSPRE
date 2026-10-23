module module_quicksort
  use module_configuration, only : RKIND
  implicit none
contains
  subroutine quicksort(array, array_size)
    implicit none
    integer,                     intent(in)    :: array_size
    real(kind=RKIND), dimension(array_size), intent(inout) :: array
    integer, dimension(array_size) :: array_index
    real(kind=RKIND), dimension(array_size) :: sorted
    integer :: i
    call sortp_1r8(array_size, array_index, array)
    sorted = (/ (array(array_index(i)), i=1, array_size) /)
    array = sorted
  end subroutine quicksort
end module module_quicksort

!---------------------------------------------------------------
! Sort a single-precision real array by index
subroutine sortp_1r4(array_size,index,value)
  integer, intent(in) :: array_size
  integer, intent(inout) :: index(array_size)
  real(4), intent(in) :: value(array_size)
  integer :: qsort_threshold = 8
  include "quicksort_inline.inc"
contains
  include "quicksort_inline_index.inc"
  logical &
  function less_than(a,b)
    integer, intent(in) :: a,b
    ! real(4), parameter :: small=1.0e-6
    ! if ( abs(value(index(a))-value(index(b))) < small ) then
    !   less_than = index(a) < index(b)
    ! else
      less_than = value(index(a)) < value(index(b))
    ! end if
  end function less_than
end subroutine sortp_1r4
!---------------------------------------------------------------
! Sort a real(kind=8) array by index
subroutine sortp_1r8(array_size,index,value)
  integer, intent(in) :: array_size
  integer, intent(inout) :: index(array_size)
  real(8), intent(in) :: value(array_size)
  integer :: qsort_threshold = 8
  include "quicksort_inline.inc"
contains
  include "quicksort_inline_index.inc"
  logical &
  function less_than(a,b)
    integer, intent(in) :: a,b
    ! real(8), parameter :: small=1.0e-6
    ! if ( abs(value(index(a))-value(index(b))) < small ) then
    !   less_than = index(a) < index(b)
    ! else
      less_than = value(index(a)) < value(index(b))
    ! end if
  end function less_than
end subroutine sortp_1r8
!---------------------------------------------------------------
! Sort an array of integers
subroutine sort_1i(array_size,i1)
  integer, intent(in) :: array_size
  integer, intent(inout), dimension(array_size) :: i1
  integer :: qsort_threshold = 8
  include "quicksort_inline.inc"
contains
  subroutine init()
  end subroutine init
  logical &
  function less_than(a,b)
    integer, intent(in) :: a,b
    if ( i1(a) == i1(b) ) then
      less_than = a < b
    else
      less_than = i1(a) < i1(b)
    end if
  end function less_than
  subroutine swap(a,b)
    integer, intent(in) :: a,b
    integer :: hold
    hold=i1(a); i1(a)=i1(b); i1(b)=hold
  end subroutine swap
! circular shift-right by one:
  subroutine rshift(left,right)
    integer, intent(in) :: left, right
    integer :: hold
    hold=i1(right); i1(left+1:right)=i1(left:right-1); i1(left)=hold
  end subroutine rshift
end subroutine sort_1i
!---------------------------------------------------------------
! Sort a string array, with any string length.
subroutine sortp_string(array_size,index,string)
  integer, intent(in) :: array_size
  integer, intent(out) :: index(array_size)
  character(len=*), intent(in) :: string(array_size)
  integer :: qsort_threshold = 8
  include "quicksort_inline.inc"
contains
  include "quicksort_inline_index.inc"
  logical &
  function less_than(a,b)
    integer, intent(in) :: a,b
    if ( string(index(a)) == string(index(b))  ) then
      less_than = ( index(a) < index(b) )
    else
      less_than = ( string(index(a)) < string(index(b)) )
    end if
  end function less_than
end subroutine sortp_string
