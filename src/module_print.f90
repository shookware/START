!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_print.f90
!> @file
!> @breif 输出模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: print
!> @breif 输出模块.
!  DESCRIPTION:
!>
!  REVISION HISTORY:
!  2017-08-02 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-08-02
!------------------------------------------------------------------------------
module mod_print

    use penf, only: R_P
    use stringifor
    implicit none

    public :: PrintComplexArray2D
    private

    contains

    !> 输出复数型的二维数组
    !! @param[in] in 流向点数
    !! @param[in] jn 法向点数
    !! @param[ln] ln 每点自由度个数
    !! @param[in] Array 被输出数组(in, jn, ln)
    !! @param[in] fn_surf 输出文件文件名
    subroutine PrintComplexArray2D(in, jn, ln, Array, fn_surf)

      implicit none
      integer, intent(in) :: in, jn, ln
      complex(R_P), intent(in) :: Array(in, jn, ln)
      type(string), intent(in) :: fn_surf
      integer :: i, j, l

      open(999, file=fn_surf%chars(), form='unformatted')
        write(999)1
        write(999)in, jn, ln*3
        write(999)((( real(Array(i, j, l)), i=1, in), j=1, jn), &
               &   ((aimag(Array(i, j, l)), i=1, in), j=1, jn), &
               &   ((abs  (Array(i, j, l)), i=1, in), j=1, jn), l=1, ln)
        close(999)

    end subroutine PrintComplexArray2D

end module mod_print
