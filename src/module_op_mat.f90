!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: mod_op_mat.f90
!> @file
!> @breif 矩阵操作模块文件_实型.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: op_mat
!> @breif 矩阵操作模块(实型).
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-03 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-03
!------------------------------------------------------------------------------
module mod_op_mat

    use penf, only: R_P
    implicit none
    private

    !> 矩阵操作类
    type, public :: op_mat_type

        private

        real(R_P), private :: Mat(5, 5) !< 实型的5*5的矩阵

        Contains

          Procedure, pass(this) :: Set => Set !< 设置矩阵的值
          procedure, pass(this) :: Get => Get !< 获得矩阵的值
          procedure, pass(this) :: Print => PrintMat !< 输出矩阵

          generic :: ASSIGNMENT(=) => Set, Get !<定义赋值操作

    end type op_mat_type

    contains

    !< 设置矩阵的值
    !! @param[in] Mat 5*5的矩阵
    subroutine Set(this, Mat)

        implicit none
        real(R_P), intent(in) :: Mat(5, 5)
        class(op_mat_type), intent(inout) :: this

        this%Mat=Mat

    end subroutine Set

    !< 获得矩阵的值
    !! @param[out] Mat 5*5的矩阵
    subroutine Get(Mat, this)

        implicit none
        real(R_P), intent(inout) :: Mat(5, 5)
        class(op_mat_type), intent(in) :: this

        Mat=this%Mat

    end subroutine Get

    !> 输出矩阵的值
    subroutine PrintMat(this)

        use mod_parameter, only: eno => ENO_BF
        implicit none
        class(op_mat_type), intent(in) :: this
        integer :: i, j

        do i=1, 4
          do j=1, 5
            write(eno, *)this%Mat(i, j), i, j
          end do
        end do

    end subroutine PrintMat

end module mod_op_mat
