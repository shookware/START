!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: /media/psf/Home/Fortran/START/src/module_tensor.f90
!> @file
!> @breif 张量模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: tensor
!> @breif 张量模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-01 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-08-01
!------------------------------------------------------------------------------
module mod_tensor

   use mod_vector
   use penf, only: R_P
   implicit none

   private

   !> 张量类型.
   type, public :: tensor_type

      private
      type(vector_type), private :: a_i(3) !<张量分量

   Contains

      generic :: Set => set_tensor, set_tensor_matrix
      generic :: Get => Get_tensor, get_tensor_matrix

      Procedure, private :: get_tensor           !<获得张量
      Procedure, private :: set_tensor           !<设置张量
      Procedure, private :: set_tensor_matrix   !从矩阵设置张量
      Procedure, private :: get_tensor_matrix   !从矩阵获得张量
      procedure :: Zeros !< 置零函数

   end type tensor_type

contains

   !> 构造函数:从矢量构造张量.
   !!
   !! @param[in] a_i 矢量\f$a_i\f$
   !! @return 一个二阶张量
   type(tensor_type) function Tensor(a_i)

      implicit none
      type(vector_type), intent(in) :: a_i(3)
      tensor%a_i = a_i

   end function Tensor

   !> 构造函数:从3*3矩阵构造二阶张量.
   !!
   !! @param[in] a 矩阵A
   !! @return 一个二阶张量
   type(tensor_type) function TensorMatrix(a)

      implicit none
      real(R_P), intent(in) :: a(3, 3)

      tensormatrix%a_i(1) = Vector(a(1, 1), a(1, 2), a(1, 3))
      tensormatrix%a_i(2) = Vector(a(2, 1), a(2, 2), a(2, 3))
      tensormatrix%a_i(3) = Vector(a(3, 1), a(3, 2), a(3, 3))

   end function TensorMatrix

   !> 用矢量设置张量.
   !!
   !! @param[in] a 矢量数组\f$a\f$
   subroutine set_tensor(this, a)

      implicit none
      type(vector_type), intent(in) :: a(3)
      class(tensor_type), intent(inout) :: this

      this%a_i = a

   end subroutine set_tensor

   !> 读取张量到矢量.
   !!
   !! @param[out] a 矢量数组\f$a\f$
   subroutine Get_tensor(this, a)

      implicit none
      type(vector_type), intent(inout) :: a(3)
      class(tensor_type), intent(in) :: this

      a = this%a_i

   end subroutine get_tensor

   !> 用3*3矩阵设置张量.
   !!
   !! @param[in] a 3*3矩阵
   subroutine set_tensor_matrix(this, a)

      implicit none
      class(tensor_type), intent(inout) :: this
      real(R_P), intent(in) :: a(3, 3)
      type(vector_type) :: a_i(3)

      call a_i(1)%Set(a(1, 1), a(1, 2), a(1, 3))
      call a_i(2)%Set(a(2, 1), a(2, 2), a(2, 3))
      call a_i(3)%Set(a(3, 1), a(3, 2), a(3, 3))

      call this%set(a_i)

   end subroutine set_tensor_matrix

   !> 读取张量到3*3矩阵.
   !!
   !! @param[out] a 3*3矩阵
   subroutine get_tensor_matrix(this, a)

      implicit none
      class(tensor_type), intent(in) :: this
      real(R_P), intent(inout) :: a(3, 3)
      type(vector_type) :: a_i(3)

      call this%Get(a_i)

      call a_i(1)%Get(a(1, 1), a(1, 2), a(1, 3))
      call a_i(2)%Get(a(2, 1), a(2, 2), a(2, 3))
      call a_i(3)%Get(a(3, 1), a(3, 2), a(3, 3))

   end subroutine get_tensor_matrix

   !> 置零函数.
subroutine Zeros(this)

    implicit none
    class(tensor_type), intent(inout) :: this

    this%a_i=VEC_NULL

end subroutine Zeros

end module mod_tensor
