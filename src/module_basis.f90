!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_basis.f90
!> @file
!> @breif 坐标系基矢量文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  Module: basis
!
!  DESCRIPTION:
!> @breif 坐标系基矢量模块.
!>
!!
!  REVISION HISTORY:
!  2017-07-26 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-07-26
!------------------------------------------------------------------------------

module mod_basis

    use mod_vector
    use penf, only: R_P
    implicit none

    private
    public :: BASIS_NULL

    real(R_P), parameter :: ONE=1.0d0 !< 1 constant\private
    real(R_P), parameter :: ZERO=0.0d0 !< 0 constant\private

    !>  坐标系基矢量类.
    !!
    type, public :: basis_type

        type(vector_type) :: e(3) !<坐标基矢量\private

    Contains

        Procedure :: SetBasis => set_basis !<设置坐标系基矢量
        procedure :: GetBasis => get_basis !<读取坐标系基矢量
        procedure :: finalize => finalize_basis !<析构函数

    end type basis_type

    !>坐标单位基矢量
    type(basis_type), parameter :: BASIS_NULL=&
        &   basis_type([VECTOR100, VECTOR010, VECTOR001])

contains

    !>  读取坐标系基矢量.
    !!
    !!  @param[out] e1 Xi方向基矢量
    !!  @param[out] e2 Eta方向基矢量
    !!  @param[out] e3 Zeta方向基矢量
    subroutine get_basis(this, e1, e2, e3)

        implicit none
        type(vector_type), intent(inout) :: e1, e2, e3
        class(basis_type), intent(in) :: this

        e1=this%e(1); e2=this%e(2); e3=this%e(3)

    end subroutine get_basis

    !>  设置坐标系基矢量.
    !!
    !!  @param[in] e1 Xi方向基矢量
    !!  @param[in] e2 Eta方向基矢量
    !!  @param[in] e3 Zeta方向基矢量
    subroutine set_basis(this, e1, e2, e3)

        implicit none
        type(vector_type), intent(in) :: e1, e2, e3
        class(basis_type), intent(inout) :: this

        this%e(1)=e1; this%e(2)=e2; this%e(3)=e3

    end subroutine set_basis

    !>  析构函数.
    !!
    Subroutine finalize_basis(this)

        Implicit None
        class(basis_type), intent(Inout) :: this

        this%e(1)=VECTOR100
        this%e(2)=VECTOR010
        this%e(3)=VECTOR001

    End Subroutine finalize_basis

end module mod_basis
