!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_lpse_dis_op_normal.f90
!> @file
!> @breif PSE扰动算子系数沿法向分布.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: lpse_dis_op_normal
!> @breif PSE扰动线性算子系数沿法向分布模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-04 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-04
!------------------------------------------------------------------------------
module mod_lpse_dis_OP_normal

    use mod_lpse_dis_OP_point
    use penf, only: R_P
    implicit none

    private

    !> PSE扰动线性算子系数沿法向分布类
    type, public :: lpse_dis_op_normal_type

        private
        integer, private :: jn !< 法向点数
        integer, private :: iloc !< 流向站位
        type(lpse_dis_op_point_type), allocatable, private :: DisOP(:) !< 扰动线性算子系数
        real(R_P), allocatable, public :: Coef(:)

    Contains

        Procedure :: Create !< 创建lpse_dis_op_normal类,分配内存
        procedure :: Set !<设置扰动线性算子系数
        procedure :: SetCoef !< 设置流向差分算子系数
        procedure :: GetPoint !< 获得某法向位置处的扰动线性算子系数
        procedure :: finalize => finalize_lpse_dis_op_normal !< 析构函数

    end type lpse_dis_op_normal_type

contains

    !> 设置流向差分算子系数
    subroutine SetCoef(this, Coef)

        use penf, only: R_P
        implicit none
        class(lpse_dis_op_normal_type),intent(inout) :: this
        real(R_P) :: Coef(:)

        if(.not. allocated(this%Coef)) allocate(this%Coef, source=Coef)

    end subroutine SetCoef


    !> 创建lpse_dis_op_normal类,分配内存
    subroutine create(this, jn)

        implicit none
        integer, intent(in) :: jn !< 法向点数
        class(lpse_dis_op_normal_type), intent(inout) :: this

        this%jn=jn
        if(.not. allocated(this%DisOP)) allocate(this%DisOP(jn))

    end subroutine create

    !> 设置扰动线性算子系数
    subroutine set(this, BFOPNorm, NormCoord, wavenum, lpsebc)

        use mod_local_normal_coordinate
        use mod_lns_OP_normal
        use mod_dis_wavenum
        use mod_lpse_bc
        use mod_difference

        implicit none
        class(lpse_dis_op_normal_type), intent(inout) :: this
        type(local_normal_coordinate_type), intent(in) :: NormCoord !< 当地站位的坐标系特征
        type(lpse_bf_op_normal_type), intent(in) :: BFOPNorm !< 基本流(LNS)算子沿法向分布
        type(dis_wavenum_lpse_type), intent(in) :: wavenum !< 扰动色散特效
        type(lpse_bc_type), intent(in) :: lpsebc(2) !< 边界条件的算子
        integer :: J

        call this%DisOP(1)%SetBC(lpsebc(1))
        do J=2, this%jn-1
            call this%DisOP(j)%Set(BFOPNorm%GetPoint(j), wavenum, NormCoord%GetPoint(j))
        enddo
        call this%DisOP(this%jn)%SetBC(lpsebc(2))

    end subroutine set

    !> 获得某法向位置的扰动线性算子系数
    function getpoint(this, j)  result(point)

        implicit none
        class(lpse_dis_op_normal_type), intent(in) :: this
        integer, intent(in) :: j !< 法向站位
        type(lpse_dis_op_point_type) :: point !< j点处的线性算子系数

        point=this%DisOP(j)

    end function getpoint

    !> 析构函数
    Subroutine finalize_lpse_dis_op_normal(this)

        Implicit None
        class(lpse_dis_op_normal_type), intent(Inout) :: this

        deallocate(this%DisOP)
        this%jn=0

    End Subroutine finalize_lpse_dis_op_normal

end module mod_lpse_dis_OP_normal
