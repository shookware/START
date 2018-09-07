!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_lns_OP_normal.f90
!> @file
!> @breif LNS算子沿法向分布文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  Module: lns_OP_normal
!
!  DESCRIPTION:
!> @breif LNS算子沿法向分布模块.
!>
!! 实际上是LNS方程的算子
!  REVISION HISTORY:
!  2017-08-09 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-08-09
!------------------------------------------------------------------------------
module mod_lns_OP_normal

    use mod_lns_OP_point
    use mod_grid
    implicit none

    private

    !> LNS基本流算子沿法向分布
    type, public :: lns_OP_normal_type

        private
        integer, private :: jn !< 法向点数
        class(lns_OP_point_type), allocatable, private :: bf_op(:) !< LNS方程系数

        Contains

          Procedure :: Create !< 创建lns_OP_normal类型,分配内存
          procedure :: Set !<设置法向各点基本流、坐标系、是否平行
          procedure :: finalize => finalize_lns_OP_normal !< 析构函数

    end type lns_OP_normal_type

    !> LPSE用算子
    type, extends(lns_OP_normal_type), public :: lpse_bf_op_normal_type

        private

        Contains

          Procedure :: GetPoint=>getpoint_lpse !<获得法向某点处的基本流算子系数

    end type lpse_bf_op_normal_type

    !> LST用算子
    type, extends(lns_OP_normal_type), public :: lst_bf_op_normal_type !<获得法向某点处的基本流算子系数

        private

        Contains

          Procedure :: GetPoint=>getpoint_lst

    end type lst_bf_op_normal_type

    contains

    subroutine create(this, jn, eqkind)

        implicit none
        class(lns_OP_normal_type), intent(inout) :: this
        integer, intent(in) :: jn
        integer, intent(in) :: eqkind
        integer, parameter :: LPSE=0
        integer, parameter :: LST=1

        this%jn=jn
        select case (eqkind)
            case (LPSE)
                if(.not. allocated(this%bf_op)) &
                &   allocate(lpse_bf_op_point_type::this%bf_op(jn))
            case (LST)
                if(.not. allocated(this%bf_op)) &
                &   allocate(lst_bf_op_point_type::this%bf_op(jn))
            case default
                stop "Please input a correct eqkind!"
        end select

    end subroutine create

    !> 设置法向各点基本流、坐标系、是否平行
    !! @param[in] BFNorm 法向基本流分布(含导数)
    !! @param[in] NormCoord 局部的坐标系
    !! @param[in] isParallel 是否考虑非平行性
    subroutine set(this, BFNorm, NormCoord, isParallel)

        use mod_local_normal_coordinate
        use mod_BF_normal

        implicit none

        class(lns_OP_normal_type), intent(inout) :: this
        type(local_normal_coordinate_type), intent(in) :: NormCoord
        type(BF_normal_type), intent(in) :: BFNorm
        logical, intent(in), optional :: isParallel
        integer :: J

        do J=1, this%jn
          call this%bf_op(j)%SetCoord(NormCoord%GetPoint(j))
          call this%bf_op(j)%SetFromFlow(BFNorm%GetPoint(j), isParallel)
        enddo
!        if( isParallel==.False.) stop

    end subroutine set

    !> 获得法向某点处的基本流算子系数
    !! @param[in] j 法向点序号
    !! @return 法向j点处的基本流算子系数
    function getpoint_lpse(this, j) result(BFOPPoint)

        implicit none
        class(lpse_bf_OP_normal_type), intent(in) :: this
        integer, intent(in) :: j
        type(lpse_BF_OP_point_type) :: BFOPPoint

        selecttype (bf_op=>this%bf_op(j))
        type is (lpse_bf_op_point_type)
            BFOPPoint=bf_op
        endselect

    end function getpoint_lpse

    !> 获得法向某点处的基本流算子系数
    !! @param[in] j 法向点序号
    !! @return 法向j点处的基本流算子系数
    function getpoint_lst(this, j) result(BFOPPoint)

        implicit none
        class(lst_bf_OP_normal_type), intent(in) :: this
        integer, intent(in) :: j
        type(lst_BF_OP_point_type) :: BFOPPoint

        selecttype (bf_op=>this%bf_op(j))
        type is (lst_bf_op_point_type)
            BFOPPoint=bf_op
        endselect

    end function getpoint_lst

    !> 析构函数
    Subroutine finalize_lns_OP_normal(this)

        Implicit None
        class(lns_OP_normal_type), intent(Inout) :: this

        this%jn=0
        deallocate(this%bf_op)

    End Subroutine finalize_lns_OP_normal

end module mod_lns_OP_normal
