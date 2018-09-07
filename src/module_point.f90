!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_context.f90
!> @file
!> @breif 坐标点类文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  Module: context
!
!  DESCRIPTION:
!> @brief 坐标点模块.
!!
!  REVISION HISTORY:
!  2015-07-26 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-07-26
!------------------------------------------------------------------------------
module mod_point

    use penf, only: R_P
    implicit none

    private

    !>  坐标点类.
    !!
    type, public :: point_type

        real(R_P), private :: x(3) !<坐标点的值

        Contains

          Procedure :: Set => setpoint !<设置坐标点
          Procedure :: Get => get_point                                    !<读坐标点
          procedure :: finalize => finalizePoint !<析构函数
          procedure :: Zeros !< 置零函数
          !procedure :: PointDistance => Point_Distance
          !Final     ::
          !generic :: operator(.distance.) => PointDistance
          generic :: operator(.abs.) => pointabs !<求点到原点的长度
          generic :: operator(+) => AddPoint !<点坐标相加
          generic :: operator(-) => MinusPoint !<点坐标相减
          generic :: operator(.addi.) => addi !<x坐标点平移
          generic :: operator(.addj.) => addj !<y坐标点平移

          procedure, private :: AddPoint !<点坐标相加
          procedure, private :: MinusPoint => minus_point !<点坐标相减
          procedure, private :: pointabs => normal_abs !<求点到原点的长度
          procedure, private :: addi => add_i !<x坐标点平移
          procedure, private :: addj => add_j !<y坐标点平移

    end type point_type


    contains

    !>  点相加.
    !!
    !!  @param[in] point1 点1
    !!  @param[in] point2 点2
    !!  @return 点1和点2的和
    function AddPoint(point1, point2) result(add)

        implicit none
        class(point_type), intent(in) :: point1
        type(point_type), intent(in) :: point2
        type(point_type) :: add

        add%x=point1%x+point2%x

    end function AddPoint

    !>  x方向上平移.
    !!
    !!  @param[in] distance x方向上平移距离
    !!  @return 平移后的点
    type(point_type) function add_i(this, distance)

        implicit none
        class(point_type), intent(in) :: this
        real(R_P), intent(in) :: distance

        add_i=this
        add_i%x(1)=add_i%x(1)+distance

    end function add_i

    !>  y方向上平移.
    !!
    !!  @param[in] distance y方向上平移距离
    !!  @return 平移后的点
    type(point_type) function add_j(this, distance)

        implicit none
        class(point_type), intent(in) :: this
        real(R_P), intent(in) :: distance

        add_j=this
        add_j%x(2)=add_j%x(2)+distance

    end function add_j

    !>  读坐标点.
    !!
    !!  @param[out] x x坐标
    !!  @param[out] y y坐标
    !!  @param[out] z z坐标
    Elemental Subroutine get_point(this, x, y, z)

        Implicit None
        Class(point_type), intent(In) :: this
        real(R_P) , intent(inout) :: x, y, z

       x=this%x(1); y=this%x(2); z=this%x(3)

    End Subroutine get_point

    !>  点相减.
    !!
    !!  @param[in] point1 点1
    !!  @param[in] point2 点2
    !!  @return 点1和点2的差
    type(point_type) function minus_point(point1, point2)

        implicit none
        class(point_type), intent(in) :: point1
        type(point_type), intent(in) :: point2
        !class(point_type), intent(inout) :: minus_point

        minus_point%x=point1%x-point2%x

    end function minus_point

    !>  点到原点的距离.
    !!  @return 点到原点的距离
    !!
    real(R_P) function normal_abs(point1)

        implicit none
        class(point_type), intent(in) :: point1

        normal_abs=sqrt(point1%x(1)**2+point1%x(2)**2+point1%x(3)**2)

    end function normal_abs

    !>  设置坐标点.
    !!
    !!  @param[in] x x坐标
    !!  @param[in] y y坐标
    !!  @param[in] z z坐标
    !!  @return 返回一个点
    Elemental Subroutine setpoint(this, x, y, z)

        Implicit None
        real(R_P), intent(In) :: x, y, z
        Class(point_type), intent(inout) :: this

        this%x(1) = x;this%x(2) = y; this%x(3) = z

    End Subroutine setpoint

    !> 置零函数.
    elemental subroutine Zeros(this)

        implicit none
        class(point_type), intent(inout) :: this

        this%x=0.0d0

    end subroutine Zeros

    !>  析构函数.
    !!
    elemental Subroutine finalizePoint(this)

        Implicit None
        class(point_type), intent(Inout) :: this

        this%x=0.0d0

    End Subroutine finalizePoint

end module mod_point
