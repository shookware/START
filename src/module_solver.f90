!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  Module: solver
!
!  DESCRIPTION:
!> @breif 稳定性求解器类.
!>
!!
!  REVISION HISTORY:
!  yyyy-mm-dd - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-09-22
!------------------------------------------------------------------------------
module mod_solver

    use mod_grid
    use mod_baseflow
    use mod_difference
    !use mod_lpse_bc
    use mod_curvature_2d
    use mod_dis
    !use mod_lpse_dis_normal
    !use mod_baseflow_org
    use penf, only: R_P

    implicit none

    public solver_type

    !> 求解器类
    type, abstract :: solver_type

        type(difference_2D_type), pointer :: Diff !< 整个二维场扰动离散的拉格朗日基函数
        type(grid2d_type), pointer :: Grid !< 二维场网格
        type(Baseflow2D_type), pointer :: BaseFlow !< 二维基本流
        type(curvature_2d_type), pointer :: Curvature !< 二维曲率
        integer :: iloc !< 当地初始流向站位编号
        integer :: Istart !< 计算域结束站位
        integer :: Iend !< 计算域结束站位
        integer :: bctype(5, 2) !< PSE方程边界条件类型

        Contains

          procedure :: Initialize !< 初始化PSE求解器
          procedure :: SetIStartEnd !< 设置起始和终止站位序号
          procedure :: SetDisDiffScheme !< 设置扰动离散格式
          procedure :: SetGrid !< 设置2D网格
          procedure :: SetBaseFlow !< 设置2D基本流
  !        procedure :: SetDBaseFlow !< 设置2D基本流的导数
          procedure :: SetCurvature !< 设置2D曲率
          procedure :: SetILoc=>set_iloc !<设置当前站位
          procedure :: SetBCType !< 设置边界条件
          procedure :: PrintILoc !< 输出流向序号

    end type solver_type

    !abstract interface
    !    subroutine av_interface(this, arrayX, arrayB) ! Lx=rhs
    !        import
    !        Implicit None
    !        class(solver_type), intent(in) :: this
    !        real*8, intent(in) :: Arrayx(:)
    !        real*8, intent(out) :: arrayB(:)
    !
    !    end subroutine av_interface
    !end interface

    contains

      subroutine Initialize(this, grid, BaseFlow, curvature, bctype)

           implicit none
           class(solver_type), intent(inout) :: this
           type(grid2d_type), intent(in) :: grid
           type(Baseflow2D_type), intent(in) :: baseflow
           type(curvature_2d_type), intent(in) :: curvature
           integer, intent(in) :: bctype(5, 2) !< 边界条件类型
      !     type(difference_2D_type), intent(in) :: diff

      !     call this%SetILoc(iloc)
           call this%SetGrid(grid)
           call this%SetBaseFlow(Baseflow)
           call this%SetCurvature(curvature)
           call this%SetBCType(BCType)
      !     this%bctype=bctype
      !     call this%SetDisDiffScheme(diff)

      end subroutine Initialize

      !> 设置起始和终止的站位序号
      subroutine SetIstartend(this, Istart, Iend)

          implicit none
          integer, intent(inout) :: Istart !< 起始流向站位序号
          integer, intent(inout) :: Iend !< 终止流向站位序号
          class(solver_type), intent(inout) :: this

          this%Istart=Istart; this%Iend=Iend
          if(this%iend==-1 .or. (this%iend>this%Grid%GetInSize()))then
            this%iend=this%Grid%GetInSize()
            iend=this%iend
          end if
          if(this%istart==-1 .or. (this%istart<1))then
            this%istart=1
            istart=1
          end if

      end subroutine SetIstartend


      !> 设置扰动离散格式
      subroutine SetDisDiffScheme(this, Diff_Dis)

           implicit none
           class(solver_type), intent(inout) :: this
           type(difference_2D_type), target, intent(in) :: Diff_Dis !< 二维场扰动离散的拉格朗日基函数

           if(.not. (associated(this%Diff, Diff_Dis))) this%Diff => Diff_Dis

      end subroutine SetDisDiffScheme

      !> 设置2D网格
      subroutine SetGrid(This, Grid)

          implicit none
          class(solver_type), intent(inout) :: this
          type(grid2d_type), target, intent(in) :: Grid !< 二维贴体网格

          if(.not. (associated(this%Grid, Grid))) this%Grid => Grid

      end subroutine SetGrid

      !> 设置2D基本流
      subroutine SetBaseFlow(this, BF)

          implicit none
          class(solver_type), intent(inout) :: this
          type(Baseflow2D_type), target, intent(in) :: BF !< 2D基本流

          if(.not. (associated(this%BaseFlow, BF))) this%BaseFlow => BF

      end subroutine SetBaseFlow

      !> 设置2D曲率
      subroutine SetCurvature(this, Curvature)

          implicit none
          class(solver_type), intent(inout) :: this
          type(curvature_2d_type), target, intent(in) :: Curvature !< 2D区率

          if(.not. (associated(this%Curvature, Curvature))) &
          &   this%Curvature => Curvature

      end subroutine SetCurvature

      !> 设置当前流向站位
      subroutine set_iloc(this, iloc)

          implicit none
          class(solver_type), intent(inout) :: this
          integer, intent(in) :: iloc !< 当前流向站位

          this%iloc=iloc

      end subroutine set_iloc

      !> 设置PSE的边界条件
      subroutine SetBCType(this, bctype)

          implicit none
          class(solver_type), intent(inout) :: this
          integer, intent(in) :: bctype(5, 2) !< 边界条件类型

#IFDEF DEBUG
          print*, BCtype
          pause
#ENDIF
          this%bctype=bctype

      end subroutine SetBCType


      !> 输出当前站位序号
      subroutine PrintILoc(this, iloc)

          implicit none
          class(solver_type), intent(in) :: this
          integer, intent(in) :: iloc

          write(*, *)"The location index now is", iloc

      end subroutine PrintILoc

end module mod_solver
