!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_context.f90
!> @file
!> @breif 计算容器类文件.
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
!> @brief 计算容器模块.
!!
!  REVISION HISTORY:
!  2015-07-25 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-07-25
!------------------------------------------------------------------------------
module mod_context

    use mod_grid
    use mod_baseflow
    use mod_dis
    use mod_difference
    use mod_curvature_2d
    use mod_dis_wavenum
    use mod_cfgio_adapter

    implicit none

    private

    !> 解算器基本容器
    !!
    !! 定义了基本的解算器类。目前包括LPSE，伴随LPSE，以及由LPSE求伴随LPSE和由伴随LPSE求LPSE。
    !!
    !! - 基本类型参数对应如下：
    !!   功能|描述
    !!   ----------|----------------------------------------
    !!   LST|线性稳定性理论分析
    !!   LPSE|计算线性LPSE
    !!   ALPSE|计算线性伴随LPSE
    !!   ALPSE_LPSE|先计算线性伴随LPSE，然后再以其为基础计算LPSE
    !!   LPSE_ALPSE|先计算线性LPSE，然后再以其为基础计算伴随LPSE
    type, public :: context_type

        private
        type(grid2d_type)         , pointer     :: grid                          !< 原网格\private
        type(grid2d_type)         , pointer     :: gridFixBody                   !< 计算网格(贴体网格)\private
        type(BaseFlow2D_type)     , pointer     :: BaseFlow                      !< 基本流\private
!        type(BaseFlow2D_Diff_type), pointer     :: dBaseFlow                     !< 基本流导数\private
        type(Difference_2D_type)  , pointer     :: diff_BF                       !<基本流差分基函数\private
        type(Difference_2D_type)  , pointer     :: diff_Dis                      !<扰动差分基函数\private
        type(Curvature_2d_type)   , pointer     :: curvature                     !<表面曲率\private
        type(Dis_type)            , allocatable :: dis                           !<扰动\private
        type(Dis_type)            , allocatable :: adjointDis                    !<伴随扰动\private
        class(*)                  , allocatable :: solver                        !<解算器类型\private

        Contains

            Procedure :: Initialize                                               !< 初始化解算器
            procedure :: Solve                                                    !< 给定求解类型进行求解

            procedure, private :: LPSE                                            !<LPSE(线性抛物化稳定性方程)解算器
            procedure, private :: ALPSE                                           !<伴随LPSE解算器
            procedure, private :: LST                                             !<ST(线性稳定性理论)解算器
            procedure, private :: ALPSE_LPSE                                      !<先用伴随LPSE求解，然后再用LPSE求解
            procedure, private :: LPSE_ALPSE                                      !<先用LPSE求解，然后再用伴随LPSE求解
            procedure, private :: NPSE

    end type context_type

    contains

    !> 伴随LPSE解算器
    subroutine ALPSE(this)

        use mod_apse
        implicit none
        class(context_type), intent(inout) :: this

        if (.not. allocated(this%AdjointDis)) allocate(this%AdjointDis)
        call this%AdjointDis%Create(this%grid)
        allocate(apse_type :: this%solver)
        select type(ALPSE => this%solver)
        type is (apse_type)
            call ALPSE%Initialize(this%Grid, this%BaseFlow, this%curvature, bctype)

            call ALPSE%SetWholeDis(this%Adjointdis)
            call ALPSE%SetDisDiffScheme(this%Diff_Dis)
            !call ALPSE%SetBCType(bctype)
            call ALPSE%SetIstartEnd(Istart, Iend)
            !! Initial Dis X=X1 in APSE
            call ALPSE%SetDisInlet_APSE(wavenum_init(1), Iend)
            !! ALPSE Advance Process
            call ALPSE%Solve()
            call ALPSE%Print(prefix)

        end select

    end subroutine ALPSE

    !> 先用伴随LPSE求解，然后再使用LPSE求解
    subroutine ALPSE_LPSE(this)

        use mod_apse_lpse
        implicit none
        class(context_type), intent(inout) :: this

        allocate(this%AdjointDis)
        call this%AdjointDis%Create(this%grid)
        allocate(this%dis)
        call this%Dis%Create(this%grid)
        allocate(apse_lpse_type :: this%solver)
        select type(APSE_LPSE => this%solver)
        type is (apse_lpse_type)
          call APSE_LPSE%Initialize(this%Grid, this%BaseFlow, &
              &                 this%curvature, bctype)
          call APSE_LPSE%SetWholeDis(this%Dis, this%adjointDis)
          call APSE_LPSE%SetDisDiffScheme(this%Diff_Dis)
          !call LPSE%SetBCType(bctype)
          call APSE_LPSE%SetIstartEnd(Istart, Iend)
          !! Initial Dis X=0 in PSE
          call APSE_LPSE%SetDisInlet(wavenum_init(1), Iloc)
          !! LPSE Advance Process
          call APSE_LPSE%Solve()
          call APSE_LPSE%Print(prefix)
        end select

   end subroutine ALPSE_LPSE

    !> 初始化解算器
    !!
    !! 给定网格、流场、曲率、差分格式基函数等基本信息，为下一步计算做准备
    !! \retval NULL
    subroutine Initialize(this)

        implicit none
        class(context_type), intent(inout) :: this
        integer, parameter                :: DIFFBF=1000, DIFFDIS=1002

        !! Create container
        allocate(this%Grid)
        call this%Grid%CreateFromPLOT3D(gridfile)
        if(.not. associated(this%GridFixBody)) &
        &   this%GridFixBody => this%Grid%TransBodyFix()
        allocate(this%Diff_BF, this%Diff_Dis)
        call this%Diff_BF%InitDiff(this%GridFixBody, 4, 4, DiffBF)
        call this%Diff_Dis%InitDiff(this%GridFixBody, 1, 4, DiffDis)
        allocate(this%BaseFlow)
        call this%BaseFlow%CreateFromPLOT3D(flowfile)
!        call this%DBaseFlow%CreateDiff(this%Baseflow, this%Diff_BF)
        allocate(this%Curvature)
        if(curvaturefile=='null') then
            call this%Curvature%CreateWithNoCurvature(this%grid%GetInSize())
        else
            call this%Curvature%CreateFromPlot3D(curvaturefile)
        endif

    end subroutine Initialize

    !> LPSE(线性抛物化稳定性方程)解算器
    subroutine LPSE(this)

        use mod_lpse
        implicit none
        class(context_type), intent(inout) :: this

        allocate(this%dis)
        call this%Dis%Create(this%grid)
        allocate(lpse_type :: this%solver)
        select type(LPSE => this%solver)
        type is (lpse_type)
            call LPSE%Initialize(this%Grid, this%BaseFlow, this%curvature, bctype)
            call LPSE%SetWholeDis(this%Dis)
            call LPSE%SetDisDiffScheme(this%Diff_Dis)
            !call LPSE%SetBCType(bctype)
            call LPSE%SetIstartEnd(Istart, Iend)
            !! Initial Dis X=0 in PSE
            call LPSE%SetDisInlet_PSE(wavenum_init(1), Istart)
            !! LPSE Advance Process
            call LPSE%Solve()
            call LPSE%Print(prefix)
        end select

    end subroutine LPSE

    !> LPSE(线性抛物化稳定性方程)解算器
    subroutine NPSE(this)

        use mod_npse
        implicit none
        class(context_type), intent(inout) :: this

        allocate(this%dis)
        call this%Dis%Create(this%grid)
        allocate(npse_type :: this%solver)
        select type(NPSE => this%solver)
        type is (npse_type)
            call NPSE%Initialize(this%Grid, this%BaseFlow, this%curvature, bctype)
            call NPSE%SetDim(mdim, ndim)
            call NPSE%SetWholeDis(this%Dis)
            call NPSE%SetDisDiffScheme(this%Diff_Dis)
            !call LPSE%SetBCType(bctype)
            call NPSE%SetIstartEnd(Istart, Iend)
            !! Initial Dis X=0 in PSE
            call NPSE%SetDisInlet_NPSE(wavenum_init, m_index, n_index, amp, istart)
            !! LPSE Advance Process
            call NPSE%Solve()
            call NPSE%Print(prefix)
        end select

    end subroutine NPSE

    !> 先用LPSE求解，再用伴随LPSE求解
    subroutine LPSE_ALPSE(this)

        use mod_lpse_apse
        implicit none
        class(context_type), intent(inout) :: this

        allocate(this%dis)
        call this%Dis%Create(this%grid)
        allocate(this%AdjointDis)
        call this%AdjointDis%Create(this%grid)
        allocate(lpse_apse_type :: this%solver)
        select type(LPSE_APSE => this%solver)
        type is (lpse_apse_type)
          call LPSE_APSE%Initialize(this%Grid, this%BaseFlow, this%curvature, bctype)
          call LPSE_APSE%SetWholeDis(this%Dis, this%adjointDis)
          call LPSE_APSE%SetDisDiffScheme(this%Diff_Dis)
          !call LPSE%SetBCType(bctype)
          call LPSE_APSE%SetIstartEnd(Istart, Iend)
          !! Initial Dis X=0 in PSE
          call LPSE_APSE%SetDisInlet(wavenum_init(1), Iloc)
          !! LPSE Advance Process
          call LPSE_APSE%Solve()
          call LPSE_APSE%Print(prefix)
        end select

    end subroutine LPSE_ALPSE

    !> LST(线性稳定性理论)解算器
    subroutine LST(this)

        use mod_lst
        implicit none
        class(context_type), intent(inout) :: this

        allocate(this%dis)
        call this%Dis%Create(this%Grid)
        allocate(lst_type :: this%solver)
        select type (LST => this%solver)
        type is (lst_type)
          call LST%Initialize(this%Grid, this%BaseFlow, this%curvature, bctype)
          call LST%SetWholeDis(this%Dis)
          call LST%SetDisDiffScheme(this%Diff_Dis)
          !call LPSE%SetBCType(bctype)
          call LST%SetIstartEnd(Istart, Iend)
          !! Initial Dis X=0 in PSE
          !! LPSE Advance Process
          call LST%Solve(iloc, wavenum_init(1))
          call LST%Print(prefix)
        end select

    end subroutine LST

    !> 给定求解类型进行求解
    !!
    !! 根据具体的求解类型，选用不同的解算器进行求解。
    !> \param[in] solver_kind 计算器类型
    !> \retval NULL
    subroutine Solve(this, solver_kind)

        implicit none
        class(context_type), intent(inout) :: this
        integer :: solver_kind
        integer, parameter :: LPSE=0, ALPSE=1, LST=2, ALPSE_LPSE=3, LPSE_ALPSE=4
        integer, parameter :: NPSE=5
        select case (solver_kind)
            case (LPSE)
                call this%LPSE()
            case (ALPSE)
                call this%ALPSE()
            case (ALPSE_LPSE)
                call this%ALPSE_LPSE()
            case (LPSE_ALPSE)
                call this%LPSE_ALPSE()
            case (NPSE)
                call this%NPSE()
            case (LST)
                call this%LST()
            case default

        end select

    end subroutine Solve

end module mod_context
