!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_pse.f90
!> @file
!> @breif PSE求解器文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: apse
!> @breif 伴随PSE求解器模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2018-06-25 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2018-06-22
!------------------------------------------------------------------------------
module mod_apse

    !use mod_grid
    !use mod_baseflow
    use mod_difference
    !use mod_lpse_bc
    !use mod_curvature_2d
    use mod_dis
    use mod_lpse_dis_normal
    use mod_baseflow_org
    use mod_lst_eqn_IR, only: lst_eqn_IR_type
    use mod_solver
    use penf, only: R_P

    use mod_dis_wavenum
    use mod_apse_eqn


    implicit none

    private

    !> APSE求解器类
    type, extends(solver_type), public :: apse_type

        private
        type(apse_eqn_type) :: apse !<伴随PSE求解器
        type(dis_type), pointer :: Dis !< 二维扰动场

        contains

        procedure :: SetDisInlet_APSE=> SetDisInlet_IR_APSE !< 设置入口扰动
        procedure :: SetWholeDis !< 设置全局求解对应的扰动

        procedure :: Print !< 输出APSE的计算结果
        procedure :: Solve=>Solve_APSE !< 推进求解APSE

        procedure, private :: finalize => finalize_lpse !<析构函数


    end type apse_type

    contains

    !> 设置求解对应的扰动
    subroutine SetWholeDis(this, Dis)

        implicit none
        class(apse_type), intent(inout) :: This
        type(dis_type), target, intent(in) :: Dis !< 二维扰动场

        if(.not. (associated(this%Dis, Dis))) this%Dis => Dis

    end subroutine SetWholeDis

    !> ALPSE求解过程.
    subroutine Solve_APSE(this)

        use mod_BF_normal
        use mod_dis_normal
        use mod_local_normal_coordinate
        use mod_difference

        implicit none
        class(apse_type), intent(inout) :: this
        type(lpse_dis_normal_type) :: DisNorm, DisNormBack
        type(BF_normal_type), save :: BFNorm
        type(local_normal_coordinate_type), save :: NormCoord
        type(Norm_Diff_Coef_type), save :: NormCoef
        real(R_P), allocatable, save :: Eta(:)
        integer :: iloc, iend
#ifdef DEBUG
        type(dis_wavenum_lpse_type) :: wavenum
#endif


        if(this%iend==-1)then
          iend=this%Grid%GetInSize()
          this%iend=iend
        end if

        if(.not. BFNorm%HasCreated()) then
          call BFNorm%Create(this%Grid%GetJnSize())
        endif

        if(.not. NormCoord%HasCreated()) then
          call NormCoord%Create(this%Grid%GetJnSize())
        endif

        if(.not. allocated(eta)) then
          allocate(Eta(this%Grid%GetJnSize()))
        endif

        call this%apse%SetSolver(method=1001)
        call this%apse%Create(this%Grid%GetJnSize())

        do iloc=This%Iend-1, this%IStart, -1

          call this%PrintILoc(iloc)
          call this%Dis%Get(iloc, DisNorm)
#ifdef DEBUG
          WaveNum= DisNorm%GetWaveNum()
          print*, Wavenum%getAlpha()
          pause
#endif

          call this%Dis%Get(iloc+1, DisNormBack )
#ifdef DEBUG
          WaveNum= DisNormBack%GetWaveNum()
          print*, Wavenum%getAlpha()
          pause
#endif
          call BFNorm%Set(iloc, this%BaseFlow, this%Diff)
          call NormCoord%SetLameFromGrid(this%grid, iloc)
          call NormCoord%SetCurvature(this%Curvature%GetPoint(iloc))
          call this%Diff%GetLocalNormDiffCoef(iloc, NormCoef)
          call this%grid%Get_iy(iloc, Eta)

          call this%apse%SolvePSE(DisNorm, iloc, DisNormBack, &
          &   BFNorm, NormCoord, NormCoef, Eta, this%bctype)

          call this%Dis%Set(iloc, DisNorm)
        end do

    end subroutine Solve_APSE

    !> 输出PSE的计算结果
    subroutine Print(this, fn_surf)

        use stringifor
        implicit none
        class(apse_type), intent(in) :: this
        type(string), intent(in) :: fn_surf !< 文件名前缀
        type(string) :: fn_lpse

        fn_lpse=fn_surf//'_APSE'

        call this%Grid%Print(fn_lpse)
        call this%Dis%Print(fn_lpse)

    end subroutine Print

    ! !> 给入口边界条件(muller法LST).
    ! !!
    ! !! 先求解入口处的LST问题,作为PSE入口边界条件
    ! subroutine SetDisInlet_muller(this, wavenum)
    !
    !     use mod_dis_wavenum, only: dis_wavenum_type, dis_wavenum_lpse_type
    !
    !     implicit none
    !     class(pse_type), intent(inout) :: this
    !     class(dis_wavenum_type), intent(in) :: wavenum !< 入口扰动色散特征
    !     type(dis_wavenum_lpse_type) :: wavenum_pse
    !     type(lpse_dis_normal_type) :: DisNorm
    !     type(lst_type) :: LST
    !     integer, parameter :: SPATIAL_LST=1, TEMPORAL_LST=2
    !     logical, parameter :: IS_PARALLEL=.false.
    !     complex(R_P) :: sigma(7)
    !
    !     call LST%pse_type()
    !     call LST%SetCurvature(this%Curvature%GetPoint(this%iloc))
    !     call LST%SetMode(SPATIAL_LST)
    !     call LST%SetEta(this%jn, this%Eta)
    !     call LST%SetBaseflow(this%BFNorm)
    !
    !     call wavenum_pse%Set(wavenum)
    !     call this%Dis%Get(This%Istart, DisNorm)
    !     call DisNorm%SetWavenum(wavenum_pse)
    !     call LST%Solve(DisNorm)
    !     call DisNorm%SetILoc(this%Istart)
    !
    !     wavenum_pse=DisNorm%GetWaveNum()
    !     sigma=wavenum_pse%getAlpha()
    !     call DisNorm%SetSigma(sigma)
    !     call this%Dis%set(this%iloc, DisNorm)
    !
    ! end subroutine SetDisInlet_muller

    ! !> 初始化伴随PSE, 给入口边界条件(反幂法求伴随LST).
    ! !!
    ! !! 先求解入口处的伴随LST问题,作为伴随PSE入口边界条件
    ! subroutine SetDisInlet_IR_APSE(this, wavenum, iloc)
    !
    !     implicit none
    !     class(apse_type), intent(inout) :: this
    !     class(dis_wavenum_type), intent(in) :: wavenum !< 初始扰动色散特征
    !     integer, intent(in) :: iloc
    !     type(dis_wavenum_lpse_type) :: wavenum_pse
    !     type(lpse_dis_normal_type) :: DisNorm
    !     type(lst_IR_type) :: LST
    !     complex(R_P) :: sigma(7)
    !     logical, parameter :: IS_PARALLEL=.false.
    !
    !     this%iloc=iloc
    !
    !     ! call this%SetIloc(this%IEnd)
    !     ! call this%SetDisDiffCoef()
    !     ! call this%SetBFOP(IS_PARALLEL)
    !     !call this%Initialize(This%Iend)
    !     call LST%Create(this%Grid%GetJnSize())
    !     call LST%Set(this%iloc, this%Grid, this%BaseFlow, &
    !                 this%Curvature, this%Diff)
    !     call LST%SetInitialWavenum(wavenum)
    !     call LST%Solve()
    !
    !     call this%Dis%Get(This%Iend, DisNorm)
    !     call LST%GetShapefunAdjoint(DisNorm)
    !     wavenum_pse=LST%GetWavenum()
    !     sigma=wavenum_pse%getAlpha()
    !     call DisNorm%SetILoc(this%IEnd)
    !     call DisNorm%SetSigma(sigma)
    !     call wavenum_pse%Set(wavenum)
    !     call DisNorm%SetWavenum(wavenum_pse)
    !
    !     call this%Dis%set(this%iloc, DisNorm)
    !     call this%Dis%Print(this%iloc, this%Grid)
    !
    ! end subroutine SetDisInlet_IR_APSE

    !> 给入口边界条件(反幂法求LST).
    !!
    !! 先求解入口处的伴随LST问题,作为APSE入口边界条件
    subroutine SetDisInlet_IR_APSE(this, wavenum, iloc_in)

        use mod_dis_wavenum, only: dis_wavenum_type, dis_wavenum_lpse_type
        use mod_BF_normal
        use mod_dis_normal
        use mod_local_normal_coordinate
        use mod_difference

        implicit none
        class(apse_type), intent(inout) :: this
        class(dis_wavenum_type), intent(in) :: wavenum !< 初始扰动色散特征
        integer, optional, intent(in) :: iloc_in
        type(dis_wavenum_lpse_type) :: wavenum_pse
        type(lpse_dis_normal_type) :: DisNorm
        type(BF_normal_type), save :: BFNorm
        type(local_normal_coordinate_type), save :: NormCoord
        type(Norm_Diff_Coef_type), save :: NormCoef
        real(R_P), allocatable, save :: Eta(:)
        logical, parameter :: ISADJOINT=.True.
        integer :: iloc

        type(lst_eqn_IR_type) :: LST

        if(present(iloc_in)) then
          this%iloc=iloc_in
        else
          this%iloc=this%Iend
        endif

        iloc=this%iloc

        if(.not. BFNorm%HasCreated()) then
          call BFNorm%Create(this%Grid%GetJnSize())
        endif

        if(.not. NormCoord%HasCreated()) then
          call NormCoord%Create(this%Grid%GetJnSize())
        endif

        if(.not. allocated(eta)) then
          allocate(Eta(this%Grid%GetJnSize()))
        endif

        call LST%Create(this%Grid%GetJnSize())

        call this%Dis%Get(iloc, DisNorm)

        call wavenum_pse%set(wavenum)
        call DisNorm%SetWavenum(wavenum_pse)

        call BFNorm%Set(iloc, this%BaseFlow, this%Diff)
        call NormCoord%SetLameFromGrid(this%grid, iloc)
        call NormCoord%SetCurvature(this%Curvature%GetPoint(iloc))
        call this%Diff%GetLocalNormDiffCoef(iloc, NormCoef)
        call this%grid%Get_iy(iloc, Eta)

        call LST%SolveLST(DisNorm, iloc, &
        &   BFNorm, NormCoord, Eta, NormCoef, ISADJOINT)


        print*, DisNorm%GetAlpha()
        print*, DisNorm%GetBeta()
        print*, DisNorm%GetOmega()
!        pause

        call this%Dis%Set(iloc, DisNorm)

        call this%Dis%Print(this%iloc, this%Grid)

    end subroutine SetDisInlet_IR_APSE

    !
    ! !> 输出当前站位序号
    ! subroutine PrintILoc(this)
    !
    !     implicit none
    !     class(apse_type), intent(in) :: this
    !
    !     write(*, *)"The location index now is", this%iloc
    !
    ! end subroutine PrintILoc

    !> 析构函数
    subroutine finalize_lpse(this)

        implicit none
        class(apse_type), intent(inout) :: this

        this%iloc=0
        this%dis=>null()

    end subroutine finalize_lpse

end module mod_apse
