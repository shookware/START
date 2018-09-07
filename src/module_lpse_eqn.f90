!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_lpse_eqn.f90
!> @file
!> @breif LPSE方程求解器文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: lpse_eqn
!> @breif LPSE方程求解器模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-05 - Initial Version
!  2018-06-25 - the wrapper is modified and the type is splitted into type_pse
!               and type_lpse. It makes that it is convinient to impliment the
!               the other method.
!> @author
!> Liu Jianxin
!> \date 2017-08-05
!------------------------------------------------------------------------------
module mod_lpse_eqn

    ! use mod_grid
    !use mod_baseflow
    use mod_local_normal_coordinate
    use mod_difference
    use mod_BF_normal
    use mod_lns_OP_normal
    use mod_lpse_dis_OP_normal
    use mod_lpse_bc
    use mod_dis_wavenum
    !use mod_curvature_2d
    !use mod_dis
    use mod_lpse_dis_normal
    use mod_pardiso_adapter
    use mod_dis_flux
    use mod_baseflow_org
    use mod_parameter, only: CPLI
    ! use mod_lst, only: lst_type
    ! use mod_lst_IR, only: lst_IR_type
!    use mod_solver
    use penf, only: R_P

    implicit none

    private

    !> LPSE求解器类
    type, public :: lpse_eqn_type

        private
        type(BF_normal_type), Pointer, private :: BFNorm !<基本流流向某站位法向分布
        type(lpse_dis_normal_type), Pointer, private :: DisNorm !<扰动流向某站位法向分布
        type(lpse_dis_normal_type), Pointer, private :: DisNormFront !<扰动流向前一站位法向分布
        type(local_normal_coordinate_type), pointer, private :: NormCoord !< 某流向展向站位的局部坐标系
        type(Norm_Diff_Coef_type), pointer, private :: NormCoef
        type(lpse_bf_op_normal_type), private :: BFOPNorm !< LNS方程在当地流向展向位置沿法向分布的系数算子
        type(lpse_dis_op_normal_type), private :: DisOPNorm !< PSE方程的线性系数算子沿法向分布
!        type(difference_2D_type), pointer, private :: Diff !< 整个二维场扰动离散的拉格朗日基函数
!        type(grid2d_type), pointer, private :: Grid !< 二维场网格
!        type(Baseflow2D_type), pointer, private :: BaseFlow !< 二维基本流
!        type(Baseflow2D_Diff_type), pointer, private :: DBaseFlow !< 二维基本流导数
!        type(curvature_2d_type), pointer, private :: Curvature !< 二维曲率
!        type(dis_type), pointer, private :: Dis !< 二维扰动场
        !type(lpse_bc_type), private :: LpseBC(2) !< LPSE边界条件算子信息
        integer, private :: iloc !< 当地流向站位编号
        integer, private :: jn !< 法向网格数
        integer, private :: bctype(5, 2) !< PSE方程边界条件类型
        type(lpse_bc_type), private :: LpseBC(2) !<PSE方程边界条件矩阵及右端项
        real(R_P), private :: Coef(2) !< 流向扰动离散的拉格朗日基函数
        real(R_P), allocatable, private :: Eta(:) !< 当地位置法向网格分布
        class(*), allocatable, private :: solver_ap !< 线性方程求解器
        integer, private :: solver_kind !< 线性代数求解器类型
        complex(R_P), allocatable, private :: nonlinearTerm(:, :) !<非线性项
        logical :: NonlinearSolver !< 标记是否计算非线性问题

        contains

        procedure :: Create !< 创建LPSE类型, 分配内存
        procedure :: SetSolver !< 设置线性代数求解器相关信息
        !procedure :: SetDisInlet=> SetDisInlet_ir !< 设置入口扰动
        procedure :: SolvePSE !< 推进求解PSE
        !procedure :: SetIStartEnd !< 设置起始和终止站位序号
        !procedure :: SetDisDiffScheme !< 设置扰动离散格式
        !procedure :: SetGrid !< 设置2D网格
        !procedure :: SetBaseFlow !< 设置2D基本流
!        procedure :: SetDBaseFlow !< 设置2D基本流的导数
        !procedure :: SetCurvature !< 设置2D曲率
        !procedure :: SetWholeDis !< 设置全局求解对应的扰动
        procedure :: SetBCType !< 设置边界条件类型
        !procedure :: Print !< 输出PSE的计算结果

        generic, private :: SetMatA => set_matA_lpse, set_matA_lpse_bsr !< 设置PSE方程离散化的线性方程组的系数矩阵\f$M\f$
        procedure, private :: SetIloc => set_iloc !< 设置当前流向站位
        procedure, private :: SetBFOP => set_BF_OP_Norm !< 设置LNS系数算子沿法向分布
!        procedure, private :: SetDisDiffCoef !<设置流向扰动离散基函数以及法向网格分布
        !procedure, private :: Initialize !< 初始化PSE
        procedure, private :: SolveWithAlpha !<LPSE求解器
!        procedure, private :: SolveStart !<在当地求解一个局部平行的LPSE问题作为LPSE初始入口条件的修正过程
        procedure, private :: finalize => finalize_lpse !<析构函数
!        procedure, private :: PrintILoc !< 输出当前站位序号
        procedure, private :: SolveDisNorm !< 求解当前站位扰动
        procedure, private :: SolveAlpha !< 求解流向波数\f$\alpha\f$
        procedure, private :: Norm_DisNorm !< 求扰动形函数范数
        procedure, private :: Norm_DxDisNorm !< 求扰动形函数的流向导数范数
        procedure, private :: SolveDisNorm_directive !< 求解当前站位扰动(直接法)
        !procedure, private :: SolveDisNorm_iterative !< 求解当前站位扰动(迭代法)
        procedure, private :: SetDisOP => set_dis_op !< 设置LPSE的线性算子系数
        procedure, private :: dsigma !< 求解LPSE物理流向复波数的形函数的修正量
        procedure, private :: ComputeSigma !< 求解LPSE的物理流向复波数
        !procedure, private :: PrintSigma => print_sigma !< 输出物理流向复波数
        procedure, private :: SetRHS => set_rhs_lpse !< 设置PSE方程离散化的线性方程组的右端项
        procedure, private :: GetX => get_X_lpse !< 获得PSE方程离散化的线性方程组的解, 并将其放入扰动形函数中
        !procedure, private :: InitializeStart !< 对作为入口边界条件的修正过程的初始化
        procedure, private :: SolveWithoutAlpha !< 求解PSE方程,但不对波数\f$\alpha\f$进行迭代
        procedure, private :: set_matA_lpse  !< 设置PSE方程离散化的线性方程组的系数矩阵\f$M\f$ (COO格式)
        procedure, private :: set_matA_lpse_bsr !< 设置PSE方程离散化的线性方程组的系数矩阵\f$M\f$(BSR格式)

    end type lpse_eqn_type

    !> 积分函数
    interface Intergal
        !> 实数积分函数
        module procedure dIntergal
        !> 复数积分函数
        module procedure zIntergal

    end interface Intergal

    contains

    !> 设置LPSE的边界条件
    subroutine SetBCType(this, bctype)

        implicit none
        class(lpse_eqn_type), intent(inout) :: this
        integer, intent(in) :: bctype(5, 2) !< 边界条件类型

        this%bctype=bctype

    end subroutine SetBCType

    !> 设置线性代数求解器相关信息
    subroutine SetSolver(this, method)

        use mod_gmres_adapter
        !use mod_gmres
        implicit none
        class(lpse_eqn_type),intent(inout) :: this
        integer, intent(in), optional :: method    !< 线性代数求解器求解方法
        integer, parameter ::  DIRECT_SOLVER=1001, ITERATIVE_SOLVER=1002
        integer, parameter :: NORMAL=0, TRANSPOSE=2
        complex(R_P) :: solutionVec(this%jn*5), RHS(this%jn*5)

        if(present(method)) then
            this%solver_kind=method
        else
            this%solver_kind=DIRECT_SOLVER
        end if

        select case (this%solver_kind)
            case (DIRECT_SOLVER)
                allocate(pardiso_adapter_type::this%solver_ap)
                !allocate(gmres_adapter_type::this%solver_ap)
            case (ITERATIVE_SOLVER)
!                allocate(gmres_adapter_type::this%solver_ap)
                 stop 'iterative solver is not contained now!'
        end select
        select type (solver_ap=>this%solver_ap)
        type is (pardiso_adapter_type)
            call solver_ap%Initialize(NORMAL)
        type is (gmres_adapter_type)
            !call solver_ap%initial(this%jn*5, SolutionVec, RHS, this)
            stop 'iterative solver is not contained now!'
        end select

    end subroutine SetSolver


    ! !> 输出PSE的计算结果
    ! subroutine Print(this, fn_surf)
    !
    !     use stringifor
    !     implicit none
    !     class(lpse_eqn_type), intent(in) :: this
    !     type(string), intent(in) :: fn_surf !< 文件名前缀
    !     type(string) :: fn_lpse
    !
    !     fn_lpse=fn_surf//'_LPSE'
    !     call this%Grid%Print(fn_lpse)
    !     call this%Dis%Print(fn_lpse)
    !     !call this%PrintSigma(fn_lpse)
    !
    ! end subroutine Print

    !> PSE推进求解
    subroutine SolvePSE(this, DisNorm, iloc, DisNormFront, &
    &               BFNorm, NormCoord, NormCoef, Eta, BCtype, &
    &               isConvergence, NonlinearTerm, IsWithAlpha)

        implicit none
        class(lpse_eqn_type), intent(inout) :: this
        type(lpse_dis_normal_type), target, intent(inout) :: DisNorm
        integer, intent(in) :: iloc
        type(lpse_dis_normal_type), target, intent(in) :: DisNormFront
        type(BF_normal_type), target, intent(in) :: BFNorm
        type(local_normal_coordinate_type), target, intent(in) :: NormCoord
        type(Norm_Diff_Coef_type), target, intent(in) :: NormCoef
        real(R_P), intent(in) :: Eta(:)
        integer, intent(in) :: BCtype(5, 2)
        complex(R_P), optional, intent(in) :: nonlinearTerm(:, :)
        logical, optional, intent(in) :: isWithAlpha
        logical, intent(out), optional :: isConvergence

        logical, parameter :: IS_PARALLEL=.false.
        integer :: times

        times=0

        call this%SetIloc(iloc)

        this%BFNorm=>BFNorm
        this%DisNormFront=>DisNormFront
        this%NormCoord=>NormCoord
        this%NormCoef=>NormCoef
        this%DisNorm=>DisNorm
#IFDEF DEBUG
        print*, 'subroutine solve lpse'
        print*, this%DisNormFront%jn, this%DisNorm%jn
        pause
#ENDIF
        this%Coef=NormCoef%Coef_di

        if(.not. allocated(this%Eta)) then
            allocate(this%Eta, source=Eta)
        else
            this%Eta=Eta
        endif


        call this%SetBCType(BCtype)

        call this%SetBFOP(IS_PARALLEL)
        call this%DisOPNorm%SetCoef(this%Coef)
        if(.not. present(nonlinearTerm)) then
          this%nonlinearTerm=0.0d0
          this%NonlinearSolver=.False.
        else
          this%nonlinearTerm=nonlinearTerm
          this%NonlinearSolver=.True.
        endif
        if(.not. present(isWithAlpha) .or. isWithAlpha) then
          call this%SolveWithAlpha(times)
        else
          call this%SolveWithoutAlpha(times)
        endif

        if(present(isConvergence)) then
          if (times>1) then
            isConvergence=.False.
          else
            isConvergence=.True.
          endif
        endif


    end subroutine SolvePSE

    ! !> 设置起始和终止的站位序号
    ! subroutine SetIstartend(this, Istart, Iend)

    !     implicit none
    !     integer, intent(inout) :: Istart !< 起始流向站位序号
    !     integer, intent(inout) :: Iend !< 终止流向站位序号
    !     class(lpse_eqn_type), intent(inout) :: this

    !     this%Istart=Istart; this%Iend=Iend
    !     if(this%iend==-1 .or. (this%iend>this%Grid%GetInSize()))then
    !       this%iend=this%Grid%GetInSize()
    !       iend=this%iend
    !     end if
    !     if(this%istart==-1 .or. (this%istart<1))then
    !       this%istart=1
    !       istart=1
    !     end if

    ! end subroutine SetIstartend

    ! !> 设置扰动离散格式
    ! subroutine SetDisDiffScheme(this, Diff_Dis)

    !      implicit none
    !      class(lpse_eqn_type), intent(inout) :: this
    !      type(difference_2D_type), target, intent(in) :: Diff_Dis !< 二维场扰动离散的拉格朗日基函数

    !      if(.not. (associated(this%Diff, Diff_Dis))) this%Diff => Diff_Dis

    ! end subroutine SetDisDiffScheme

    ! !> 设置2D网格
    ! subroutine SetGrid(This, Grid)

    !     implicit none
    !     type(grid2d_type), target, intent(in) :: Grid !< 二维贴体网格
    !     class(lpse_eqn_type), intent(inout) :: this

    !     if(.not. (associated(this%Grid, Grid))) this%Grid => Grid

    ! end subroutine SetGrid

    ! !> 设置2D基本流
    ! subroutine SetBaseFlow(this, BF)

    !     implicit none
    !     type(Baseflow2D_type), target, intent(in) :: BF !< 2D基本流
    !     class(lpse_eqn_type), intent(inout) :: this

    !     if(.not. (associated(this%BaseFlow, BF))) this%BaseFlow => BF

    ! end subroutine SetBaseFlow

    !!> 设置2D基本流的导数
    !subroutine SetDBaseFlow(this, DBF)
    !
    !    implicit none
    !    type(Baseflow2D_Diff_type), target, intent(in) :: DBF !<2D基本流导数
    !    class(lpse_eqn_type), intent(inout) :: this
    !
    !    if(.not. (associated(this%DBaseFlow, DBF))) this%DBaseFlow => DBF
    !
    !end subroutine SetDBaseFlow

    ! !> 设置2D曲率
    ! subroutine SetCurvature(this, Curvature)

    !     implicit none
    !     type(curvature_2d_type), target, intent(in) :: Curvature !< 2D区率
    !     class(lpse_eqn_type), intent(inout) :: this

    !     if(.not. (associated(this%Curvature, Curvature))) &
    !     &   this%Curvature => Curvature

    ! end subroutine SetCurvature

    ! !> 设置求解对应的扰动
    ! subroutine SetWholeDis(this, Dis)

    !     implicit none
    !     type(dis_type), target, intent(in) :: Dis !< 二维扰动场
    !     class(lpse_eqn_type), intent(inout) :: This

    !     if(.not. (associated(this%Dis, Dis))) this%Dis => Dis

    ! end subroutine SetWholeDis

    !> 创建LPSE类型, 分配内存
    subroutine create(this, jn)

        implicit none
        class(lpse_eqn_type), intent(inout) :: this
        integer, intent(in) :: jn !< 法向点数
        integer, parameter :: NORMAL=0, TRANSPOSE=2
        integer :: eqkind

        this%jn=jn
        eqkind=0
        call this%BFOPNorm%Create(this%jn, eqkind)
        call this%DisOPNorm%Create(this%jn)
        allocate(this%nonlinearTerm(5, this%jn))
        this%NonlinearSolver=.False.

    end subroutine create
    !
    ! !> 初始化PSE.
    ! !!
    ! !! 利用设置条件完成PSE初始化
    ! subroutine Initialize(this)
    !
    !     implicit none
    !     class(lpse_eqn_type), intent(inout) :: this
    !
    !
    !
    ! end subroutine Initialize

!     !> 给入口边界条件(muller法LST).
!     !!
!     !! 先求解入口处的LST问题,作为PSE入口边界条件
!     subroutine SetDisInlet_muller(this, wavenum)
!
!         implicit none
!         class(lpse_eqn_type), intent(inout) :: this
!         class(dis_wavenum_type), intent(in) :: wavenum !< 入口扰动色散特征
!         type(dis_wavenum_lpse_type) :: wavenum_pse
!         type(lpse_dis_normal_type) :: DisNorm
!         type(lst_type) :: LST
!         integer, parameter :: SPATIAL_LST=1, TEMPORAL_LST=2
!         logical, parameter :: IS_PARALLEL=.false.
!         complex(R_P) :: sigma(7)
!
!         call This%Initialize(this%Istart)
!         call LST%Initialize()
!         call LST%SetCurvature(this%Curvature%GetPoint(this%iloc))
!         call LST%SetMode(SPATIAL_LST)
!         call LST%SetEta(this%jn, this%Eta)
!         call LST%SetBaseflow(this%BFNorm)
!
!         call wavenum_pse%Set(wavenum)
!         call this%Dis%Get(This%Istart, DisNorm)
!         call DisNorm%SetWavenum(wavenum_pse)
!         call LST%Solve(DisNorm)
!         call DisNorm%SetILoc(this%Istart)
!
!         wavenum_pse=DisNorm%GetWaveNum()
!         sigma=wavenum_pse%getAlpha()
!         call DisNorm%SetSigma(sigma)
!         call this%Dis%set(this%iloc, DisNorm)
! !        call this%Dis%Print(this%iloc, this%Grid)
! !        call this%InitializeStart()
!
! end subroutine SetDisInlet_muller
!
!     !> 给入口边界条件(反幂法求LST).
!     !!
!     !! 先求解入口处的LST问题,作为PSE入口边界条件
!     subroutine SetDisInlet_IR(this, wavenum)
!
!         implicit none
!         class(lpse_eqn_type), intent(inout) :: this
!         class(dis_wavenum_type), intent(in) :: wavenum !< 初始扰动色散特征
!         type(dis_wavenum_lpse_type) :: wavenum_pse
!         type(lpse_dis_normal_type) :: DisNorm
!         type(lst_IR_type) :: LST
!         complex(R_P) :: sigma(7)
!
!         call this%Initialize(this%istart)
!         call LST%Create(this%jn)
!         call LST%Set(this%iloc, this%Grid, this%BaseFlow, &
!                     this%Curvature, this%Diff)
!         call LST%SetInitialWavenum(wavenum)
!         call LST%Solve()
!
!         call this%Dis%Get(This%Istart, DisNorm)
!         call LST%GetShapefun(DisNorm)
!         wavenum_pse=LST%GetWavenum()
!         sigma=wavenum_pse%getAlpha()
!         call DisNorm%SetILoc(this%Istart)
!         call DisNorm%SetSigma(sigma)
!         call wavenum_pse%Set(wavenum)
!         call DisNorm%SetWavenum(wavenum_pse)
!
!         call this%Dis%set(this%iloc, DisNorm)
!         call this%Dis%Print(this%iloc, this%Grid)
!         !call this%InitializeStart()
!
!     end subroutine SetDisInlet_IR

    ! !> 对作为入口边界条件的修正过程的初始化
    ! subroutine InitializeStart(this)
    !
    !     implicit none
    !     class(lpse_eqn_type), intent(inout) :: this
    !     logical, parameter :: IS_PARALLEL=.true.
    !
    !     call this%SetBFOP(IS_PARALLEL)
    !     call this%SetDisDiffCoef()
    !     call this%SolveStart()
    !
    ! end subroutine InitializeStart

    !> 设置当前流向站位
    subroutine set_iloc(this, iloc)

        implicit none
        integer, intent(in) :: iloc !< 当前流向站位
        class(lpse_eqn_type), intent(inout) :: this

        this%iloc=iloc

    end subroutine set_iloc

    !> 设置LPSE边界条件
    subroutine set_lpse_bctype(this, bctype)

        implicit none
        class(lpse_eqn_type), intent(inout) :: this
        integer, intent(in) :: bctype(5, 2) !< 边界条件类型

        this%bctype=bctype

    end subroutine set_lpse_bctype

    !> 设置LNS系数算子沿法向分布
    subroutine set_BF_OP_Norm(this, isParallel)

        implicit none
!        type(grid2d_type), intent(in) :: grid
!        type(Baseflow2D_type), intent(in) :: Flow
!        type(difference_2D_type), intent(in) :: Diff
!        type(curvature_2d_type), intent(in) :: curvature
        class(lpse_eqn_type), intent(inout) :: this
        logical, intent(in) :: isParallel !< 流场是否考虑非平行性
!        type(Baseflow2D_Diff_type), intent(in) :: DiffFlow

        call this%BFOPNorm%Set(this%BFNorm, this%NormCoord, isParallel)

    end subroutine set_BF_OP_Norm

    !> 设置LPSE的线性算子系数
    subroutine set_dis_op(this, WaveNum)

        implicit none
        class(lpse_eqn_type), intent(inout) :: this
        type(dis_wavenum_lpse_type), intent(in) :: WaveNum !< 扰动色散特征
#IFDEF DEBUG
        ! print*, 'set_lpsebc'
        ! print*, this%bctype(:,1)
#ENDIF
        call this%LpseBC(1)%set(this%bctype(:, 1), this%BFOPNorm%GetPoint(1), wavenum, this%NormCoord%GetPoint(1))
#IFDEF DEBUG
        ! print*, 'set_lpsebc_wall'
#ENDIF
        call this%LpseBC(2)%set(this%bctype(:, 2), this%BFOPNorm%GetPoint(this%jn), wavenum, this%NormCoord%GetPoint(this%jn))
#IFDEF DEBUG
        ! print*, 'set_lpsebc_farfield'
        ! print*, this%BCtype(:,1)
        ! print*, this%BCtype(:,2)
#ENDIF
        call this%DisOPNorm%Set(this%BFOPNorm, this%NormCoord, wavenum, this%lpsebc)

    end subroutine set_dis_op

    !> 设置流向扰动离散基函数和法向网格分布
!     subroutine SetDisDiffCoef(this)
!
!         implicit none
!         class(lpse_eqn_type), intent(inout) :: this
! !        type(difference_2D_type), pointer, intent(in) :: Diff
!         real(R_P) :: Coef(2)
!
!         call this%DisOPNorm%SetCoef(this%Coef)
!
!         if( .not. allocated(this%eta)) &
!             &   allocate(this%Eta(this%jn))
!         !call this%Diff%GetEta(this%iloc, this%Eta)
!
!     end subroutine SetDisDiffCoef

    !> LPSE求解器
    subroutine SolveWithAlpha(this, times)

        use mod_parameter, only: EPS, EPS_REL
        implicit none
        class(lpse_eqn_type), intent(inout) :: this
        type(dis_wavenum_lpse_type) :: WaveNum
        type(lpse_dis_normal_type) :: DxDisNorm!, DisNormFront, DisNorm
        complex(R_P) :: AlphaOld, AlphaNew, AlphaFront
        complex(R_P) :: Dx_alpha
        integer, intent(out) :: times

        times=0

        !call This%Dis%Get(this%iloc-1, DisNormFront)
          !DisNormFront=this%DisNormFront
          !call this%PrintILoc()
          !!用iloc-1位置当DisNorm初值
          ! print*, 'alpha'
          ! print*, this%DisNorm%getAlpha()
          ! print*, this%DisNormFront%GetAlpha()
          ! pause

          !DisNorm=DisNormFront
          call this%DisNorm%SetILoc(this%iloc)!; call DisNorm%SetDiffIloc(this%iloc)
          WaveNum= this%DisNorm%GetWaveNum()


          AlphaNew=wavenum%getAlpha()
          AlphaFront=this%DisNormFront%GetAlpha()
          AlphaOld=0.0d0
          Dx_alpha=this%Coef(2)*AlphaNew+this%Coef(1)*AlphaFront

#IFDEF DEBUG
          write(*, *)this%Coef
          write(*, *) "AlphaOld=", (AlphaOld)
          write(*, *) "AlphaNew=", (AlphaNew)
          write(*, *) "AlphaFront=", (AlphaFront)
          print*, this%Coef(2)*AlphaNew+this%Coef(1)*AlphaFront
          pause
#ENDIF

          do while ( abs(AlphaNew-AlphaOld)/abs(AlphaNew) >=EPS_REL)
              call WaveNum%SetDxAlpha(Dx_alpha)
#IFDEF DEBUG
              print*, 'Dx_alpha=', Dx_Alpha
              pause
#ENDIF
              call this%SetDisOP(WaveNum)
              call this%SolveDisNorm(this%DisNorm, this%DisNormFront)
              call this%SolveAlpha(this%DisNorm, this%DisNormFront, DxDisNorm)
              WaveNum= this%DisNorm%GetWaveNum()
              AlphaOld=AlphaNew
              AlphaNew=WaveNum%getAlpha()
              Dx_alpha=WaveNum%GetDxAlpha()
!if(maxval(this%nonlinearTerm))
!              if(.not. this%NonlinearSolver) then
                write(*, *) "AlphaOld=", (AlphaOld)
                write(*, *) "AlphaNew=", (AlphaNew)
!              endif
#IFDEF DEBUG
              write(*, *) "AlphaFront=", (AlphaFront)
#ENDIF
!              if(.not. this%NonlinearSolver) then
                write(*, *) "AlphaErr=", abs(AlphaNew-AlphaOld)
!              endif
               times=times+1
          end do

!          print*, times
!          pause
!          if(this%NonlinearSolver) then
            write(*, *) "AlphaErr=", abs(AlphaNew-AlphaOld)
            write(*, *) "AlphaNew=", AlphaNew
!          endif

          !!算增长率
          call this%ComputeSigma(this%DisNorm, DxDisNorm)


          !this%DisNorm=DisNorm

    end subroutine SolveWithAlpha

    !>  求解LPSE方程,但不对波数\f$\alpha\f$进行迭代
    subroutine SolveWithoutAlpha(this, times)

        use mod_parameter, only: EPS
        implicit none
        class(lpse_eqn_type), intent(inout) :: this
        type(dis_wavenum_lpse_type) :: WaveNum !< 扰动色散特征
        type(lpse_dis_normal_type) :: DxDisNorm, DisNormFront, DisNorm
        type(lpse_dis_normal_type) :: DisNorm1
        complex(R_P) :: Dx_alpha, AlphaFront, AlphaNew
        integer, intent(out) :: times

        !DisNormFront=this%DisNormFront
        ! wavenum=DisNormFront%GetWaveNum()
        AlphaFront=this%disNormFront%getAlpha()

        !call this%PrintILoc()

        Wavenum=this%DisNorm%GetWaveNum()
        alphaNew=wavenum%getAlpha()
        Dx_alpha=this%Coef(2)*AlphaNew+this%Coef(1)*AlphaFront
        call WaveNum%SetDxAlpha(Dx_alpha)

        !!用iloc-1位置当DisNorm初值
        !DisNorm=DisNormFront

#IFDEF DEBUG
        print*, this%DisNorm%jn, this%disNormFront%jn
        print*, DisNorm%jn, DisNormFront%jn
        pause 'solvewithalpha subroutine'
#ENDIF

        call this%DisNorm%SetILoc(this%iloc)!; call DisNorm%SetDiffIloc(this%iloc)

        call this%SetDisOP(WaveNum)
        call this%SolveDisNorm(this%DisNorm, this%DisNormFront)
        call this%DisNorm%SetWavenum(WaveNum)
        DisNorm1=(this%Coef(2).mx. this%DisNorm)
        !DisNorm1=(this%Coef(1).mx.DisNormFront)
        DxDisNorm=(this%Coef(1).mx.this%DisNormFront) .add. DisNorm1
        !write(*, *), WaveNum%getAlpha()
        !!算增长率
        call this%ComputeSigma(this%DisNorm, DxDisNorm)

        times=1

        !this%DisNorm=DisNorm

    end subroutine SolveWithoutAlpha

    ! !> 在当地求解一个局部平行的LPSE问题作为LPSE初始入口条件的修正过程
    ! subroutine SolveStart(this)
    !
    !     use mod_parameter, only: EPS, EPS_REL
    !     implicit none
    !     class(lpse_eqn_type), intent(inout) :: this
    !     type(dis_wavenum_lpse_type) :: WaveNum
    !     type(lpse_dis_normal_type) :: DisNormFront, DisNorm, DxDisNorm
    !     complex(R_P) :: AlphaOld, AlphaNew, AlphaFront
    !     complex(R_P) :: Dx_alpha
    !
    !     DisNormFront=this%DisNormFront
    !     !!用iloc-1位置当DisNorm初值
    !     DisNorm=DisNormFront
    !     call DisNorm%SetILoc(this%iloc)!; call DisNorm%SetDiffIloc(this%iloc)
    !     WaveNum= DisNorm%GetWaveNum()
    !
    !     AlphaNew=WaveNum%getAlpha()
    !     AlphaFront=AlphaNew
    !     AlphaOld=0.0d0
    !     Dx_alpha=this%Coef(2)*AlphaNew+this%Coef(1)*AlphaFront
    !
    !     do while ( abs(AlphaNew-AlphaOld)/abs(AlphaOld) >=EPS_REL)
    !         call WaveNum%SetDxAlpha(Dx_alpha)
    !         call this%SetDisOP(WaveNum)
    !        ! call this%CheckDisNorm(DisNormBack)
    !         call this%SolveDisNorm(DisNorm, DisNormFront)
    !         call this%SolveAlpha(DisNorm, DisNormFront, DxDisNorm)
    !         WaveNum= DisNorm%GetWaveNum()
    !         AlphaOld=AlphaNew
    !         AlphaNew=WaveNum%getAlpha()
    !         Dx_alpha=WaveNum%GetDxAlpha()
    !         DisNormFront=DisNorm
    !      end do
    !
    !     !!算增长率
    !     call this%ComputeSigma(DisNorm, DxDisNorm)
    !
    !     call This%Dis%Set(this%iloc, DisNorm)
    !
    ! end subroutine SolveStart

    !> 求解LPSE的物理流向复波数
    subroutine ComputeSigma(this, DisNorm, DxDisNorm)

        implicit none
        class(lpse_eqn_type), intent(in) :: this
        type(lpse_dis_normal_type), intent(inout) :: DisNorm !< 当前站位扰动沿法向分布
        type(lpse_dis_normal_type), intent(inout) :: DxDisNorm !< 当前站位扰动的流向导数沿法向分布
        complex(R_P) :: sigma(7)
        type(dis_wavenum_lpse_type) :: wavenum
        integer :: j_m, j_r, j_u, j_v, j_T, j_w, j
        complex(R_P) :: Alpha
        real(R_P) :: NormMax(7), Norm(7)
        complex(R_P) :: rho, u, v, w, T
        real(R_P) :: Rho0, U0, V0, W0, T0
        type(dis_flux_ij_type) :: FluxPoint
        type(BF_flux_org_ij_type) :: BFPoint
#IFDEF DEBUG
        print*, "subroutine computesigma"
        print*, DisNorm%jn, DxDisNorm%jn
        pause
#ENDIF
        wavenum=DisNorm%GetWaveNum()
        Alpha=wavenum%GetAlpha()
        j_m=1; j_r=1; j_u=1; j_v=1; j_T=1; j_w=1
        NormMax=-10.0d0
        do j=1, this%jn
          call DisNorm%Get(j, FluxPoint)
          call FluxPoint%get(rho, u, v, w,  T)
          BFPoint=this%BFNorm%GetBFPoint(j)
          !BFPoint=this%BaseFlow%GetPoint(this%iloc, j)
          call BFPoint%Get(rho0, U0, V0, W0, T0)
          Norm(1)=abs(rho0*u+U0*rho)
          Norm(2)=abs(rho)
          Norm(3)=abs(u)
          Norm(4)=abs(v)
          Norm(5)=abs(T)
          Norm(7)=abs(w)

          if(Norm(1)>NormMax(1)) then
              j_m=j; NormMax(1)=Norm(1)
          endif
          if(Norm(2)>=NormMax(2)) then
              j_r=j; NormMax(2)=Norm(2)
          endif
          if(Norm(3)>=NormMax(3)) then
              j_u=j; NormMax(3)=Norm(3)
          endif
          if(Norm(4)>=NormMax(4)) then
              j_v=j; NormMax(4)=Norm(4)
          endif
          if(Norm(5)>=NormMax(5)) then
              j_T=j; NormMax(5)=Norm(5)
          endif
          if(Norm(7)>=NormMax(7)) then
              j_w=j; NormMax(7)=Norm(7)
          endif
        end do
        sigma(1)=alpha-CPLI*this%dsigma(DisNorm, DxDisNorm, j_m, 'm')
        sigma(2)=alpha-CPLI*this%dsigma(DisNorm, DxDisNorm, j_r, 'r')
        sigma(3)=alpha-CPLI*this%dsigma(DisNorm, DxDisNorm, j_u, 'u')
        sigma(4)=alpha-CPLI*this%dsigma(DisNorm, DxDisNorm, j_v, 'v')
        sigma(5)=alpha-CPLI*this%dsigma(DisNorm, DxDisNorm, j_T, 'T')
        sigma(6)=Alpha-CPLI*this%dsigma(DisNorm, DxDisNorm, 1,   'E')
        sigma(7)=Alpha-CPLI*this%dsigma(DisNorm, DxDisNorm, j_w, 'w')

        call DisNorm%SetSigma(Sigma)

    end subroutine ComputeSigma

    !> 求解LPSE物理流向复波数的形函数的修正量
    complex(R_P) function dsigma(this, DisNorm, DxDisNorm, jloc, flag)

        implicit none
        class(lpse_eqn_type), intent(in) :: this
        type(lpse_dis_normal_type), intent(in) :: DisNorm !< 当前站位扰动
        type(lpse_dis_normal_type), intent(in) :: DxDisNorm !< 当前站位扰动沿流向导数
        integer, intent(in) :: jloc
        character, intent(in) :: Flag
        complex(R_P) :: rho, u, v, w, T
        complex(R_P) :: dx_rho, dx_u, dx_v, dx_w, dx_T
        real(R_P) :: rho0, u0, v0, w0, T0
        real(R_P) :: dx_rho0, dx_u0, dx_v0, dx_w0, dx_T0
        type(dis_flux_ij_type) :: FluxPoint
        type(bf_point_type) :: BFPoint
        type(BF_flux_org_ij_type) :: BFFluxPoint, dx_BFFluxPoint, tmp
        complex(R_P) :: energy, energy_dx
        integer, parameter :: method=4
#IFDEF DEBUG
        print*, 'function dsigma'
        print*, disNorm%jn, DxDisNorm%jn, jloc
        pause
#ENDIF

        call DisNorm%Get(jloc, FluxPoint)
        call FluxPoint%get(rho, u, v, w, T)
        call DxDisNorm%Get(Jloc, FluxPoint)
        call FluxPoint%get(dx_rho, dx_u, dx_v, dx_w, dx_T)

        if(abs(rho)<1e-12 .and. abs(u)<1e-12 .and. abs(v)<1e-12 .and. &
        & abs(w)<1e-12 .and. abs(T)< 1e-12) then
          dsigma=0.0d0
        else

          select case (Flag)
              case ('m')
                !BFPoint=this%BaseFlow%GetPoint(this%iloc, jloc)
                BFPoint=this%BFNorm%GetPoint(jloc)
                call BFPoint%Get(BFFluxPoint, dx_BFFluxPoint, tmp, tmp)
                call BFFluxPoint%Get(rho0, U0, V0, W0, T0)
                !BFPoint=this%DBaseFlow%GetPartial_Xi(this%iloc, jloc)
                call Dx_BFFluxPoint%Get(dx_rho0, dx_U0, dx_V0, dx_W0, dx_T0)
                ! print*, jloc
                ! print*, rho0, U0
                ! print*, rho, u
                ! pause
                dsigma=(dx_rho0*u+rho0*dx_u+u0*dx_rho+dx_u0*rho)/(rho0*u+U0*rho)
              case ('r')
                dsigma=dx_rho/rho
              case ('u')
                dsigma=dx_u/u
              case ('v')
                dsigma=dx_v/v
              case ('w')
                if(abs(w)==0.0d0) then
                  dsigma=0.0d0
                else
                  dsigma=dx_w/w
                endif
              case ('E')
                  energy_dx=this%Norm_DxDisNorm(DisNorm, DxDisNorm, method)
                  energy=this%Norm_DisNorm(DisNorm, method)
                  dsigma=energy_dx/energy
              case ('T')
                dsigma=dx_T/T
            end select
          endif
    end function dsigma

    ! !> 输出物理流向复波数
    ! subroutine print_sigma(this, fn_surf)

    !     use stringifor
    !     implicit none
    !     class(lpse_eqn_type), intent(in) :: this
    !     type(string), intent(in) :: fn_surf !< 文件名前缀
    !     integer :: i, l, j
    !     complex(R_P) :: sigma(7, This%Istart:This%Iend)
    !     real(R_P) :: GrowthRate(7, This%Istart:This%Iend)
    !     real(R_P) :: Nfact(7, This%Istart:This%Iend)
    !     real(R_P) :: Amp(7, This%Istart:This%Iend)
    !     real(R_P) :: Amp_pse(7, This%Istart:This%Iend)
    !     real(R_P) :: Amp_local(7)
    !     type(lpse_dis_normal_type) :: DisNorm
    !     complex(R_P) :: rho, u, v, w, T
    !     type(BF_flux_org_ij_type) :: BFFlux(1, 1:this%jn)
    !     type(dis_flux_ij_type) :: Flux(this%jn)
    !     real(R_P) :: Rho0, U0, V0, W0, T0
    !     integer :: in
    !     real(R_P), allocatable :: xx(:)

    !     In=this%Grid%GetInSize()
    !     if(.not. allocated(xx)) allocate(xx(in))
    !     call this%Grid%Get_jx(1, xx)
    !     do i=This%istart, this%Iend
    !       call this%Dis%Get(i, DisNorm)
    !       call DisNorm%GetSigma(sigma(:, i))
    !     end do
    !     GrowthRate=-aimag(sigma)

    !     Amp(:, This%Istart)=1.0d0
    !     Nfact(:, This%Istart)=0.0d0
    !     amp_pse(:, This%Istart)=1.0d0
    !     do i=This%Istart+1, This%Iend
    !       do l=1, 7
    !         Nfact(l, i)=Nfact(l, i-1)+(GrowthRate(l, i)+GrowthRate(l, i-1))*0.5d0 &
    !                                   *(xx(i)-xx(i-1))
    !         Amp(l, i)=exp(Nfact(l, i))
    !       enddo
    !       call this%Dis%Get(i, DisNorm)
    !       DisNorm=Amp(6, i) .mx. DisNorm
    !       call this%BaseFlow%GetPart(i, i, 1, this%jn, BFFlux(:, :))
    !       call DisNorm%GetFlux(Flux)
    !       Amp_pse(:, i)=0.0d0
    !       do j=1, this%jn
    !         call BFFlux(1, j)%get(rho0, u0, v0, w0, T0)
    !         call Flux(j)%get(rho, u, v, w, T)
    !         Amp_local(1)=abs(rho*u0+rho0+u)
    !         Amp_local(2)=abs(rho)
    !         Amp_local(3)=abs(u)
    !         Amp_local(4)=abs(v)
    !         Amp_local(5)=abs(T)
    !         Amp_local(6)=0.0d0
    !         Amp_local(7)=abs(w)
    !         do l=1, 7
    !           if(amp_local(l)>=Amp_pse(l, i)) Amp_pse(l, i)=amp_local(l)
    !         enddo
    !       enddo
    !       amp_pse(6, i)=amp(6, i)
    !     end do

    !     open(997, file=fn_surf//'_Alf.plt', form='formatted')
    !     write(997, *)"variables=x, N_u, a_r, s_u, A_u"
    !     do i=This%Istart, this%Iend
    !       write(997, '(5ES20.8)')xx(i), Nfact(3, i), real(sigma(6, i)), -aimag(sigma(3, i)), Amp(3, i)
    !     enddo
    !     close(997)

    !     open(996, file=fn_surf//'_GrowthRate.plt', form='formatted')
    !     write(996, *)"variables=x, s_rhou, s_rho, s_u, s_v, s_T, s_E, s_w"
    !     do i=This%Istart, this%Iend
    !       write(996, '(8ES20.8)')xx(i), (GrowthRate(l, i), l=1, 7)
    !     enddo
    !     close(996)

    !     open(995, file=fn_surf//'_Amp.plt', form='formatted')
    !     write(995, *)"variables=x, A_rhou, A_rho, A_u, A_v, A_T, A_E, A_w"
    !     do i=This%Istart, this%Iend
    !       write(995, '(8ES20.8)')xx(i), (Amp(l, i), l=1, 7)
    !     enddo
    !     close(995)

    !     open(994, file=fn_surf//'_Nfact.plt', form='formatted')
    !     write(994, *)"variables=x, N_rhou, N_rho, N_u, N_v, N_T, N_E, N_w"
    !     do i=This%Istart, this%Iend
    !       write(994, '(8ES20.8)')xx(i), (Nfact(l, i), l=1, 7)
    !     enddo
    !     close(994)

    !     open(993, file=fn_surf//'_Amp_pse.plt', form='formatted')
    !     write(993, *)"variables=x, A_rhou, A_rho, A_u, A_v, A_T, A_E, A_w"
    !     do i=This%Istart, this%Iend
    !       write(993, '(8ES20.8)')xx(i), (Amp_pse(l, i), l=1, 7)
    !     enddo
    !     close(993)

    ! end subroutine print_sigma

    !> 求解当前站位扰动
    subroutine SolveDisNorm(this, DisNorm, DisNormFront)

        implicit none
        class(lpse_eqn_type), intent(inout) :: this
        type(lpse_dis_normal_type), intent(inout) :: DisNorm !< 当前站位饶松
        type(lpse_dis_normal_type), intent(in) :: DisNormFront !< 前一站位扰动
        integer :: Solve_method
        integer, parameter :: DIRECTIVESOLVER=1001, ITERATIVESOLVER=1002

        Solve_method=this%solver_kind

        select case (Solve_method)
            case (DIRECTIVESOLVER)
              call this%SolveDisNorm_directive(DisNorm, DisNormFront)
            case (ITERATIVESOLVER)
  !            call this%SolveDisNorm_iterative(DisNorm, DisNormFront)
              stop 'Iterative Solver is not able to be used now!'
            case default
              stop "Please input a correct PSE solver type!"
        end select

    end subroutine SolveDisNorm

    ! !> 求解当前站位扰动(迭代法)
    ! subroutine SolveDisNorm_iterative(this, DisNorm, DisNormFront)
    !
    !     use mod_sparse_matrix
    !     use mod_dis_flux
    !     use mod_gmres_adapter
    !     use mod_gmres
    !     implicit none
    !     class(lpse_eqn_type), intent(inout) :: this
    !     type(lpse_dis_normal_type), intent(inout) :: DisNorm !< 当前站位扰动
    !     type(lpse_dis_normal_type), intent(in) :: DisNormFront !< 前一站位扰动
    !     type(lpse_dis_normal_type) :: RHS_DisNorm
    !     type(dis_wavenum_lpse_type) :: wavenum
    !     type(mat_coo_type), save :: MatAcoo
    !     complex(R_P) :: RHS(this%jn*5), SolutionVec(this%jn*5)
    !     real(R_P) :: arrayX(this%jn*5*2), ArrayB(this%jn*5*2), res(this%jn*5*2)
    !     integer :: iloc
    !     integer :: J
    !
    !     RHS_DisNorm=(-this%coef(1)) .mx. DisNormFront
    !     wavenum=DisNorm%GetWaveNum()
    !     select type (solver_ap=>this%solver_ap)
    !     type is (gmres_adapter_type)
    !         call this%SetRHS(RHS_DisNorm, RHS)
    !         call DisNorm%getflux_array(SolutionVec)
    !         !call solver_ap%initial(this%jn*5, SolutionVec, RHS, this)
    !         call solver_ap%Set(SolutionVec, RHS)
    !         if( .not. MatAcoo%IsCreate()) then
    !           call MatAcoo%Create(this%jn*5, this%jn*5*5*5-2*2*2*5*5)
    !         else
    !           call MatAcoo%zeros()
    !         endif
    !         call this%SetMatA(MatAcoo)
    !         call solver_ap%SetMatA(MatAcoo)
    !         call solver_ap%SetOP(AV)
    !         call solver_ap%solve()
    !         call solver_ap%getX(RHS)
    !     end select
    !     call this%GetX(RHS, DisNorm)
    !     call DisNorm%SetWavenum(wavenum)
    !
    ! end subroutine SolveDisNorm_iterative

!     subroutine av(self, arrayX, arrayB)
!
!         use mod_gmres
!         implicit none
!         class(solver_type), intent(in) :: self
!         real(R_P), intent(in) :: ArrayX(:)
!         real(R_P), intent(out) :: ArrayB(:)
!         complex(R_P), allocatable :: X(:)
!         complex(R_P), allocatable :: B(:)
!         type(lpse_dis_normal_type), save :: DisNorm
!         type(lpse_dis_normal_type), save :: DisNormB
!
!         allocate(X(size(ArrayX)/2))
!         allocate(B(size(ArrayB)/2))
!         call array_complex(ArrayX, X)
!         select type (this=>self)
!         type is(lpse_eqn_type)
!             call DisNorm%set(X, this%iloc)
!             call LPSE_AV(this, DisNorm, DisNormB)
!         end select
!         call DisNormB%getflux_array(B)
!         call array_real(B, ArrayB)
!         deallocate(X)
!         deallocate(B)
!     end subroutine av
!
!     subroutine LPSE_AV(this, DisNorm, DisNorm_RHS)
!
!         use mod_lpse_dis_OP_point
!         implicit none
!         class(lpse_eqn_type), intent(in), target :: this
!         type(lpse_dis_normal_type), intent(in) :: DisNorm
!         type(lpse_dis_normal_type), intent(out) :: DisNorm_RHS
!         type(dis_flux_ij_type) :: DisFlux(this%jn)
!         type(dis_flux_ij_type) :: DyDisFlux(this%jn)
!         type(dis_flux_ij_type) :: DyyDisFlux(this%jn)
!         type(dis_flux_ij_type) :: DisRHS(this%jn)
!         integer :: j
!         real(R_P), pointer :: Coef_dx(:)
!         complex(R_P), dimension(5, 5) :: A, B, D, Vyy
!         type(lpse_dis_op_point_type) :: DisOP
!
!         call DisNorm%GetFlux(DisFlux)
!         call this%Diff%Partial_Eta_DisFluxIJ(this%iloc, DisFlux, DyDisFlux)
! !        call this%Diff%Partial_Eta2_DisFluxIJ(this%iloc, this%iloc, DisFlux, DyyDisFlux)
!         call this%Diff%Partial_Eta2_DisFluxIJ(this%iloc, DisFlux, DyyDisFlux)
!         Coef_dx=>this%Coef
!
!         do j=1, this%jn
!             DisOP=this%DisOPNorm%GetPoint(j)
!             call DisOP%Get(A, B, D, Vyy)
!             DisRHS(j)=((A*Coef_dx(size(Coef_dx)) + D) .mx. DisFlux(j)) + &
!                     & (B .mx. DyDisFlux(j)) - &
!                     & (Vyy .mx. DyyDisFlux(j))
!         end do
!         call DisNorm_RHS%Set(DisRHS, this%iloc)
!
!     end subroutine LPSE_AV
!
!     subroutine mv_inv(self, arrayX, arrayB)
!
!         use mod_gmres
!         implicit none
!         class(solver_type), intent(in) :: self
!         real(R_P), intent(in) :: ArrayX(:)
!         real(R_P), intent(out) :: ArrayB(:)
!         complex(R_P), allocatable :: X(:)
!         complex(R_P), allocatable :: B(:)
!         type(lpse_dis_normal_type), save :: DisNorm
!         type(lpse_dis_normal_type), save :: DisNormB
!
!         allocate(X(size(ArrayX)/2))
!         allocate(B(size(ArrayB)/2))
!         call array_complex(ArrayX, X)
!         select type (this=>self)
!         type is(lpse_eqn_type)
!             call DisNorm%set(X, this%iloc)
!             call LPSE_MV_inv(this, DisNorm, DisNormB)
!         end select
!         call DisNormB%getflux_array(B)
!         call array_real(B, ArrayB)
!         deallocate(X)
!         deallocate(B)
!
!     end subroutine mv_inv
!
!     subroutine mv(self, arrayB, arrayX)
!
!         use mod_gmres
!         implicit none
!         class(solver_type), intent(in) :: self
!         real(R_P), intent(in) :: ArrayB(:)
!         real(R_P), intent(out) :: ArrayX(:)
!         complex(R_P), allocatable :: X(:)
!         complex(R_P), allocatable :: B(:)
!         type(lpse_dis_normal_type), save :: DisNorm
!         type(lpse_dis_normal_type), save :: DisNormB
!
!         allocate(X(size(ArrayX)/2))
!         allocate(B(size(ArrayB)/2))
!         call array_complex(ArrayB, B)
!         select type (this=>self)
!         type is(lpse_eqn_type)
!             call DisNormB%set(B, this%iloc)
!             call LPSE_MV(this, DisNorm, DisNormB)
!         end select
!         call DisNorm%getflux_array(X)
!         call array_real(X, ArrayX)
!         deallocate(X)
!         deallocate(B)
!
!     end subroutine mv
!
!     subroutine LPSE_mv(this, DisNorm, DisNorm_RHS)
!
!         use mod_lpse_dis_OP_point
!         implicit none
!         class(lpse_eqn_type),intent(in), target :: this
!         type(lpse_dis_normal_type), intent(out) :: DisNorm
!         type(lpse_dis_normal_type), intent(in) :: DisNorm_RHS
!         type(dis_flux_ij_type) :: Flux(this%jn), Flux_w(this%jn)
!         type(dis_flux_ij_type) :: Flux_tmp
!         integer :: j, l
!         real(R_P), pointer :: Coef_dx(:)
!         complex(R_P), dimension(5, 5) :: A, B, D, Vyy
!         type(lpse_dis_op_point_type) :: DisOP
!         integer, parameter :: XI1 = 1, ETA1 = 2, ETA2 = 3
!         real(R_P) :: Coef_dy(5), Coef_dyy(5), Coef_dx1
!         complex(R_P), dimension(5,5) :: mU, mD, mL
!         real(R_P) :: dy, idy, idy2
!         real(R_P) :: rhox, rhoy, rhoz
!         real(R_P) :: mI(5,5)
!
!         Coef_dx=>this%Coef
!         call DisNorm_RHS%GetFlux(Flux)
!         do j=1, this%jn
! !            call flux(j)%get(fluxArray(:,j))
!             Flux_w(j)=DIS_FLUXIJ_NULL
!         enddo
!
!         do l=1, 5
!             mI(l, l)=1.0d0
!         end do
!         l=2000
!
!         ! (L+D)X=w
!         do j=1, this%jn
!             call this%Diff%GetETACoef(this%iloc, j, ETA1, Coef_dy)
!             call this%Diff%GetETACoef(this%iloc, j, ETA2, Coef_dyy)
!             DisOP=this%DisOPNorm%GetPoint(j)
!             call DisOP%GetSpt(rhox, rhoy, rhoz)
!             call DisOP%Get(A, B, D, Vyy)
!             if( j==1) then
!                 flux_tmp= DIS_FLUXIJ_NULL
!                 dy=(this%eta(j+1)-this%eta(j))
!                 idy=1.0d0/dy
!                 idy2=idy**2
!                 mD=(Coef_dx(size(Coef_dx))* A + D + (-1.0d0*idy)* B)
!                 Flux_tmp=Flux(j)
!                 Flux_w(j)=mD .div. Flux_tmp
!             else if ( j==2 ) then
!                 flux_tmp= DIS_FLUXIJ_NULL
!                 dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                 idy=1.0d0/dy
!                 idy2=idy**2
!                 mL=(-0.5d0*idy) * B - (1.d0*idy2) * Vyy-0.5d0*(1.0d0)*idy*rhoy*mI
!                 Flux_tmp=(mL .mx. Flux_w(j-1))+Flux_tmp
!                 mD=(Coef_dx(size(Coef_dx))* A + D - (-2.0d0*idy2)* Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!                 flux_tmp=Flux(j)-Flux_tmp
!                 Flux_w(j)=mD .div. Flux_tmp
!             else if ( j==this%jn-1) then
!                 Flux_tmp = DIS_FLUXIJ_NULL
!                 dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                 idy=1.0d0/dy
!                 idy2=idy**2
!                 mL=(-0.5d0*idy) * B - (1.d0*idy2) * Vyy-0.5d0*(1.0d0)*idy*rhoy*mI
!                 Flux_tmp=Flux_tmp+ (mL .mx. Flux_w(j-1))
!                 mD = (Coef_dx(size(Coef_dx))* A + D - (-2.0d0*idy2)* Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!                 Flux_tmp=flux(j)-Flux_tmp
!                 Flux_w(j)=mD .div. Flux_tmp
!             else if (j ==this%jn) then
!                 Flux_tmp = DIS_FLUXIJ_NULL
!                 dy=(this%eta(j)-this%eta(j-1))
!                 idy=1.0d0/dy
!                 idy2=idy**2
!                 mL=(-1.0d0*idy) * B
!                 Flux_tmp=Flux_tmp+ (mL .mx. Flux_w(j-1))
!                 Flux_tmp=Flux(j)-Flux_tmp
!                 mD = (Coef_dx(size(Coef_dx))* A + D + (1.0d0*idy)* B)
!                 Flux_w(j)=mD .div. Flux_tmp
!             else
!                 Flux_tmp = DIS_FLUXIJ_NULL
!                 dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                 idy=1.0d0/dy
!                 idy2=idy**2
!                 mL=(-0.5d0*idy) * B - (1.d0*idy2) * Vyy-0.5d0*(1.0d0)*idy*rhoy*mI
!                 Flux_tmp=Flux_tmp + (mL .mx. Flux_w(j-1))
!                 Flux_tmp=Flux(j)-Flux_tmp
!                 mD = (Coef_dx(size(Coef_dx))* A + D - (-2.0d0*idy2)* Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!                 Flux_w(j)=mD .div. Flux_tmp
!             end if
!         enddo
!
!        ! D^(-1)w=y
!        do j=1, this%jn
!            call this%Diff%GetETACoef(this%iloc, j, ETA1, Coef_dy)
!            call this%Diff%GetETACoef(this%iloc, j, ETA2, Coef_dyy)
!            DisOP=this%DisOPNorm%GetPoint(j)
!            call DisOP%GetSpt(rhox, rhoy, rhoz)
!            call DisOP%Get(A, B, D, Vyy)
!            if(j==1)then
!                dy=(this%eta(j+1)-this%eta(j))
!                idy=1.0d0/dy
!                idy2=idy**2
!                mD = (Coef_dx(size(Coef_dx))* A + D +(-1.0d0)*idy* B)
!            else if( j==2 )then
!                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                idy=1.0d0/dy
!                idy2=idy**2
!                mD = (Coef_dx(size(Coef_dx))* A + D -(-2.0d0)*idy2* Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!            else if( j==this%jn-1 )then
!                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                idy=1.0d0/dy
!                idy2=idy**2
!                mD = (Coef_dx(size(Coef_dx))* A + D -(-2.0d0)*idy2* Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!            else if( j==this%jn   )then
!                dy=(this%eta(j)-this%eta(j-1))
!                idy=1.0d0/dy
!                idy2=idy**2
!                mD = (Coef_dx(size(Coef_dx))* A + D +(1.0d0)*idy* B)
!            else
!                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                idy=1.0d0/dy
!                idy2=idy**2
!                mD = (Coef_dx(size(Coef_dx))* A + D -(-2.0d0)*idy2* Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!            end if
!            Flux_w(j)=mD .mx. Flux_w(j)
!            Flux(j)=DIS_FLUXIJ_NULL
!        end do
!
!         ! (U+D)y=b
!         do j=this%jn, 2, -1
!             call this%Diff%GetETACoef(this%iloc, j, ETA1, Coef_dy)
!             call this%Diff%GetETACoef(this%iloc, j, ETA2, Coef_dyy)
!             DisOP=this%DisOPNorm%GetPoint(j)
!             call DisOP%GetSpt(rhox, rhoy, rhoz)
!             call DisOP%Get(A, B, D, Vyy)
!             if(j==1) then
!                 dy=(this%eta(j+1)-this%eta(j))
!                 idy=1.0d0/dy
!                 idy2=idy**2
!                 Flux_tmp = DIS_FLUXIJ_NULL
!                 mU = (1.0d0*idy) * B
!                 Flux_tmp=Flux_tmp+(mU .mx. Flux(j+1))
!                 Flux_tmp=Flux_w(j)-Flux_tmp
!                 mD = (Coef_dx(size(Coef_dx))* A + D + (-1.0d0*idy)* B)
!                 Flux(j)=mD .div. Flux_tmp
!             else if(j==2)then
!                 dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                 idy=1.0d0/dy
!                 idy2=idy**2
!                 Flux_tmp = DIS_FLUXIJ_NULL
!                 mU = (0.5d0*idy) * B - (1.0d0*idy2) * Vyy - 0.5d0*(1.0d0)*idy*rhoy*mI
!                 Flux_tmp=Flux_tmp+(mU .mx. Flux(j+1))
!                 Flux_tmp=Flux_w(j)-Flux_tmp
!                 mD = (Coef_dx(size(Coef_dx))* A + D - (-2.0d0)*idy2* Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!                 Flux(j)=mD .div. Flux_tmp
!             else if(j==this%jn-1) then
!                 dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                 idy=1.0d0/dy
!                 idy2=idy**2
!                 Flux_tmp = DIS_FLUXIJ_NULL
!                 mU = (0.5d0*idy) * B - (1.0d0*idy2) * Vyy - 0.5d0*(1.0d0)*idy*rhoy*mI
!                 Flux_tmp=Flux_tmp+(mU .mx. Flux(j+1))
!                 Flux_tmp = Flux_w(j)-Flux_tmp
!                 mD = (Coef_dx(size(Coef_dx))* A + D - (-2.0d0)*idy2* Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!                 Flux(j)=mD .div. Flux_tmp
!             else if(j==this%jn) then
!                 dy=(this%eta(j)-this%eta(j-1))
!                 idy=1.0d0/dy
!                 idy2=idy**2
!                 Flux_tmp = DIS_FLUXIJ_NULL
!                 mD = (Coef_dx(size(Coef_dx))* A + D + (1.0d0)*idy* B)
!                 Flux_tmp = Flux_w(j)
!                 Flux(j)=mD .div. Flux_tmp
!             else
!                 dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                 idy=1.0d0/dy
!                 idy2=idy**2
!                 Flux_tmp = DIS_FLUXIJ_NULL
!                 mU = (0.5d0*idy) * B - (1.0d0*idy2) * Vyy - 0.5d0*(1.0d0)*idy*rhoy*mI
!                 Flux_tmp=Flux_tmp+(mU .mx. Flux(j+1))
!                 Flux_tmp=Flux_w(j)-Flux_tmp
!                 mD = (Coef_dx(size(Coef_dx))* A + D - (-2.0d0)*idy2* Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!                 Flux(j)=mD .div. Flux_tmp
!             end if
!         end do
!
!         call DisNorm%Set(Flux, this%iloc)
!
!     end subroutine LPSE_mv

!     subroutine LPSE_MV_inv(this, DisNorm, DisNorm_RHS)
!
!         use mod_lpse_dis_OP_point
!         implicit none
!         class(lpse_eqn_type), intent(in), target :: this
!         type(lpse_dis_normal_type), intent(in) :: DisNorm
!         type(lpse_dis_normal_type), intent(out) :: DisNorm_RHS
!         type(dis_flux_ij_type) :: Flux(this%jn), Flux_w(this%jn)
!         type(dis_flux_ij_type) :: Flux_tmp, flux_tmp1
!         type(dis_flux_ij_type) :: fluxb(this%jn), fluxx(this%jn)
!
! !        type(dis_flux_ij_type) :: DyFlux(this%jn)
! !        type(dis_flux_ij_type) :: DyyFlux(this%jn)
! !        type(dis_flux_ij_type) :: DisRHS(this%jn)
!         integer :: j, l
!         real(R_P), pointer :: Coef_dx(:)
!         complex(R_P), dimension(5, 5) :: A, B, D, Vyy
!         type(lpse_dis_op_point_type) :: DisOP
!         integer, parameter :: XI1 = 1, ETA1 = 2, ETA2 = 3
!         real(R_P) :: Coef_dy(5), Coef_dyy(5), Coef_dx1
!         complex(R_P), dimension(5,5) :: mU, mD, mL
!         real(R_P) :: omega
!         real(R_P) :: rhox, rhoy, rhoz
!         real(R_P) :: mI(5,5), dy, idy, idy2
!
!         mI=0.0d0
!         do l=1,5
!             mI(l,l)=1.0d0
!         enddo
!
!         omega=1.0d0
!  !       complex(R_P), dimension(5,5) :: fluxArray
!
!         Coef_dx=>this%Coef
!         call DisNorm%GetFlux(Flux)
!         do j=1, this%jn
! !            call flux(j)%get(fluxArray(:,j))
!             Flux_w(j)=DIS_FLUXIJ_NULL
!         enddo
!
!         l=20000
!
! !        ! (U+D)y=b
! !        do j=this%jn, 2, -1
! !            call this%Diff%GetETACoef(this%iloc, j, ETA1, Coef_dy)
! !            call this%Diff%GetETACoef(this%iloc, j, ETA2, Coef_dyy)
! !            DisOP=this%DisOPNorm%GetPoint(j)
! !            call DisOP%GetSpt(rhox, rhoy, rhoz)
! !            call DisOP%Get(A, B, D, Vyy)
! !            if(j==1) then
! !                dy=(this%eta(j+1)-this%eta(j))
! !                idy=1.0d0/dy
! !                idy2=idy**2
! !                mD = (Coef_dx(size(Coef_dx))* A + D + (-1.0d0*idy)* B)
! !                Flux_w(j)=mD .mx. Flux(j)
! !                mU = (1.0d0*idy) * B
! !                Flux_w(j)=Flux_w(j)+(mU .mx. Flux(j+1))
! !            else if(j==2)then
! !                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
! !                idy=1.0d0/dy
! !                idy2=idy**2
! !                mD = (Coef_dx(size(Coef_dx))* A + D - (-2.0d0*idy2)* Vyy)-0.5d0*(-2.0d0*idy)*rhoy*mI
! !                Flux_w(j)=mD .mx. Flux(j)
! !                mU = (0.5d0*idy) * B - (1.0d0*idy2) * Vyy-0.5d0*1.0d0*idy*rhoy*mI
! !                Flux_w(j)=Flux_w(j)+(mU .mx. Flux(j+1))
! !            else if(j==this%jn-1) then
! !                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
! !                idy=1.0d0/dy
! !                idy2=idy**2
! !                mD = (Coef_dx(size(Coef_dx))* A + D - (-2.0d0*idy2)* Vyy)-0.5d0*(-2.0d0*idy)*rhoy*mI
! !                Flux_w(j)=mD .mx. Flux(j)
! !                mU = (0.5d0*idy) * B - (1.0d0*idy2) * Vyy-0.5d0*1.0d0*idy*rhoy*mI
! !                Flux_w(j)=Flux_w(j)+(mU .mx. Flux(j+1))
! !            else if(j==this%jn) then
! !                dy=(this%eta(j)-this%eta(j-1))
! !                idy=1.0d0/dy
! !                idy2=idy**2
! !                mD = (Coef_dx(size(Coef_dx))* A + D + (1.0d0*idy)* B )
! !                Flux_w(j)=mD .mx. Flux(j)
! !            else
! !                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
! !                idy=1.0d0/dy
! !                idy2=idy**2
! !                mD = (Coef_dx(size(Coef_dx))* A + D - (-2.0d0*idy2)* Vyy)-0.5d0*(-2.0d0*idy)*rhoy*mI
! !                Flux_w(j)=mD .mx. Flux(j)
! !                mU = (0.5d0*idy) * B - (1.0d0*idy2) * Vyy-0.5d0*1.0d0*idy*rhoy*mI
! !                Flux_w(j)=Flux_w(j)+(mU .mx. Flux(j+1))
! !            end if
! !        end do
! !
!
! do j=1, this%jn
!     Flux_w(j)=Flux(j)
! end do
!        ! D^(-1)w=y
!        do j=1, this%jn
! !           call this%Diff%GetETACoef(this%iloc, j, ETA1, Coef_dy)
! !           call this%Diff%GetETACoef(this%iloc, j, ETA2, Coef_dyy)
!            DisOP=this%DisOPNorm%GetPoint(j)
!            call DisOP%Get(A, B, D, Vyy)
!            call DisOP%GetSpt(rhox, rhoy, rhoz)
!            if(j==1)then
!                dy=(this%eta(j+1)-this%eta(j))
!                idy=1.0d0/dy
!                idy2=idy**2
!                mD =(Coef_dx(size(Coef_dx))* A+D-1.0*idy*B)
!            else if( j==2 )then
!                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                idy=1.0d0/dy
!                idy2=idy**2
!                mD = (Coef_dx(size(Coef_dx))*A+D-(-2.0d0)*idy2*Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!            else if( j==this%jn-1 )then
!                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                idy=1.0d0/dy
!                idy2=idy**2
!                mD = (Coef_dx(size(Coef_dx))*A+D-(-2.0d0)*idy2*Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!            else if( j==this%jn   )then
!                dy=(this%eta(j)-this%eta(j-1))
!                idy=1.0d0/dy
!                idy2=idy**2
!                mD =(Coef_dx(size(Coef_dx))*A+D+1.0*idy*B)
!            else
!                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
!                idy=1.0d0/dy
!                idy2=idy**2
!                mD = (Coef_dx(size(Coef_dx))*A+D-(-2.0d0)*idy2*Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
!            end if
!            Flux_w(j)=mD .div. Flux_w(j)
!            Flux(j)=DIS_FLUXIJ_NULL
!        end do
! do j=1, this%jn
!     Flux(j)=Flux_w(j)
! end do
!
! !        ! (omega*L+D)X=w
! !        do j=1, this%jn
! !            call this%Diff%GetETACoef(this%iloc, j, ETA1, Coef_dy)
! !            call this%Diff%GetETACoef(this%iloc, j, ETA2, Coef_dyy)
! !            DisOP=this%DisOPNorm%GetPoint(j)
! !            call DisOP%GetSpt(rhox, rhoy, rhoz)
! !            call DisOP%Get(A, B, D, Vyy)
! !            if( j==1) then
! !                dy=(this%eta(j+1)-this%eta(j))
! !                idy=1.0d0/dy
! !                idy2=idy**2
! !                mD =(Coef_dx(size(Coef_dx))* A+D-1.0*idy*B)
! !                Flux(j)=mD .mx. Flux_w(j)
! !
! !            else if ( j==2 ) then
! !                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
! !                idy=1.0d0/dy
! !                idy2=idy**2
! !                mD = (Coef_dx(size(Coef_dx))*A+D-(-2.0d0)*idy2*Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
! !
! !                Flux(j)=mD .mx. Flux_w(j)
! !                mL=-0.5d0*idy*B-(1.0d0)*idy2*Vyy-0.5d0*(1.0d0)*idy*rhoy*mI
! !                Flux(j)=Flux(j)+ omega* ( mL .mx. Flux_w(j-1))
! !
! !            else if ( j==this%jn-1) then
! !                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
! !                idy=1.0d0/dy
! !                idy2=idy**2
! !                mD = (Coef_dx(size(Coef_dx))*A+D-(-2.0d0)*idy2*Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
! !
! !                Flux(j)=mD .mx. Flux_w(j)
! !                mL=-0.5d0*idy*B-(1.0d0)*idy2*Vyy-0.5d0*(1.0d0)*idy*rhoy*mI
! !                Flux(j)=Flux(j)+ omega* ( mL .mx. Flux_w(j-1))
! !
! !            else if (j ==this%jn) then
! !                dy=(this%eta(j)-this%eta(j-1))
! !                idy=1.0d0/dy
! !                idy2=idy**2
! !                mD = (Coef_dx(size(Coef_dx))* A + D + 1.0d0*idy*B)
! !                Flux(j)=mD .mx. Flux_w(j)
! !                mL=(-1.0d0)*idy*B
! !                Flux(j)=Flux(j)+omega*(mL .mx. Flux_w(j-1))
! !            else
! !
! !                dy=(this%eta(j+1)-this%eta(j-1))*0.5d0
! !                idy=1.0d0/dy
! !                idy2=idy**2
! !                mD = (Coef_dx(size(Coef_dx))*A+D-(-2.0d0)*idy2*Vyy)-0.5d0*(-2.0d0)*idy*rhoy*mI
! !
! !                Flux(j)=mD .mx. Flux_w(j)
! !                mL=-0.5d0*idy*B-(1.0d0)*idy2*Vyy-0.5d0*(1.0d0)*idy*rhoy*mI
! !                Flux(j)=Flux(j)+ omega* ( mL .mx. Flux_w(j-1))
! !
! !            end if
! !        enddo
!
!         call DisNorm_RHS%Set(Flux, this%iloc)
!
!     end subroutine LPSE_MV_inv

    !> 求解当前站位扰动(直接法)
    subroutine SolveDisNorm_directive(this, DisNorm, DisNormFront)

        use mod_sparse_matrix
        use mod_dis_flux
        implicit none
        class(lpse_eqn_type), intent(inout) :: this
        type(lpse_dis_normal_type), intent(inout) :: DisNorm !< 当前站位扰动
        type(lpse_dis_normal_type), intent(in) :: DisNormFront !< 前一站位扰动
        type(lpse_dis_normal_type) :: RHS_DisNorm
        type(dis_wavenum_lpse_type) :: wavenum
!        type(mat_coo_type), save :: MatAcoo
        type(mat_bsr_type), save :: MatAbsr
        complex(R_P) :: RHS(this%jn*5)
        integer :: iloc
        integer :: J

        RHS_DisNorm=(-this%coef(1)) .mx. DisNormFront

        wavenum=DisNorm%GetWaveNum()
        select type (solver_ap=>this%solver_ap)
        type is (pardiso_adapter_type)
            call solver_ap%SetNDim(this%jn, 5)
            call solver_ap%SetNNZ(this%jn*5-3*2, 5)
            !if( .not. MatAcoo%IsCreate()) then
            !  call MatAcoo%Create(this%jn*5, this%jn*5*5*5-2*2*2*5*5)
            !else
            !  call MatAcoo%zeros()
            !endif
            if( .not. MatAbsr%IsCreate()) then
              call MatAbsr%Create(this%jn, this%jn*5-3*2, 5)
            else
              call MatAbsr%zeros()
            endif

!            call this%SetMatA(MatAcoo)
            call this%SetMatA(MatAbsr)
!            call solver_ap%SetMatA(MatAcoo)
            call solver_ap%SetMatA(MatAbsr)
            call this%SetRHS(RHS_DisNorm, RHS)
            !call this%solver_ap%Solve(RHS_DisNorm, DisNorm, this%DisOPNorm)
            call solver_ap%Solve(RHS)
        end select
        call this%GetX(RHS, DisNorm)
        call DisNorm%SetWavenum(wavenum)

    end subroutine SolveDisNorm_directive

    !>设置PSE方程离散化的线性方程组的系数矩阵\f$M\f$
    subroutine set_matA_lpse(this, MatACoo)

        use mod_lpse_dis_OP_point
        use mod_sparse_matrix
        implicit none
        class(lpse_eqn_type), target, intent(in) :: this
        type(mat_coo_type), intent(inout) :: MatAcoo !< 方程组系数矩阵
        integer :: iloc, jn
        !type(difference_2D_type), pointer :: DisDiff
        type(lpse_dis_op_normal_type), pointer :: DisOPNorm
        type(lpse_dis_op_point_type) :: DisOP
        complex(R_P), dimension(5, 5) :: A, B, D, Vyy
        real(R_P) :: Coef_dy(5), Coef_dyy(5), Coef_dx1
        real(R_P), pointer :: Coef_dx(:)
        complex(R_P), dimension(-2:2, 5, 5) :: tmpIJ
        integer :: j, l
        integer, parameter :: XI1 = 1, ETA1 = 2, ETA2 = 3

        iloc = this%iloc
        jn = this%jn
        !if(.not. (associated(DisDiff, this%Diff))) DisDiff => this%Diff
        if(.not. (associated(DisOPNorm, this%DisOPNorm))) &
        &   DisOPNorm => this%DisOPNorm
        if(.not. (associated(Coef_dx, this%Coef))) Coef_dx => this%Coef
        Coef_dx1=Coef_dx(size(Coef_dx))

!        do j=1, 1
         j=1
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, B, D, Vyy)
          tmpIJ=0.0d0
          tmpIJ(0, :, :)=A*Coef_dx1+D+B*Coef_dy(1)-Vyy*Coef_dyy(1)
          tmpIJ(1, :, :)=             B*Coef_dy(2)-Vyy*Coef_dyy(2)
          tmpIJ(2, :, :)=             B*Coef_dy(3)-Vyy*Coef_dyy(3)
          do l=0, 2
             call MatAcoo%set(5, j, l+1, tmpIJ(l, :, :))
          end do
!        end do

!        do j=1, 2
         j=2
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, B, D, Vyy)
          tmpIJ=0.0d0
          tmpIJ(0, :, :)=             B*Coef_dy(1)-Vyy*Coef_dyy(1)
          tmpIJ(1, :, :)=A*Coef_dx1+D+B*Coef_dy(2)-Vyy*Coef_dyy(2)
          tmpIJ(2, :, :)=             B*Coef_dy(3)-Vyy*Coef_dyy(3)
          do l=0, 2
             call MatAcoo%set(5, j, l+1, tmpIJ(l, :, :))
          end do

        do j=3, jn-2
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, B, D, Vyy)
          tmpIJ=0.0d0
          tmpIJ(-2, :, :)=B*Coef_dy(1)-Vyy*Coef_dyy(1)
          tmpIJ(-1, :, :)=B*Coef_dy(2)-Vyy*Coef_dyy(2)
          tmpIJ( 0, :, :)=B*Coef_dy(3)-Vyy*Coef_dyy(3)+A*Coef_dx1+D
          tmpIJ( 1, :, :)=B*Coef_dy(4)-Vyy*Coef_dyy(4)
          tmpIJ( 2, :, :)=B*Coef_dy(5)-Vyy*Coef_dyy(5)
          do l=-2, 2
            call MatAcoo%set(5, j, j+l, tmpIJ(l, :, :))
          enddo
        enddo

!        do j=jn-1, jn
        j=jn-1
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, B, D, Vyy)
          tmpIJ=0.0d0
          tmpIJ(-2, :, :)=B*Coef_dy(3)-Vyy*Coef_dyy(3)
          tmpIJ(-1, :, :)=B*Coef_dy(4)-Vyy*Coef_dyy(4)+A*Coef_dx1+D
          tmpIJ( 0, :, :)=B*Coef_dy(5)-Vyy*Coef_dyy(5)
          do l=-2, 0
            call MatAcoo%set(5, j, jn+l, tmpIJ(l, :, :))
          enddo
!        enddo
        j=jn
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, B, D, Vyy)
          tmpIJ=0.0d0
          tmpIJ(-2, :, :)=B*Coef_dy(3)-Vyy*Coef_dyy(3)
          tmpIJ(-1, :, :)=B*Coef_dy(4)-Vyy*Coef_dyy(4)
          tmpIJ( 0, :, :)=B*Coef_dy(5)-Vyy*Coef_dyy(5)+A*Coef_dx1+D
          do l=-2, 0
            call MatAcoo%set(5, j, jn+l, tmpIJ(l, :, :))
          enddo

    end subroutine set_matA_lpse

    !>设置PSE方程离散化的线性方程组的系数矩阵\f$M\f$
    subroutine set_matA_lpse_bsr(this, MatAbsr)

        use mod_lpse_dis_OP_point
        use mod_sparse_matrix
        implicit none
        class(lpse_eqn_type), target, intent(in) :: this
        type(mat_bsr_type), intent(inout) :: MatAbsr !< 方程组系数矩阵
        integer :: iloc, jn
        !type(difference_2D_type), pointer :: DisDiff
        type(lpse_dis_op_normal_type), pointer :: DisOPNorm
        type(lpse_dis_op_point_type) :: DisOP
        complex(R_P), dimension(5, 5) :: A, B, D, Vyy
        real(R_P) :: Coef_dy(5), Coef_dyy(5), Coef_dx1
        real(R_P), pointer :: Coef_dx(:)
        complex(R_P), dimension(-3:3, 5, 5) :: tmpIJ !<检查错误用
!        complex(R_P), dimension(-2:2, 5, 5) :: tmpIJ
        integer :: j, l
        integer, parameter :: XI1 = 1, ETA1 = 2, ETA2 = 3

        iloc = this%iloc
        jn = this%jn
        !if(.not. (associated(DisDiff, this%Diff))) DisDiff => this%Diff
        if(.not. (associated(DisOPNorm, this%DisOPNorm))) &
        &   DisOPNorm => this%DisOPNorm
        if(.not. (associated(Coef_dx, this%Coef))) Coef_dx => this%Coef
        Coef_dx1=Coef_dx(size(Coef_dx))

!        do j=1, 1
         j=1
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, B, D, Vyy)
          tmpIJ=0.0d0
          tmpIJ(0, :, :)=A*Coef_dx1+D+B*Coef_dy(1)-Vyy*Coef_dyy(1)
          tmpIJ(1, :, :)=             B*Coef_dy(2)-Vyy*Coef_dyy(2)
          tmpIJ(2, :, :)=             B*Coef_dy(3)-Vyy*Coef_dyy(3)

          do l=0, 2
             !call MatAcoo%set(5, j, l+1, tmpIJ(l, :, :))
             call MatAbsr%set(tmpIJ(l, :, :), j, l+1)
          end do

!
!           block
!             use mod_debug
!             integer :: m1, m2
!             if(mm==0 .and. nn==0)then
!           ! print*, Coef_dy(1:5)
!           ! print*, Coef_dyy(1:5)
!           ! pause 'coef'
!           open(117, file='check_matrix.plt')
!           write(117, *)"variables='j', 'l', 'm12', 'value'"
!           write(117, *)'Zone i=',25,'j=', 5,'k=',this%jn
!
!           do l=-2, 2
!             do m1=1, 5
!               do m2=1, 5
!                 write(117, '(3I10, E20.7)')j, l, m1*10+m2, real(tmpij(l, m1, m2))
!               enddo
!             enddo
!           enddo
!           endif
!
! !        end do

!        do j=1, 2
         j=2
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, B, D, Vyy)
          tmpIJ=0.0d0

          tmpIJ(-1, :, :)=             B*Coef_dy(1)-Vyy*Coef_dyy(1)
          tmpIJ( 0, :, :)=A*Coef_dx1+D+B*Coef_dy(2)-Vyy*Coef_dyy(2)
          tmpIJ( 1, :, :)=             B*Coef_dy(3)-Vyy*Coef_dyy(3)
          tmpIJ( 2, :, :)=             B*Coef_dy(4)-Vyy*Coef_dyy(4)

          do l=-1, 2
             !call MatAcoo%set(5, j, l+1, tmpIJ(l, :, :))
             call MatAbsr%set(tmpIJ(l, :, :), j, l+j)
          end do

          ! if(mm==0 .and. nn==0)then
          !   ! print*, Coef_dy(1:5)
          !   ! print*, Coef_dyy(1:5)
          !   ! pause 'coef'
          ! do l=-2, 2
          !   do m1=1, 5
          !     do m2=1, 5
          !       write(117, '(3I10, E20.7)')j, l, m1*10+m2, real(tmpij(l, m1, m2))
          !     enddo
          !   enddo
          ! enddo
          ! endif


        do j=3, jn-2
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, B, D, Vyy)
          tmpIJ=0.0d0
          tmpIJ(-2, :, :)=B*Coef_dy(1)-Vyy*Coef_dyy(1)
          tmpIJ(-1, :, :)=B*Coef_dy(2)-Vyy*Coef_dyy(2)
          tmpIJ( 0, :, :)=B*Coef_dy(3)-Vyy*Coef_dyy(3)+A*Coef_dx1+D
          tmpIJ( 1, :, :)=B*Coef_dy(4)-Vyy*Coef_dyy(4)
          tmpIJ( 2, :, :)=B*Coef_dy(5)-Vyy*Coef_dyy(5)
          do l=-2, 2
            !call MatAcoo%set(5, j, j+l, tmpIJ(l, :, :))
            call MatAbsr%set(tmpIJ(l, :, :), j, j+l)
          enddo
        !   if(mm==0 .and. nn==0)then
        !     if(j==30)then
        !     print*, A(3, 1)
        !     print*, B(3, 1)
        !     print*, D(3, 1)
        !     print*, Vyy(3, 1)
        !     print*, Coef_dx1*A(3, 1)
        !     print*, Coef_dy(3)*B(3, 1)
        !     print*, Coef_dyy(3)*Vyy(3, 1)
        !     print*, tmpij(0, 3, 1)
        !     pause 'check matrix'
        !     endif
        !     ! print*, Coef_dy(1:5)
        !     ! print*, Coef_dyy(1:5)
        !     ! pause 'coef'
        !   do l=-2, 2
        !     do m1=1, 5
        !       do m2=1, 5
        !         write(117, '(3I10, E20.7)')j, l, m1*10+m2, real(tmpij(l, m1, m2))
        !       enddo
        !     enddo
        !   enddo
        !   endif
        !
         enddo

!        do j=jn-1, jn
        j=jn-1
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, B, D, Vyy)
          tmpIJ=0.0d0
          tmpIJ(-2, :, :)=B*Coef_dy(2)-Vyy*Coef_dyy(2)
          tmpIJ(-1, :, :)=B*Coef_dy(3)-Vyy*Coef_dyy(3)
          tmpIJ( 0, :, :)=B*Coef_dy(4)-Vyy*Coef_dyy(4)+A*Coef_dx1+D
          tmpIJ( 1, :, :)=B*Coef_dy(5)-Vyy*Coef_dyy(5)
          do l=-2, 1
!            call MatAcoo%set(5, j, jn+l, tmpIJ(l, :, :))
            call MatAbsr%set(tmpIJ(l, :, :), j, j+l)
          enddo
! !        enddo
! if(mm==0 .and. nn==0)then
!   ! print*, Coef_dy(1:5)
!   ! print*, Coef_dyy(1:5)
!   ! pause 'coef'
! do l=-2, 2
!   do m1=1, 5
!     do m2=1, 5
!       write(117, '(3I10, E20.7)')j, l, m1*10+m2, real(tmpij(l, m1, m2))
!     enddo
!   enddo
! enddo
! endif
        j=jn
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, B, D, Vyy)
          tmpIJ=0.0d0
          tmpIJ(-2, :, :)=B*Coef_dy(3)-Vyy*Coef_dyy(3)
          tmpIJ(-1, :, :)=B*Coef_dy(4)-Vyy*Coef_dyy(4)
          tmpIJ( 0, :, :)=B*Coef_dy(5)-Vyy*Coef_dyy(5)+A*Coef_dx1+D

          do l=-2, 0
            !call MatAcoo%set(5, j, jn+l, tmpIJ(l, :, :))
            call MatAbsr%Set(tmpIJ(l, :, :), j, jn+l)
          enddo

        !   if(mm==0 .and. nn==0)then
        !     ! print*, Coef_dy(1:5)
        !     ! print*, Coef_dyy(1:5)
        !     ! pause 'coef'
        !   do l=-2, 2
        !     do m1=1, 5
        !       do m2=1, 5
        !         write(117, '(3I10, E20.7)')j, l, m1*10+m2, real(tmpij(l, m1, m2))
        !       enddo
        !     enddo
        !   enddo
        !
        !   close(117)
        !   !stop
        ! endif
        ! endblock

     end subroutine set_matA_lpse_bsr

    !> 设置LPSE方程离散化的线性方程组的右端项
    subroutine set_rhs_lpse(this, DisNorm, RHS)

        use mod_lpse_dis_normal
        use mod_dis_flux
        use mod_vector_cmplx
        use mod_lpse_dis_OP_normal
        use mod_lpse_dis_OP_point
        implicit none
        class(lpse_eqn_type), intent(in) :: this
        type(lpse_dis_normal_type), intent(in) :: DisNorm !< 扰动
        complex(R_P), intent(out) :: RHS(this%jn*5) !< LPSE离散化线性方程组右端项
!        type(lpse_dis_op_normal_type), intent(in) :: DisOPNorm
        type(lpse_dis_op_point_type) :: DisOP
        type(dis_flux_ij_type) :: FluxPoint
        integer :: j, l
        complex(R_P) :: rho, T, U, V, W
        type(vector_cmplx_type) :: Vel
        complex(R_P) :: A(5, 5)
        complex(R_P) :: rhsbc(5)

        ! block

        do j=1, this%jn
          call DisNorm%Get(j, FluxPoint)
          DisOP= this%DisOPNorm%GetPoint(j)
          call DisOP%Get(A, 'A')
          FluxPoint=A .mx. FluxPoint

          call FluxPoint%get(Rho, Vel, T)
          call Vel%Get(U, V, W)

          rhs((j-1)*5+1)=Rho+this%nonlinearTerm(1, j)
          rhs((j-1)*5+2)=U+this%nonlinearTerm(2, j)
          rhs((j-1)*5+3)=V+this%nonlinearTerm(3, j)
          rhs((j-1)*5+4)=W+this%nonlinearTerm(4, j)
          rhs((j-1)*5+5)=T+this%nonlinearTerm(5, j)
          !
          ! block
          !   use mod_debug
          ! if(j==301 .and. mm==0 .and. nn==0) then
          !   write(*, *) 'RHS, nonlinear..'
          !   write(*, '(4E20.7)') rho, this%nonlinearterm(1, j)
          !   write(*, '(4E20.7)') u, this%nonlinearterm(2, j)
          !   write(*, '(4E20.7)') v, this%nonlinearterm(3, j)
          !   write(*, '(4E20.7)') w, this%nonlinearterm(4, j)
          !   write(*, '(4E20.7)') T, this%nonlinearterm(5, j)
          !   pause 'j301 00'
          !
          ! endif
          ! endblock

        enddo

        !! update boundary condition in the RHS

        j=1
        call this%LpseBC(1)%Get(rhsbc)
        do l=1, 5
          if(abs(rhsbc(l))>=1.0d10) then
            rhs(l)=rhs(l)
          else
            rhs(l)=rhsbc(l)
          endif
        enddo

        j=this%jn
        call this%LpseBC(2)%Get(rhsbc)
        do l=1, 5
          if(abs(rhsbc(l))>=1.0d10) then
            rhs((j-1)*5+l)=rhs((j-1)*5+l)
          else
            rhs((j-1)*5+l)=rhsbc(l)
          endif
        enddo

      !   block
      !     use mod_debug
      !     integer :: l
      !
      !
      !     if(nn==0 .and. mm==0)then
      !       write(*, *)'start check rhs 0 0'
      !       do j=this%jn, 1, -1
      !         print*, 'location', j-1
      !         do l=1, 5
      !           write(277+l, '(I10, 2E20.7)')j, rhs((j-1)*5+l)
      !         enddo
      !         !pause
      !       enddo
      !
      !       pause 'rhs done'
      !     endif
      !   endblock
      !
      ! !endblock

    end subroutine set_rhs_lpse

    !> 得PSE方程离散化的线性方程组的解, 并将其放入扰动形函数中
    subroutine get_X_lpse(this, ArrayX, X)

        use mod_lpse_dis_normal
        use mod_dis_flux
        use mod_vector_cmplx
        implicit none
        type(lpse_dis_normal_type), intent(inout) :: X !< 扰动形函数
        complex(R_P), intent(in) :: ArrayX(:) !< 线性方程组数组形式解
        class(lpse_eqn_type), intent(in) :: this
        type(dis_flux_ij_type) :: FluxPoint
        integer :: j
        complex(R_P) :: rho, T, U, V, W
        type(vector_cmplx_type) :: Vel

        call X%Create(this%jn)

        ! block
        !   use mod_debug

        do j=1, this%jn

          Rho=ArrayX((j-1)*5+1)
          U  =ArrayX((j-1)*5+2)
          V  =ArrayX((j-1)*5+3)
          W  =ArrayX((j-1)*5+4)
          T  =ArrayX((j-1)*5+5)
          ! if(mm==0 .and. nn==0)then
          !   print*,j
          !   print*, rho
          !   print*, u
          !   print*, v
          !   print*, w
          !   print*, T
          !   pause
          !
          ! endif
          call Vel%Set(U, V, W)
          call FluxPoint%set(rho, Vel, T)
          call X%SetFluxPoint(j, FluxPoint)
        enddo

        ! block
        ! use mod_debug
        ! integer :: j, l
        ! if(nn==0 .and. mm==0)then
        !   do j=1, this%jn
        !     write(*, *)'start check X 0 0', j-1
        !     do l=1, 5
        !       write(*, '(I10, 2E20.7)')j, ArrayX((j-1)*5+l)
        !       write(377+l, '(I10, 2E20.7)')j, ArrayX((j-1)*5+l)
        !     enddo
        !     !pause
        !   enddo
        !   pause 'X done'
        !   stop
        ! endif
        ! endblock
        !

      ! endblock
    end subroutine get_X_lpse

    !> 求解流向波数\f$\alpha\f$
    subroutine SolveAlpha(this, DisNorm, DisNormFront, DxDisNorm) !! 计算alpha和dxalpha

        use mod_parameter, only: CPLI
        implicit none
        class(lpse_eqn_type), intent(in) :: this
        type(lpse_dis_normal_type), intent(inout) :: DisNorm !< 当前站位扰动
        type(lpse_dis_normal_type), intent(in) :: DisNormFront !< 上一站位扰动
        type(lpse_dis_normal_type), intent(out) :: DxDisNorm !<当前站位扰动沿流向导数
        type(lpse_dis_normal_type) :: DisNorm1, DisNorm2
        type(dis_wavenum_lpse_type) :: Wavenum, WavenumFront
        complex(R_P) :: AlphaOld, AlphaNew, AlphaFront, DxAlpha
        integer :: method
        integer, parameter :: ENERGEY_Method=4
        complex(R_P) :: energy_dx, energy

        method=ENERGEY_Method

        WaveNum= DisNorm%GetWaveNum()
        AlphaOld=WaveNum%getAlpha()
        WavenumFront=DisNormFront%GetWaveNum()
        AlphaFront=WavenumFront%getAlpha()

        DisNorm2=(this%Coef(2).mx.DisNorm)
        DisNorm1=(this%Coef(1).mx.DisNormFront)

        DxDisNorm=DisNorm1 .add. DisNorm2
        energy_dx=this%Norm_DxDisNorm(DisNorm, DxDisNorm, method)
        energy=this%Norm_DisNorm(DisNorm, method)
        AlphaNew=AlphaOld-CPLI* energy_dx/ energy
        DxAlpha=this%Coef(1)*AlphaFront+This%Coef(2)*AlphaNew
        call Wavenum%SetAlpha(AlphaNew)
        call Wavenum%SetDxAlpha(DxAlpha)
        call DisNorm%SetWavenum(wavenum)

    end subroutine

    !> 求扰动形函数范数
    function Norm_DisNorm(this, DisNorm, method) result(Norm)

        use mod_dis_flux
        use mod_vector_cmplx
        use mod_baseflow_org
        implicit none
        class(lpse_eqn_type), intent(in) :: this
        type(lpse_dis_normal_type), intent(in) :: DisNorm !< 扰动
        integer, intent(in) :: method !< 范数标准
        complex(R_P) :: Norm !< 扰动范数
        complex(R_P) :: InputVar(this%jn)
        integer :: J
        type(dis_flux_ij_type) :: Flux(this%jn)
        type(vector_cmplx_type) :: Vel
        integer, parameter :: ENERGEY_Method=4

        select case (method)
            case (ENERGEY_Method)
              call DisNorm%GetFlux(Flux)
              do j=1, this%jn
                Vel=Flux(j)%GetVel()
                InputVar(j)=Vel%SelfCmplxDot()
              enddo
              Norm=Intergal(this%Jn, InputVar, This%Eta)
            case default
                stop "Please select a correct method for Normalizing the Disturbance!"
        end select

    end function Norm_DisNorm

    !> 求扰动形函数的流向导数的范数
    function Norm_DxDisNorm(this, DisNorm, DxDisNorm, method) result(Norm)

        use mod_dis_flux
        use mod_vector_cmplx
        use mod_baseflow_org
        implicit none
        class(lpse_eqn_type), intent(in) :: this
        type(lpse_dis_normal_type), intent(in) :: DisNorm !< 扰动
        type(lpse_dis_normal_type), intent(in) :: DxDisNorm !< 扰动沿流向导数
        integer, intent(in) :: method !< 归一化范数准则
        complex(R_P) :: Norm !< 归一化后的扰动流向导数范数
        integer :: j
        complex(R_P) :: InputVar(this%jn)
        type(dis_flux_ij_type) :: Flux(this%jn), DxFlux(this%jn)
        type(vector_cmplx_type) :: Vel, DxVel
        integer, parameter :: ENERGEY_Method=4

        select case (method)
            case (ENERGEY_Method)
              call DisNorm%GetFlux(Flux)
              call DxDisNorm%GetFlux(DxFlux)
              do j=1, this%jn
                Vel=Flux(j)%GetVel(); DxVel=DxFlux(j)%GetVel()
                InputVar(j)=Vel .dot. DxVel
              enddo
              Norm=Intergal(this%Jn, InputVar, This%Eta)
            case default
                stop "Please select a correct method for Normalizing the Disturbance!"
        end select

    end function Norm_DxDisNorm

    !> 复数积分函数
    function zIntergal(n, Var, Y) Result(IntValue)

        implicit none

        integer, intent(in) :: n !< 离散点数量
        complex(R_P), intent(in) :: Var(n) !< 被积函数的函数值
        real(R_P), intent(in) :: Y(n) !< 被积函数函数值对应的自变量
        complex(R_P) :: IntValue !< 积分值
        integer :: i

        IntValue=0.0d0
        do i=1, n-1
          IntValue=IntValue+0.5d0*(Var(i)+Var(i+1))*(Y(i+1)-Y(i))
        end do

    end function zIntergal

    !> 实数积分函数
    function dIntergal(n, Var, Y) Result(IntValue)

        implicit none

        integer, intent(in) :: n !< 离散点数量
        real(R_P), intent(in) :: Var(n) !< 被积函数的函数值
        real(R_P), intent(in) :: Y(n) !< 被积函数函数值对应的自变量
        real(R_P) :: IntValue !< 积分值
        integer :: i

        IntValue=0.0d0
        do i=1, n-1
          IntValue=IntValue+0.5d0*(Var(i)+Var(i+1))*(Y(i+1)-Y(i))
        end do

    end function dIntergal

    !> 析构函数
    subroutine finalize_lpse(this)

        implicit none
        class(lpse_eqn_type), intent(inout) :: this

        this%jn=0
        call this%BFOPNorm%finalize
        call this%DisOPNorm%finalize
        call this%BFNorm%finalize
        call this%NormCoord%finalize
        this%iloc=0
        select type (solver_ap=>this%solver_ap)
        type is (pardiso_adapter_type)
            call solver_ap%Finalize
        end select

    end subroutine finalize_lpse

end module mod_lpse_eqn
