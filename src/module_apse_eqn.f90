!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_apse.f90
!> @file
!> @breif 伴随LPSE求解器文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: module_apse
!> @breif 伴随LPSE求解器模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2018-06-26 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2018-06-26
!------------------------------------------------------------------------------

module mod_apse_eqn

    use mod_local_normal_coordinate
    use mod_difference
    use mod_BF_normal
    use mod_lns_OP_normal
    use mod_lpse_dis_OP_normal
    use mod_lpse_bc
    use mod_dis_wavenum
    use mod_lpse_dis_normal
    use mod_pardiso_adapter
    use mod_dis_flux
    use mod_baseflow_org
    use mod_parameter, only: CPLI
    use penf, only: R_P

    implicit none

    private

    !> 伴随LPSE求解器类
    type, public :: apse_eqn_type

        private
        type(BF_normal_type), Pointer, private :: BFNorm !<基本流流向某站位法向分布
        type(lpse_dis_normal_type), Pointer, private :: DisNorm !<扰动流向某站位法向分布
        type(lpse_dis_normal_type), Pointer, private :: DisNormBack !<扰动流向后一站位法向分布
        type(local_normal_coordinate_type), pointer, private :: NormCoord !< 某流向展向站位的局部坐标系
        type(Norm_Diff_Coef_type), pointer, private :: NormCoef
        type(lpse_bf_op_normal_type), private :: BFOPNorm !< LNS方程在当地流向展向位置沿法向分布的系数算子
        type(lpse_dis_op_normal_type), private :: DisOPNorm !< PSE方程的线性系数算子沿法向分布
        integer, private :: iloc !< 当地流向站位编号
        integer, private :: jn !< 法向网格数
        integer, private :: bctype(5, 2) !< PSE方程边界条件类型
        real(R_P), private :: Coef(2) !< 流向扰动离散的拉格朗日基函数
        real(R_P), allocatable, private :: Eta(:) !< 当地位置法向网格分布
        class(*), allocatable, private :: solver_ap !< 线性方程求解器
        integer, private :: solver_kind !< 线性代数求解器类型

        contains

        procedure :: Create !< 创建LPSE类型, 分配内存
        procedure :: SetSolver !< 设置线性系统求解器类型
        !procedure :: SetDisInlet=>SetDisInlet_IR !<用LST给定伴随PSE入口
        procedure :: SolvePSE !< 推进求解伴随PSE

        procedure :: SetBCType  !< 设置边界条件类型

        generic, private :: SetMatA => set_matA_lpse, set_matA_lpse_bsr !<设置伴随PSE方程离散化的线性方程组的系数矩阵\f$M\f$

        procedure, private :: SetIloc => set_iloc !< 设置当前流向站位
        procedure, private :: SetBFOP => set_BF_OP_Norm !< 设置LNS系数算子沿法向分布
        procedure, private :: SolveWithAlpha !<伴随LPSE求解器
        procedure, private :: finalize => finalize_lpse !<析构函数
        !procedure, private :: PrintILoc !< 输出当前站位序号
        procedure, private :: SolveDisNorm !< 求解当前站位扰动
        procedure, private :: SolveAlpha !< 求解流向波数\f$\alpha\f$
        procedure, private :: Norm_DisNorm !< 求扰动形函数范数
        procedure, private :: Norm_DxDisNorm !< 求扰动形函数的流向导数范数
        procedure, private :: SolveDisNorm_directive !< 求解当前站位扰动(直接法)

        procedure, private :: SetDisOP => set_dis_op !< 设置伴随LPSE的线性算子系数
        procedure, private :: set_matA_lpse !< 设置伴随PSE方程离散化的线性方程组的系数矩阵\f$M\f$(COO)
        procedure, private :: set_matA_lpse_bsr !< 设置伴随PSE方程离散化的线性方程组的系数矩阵\f$M\f$(BSR)
        procedure, private :: SetRHS => set_rhs_lpse !< 设置伴随PSE方程离散化的线性方程组的右端项
        procedure, private :: GetX => get_X_lpse !< 获得伴随PSE方程离散化的线性方程组的解, 并将其放入扰动形函数中
        procedure, private :: SolveWithoutAlpha !< 求解伴随PSE方程,但不对波数\f$\alpha\f$进行迭代

    end type apse_eqn_type

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
        class(apse_eqn_type), intent(inout) :: this
        integer, intent(in) :: bctype(5, 2) !< 边界条件类型

        this%bctype=bctype

    end subroutine SetBCType
    !> 设置线性代数求解器相关信息
    subroutine SetSolver(this, method)

        implicit none
        class(apse_eqn_type),intent(inout) :: this
        integer, intent(in), optional :: method    !< 线性代数求解器求解方法
        integer, parameter ::  DIRECT_SOLVER=1001, ITERATIVE_SOLVER=1002
        integer, parameter :: NORMAL=0, TRANSPOSE=2

        if(present(method)) then
            this%solver_kind=method
        else
            this%solver_kind=DIRECT_SOLVER
        end if

        select case (this%solver_kind)
            case (DIRECT_SOLVER)
                allocate(pardiso_adapter_type::this%solver_ap)
            case (ITERATIVE_SOLVER)
!                allocate(pardiso_adapter_type::solver_ap)
        end select
        select type (solver_ap=>this%solver_ap)
        type is (pardiso_adapter_type)
             call solver_ap%Initialize(NORMAL)
        end select

    end subroutine SetSolver



    !> PSE推进求解
    subroutine SolvePSE(this, DisNorm, iloc, DisNormBack, &
    &               BFNorm, NormCoord, NormCoef, Eta, BCtype, IsWithAlpha)

        implicit none
        class(apse_eqn_type), intent(inout) :: this
        type(lpse_dis_normal_type), target, intent(inout) :: DisNorm
        integer, intent(in) :: iloc
        type(lpse_dis_normal_type), target, intent(in) :: DisNormBack
        type(BF_normal_type), target, intent(in) :: BFNorm
        type(local_normal_coordinate_type), target, intent(in) :: NormCoord
        type(Norm_Diff_Coef_type), target, intent(in) :: NormCoef
        real(R_P), intent(in) :: Eta(:)
        integer, intent(in) :: BCtype(5, 2)
        logical, optional, intent(in) :: isWithAlpha

        logical, parameter :: IS_PARALLEL=.false.

        call this%SetIloc(iloc)

        this%BFNorm=>BFNorm
        this%DisNormBack=>DisNormBack
        this%NormCoord=>NormCoord
        this%NormCoef=>NormCoef
        this%DisNorm=>DisNorm

        this%Coef=NormCoef%Coef_di

        this%Eta=Eta

        call this%SetBCType(BCtype)

        call this%SetBFOP(IS_PARALLEL)
        call this%DisOPNorm%SetCoef(this%Coef)
        if(.not. present(isWithAlpha) .or. isWithAlpha) then
          call this%SolveWithAlpha()
        else
          call this%SolveWithoutAlpha()
        endif


    end subroutine SolvePSE




    !> 创建LPSE类型, 分配内存
    subroutine create(this, jn)

        implicit none
        class(apse_eqn_type), intent(inout) :: this
        integer, intent(in) :: jn !< 法向点数
        integer, parameter :: NORMAL=0, TRANSPOSE=2
        integer :: eqkind

        this%jn=jn
        eqkind=0
        call this%BFOPNorm%Create(this%jn, eqkind)
        call this%DisOPNorm%Create(this%jn)
        allocate(this%eta(this%jn))

    end subroutine create




    !
     !> 设置当前流向站位
     subroutine set_iloc(this, iloc)

         implicit none
         integer, intent(in) :: iloc
         class(apse_eqn_type), intent(inout) :: this

         this%iloc=iloc

     end subroutine set_iloc

    !> 设置LNS系数算子沿法向分布
    subroutine set_BF_OP_Norm(this, isParallel)

        implicit none

        class(apse_eqn_type), intent(inout) :: this
        logical, intent(in) :: isParallel !< 流场是否考虑非平行性




        call this%BFOPNorm%Set(this%BFNorm, this%NormCoord, isParallel)

    end subroutine set_BF_OP_Norm

    !> 设置伴随LPSE的线性算子系数
    subroutine set_dis_op(this, WaveNum)

        implicit none
        class(apse_eqn_type), intent(inout) :: this
        type(dis_wavenum_lpse_type), intent(in) :: WaveNum !< 扰动色散特征
        type(lpse_bc_type) :: LpseBC(2)

        call LpseBC(1)%set(this%bctype(:, 1), this%BFOPNorm%GetPoint(1), wavenum, this%NormCoord%GetPoint(1))
        call LpseBC(2)%set(this%bctype(:, 2), this%BFOPNorm%GetPoint(this%jn), wavenum, this%NormCoord%GetPoint(this%jn))
        call this%DisOPNorm%Set(this%BFOPNorm, this%NormCoord, wavenum, lpsebc)

    end subroutine set_dis_op



    !> 伴随LPSE求解器
    subroutine SolveWithAlpha(this)

        use mod_parameter, only: EPS, EPS_REL
        implicit none
        class(apse_eqn_type), intent(inout) :: this
        type(dis_wavenum_lpse_type) :: WaveNum
        type(lpse_dis_normal_type) :: DisNormBack, DisNorm, DxDisNorm
        complex(R_P) :: AlphaOld, AlphaNew, AlphaBack
        complex(R_P) :: Dx_alpha
        logical :: LoopFlag

        !call This%Dis%Get(this%iloc+1, DisNormBack)
	      DisNormBack=this%DisNormBack
        !call this%PrintILoc()
        !!用iloc-1位置当DisNorm初值
        DisNorm=DisNormBack!this%DisNorm
        call DisNorm%SetILoc(this%iloc)
        WaveNum= DisNorm%GetWaveNum()

        ! print*, DisNorm%GetAlpha()
        ! print*, DisNormBack%GetAlpha()
        ! pause

        AlphaNew=WaveNum%getAlpha()
        AlphaBack=this%DisNormBack%GetAlpha()
        AlphaOld=0.0d0
        Dx_alpha=this%Coef(1)*AlphaNew+this%Coef(2)*AlphaBack
        LoopFlag=.True.
        do while (LoopFlag)
            call WaveNum%SetDxAlpha(Dx_alpha)
            call this%SetDisOP(WaveNum)
            call this%SolveDisNorm(DisNorm, DisNormBack)
            call this%SolveAlpha(DisNorm, DisNormBack, DxDisNorm)
            WaveNum= DisNorm%GetWaveNum()
            AlphaOld=AlphaNew
            AlphaNew=WaveNum%getAlpha()
            Dx_alpha=WaveNum%GetDxAlpha()

              write(*, *) "AlphaOld=", (AlphaOld)
              write(*, *) "AlphaNew=", (AlphaNew)
              write(*, *) "AlphaErr=", abs(AlphaNew-AlphaOld)

            LoopFlag=abs(AlphaNew-AlphaOld)/abs(AlphaNew) >=EPS_REL

        end do
        write(*, *) "AlphaErr=", abs(AlphaNew-AlphaOld)
        write(*, *)alphaNew

          !!算增长率

          this%DisNorm=DisNorm

    end subroutine SolveWithAlpha

    !>  求解伴随LPSE方程,但不对波数\f$\alpha\f$进行迭代
    subroutine SolveWithoutAlpha(this)

        use mod_parameter, only: EPS
        implicit none
        class(apse_eqn_type), intent(inout) :: this
        type(dis_wavenum_lpse_type):: WaveNum !< 扰动色散特征
        type(lpse_dis_normal_type) :: DisNormBack, DisNorm, DxDisNorm
        type(lpse_dis_normal_type) :: DisNorm1
        complex(R_P) :: Dx_alpha, alphaNew, AlphaBack

        DisNormBack=this%DisNormBack
        !call this%PrintILoc()
        wavenum=DisNormBack%GetWaveNum()
        AlphaBack=wavenum%getAlpha()

        Wavenum=this%DisNorm%GetWaveNum()
        alphanew=wavenum%getAlpha()
        !!用iloc-1位置当DisNorm初值
        DisNorm=DisNormBack
        call DisNorm%SetILoc(this%iloc)!; call DisNorm%SetDiffIloc(this%iloc)

        Dx_alpha=this%Coef(1)*AlphaNew+this%Coef(2)*AlphaBack
        call Wavenum%SetDxAlpha(Dx_alpha)
        call this%SetDisOP(WaveNum)
        call this%SolveDisNorm(DisNorm, DisNormBack)
        call DisNorm%SetWavenum(WaveNum)
        DisNorm1=(this%Coef(1).mx.DisNorm)
        !DisNorm1=(this%Coef(1).mx.DisNormFront)
        DxDisNorm=(this%Coef(2).mx.DisNormBack) .add. DisNorm1
        !write(*, *), WaveNum%getAlpha()
        !!算增长率
     !   call this%ComputeSigma(DisNorm, DxDisNorm)
      !  call This%Dis%Set(this%iloc, DisNorm)
        this%DisNorm=DisNorm

    end subroutine SolveWithoutAlpha


    !
    !
    !
    ! !> 输出当前站位序号
    ! subroutine PrintILoc(this)
    !
    !     implicit none
    !     class(apse_eqn_type), intent(in) :: this
    !
    !     write(*, *)"The location index now is", this%iloc
    !
    ! end subroutine PrintILoc

    !> 求解当前站位扰动
    subroutine SolveDisNorm(this, DisNorm, DisNormBack)

        implicit none
        class(apse_eqn_type), intent(inout) :: this
        type(lpse_dis_normal_type), intent(inout) :: DisNorm !< 当前站位扰动
        type(lpse_dis_normal_type), intent(in) :: DisNormBack !< 后一站位扰动
        integer :: Solve_method
        integer, parameter :: DIRECTIVESOLVER=1001, ITERATIVESOLVER=1002

        Solve_method=this%solver_kind

        select case (Solve_method)
            case (DIRECTIVESOLVER)
              call this%SolveDisNorm_directive(DisNorm, DisNormBack)
            case (ITERATIVESOLVER)
              !call this%SolveDisNorm_interative(DisNorm, DisNormBack)
            case default
              stop "Please input a correct PSE solver type!"
        end select

    end subroutine SolveDisNorm

!     !> 求解当前站位扰动(迭代法)
!     subroutine SolveDisNorm_interative(this, DisNorm, DisNormBack)
!
!         implicit none
!         class(apse_eqn_type), intent(in) :: this
!         type(lpse_dis_normal_type), intent(inout) :: DisNorm !< 当前站位扰动
!         type(lpse_dis_normal_type), intent(in) :: DisNormBack !< 后一站位扰动
!         real(R_P), dimension(this%jn*5) :: RHS
!         real(R_P), dimension(this%jn*5) :: SolutionVec
!
! !        call DisNorm%UpdateDiff(this%Diff)
!
!
!     end subroutine SolveDisNorm_interative

    !> 求解当前站位扰动(直接法)
    subroutine SolveDisNorm_directive(this, DisNorm, DisNormBack)

        use mod_sparse_matrix
        implicit none
        class(apse_eqn_type), intent(inout) :: this
        type(lpse_dis_normal_type), intent(inout) :: DisNorm !< 当前站位扰动
        type(lpse_dis_normal_type), intent(in) :: DisNormBack !< 后一站位扰动
        type(lpse_dis_normal_type) :: RHS_DisNorm
        type(dis_wavenum_lpse_type) :: wavenum
        type(mat_bsr_type), save :: MatAbsr
        complex(R_P) :: RHS(this%jn*5)

        RHS_DisNorm=(this%coef(2)) .mx. DisNormBack

        wavenum=DisNorm%GetWaveNum()
        select type (solver_ap=> this%solver_ap)
        type is (pardiso_adapter_type)
            call solver_ap%SetNDim(this%jn, 5)
            call solver_ap%SetNNZ(this%jn*5-3*2, 5)

            !if( .not. MatAcoo%ISCreate()) then
            !  call MatAcoo%Create(this%jn*5, this%jn*5*5*5-2*2*2*5*5)
            !else
            !  call MatAcoo%zeros()
            !endif
            !call this%SetMatA(MatAcoo)
            !call solver_ap%SetMatA(MatAcoo)
            if( .not. MatAbsr%IsCreate()) then
              call MatAbsr%Create(this%jn, this%jn*5-3*2, 5)
            else
              call MatAbsr%zeros()
            endif
            call this%SetMatA(MatAbsr)
            call solver_ap%SetMatA(MatAbsr)

            call this%SetRHS(RHS, RHS_DisNorm)
            !call this%solver_ap%Solve(RHS_DisNorm, DisNorm, this%DisOPNorm)
            call solver_ap%Solve(RHS)
        end select
        call this%GetX(RHS, DisNorm)
        call DisNorm%SetWavenum(wavenum)

    end subroutine SolveDisNorm_directive



    !>设置伴随PSE方程离散化的线性方程组的系数矩阵\f$M\f$(COO)
    subroutine set_matA_lpse(this, MatACoo)

        use mod_lpse_dis_OP_point
        use mod_sparse_matrix
        implicit none
        class(apse_eqn_type), target, intent(in) :: this
        type(mat_coo_type), intent(inout) :: MatAcoo !< 方程组系数矩阵
        integer :: iloc, jn
        type(difference_2D_type), pointer :: DisDiff
        type(lpse_dis_op_normal_type), pointer :: DisOPNorm
        type(lpse_dis_op_point_type) :: DisOP
        complex(R_P), dimension(5, 5) :: A, D, DxA
        complex(R_P), dimension(5, 5, this%jn) :: B, Vyy
        real(R_P) :: Coef_dy(5), Coef_dyy(5), Coef_dx1
        real(R_P), pointer :: Coef_dx(:)
        complex(R_P), dimension(-2:2, 5, 5) :: tmpIJ
        integer :: j, l
        integer, parameter :: XI1 = 1, ETA1 = 2, ETA2 = 3

        iloc = this%iloc
        jn = this%jn
        !if (.not. (associated(DisDiff, this%Diff))) DisDiff => this%Diff
        DisOPNorm => this%DisOPNorm
        Coef_dx => this%Coef
        Coef_dx1=-Coef_dx(1)

        do j=1, jn
            DisOP=DisOPNorm%GetPoint(j)
            call DisOP%Get(B(:, :, j), "B")
            call DisOP%Get(Vyy(:, :, j), "V")
        end do

!        do j=1, 1
         j=1
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, 'A')
          call DisOP%Get(D, 'D')
          call DisOP%get(DxA, 'X')
          tmpIJ(0, :, :)= transpose(A)*Coef_dx1+transpose(D) &
                       & -transpose(B(:, :, j  ))*Coef_dy(1) &
                       & -transpose(Vyy(:, :, j  ))*Coef_dyy(1) &
                       & -transpose(DxA)
          tmpIJ(1, :, :)=-transpose(B(:, :, j+1))*Coef_dy(2) &
                       & -transpose(Vyy(:, :, j+1))*Coef_dyy(2)
          tmpIJ(2, :, :)=-transpose(B(:, :, j+2))*Coef_dy(3)- &
                       &  transpose(Vyy(:, :, j+2))*Coef_dyy(3)
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
          call DisOP%Get(A, 'A')
          call DisOP%Get(D, 'D')
          call DisOP%get(DxA, 'X')
          tmpIJ(0, :, :)=-transpose(B(:, :, j-1))*Coef_dy(1)-transpose(Vyy(:, :, j-1))*Coef_dyy(1)
          tmpIJ(1, :, :)= transpose(A)*Coef_dx1+transpose(D) &
                      &  -transpose(B(:, :, j  ))*Coef_dy(2)-transpose(Vyy(:, :, j  ))*Coef_dyy(2)-transpose(DxA)
          tmpIJ(2, :, :)=-transpose(B(:, :, j+1))*Coef_dy(3)-transpose(Vyy(:, :, j+1))*Coef_dyy(3)
          do l=0, 2
             call MatAcoo%set(5, j, l+1, tmpIJ(l, :, :))
          end do

        do j=3, jn-2
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, 'A')
          call DisOP%Get(D, 'D')
          call DisOP%get(DxA, 'X')
          tmpIJ(-2, :, :)=-transpose(B(:, :, j-2))*Coef_dy(1) &
                        & -transpose(Vyy(:, :, j-2))*Coef_dyy(1)
          tmpIJ(-1, :, :)=-transpose(B(:, :, j-1))*Coef_dy(2) &
                        & -transpose(Vyy(:, :, j-1))*Coef_dyy(2)
          tmpIJ( 0, :, :)=-transpose(B(:, :, j  ))*Coef_dy(3) &
                        & -transpose(Vyy(:, :, j  ))*Coef_dyy(3)+transpose(A)*Coef_dx1+transpose(D)-transpose(DxA)
          tmpIJ( 1, :, :)=-transpose(B(:, :, j+1))*Coef_dy(4)-transpose(Vyy(:, :, j+1))*Coef_dyy(4)
          tmpIJ( 2, :, :)=-transpose(B(:, :, j+2))*Coef_dy(5)-transpose(Vyy(:, :, j+2))*Coef_dyy(5)
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
          call DisOP%Get(A, 'A')
          call DisOP%Get(D, 'D')
          call DisOP%get(DxA, 'X')
          tmpIJ(-2, :, :)=-transpose(B(:, :, j-1))*Coef_dy(3)-transpose(Vyy(:, :, j-1))*Coef_dyy(3)
          tmpIJ(-1, :, :)=-transpose(B(:, :, j  ))*Coef_dy(4)-transpose(Vyy(:, :, j  ))*Coef_dyy(4)+ &
          &                transpose(A)*Coef_dx1+transpose(D)-transpose(DxA)
          tmpIJ( 0, :, :)=-transpose(B(:, :, j+1))*Coef_dy(5)-transpose(Vyy(:, :, j+1))*Coef_dyy(5)
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
          call DisOP%Get(A, 'A')
          call DisOP%Get(D, 'D')
          call DisOP%get(DxA, 'X')
          tmpIJ(-2, :, :)=-transpose(B(:, :, j-2))*Coef_dy(3)-transpose(Vyy(:, :, j-2))*Coef_dyy(3)
          tmpIJ(-1, :, :)=-transpose(B(:, :, j-1))*Coef_dy(4)-transpose(Vyy(:, :, j-1))*Coef_dyy(4)
          tmpIJ( 0, :, :)=-transpose(B(:, :, j  ))*Coef_dy(5)-transpose(Vyy(:, :, j  ))*Coef_dyy(5) &
                    &     +transpose(A)*Coef_dx1+transpose(D)-transpose(DxA)
          do l=-2, 0
            call MatAcoo%set(5, j, jn+l, tmpIJ(l, :, :))
          enddo

          DisOPNorm=>null()
          Coef_dx=>null()

     end subroutine set_matA_lpse

    !>设置伴随PSE方程离散化的线性方程组的系数矩阵\f$M\f$(BSR)
    subroutine set_matA_lpse_bsr(this, MatAbsr)

        use mod_lpse_dis_OP_point
        use mod_sparse_matrix
        implicit none
        class(apse_eqn_type), target, intent(in) :: this
        type(mat_bsr_type), intent(inout) :: MatAbsr !< 方程组系数矩阵
        integer :: iloc, jn
        type(difference_2D_type), pointer :: DisDiff
        type(lpse_dis_op_normal_type), pointer :: DisOPNorm
        type(lpse_dis_op_point_type) :: DisOP
        complex(R_P), dimension(5, 5) :: A, D, DxA
        complex(R_P), dimension(5, 5, this%jn) :: B, Vyy
        real(R_P) :: Coef_dy(5), Coef_dyy(5), Coef_dx1
        real(R_P), pointer :: Coef_dx(:)
        complex(R_P), dimension(-2:2, 5, 5) :: tmpIJ
        integer :: j, l
        integer, parameter :: XI1 = 1, ETA1 = 2, ETA2 = 3

        iloc = this%iloc
        jn = this%jn
        !if (.not. (associated(DisDiff, this%Diff))) DisDiff => this%Diff
        DisOPNorm => this%DisOPNorm
        Coef_dx => this%Coef
        Coef_dx1=-Coef_dx(1)

        do j=1, jn
            DisOP=DisOPNorm%GetPoint(j)
            call DisOP%Get(B(:, :, j), "B")
            call DisOP%Get(Vyy(:, :, j), "V")
        end do

!        do j=1, 1
         j=1
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, 'A')
          call DisOP%Get(D, 'D')
          call DisOP%get(DxA, 'X')
          tmpIJ(0, :, :)= transpose(A)*Coef_dx1+transpose(D) &
                       & -transpose(B(:, :, j  ))*Coef_dy(1) &
                       & -transpose(Vyy(:, :, j  ))*Coef_dyy(1) &
                       & -transpose(DxA)
          tmpIJ(1, :, :)=-transpose(B(:, :, j+1))*Coef_dy(2) &
                       & -transpose(Vyy(:, :, j+1))*Coef_dyy(2)
          tmpIJ(2, :, :)=-transpose(B(:, :, j+2))*Coef_dy(3)- &
                       &  transpose(Vyy(:, :, j+2))*Coef_dyy(3)
          do l=0, 2
             call MatAbsr%Set(tmpIJ(l,:,:), j, l+1)
          end do
!        end do

!        do j=1, 2
         j=2
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, 'A')
          call DisOP%Get(D, 'D')
          call DisOP%get(DxA, 'X')
          tmpIJ(-1, :, :)=-transpose(B(:, :, j-1))*Coef_dy(1)- &
          &               transpose(Vyy(:, :, j-1))*Coef_dyy(1)
          tmpIJ(0, :, :)= transpose(A)*Coef_dx1+transpose(D) &
                      &  -transpose(B(:, :, j  ))*Coef_dy(2)- &
                      &   transpose(Vyy(:, :, j  ))*Coef_dyy(2)-transpose(DxA)
          tmpIJ(1, :, :)=-transpose(B(:, :, j+1))*Coef_dy(3)- &
                      &   transpose(Vyy(:, :, j+1))*Coef_dyy(3)
          tmpIJ(2, :, :)=-transpose(B(:, :, j+2))*Coef_dy(4)- &
                      &   transpose(Vyy(:, :, j+2))*Coef_dyy(4)
          do l=-1, 2
             call MatAbsr%Set(tmpIJ(l, :, :), j, l+j)
          end do

        do j=3, jn-2
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, 'A')
          call DisOP%Get(D, 'D')
          call DisOP%get(DxA, 'X')
          tmpIJ(-2, :, :)=-transpose(B(:, :, j-2))*Coef_dy(1) &
                        & -transpose(Vyy(:, :, j-2))*Coef_dyy(1)
          tmpIJ(-1, :, :)=-transpose(B(:, :, j-1))*Coef_dy(2) &
                        & -transpose(Vyy(:, :, j-1))*Coef_dyy(2)
          tmpIJ( 0, :, :)=-transpose(B(:, :, j  ))*Coef_dy(3) &
                        & -transpose(Vyy(:, :, j  ))*Coef_dyy(3)+ &
                        &  transpose(A)*Coef_dx1+transpose(D)-transpose(DxA)
          tmpIJ( 1, :, :)=-transpose(B(:, :, j+1))*Coef_dy(4) &
                        & -transpose(Vyy(:, :, j+1))*Coef_dyy(4)
          tmpIJ( 2, :, :)=-transpose(B(:, :, j+2))*Coef_dy(5) &
                        & -transpose(Vyy(:, :, j+2))*Coef_dyy(5)
          do l=-2, 2
            call MatAbsr%Set(tmpIJ(l, :, :), j, j+l)
          enddo
        enddo

!        do j=jn-1, jn
        j=jn-1
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, 'A')
          call DisOP%Get(D, 'D')
          call DisOP%get(DxA, 'X')
          tmpIJ(-2, :, :)=-transpose(B(:, :, j-2))*Coef_dy(2) &
          &               -transpose(Vyy(:, :, j-2))*Coef_dyy(2)
          tmpIJ(-1, :, :)=-transpose(B(:, :, j-1))*Coef_dy(3) &
          &               -transpose(Vyy(:, :, j-1))*Coef_dyy(3)
          tmpIJ( 0, :, :)=-transpose(B(:, :, j  ))*Coef_dy(4) &
          &               -transpose(Vyy(:, :, j  ))*Coef_dyy(4) &
          &               +transpose(A)*Coef_dx1+transpose(D)-transpose(DxA)
          tmpIJ( 1, :, :)=-transpose(B(:, :, j+1))*Coef_dy(5) &
          &               -transpose(Vyy(:, :, j+1))*Coef_dyy(5)
          do l=-2, 1
            call MatAbsr%Set(tmpIJ(l, :, :), j, j+l)
          enddo
!        enddo
        j=jn
          !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
          !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
          Coef_dy=this%NormCoef%Coef_dj(:, j)
          Coef_dyy=this%NormCoef%Coef_djj(:, j)
          DisOP=DisOPNorm%GetPoint(j)
          call DisOP%Get(A, 'A')
          call DisOP%Get(D, 'D')
          call DisOP%get(DxA, 'X')
          tmpIJ(-2, :, :)=-transpose(B(:, :, j-2))*Coef_dy(3) &
                      &   -transpose(Vyy(:, :, j-2))*Coef_dyy(3)
          tmpIJ(-1, :, :)=-transpose(B(:, :, j-1))*Coef_dy(4) &
                      &   -transpose(Vyy(:, :, j-1))*Coef_dyy(4)
          tmpIJ( 0, :, :)=-transpose(B(:, :, j  ))*Coef_dy(5) &
                      &   -transpose(Vyy(:, :, j  ))*Coef_dyy(5) &
                      &   +transpose(A)*Coef_dx1+transpose(D)-transpose(DxA)
          do l=-2, 0
            call MatAbsr%set(tmpIJ(l, :, :), j, jn+l)
          enddo

          DisOPNorm=>null()
          Coef_dx=>null()

     end subroutine set_matA_lpse_bsr

    !> 设置伴随LPSE方程离散化的线性方程组的右端项
    subroutine set_rhs_lpse(this, RHS, DisNorm)

        use mod_lpse_dis_normal
        use mod_dis_flux
        use mod_vector_cmplx
        use mod_lpse_dis_OP_normal
        use mod_lpse_dis_OP_point
        implicit none
        class(apse_eqn_type), intent(in) :: this
        complex(R_P), intent(out) :: RHS(this%jn*5) !< 伴随LPSE离散化线性方程组右端项
!        type(lpse_dis_op_normal_type), intent(in) :: DisOPNorm
        type(lpse_dis_normal_type), intent(in) :: DisNorm !< 扰动
        type(lpse_dis_op_point_type) :: DisOP
        type(dis_flux_ij_type) :: FluxPoint
        integer :: j
        complex(R_P) :: rho, T, U, V, W
        type(vector_cmplx_type) :: Vel
        complex(R_P) :: A(5, 5)

        do j=1, this%jn
          call DisNorm%Get(j, FluxPoint)
          DisOP= this%DisOPNorm%GetPoint(j)
          call DisOP%Get(A, 'A')
          FluxPoint=transpose(A) .mx. FluxPoint

          call FluxPoint%get(Rho, Vel, T)
          call Vel%Get(U, V, W)

          rhs((j-1)*5+1)=Rho
          rhs((j-1)*5+2)=U
          rhs((j-1)*5+3)=V
          rhs((j-1)*5+4)=W
          rhs((j-1)*5+5)=T
        enddo

    end subroutine set_rhs_lpse

    !> 求得伴随PSE方程离散化的线性方程组的解, 并将其放入扰动形函数中
    subroutine get_X_lpse(this, ArrayX, X)

        use mod_lpse_dis_normal
        use mod_dis_flux
        use mod_vector_cmplx
        implicit none
        type(lpse_dis_normal_type), intent(inout) :: X !< 扰动形函数
        complex(R_P), intent(in) :: ArrayX(:) !< 线性方程组数组形式解
        class(apse_eqn_type), intent(inout) :: this
        type(dis_flux_ij_type) :: FluxPoint
        integer :: j
        complex(R_P) :: rho, T, U, V, W
        type(vector_cmplx_type) :: Vel

        call X%Create(this%jn)
        do j=1, this%jn
          Rho=ArrayX((j-1)*5+1)
          U  =ArrayX((j-1)*5+2)
          V  =ArrayX((j-1)*5+3)
          W  =ArrayX((j-1)*5+4)
          T  =ArrayX((j-1)*5+5)
          call Vel%Set(U, V, W)
          call FluxPoint%set(rho, Vel, T)
          call X%SetFluxPoint(j, FluxPoint)
        enddo

    end subroutine get_X_lpse

    !> 求解流向波数\f$\alpha\f$
    subroutine SolveAlpha(this, DisNorm, DisNormBack, DxDisNorm) !! 计算alpha和dxalpha

        use mod_parameter, only: CPLI
        implicit none
        class(apse_eqn_type), intent(in) :: this
        type(lpse_dis_normal_type), intent(inout) :: DisNorm !< 当前站位扰动
        type(lpse_dis_normal_type), intent(in) :: DisNormBack !< 后一站位扰动
        type(lpse_dis_normal_type), intent(out) :: DxDisNorm !<当前站位扰动沿流向导数
        type(lpse_dis_normal_type) :: DisNorm1, DisNorm2
        type(dis_wavenum_lpse_type) :: Wavenum, WavenumBack
        complex(R_P) :: AlphaOld, AlphaNew, AlphaBack, DxAlpha
        integer :: method
        integer, parameter :: ENERGEY_Method=4
        complex(R_P) :: energy_dx, energy

        method=ENERGEY_Method

        WaveNum= DisNorm%GetWaveNum()
        AlphaOld=WaveNum%getAlpha()
        WavenumBack=DisNormBack%GetWaveNum()
        AlphaBack=WavenumBack%getAlpha()
        DisNorm2=(this%Coef(1).mx.DisNorm)
        DisNorm1=(this%Coef(2).mx.DisNormBack)
        DxDisNorm=DisNorm1 .add. DisNorm2
        energy_dx=this%Norm_DxDisNorm(DisNorm, DxDisNorm, method)
        energy=this%Norm_DisNorm(DisNorm, method)
        AlphaNew=AlphaOld+CPLI* energy_dx/ energy
        DxAlpha=This%Coef(1)*AlphaNew+this%Coef(2)*AlphaBack
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
        class(apse_eqn_type), intent(in) :: this
        type(lpse_dis_normal_type), intent(in) :: DisNorm !< 扰动
        integer, intent(in) :: method !< 范数标准
        integer, parameter :: ENERGEY_Method=4
        complex(R_P) :: Norm !< 扰动范数
        complex(R_P) :: InputVar(this%jn)
        integer :: J
        type(dis_flux_ij_type) :: Flux(this%jn)
        type(vector_cmplx_type) :: Vel


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
        class(apse_eqn_type), intent(in) :: this
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
        class(apse_eqn_type), intent(inout) :: this

        this%jn=0
        call this%BFOPNorm%finalize
        call this%DisOPNorm%finalize
        call this%BFNorm%finalize
        call this%NormCoord%finalize
        this%iloc=0
        select type(solver_ap=> this%solver_ap)
        type is (pardiso_adapter_type)
            call solver_ap%Finalize
        end select

    end subroutine finalize_lpse

end module mod_apse_eqn
