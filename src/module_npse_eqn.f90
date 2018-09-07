!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_npse_eqn.f90
!> @file
!> @breif LPSE求解器文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: npse
!> @breif NPSE方程求解器模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2018-06-26 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2018-06-26
!------------------------------------------------------------------------------
module mod_npse_eqn

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
    use mod_lpse_eqn
    use mod_nonlinear

    ! use mod_lst, only: lst_type
    ! use mod_lst_IR, only: lst_IR_type
!    use mod_solver
    use penf, only: R_P

    implicit none

    private

    !> LPSE求解器类
    type, public :: npse_eqn_type

        private
        type(BF_normal_type), Pointer, private :: BFNorm !<基本流流向某站位法向分布
        type(lpse_dis_normal_type), Pointer, private :: DisNorm(:, :) !<扰动流向某站位法向分布
        type(lpse_dis_normal_type), Pointer, private :: DisNormFront(:, :) !<扰动流向前一站位法向分布
        type(local_normal_coordinate_type), pointer, private :: NormCoord !< 某流向展向站位的局部坐标系
        type(Norm_Diff_Coef_type), pointer, private :: NormCoef
        ! type(lpse_bf_op_normal_type), private :: BFOPNorm !< LNS方程在当地流向展向位置沿法向分布的系数算子
        ! type(lpse_dis_op_normal_type), private :: DisOPNorm !< PSE方程的线性系数算子沿法向分布
        integer, private :: iloc !< 当地流向站位编号
        integer, private :: jn !< 法向网格数
        integer, private :: mdim !< 时间方向离散谱
        integer, private :: ndim !< 空间展向离散谱
        integer, private, allocatable :: bctype(:, :, :, :) !< PSE方程边界条件类型
        real(R_P), private :: Coef(2) !< 流向扰动离散的拉格朗日基函数
        real(R_P), allocatable, private :: Eta(:) !< 当地位置法向网格分布
        logical, allocatable, private :: isWithAlpha(:, :) !< 是否迭代Alpha
        type(lpse_eqn_type), private :: LPSE !< LPSE求解器
        complex(R_P), allocatable, private :: nonlinearTerm(:, :, :, :) !< 非线性项
        complex(R_P), allocatable, private :: nonlinearTermOld(:, :, :, :) !<前次非线性项
        complex(R_P), allocatable, private :: IntergalAlf(:, :, :)
        complex(R_P), allocatable, private :: AlfOld(:, :)
        complex(R_P), allocatable, private :: AlfNew(:, :)
        type(nonlinear_type), private :: nonlinearSolver
        logical, private :: isCheckAmp !<是否已经检查过幅值确认是否迭代

        contains

        procedure :: Create !< 创建LPSE类型, 分配内存
        !procedure :: SetSolver !< 设置线性代数求解器相关信息
        procedure :: SolveNPSE !< 推进求解PSE
        procedure :: SetBCType !< 设置边界条件类型

        !generic, private :: SetMatA => set_matA_lpse, set_matA_lpse_bsr !< 设置PSE方程离散化的线性方程组的系数矩阵\f$M\f$

        procedure, private :: SetIloc => set_iloc !< 设置当前流向站位
        !procedure, private :: SetBFOP => set_BF_OP_Norm !< 设置LNS系数算子沿法向分布
        !procedure, private :: SolveWithAlpha !<LPSE求解器
        procedure, private :: PrintILoc !< 输出当前站位序号
        procedure, private :: isConverge !< 判断是否收敛
        procedure, private :: ComputeNonlinearTerm !< 计算非线性项
        procedure, private :: ComputeIntergalAlf !<计算alpha历史积分值
        procedure, private :: CheckAmplitude !< 检查扰动幅值
        procedure, private :: SetAlpha !< 设置扰动波数
        procedure, private :: GetAlf !<设置站位的波数

    end type npse_eqn_type

    contains

    !> 设置LPSE的边界条件
    subroutine SetBCType(this, bctype)

        implicit none
        class(npse_eqn_type), intent(inout) :: this
        integer, intent(in) :: bctype(5, 2) !< 边界条件类型
        integer :: mCount, nCount

        do nCount=-this%ndim, this%ndim
          do mCount=0, this%mdim
            this%bctype(:, :, mCount, nCount)=bctype
          enddo
        enddo

        ! 对平均流修正单独设置边界条件

        this%bctype(3, 2, 0, 0)=3001

    end subroutine SetBCType

    !> PSE推进求解
    subroutine SolveNPSE(this, DisNorm, iloc, DisNormFront, &
    &               BFNorm, NormCoord, NormCoef, Eta, BCtype, isWithAlpha)

        use mod_cfgio_adapter, only: isCrossflow, &
                                &    Coef_Nonlinear_set=>Coef_Nonlinear, &
                                &    Coef_Non_sub_loc
        implicit none
        class(npse_eqn_type), intent(inout) :: this
        type(lpse_dis_normal_type), target, intent(inout) :: DisNorm(0:, -this%ndim:)
        integer, intent(in) :: iloc
        type(lpse_dis_normal_type), target, intent(in) :: DisNormFront(0:, -this%ndim:)
        type(BF_normal_type), target, intent(in) :: BFNorm
        type(local_normal_coordinate_type), target, intent(in) :: NormCoord
        type(Norm_Diff_Coef_type), target, intent(in) :: NormCoef
        real(R_P), intent(in) :: Eta(:)
        logical, intent(in), OPTIONAL :: isWithAlpha(:,:)
        integer, intent(in) :: BCtype(5, 2)
        logical :: isLPSEConvergence(0:this%mdim, -this%ndim:this%ndim)

        integer :: mCount, nCount

        logical, parameter :: IS_PARALLEL=.false.
        integer :: iteriaNum
        real(R_P) :: Coef_Nonlinear

        isLPSEConvergence=.False.

        call this%SetIloc(iloc)

        this%BFNorm=>BFNorm
        this%DisNormFront=>DisNormFront
        this%NormCoord=>NormCoord
        this%NormCoef=>NormCoef
        this%DisNorm=>DisNorm

        this%Coef=NormCoef%Coef_di

        this%Eta=Eta

        call this%SetBCType(BCtype)
#IFDEF DEBUG
          print*, size(this%DisNorm, dim=1)
          print*, size(this%DisNorm, dim=2)

          print*, this%DisNorm(0, 1)%GetAlpha()

          print*, 1; pause
#ENDIF

        if(present(isWithAlpha)) then
          this%isWithAlpha=isWithAlpha
          this%isCheckAmp=.True.
        endif
        if (.not. this%isCheckAmp) call this%CheckAmplitude()
        this%AlfOld=this%GetAlf(this%DisNormFront)

        call this%SetAlpha(isCrossFlow)
        this%AlfNew=this%GetAlf(this%DisNorm)
        this%nonlinearTermOld=0.0d0
        call this%ComputeIntergalAlf()
        call this%ComputeNonlinearTerm()

        !do while( .not. this%isConverge())
        !isNotConv=.False.

        iteriaNum=0

        !print*, all(isLPSEConvergence)
        do while( .not. all(isLPSEConvergence) .and. iteriaNum<=200)

          call this%ComputeIntergalAlf()
          this%nonlinearTermOld=this%nonlinearTerm
          call this%ComputeNonlinearTerm()

#IFDEF DEBUG


          print*, this%iloc
#ENDIF
          do nCount=-this%ndim, this%ndim
            do mCount=0, this%mdim

              !print*, "m_index=", mCount, "n_index=", nCount
              !if(nCount ==0 .and. mCount==0) cycle
              if(this%isWithAlpha(mCount, nCount)) &
              &   write(*, *)'m=',mCount,',n=', nCount


#IFDEF DEBUG
              ! if(mCount==0 .and. nCount==0) then
              !       print*, '....................................'
              !       print*, 'check DisNorm'
              !       print*, mcount, nCount
              !       print*, this%DisNorm(0, 0)%jn
              !       print*, this%DisNorm(0, 0)%GetAlpha()
              !       print*, this%DisNorm(0, 0)%GetBeta()
              !       print*, this%DisNorm(0, 0)%GetOmega()
              !       print*, this%DisNorm(0, 0)%Uamp()
              !       print*, this%DisNorm(0, 0)%Tamp()
              !       pause
              ! endif
              !print*, this%DisNorm(mCount, nCount)%jn
#ENDIF
              !
              ! print*, mCount, nCount, this%isWithAlpha(mCount, nCount)

              ! pause
#IFDEF DEBUG
block

  use mod_debug

  nn=nCount
  mm=mcount
  ! print*, 'print bc', mm, nn
  ! print*, this%bctype(:, :, mm, nn)
  ! pause


endblock
#ENDIF
              if(iloc< Coef_Non_sub_loc) then
                Coef_Nonlinear=1.0d0
              else
                Coef_Nonlinear=Coef_NOnlinear_set
              endif


              call this%LPSE%SolvePSE(this%DisNorm(mCount, nCount), &
              & iloc, this%DisNormFront(mCount, nCount), this%BFNorm, &
              & this%NormCoord, this%NormCoef, this%Eta, &
              & this%bctype(:, :, mCount, nCount), &
              & nonlinearTerm=Coef_Nonlinear*this%nonlinearTerm(:, :, mCount, nCount)+ &
                   (1.0d0-Coef_Nonlinear)*this%nonlinearTermOld(:, :, mCount, nCount), &
              & isWithAlpha=this%IsWithAlpha(mCount, nCount), &
              & isConvergence=isLPSEConvergence(mCount, nCount))

              ! if(mCount==0 .and. nCount==0) then
              !       print*, '....................................'
              !       print*, 'check DisNorm'
              !       print*, mcount, nCount
              !       print*, this%DisNorm(0, 0)%jn
              !       print*, this%DisNorm(0, 0)%GetAlpha()
              !       print*, this%DisNorm(0, 0)%GetBeta()
              !       print*, this%DisNorm(0, 0)%GetOmega()
              !       print*, this%DisNorm(0, 0)%Uamp()
              !       print*, this%DisNorm(0, 0)%Tamp()
              !       pause
              ! endif

            end do
          end do


#IFDEF DEBUG
          print*, this%isConverge()
          pause
#ENDIF
          call this%SetAlpha(isCrossFlow)
          this%AlfNew=this%GetAlf(this%DisNorm)

          !print*, all(isLPSEConvergence)
          iteriaNum=iteriaNum+1

        enddo
#IFDEF DEBUG
        print*, this%isConverge()
        pause
#ENDIF
        this%AlfOld=this%AlfNew
        this%IntergalAlf(1, :, :)=this%IntergalAlf(2, :, :)

    end subroutine SolveNPSE

    function GetAlf(this, DisNorm) result(Alf)
      implicit none
      class(npse_eqn_type), intent(inout) :: this
      class(lpse_dis_normal_type) :: DisNorm(0:, -this%ndim:)
      complex(R_P) :: Alf(0:this%mdim, -this%ndim:this%ndim)
      integer :: mCount, nCount

      do nCount=-this%ndim, this%ndim
        do mCount=0, this%mdim
          Alf(mCount, nCount)= &
          & DisNorm(mCount, nCount)%GetAlpha()
        enddo
      enddo

    end function GetAlf

    subroutine SetAlpha(this, isCrossflow)
      implicit none
      class(npse_eqn_type), intent(inout) :: this
      logical, intent(in) :: isCrossflow

      if(isCrossflow) then
        call SetAlpha_crossflow(this)
      else
        call SetAlpha_streamwise(this)
      endif

    end subroutine SetAlpha

    !> 设置扰动波数(流向不稳定性)
    subroutine SetAlpha_streamwise(this)

      use mod_cfgio_adapter, only: omega0, beta0
      implicit none
      class(npse_eqn_type), intent(inout) :: this
      integer :: mCount, nCount
      real(R_P) :: PhaseVel
      complex(R_P) :: Omega, Alpha, Beta

#IFDEF DEBUG
      ! pause 140
      ! print*, this%DisNorm(0, 1)%GetAlpha()
      ! pause 141
#ENDIF
      PhaseVel=0.0d0
      outer1: do nCount=-this%ndim, this%ndim
        do mCount=1, this%mdim
          if(this%isWithAlpha(mCount, nCount)) then
            Alpha= this%DisNorm(mCount, nCount)%getAlpha()
            Omega= this%DisNorm(mCount, nCount)%getOmega()
            PhaseVel=real(Omega)/real(Alpha)
            print*, 'phase velocity is', PhaseVel
            exit outer1
          endif
        enddo
      enddo outer1

      do nCount=-this%ndim, this%ndim
        do mCount=0, this%mdim
          omega=mCount*omega0
          beta=nCount*beta0

          alpha=Omega/PhaseVel

          call this%DisNorm(mCount, nCount)%SetOmega(omega)
          call this%DisNorm(mCount, nCount)%SetBeta(beta)
          if(.not. this%isWithAlpha(mCount, nCount)) &
          &  call this%DisNorm(mCount, nCount)%SetAlpha(alpha)
        enddo
      enddo

    end subroutine SetAlpha_streamwise

    !> 设置扰动波数(横流不稳定性)
    subroutine SetAlpha_crossflow(this)

      use mod_cfgio_adapter, only: omega0, beta0
      implicit none
      class(npse_eqn_type), intent(in) :: this
      integer :: mCount, nCount
      real(R_P) :: PhaseAngle
      complex(R_P) :: Omega, Alpha, Beta

      PhaseAngle=0.0d0
      outer2: do nCount=-this%ndim, this%ndim
        do mCount=0, this%mdim
          if(this%isWithAlpha(mCount, nCount)) then
            Alpha= this%DisNorm(mCount, nCount)%getAlpha()
            Beta= this%DisNorm(mCount, nCount)%getBeta()
            PhaseAngle=real(Alpha)/real(Beta)
            exit outer2
          endif
        enddo
      enddo outer2

      do nCount=-this%ndim, this%ndim
        do mCount=0, this%mdim
          omega=mCount*omega0
          beta=nCount*beta0
          alpha=PhaseAngle*beta
          call this%DisNorm(mCount, nCount)%SetOmega(omega)
          call this%DisNorm(mCount, nCount)%SetBeta(beta)
          if(.not. this%isWithAlpha(mCount, nCount)) &
          &  call this%DisNorm(mCount, nCount)%SetAlpha(alpha)
        enddo
      enddo

    end subroutine SetAlpha_crossflow

    !> 检查扰动幅值，确定是否需要迭代
    subroutine CheckAmplitude(this)

      use mod_parameter, only: CriticalT, CriticalU, CPLI
      implicit none
      class(npse_eqn_type), intent(inout) :: this

      integer :: mCount, nCount
      real(R_P) :: Amp(0:this%mdim, -this%ndim:this%ndim)

      Amp=exp(CPLI*this%IntergalAlf(1, :, :))
#IFDEF DEBUG
      ! pause  555
#ENDIF

      do nCount=-this%ndim, this%ndim
        do mCount=0, this%mdim
          if(this%DisNorm(mCount, nCount)%TAmp()>=CriticalT) then
            this%isWithAlpha(mCount, nCount)=.True.
          else
            this%isWithAlpha(mCount, nCount)=.False.
          endif
          if(this%DisNorm(mCount, nCount)%UAmp()>=CriticalU) then
            this%isWithAlpha(mCount, nCount)=.True.
          else
            this%isWithAlpha(mCount, nCount)=.False.
          endif
        enddo
      enddo

      this%isWithAlpha(0, 0)=.False.
#IFDEF DEBUG
      ! pause  666
#ENDIF

      this%isCheckAmp=.True.

    end subroutine CheckAmplitude

    !> 计算非线性项
    subroutine ComputeNonlinearTerm(this)

      implicit none
      class(npse_eqn_type), intent(inout) :: this
#IFDEF DEBUG
      print*, 'subroutine compute nonlinear term..'
      print*, this%DisNorm%jn, this%disNormfront%jn
#ENDIF
      call this%nonlinearSolver%SolveNonlinearTerm(this%DisNorm, this%DisNormFront, &
      &             this%BFNorm, this%NormCoord, this%NormCoef, &
      &             this%IntergalAlf(2, :, :), this%AlfNew)


#IFDEF DEBUG

      write(*, *) '............'
      write(*, *) 'Check nonlinearTerm'
      block
        integer :: m, n, l, j
        print*, this%mdim, this%ndim
        do m=0, this%mdim, 1
          do n=0, 0
            do l=1, 5
            !-this%ndim, this%ndim
            write(*, *)m, l, maxval(abs(this%nonlinearTerm(l, :, m, n))),maxloc(abs(this%nonlinearTerm(l, :, m, n)))
            enddo
          enddo
        enddo

      ! do j=1, this%jn
      !
      !   write(*, *)j, abs(this%nonlinearTerm(5, j, 0, 0))
      !
      ! enddo


      pause
      endblock

#ENDIF

    end subroutine ComputeNonlinearTerm

    !>计算$\f\alpha$\f的历史积分值
    subroutine ComputeIntergalAlf(this)
      implicit none
      class(npse_eqn_type), intent(inout) :: this
      integer :: mCount, nCount
      complex(R_P) :: Alf(0:this%mdim, -this%ndim:this%ndim)

      Alf=this%DisNorm%GetAlpha()

          this%IntergalAlf(2, :, :)= &
        & this%IntergalAlf(1, :, :)+ &
        & 0.5d0*(this%AlfOld+Alf)*1.0d0/abs(this%coef(2))
#IFDEF DEBUG
        print*, 'Alf'
        print*, Alf
        print*, 'intergalalf'

        !print*, this%IntergalAlf(1,:,:)
        print*, this%IntergalAlf(2,:,:)
#ENDIF
    end subroutine ComputeIntergalAlf

    !> 判断是否收敛
    logical function isConverge(this)

      use mod_parameter, only: EPS_REL, EPS
      implicit none
      class(npse_eqn_type), intent(in) :: this
      real(R_P), save :: maxNonlinearTerm
      real(R_P), save :: maxNonlinearTermOld
      integer, save :: iterationNum=0

      isConverge=.False.
      iterationNum=iterationNum+1
      maxNonlinearTerm=maxval(abs(this%nonlinearTerm))
      maxNonlinearTermOld=maxval(abs(this%nonlinearTermOld))
#IFDEF DEBUG
      print*, this%nonlinearTerm(1, 1, 1, 1)

      print*, 'maxNonlinearTermValue=', maxNonlinearTerm
      print*, 'maxNonlinearTermValueOld=', maxNonlinearTermOld
      print*, abs(maxNonlinearTerm-maxNonlinearTermOld)/maxNonlinearTerm
      pause
#ENDIF

      if(abs(maxNonlinearTerm-maxNonlinearTermOld)/maxNonlinearTerm<=EPS_REL &
      &   .or. maxNonlinearTerm<=EPS*1E-4) then
        isConverge=.True.
        iterationNum=0
      endif
      if(iterationNum>=40) then
        isConverge=.True.
        iterationNum=0
      endif
      !

    end function isConverge

    !> 创建NPSE类型, 分配内存
    subroutine Create(this, mdim, ndim, jn)

        implicit none
        class(npse_eqn_type), intent(inout) :: this
        integer, intent(in) :: mdim
        integer, intent(in) :: ndim
        integer, intent(in) :: jn !< 法向点数
        integer, parameter :: NORMAL=0, TRANSPOSE=2
        integer :: eqkind

        print*, 'NPSE is creating.....'

        this%jn=jn
        print*, 'jn=', jn
        eqkind=0
        !call this%BFOPNorm%Create(this%jn, eqkind)
        !call this%DisOPNorm%Create(this%jn)
        allocate(this%Eta(this%jn))

        this%mdim=mdim
        this%ndim=ndim

        print*, 'm=', mdim, 'n=', ndim

        if(.not. allocated(this%isWithAlpha)) then
          allocate(this%isWithAlpha(0:this%mdim, -this%ndim:this%ndim))
          allocate(this%bctype(5,2,0:this%mdim, -this%ndim:this%ndim))
        endif

        print*, 'LPSE solver is creating....'
        call this%LPSE%SetSolver(method=1001)
        call this%LPSE%Create(jn)
        if(.not. allocated(this%nonlinearTerm)) then
          allocate(this%nonlinearTerm(5, jn, 0:this%mdim, -this%ndim:this%ndim))
          allocate(this%nonlinearTermOld(5, jn, 0:this%mdim, -this%ndim:this%ndim))
        endif
        this%nonlinearterm=0.0d0
        this%nonlinearTermOld=0.0d0
        if(.not. allocated(this%IntergalAlf)) &
        &   allocate(this%IntergalAlf(2, 0:this%mdim, -this%ndim:this%ndim))
        this%IntergalAlf=0.0d0
        if(.not. allocated(this%AlfOld)) &
        &   allocate(this%AlfOld(0:this%mdim, -this%ndim:this%ndim))
        this%AlfOld=0.0d0
        if(.not. allocated(this%AlfNew)) &
        &   allocate(this%AlfNew(0:this%mdim, -this%ndim:this%ndim))
        this%AlfNew=0.0d0
#IFDEF DEBUG
        print*, this%nonlinearSolver%isCreate()
        pause
#ENDIF
        if(.not. this%nonlinearSolver%isCreate()) &
        & call this%nonlinearSolver%Create(this%jn, this%mdim, this%ndim, this%nonlinearTerm)

        this%isCheckAmp=.False.

    end subroutine Create

    !> 设置当前流向站位
    subroutine set_iloc(this, iloc)

        implicit none
        integer, intent(in) :: iloc !< 当前流向站位
        class(npse_eqn_type), intent(inout) :: this

        this%iloc=iloc

    end subroutine set_iloc

    !> 输出当前站位序号
    subroutine PrintILoc(this)

        implicit none
        class(npse_eqn_type), intent(in) :: this

        write(*, *)"The location index now is", this%iloc

    end subroutine PrintILoc

end module mod_npse_eqn
