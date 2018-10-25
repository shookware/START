!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_nonlinear_term.f90
!> @file
!> @breif 非线性项计算.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: mod_nonlinear
!> @breif 非线性项模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2018-07-08 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2018-07-08
!------------------------------------------------------------------------------
module mod_nonlinear

  use penf, only: R_P
  use mod_BF_normal
  use mod_lpse_dis_normal
  use mod_local_normal_coordinate
  use mod_difference
  use mod_spectral_method
  use mod_cfgio_adapter, only: omega0, beta0
  implicit none

  public :: nonlinear_type
  private

  !> 各种临时扰动变量(谱空间)
  type dis_spt_type

    complex(R_P), allocatable :: rho(:, :, :)
    complex(R_P), allocatable :: u(:, :, :)
    complex(R_P), allocatable :: v(:, :, :)
    complex(R_P), allocatable :: w(:, :, :)
    complex(R_P), allocatable :: T(:, :, :)

    complex(R_P), allocatable :: rhox(:, :, :)
    complex(R_P), allocatable :: ux(:, :, :)
    complex(R_P), allocatable :: vx(:, :, :)
    complex(R_P), allocatable :: wx(:, :, :)
    complex(R_P), allocatable :: Tx(:, :, :)

    complex(R_P), allocatable :: rhoy(:, :, :)
    complex(R_P), allocatable :: uy(:, :, :)
    complex(R_P), allocatable :: vy(:, :, :)
    complex(R_P), allocatable :: wy(:, :, :)
    complex(R_P), allocatable :: Ty(:, :, :)

    complex(R_P), allocatable :: rhoz(:, :, :)
    complex(R_P), allocatable :: uz(:, :, :)
    complex(R_P), allocatable :: vz(:, :, :)
    complex(R_P), allocatable :: wz(:, :, :)
    complex(R_P), allocatable :: Tz(:, :, :)

    complex(R_P), allocatable :: rhot(:, :, :)
    complex(R_P), allocatable :: ut(:, :, :)
    complex(R_P), allocatable :: vt(:, :, :)
    complex(R_P), allocatable :: wt(:, :, :)
    complex(R_P), allocatable :: Tt(:, :, :)

    complex(R_P), allocatable :: uxx(:, :, :)
    complex(R_P), allocatable :: vxx(:, :, :)
    complex(R_P), allocatable :: wxx(:, :, :)
    complex(R_P), allocatable :: Txx(:, :, :)

    complex(R_P), allocatable :: uyy(:, :, :)
    complex(R_P), allocatable :: vyy(:, :, :)
    complex(R_P), allocatable :: wyy(:, :, :)
    complex(R_P), allocatable :: Tyy(:, :, :)

    complex(R_P), allocatable :: uzz(:, :, :)
    complex(R_P), allocatable :: vzz(:, :, :)
    complex(R_P), allocatable :: wzz(:, :, :)
    complex(R_P), allocatable :: Tzz(:, :, :)

    complex(R_P), allocatable :: uxy(:, :, :)
    complex(R_P), allocatable :: vxy(:, :, :)
    complex(R_P), allocatable :: uxz(:, :, :)
    complex(R_P), allocatable :: wxz(:, :, :)
    complex(R_P), allocatable :: vyz(:, :, :)
    complex(R_P), allocatable :: wyz(:, :, :)

  contains

    procedure :: Create=>Create_dis_spt
    procedure :: Finalize=> Finalize_dis_spt

  end type dis_spt_type

  !> 各种临时扰动变量(谱空间)
  type dis_phy_type

    real(R_P), allocatable :: rho(:, :)
    real(R_P), allocatable :: u(:, :)
    real(R_P), allocatable :: v(:, :)
    real(R_P), allocatable :: w(:, :)
    real(R_P), allocatable :: T(:, :)

    real(R_P), allocatable :: rhox(:, :)
    real(R_P), allocatable :: ux(:, :)
    real(R_P), allocatable :: vx(:, :)
    real(R_P), allocatable :: wx(:, :)
    real(R_P), allocatable :: Tx(:, :)

    real(R_P), allocatable :: rhoy(:, :)
    real(R_P), allocatable :: uy(:, :)
    real(R_P), allocatable :: vy(:, :)
    real(R_P), allocatable :: wy(:, :)
    real(R_P), allocatable :: Ty(:, :)

    real(R_P), allocatable :: rhoz(:, :)
    real(R_P), allocatable :: uz(:, :)
    real(R_P), allocatable :: vz(:, :)
    real(R_P), allocatable :: wz(:, :)
    real(R_P), allocatable :: Tz(:, :)

    real(R_P), allocatable :: rhot(:, :)
    real(R_P), allocatable :: ut(:, :)
    real(R_P), allocatable :: vt(:, :)
    real(R_P), allocatable :: wt(:, :)
    real(R_P), allocatable :: Tt(:, :)

    real(R_P), allocatable :: uxx(:, :)
    real(R_P), allocatable :: vxx(:, :)
    real(R_P), allocatable :: wxx(:, :)
    real(R_P), allocatable :: Txx(:, :)

    real(R_P), allocatable :: uyy(:, :)
    real(R_P), allocatable :: vyy(:, :)
    real(R_P), allocatable :: wyy(:, :)
    real(R_P), allocatable :: Tyy(:, :)

    real(R_P), allocatable :: uzz(:, :)
    real(R_P), allocatable :: vzz(:, :)
    real(R_P), allocatable :: wzz(:, :)
    real(R_P), allocatable :: Tzz(:, :)

    real(R_P), allocatable :: uxy(:, :)
    real(R_P), allocatable :: vxy(:, :)
    real(R_P), allocatable :: uxz(:, :)
    real(R_P), allocatable :: wxz(:, :)
    real(R_P), allocatable :: vyz(:, :)
    real(R_P), allocatable :: wyz(:, :)
    real(R_P), allocatable :: FTerm(:, :, :)
  contains
    procedure :: Create=>Create_dis_phy
    procedure :: Finalize=>Finalize_dis_phy

  end type dis_phy_type

  !> 非线性类型
  type nonlinear_type
    private
    integer, private :: jn !<法向点数
    integer, private :: mdim !<时间离散谱宽度
    integer, private :: ndim !<展向离散谱宽度
    type(BF_normal_type), Pointer, private :: BFNorm !<基本流流向某站位法向分布
    type(lpse_dis_normal_type), Pointer, private :: DisNorm(:, :) !<扰动流向某站位法向分布
    type(lpse_dis_normal_type), Pointer, private :: DisNormFront(:, :) !<扰动流向前站位位法向分布
    type(local_normal_coordinate_type), pointer, private :: NormCoord !< 某流向展向站位的局部坐标系
    type(Norm_Diff_Coef_type), pointer, private :: NormCoef !< 局部差分系数
    complex(R_P), allocatable, private :: IntergalAlf(:, :)
    complex(R_P), allocatable, private :: Alf(:, :)

    complex(R_P), private, pointer :: nonlinearTerm(:, :, :, :) !<非线性项
    type(spectral_method_type), private :: SPT
    type(dis_spt_type), private :: dis_spt
    type(dis_phy_type), private :: dis_phy
    !complex(R_P)

  contains
    procedure :: Create !<创建非线性项
    procedure :: SolveNonlinearTerm !<求解非线性项
    procedure :: Finalize !< 析构函数
    procedure :: IsCreate !< 判断是否创建该非线性项

    generic, private :: ComputeF => computeF1

    procedure, private :: Solve !< 求解过程
    procedure, private :: spectral_to_physical !<谱空间变换到物理空间
    procedure, private :: physical_to_spectral !<物理空间变换到谱空间
    procedure, private :: ComputeF1 !<计算非线性项
    procedure, private :: ComputeF2 !<计算非线性项

  end type nonlinear_type

  contains

    subroutine Create_dis_spt(this, m, n, jn)
      implicit none
      class(dis_spt_type), intent(inout) :: this
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, intent(in) ::jn

      allocate(this%rho(0:m, -n:n, jn))
      allocate(this%u(0:m, -n:n, jn))
      allocate(this%v(0:m, -n:n, jn))
      allocate(this%w(0:m, -n:n, jn))
      allocate(this%T(0:m, -n:n, jn))

      allocate(this%rhox(0:m, -n:n, jn))
      allocate(this%ux(0:m, -n:n, jn))
      allocate(this%vx(0:m, -n:n, jn))
      allocate(this%wx(0:m, -n:n, jn))
      allocate(this%Tx(0:m, -n:n, jn))

      allocate(this%rhoy(0:m, -n:n, jn))
      allocate(this%uy(0:m, -n:n, jn))
      allocate(this%vy(0:m, -n:n, jn))
      allocate(this%wy(0:m, -n:n, jn))
      allocate(this%Ty(0:m, -n:n, jn))

      allocate(this%rhoz(0:m, -n:n, jn))
      allocate(this%uz(0:m, -n:n, jn))
      allocate(this%vz(0:m, -n:n, jn))
      allocate(this%wz(0:m, -n:n, jn))
      allocate(this%Tz(0:m, -n:n, jn))

      allocate(this%rhot(0:m, -n:n, jn))
      allocate(this%ut(0:m, -n:n, jn))
      allocate(this%vt(0:m, -n:n, jn))
      allocate(this%wt(0:m, -n:n, jn))
      allocate(this%Tt(0:m, -n:n, jn))

      allocate(this%uxx(0:m, -n:n, jn))
      allocate(this%vxx(0:m, -n:n, jn))
      allocate(this%wxx(0:m, -n:n, jn))
      allocate(this%Txx(0:m, -n:n, jn))

      allocate(this%uyy(0:m, -n:n, jn))
      allocate(this%vyy(0:m, -n:n, jn))
      allocate(this%wyy(0:m, -n:n, jn))
      allocate(this%Tyy(0:m, -n:n, jn))

      allocate(this%uzz(0:m, -n:n, jn))
      allocate(this%vzz(0:m, -n:n, jn))
      allocate(this%wzz(0:m, -n:n, jn))
      allocate(this%Tzz(0:m, -n:n, jn))

      allocate(this%uxy(0:m, -n:n, jn))
      allocate(this%vxy(0:m, -n:n, jn))
      allocate(this%uxz(0:m, -n:n, jn))
      allocate(this%wxz(0:m, -n:n, jn))
      allocate(this%vyz(0:m, -n:n, jn))
      allocate(this%wyz(0:m, -n:n, jn))

    end subroutine Create_dis_spt

    subroutine Create_dis_phy(this, m, n)
      implicit none
      class(dis_phy_type), intent(inout) :: this
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer :: msp, nsp

      msp=(m+1)*8-1
      nsp=(n+1)*8-1


      allocate(this%rho(0:msp, 0:nsp))
      allocate(this%u(0:msp, 0:nsp))
      allocate(this%v(0:msp, 0:nsp))
      allocate(this%w(0:msp, 0:nsp))
      allocate(this%T(0:msp, 0:nsp))

      allocate(this%rhox(0:msp, 0:nsp))
      allocate(this%ux(0:msp, 0:nsp))
      allocate(this%vx(0:msp, 0:nsp))
      allocate(this%wx(0:msp, 0:nsp))
      allocate(this%Tx(0:msp, 0:nsp))

      allocate(this%rhoy(0:msp, 0:nsp))
      allocate(this%uy(0:msp, 0:nsp))
      allocate(this%vy(0:msp, 0:nsp))
      allocate(this%wy(0:msp, 0:nsp))
      allocate(this%Ty(0:msp, 0:nsp))

      allocate(this%rhoz(0:msp, 0:nsp))
      allocate(this%uz(0:msp, 0:nsp))
      allocate(this%vz(0:msp, 0:nsp))
      allocate(this%wz(0:msp, 0:nsp))
      allocate(this%Tz(0:msp, 0:nsp))

      allocate(this%rhot(0:msp, 0:nsp))
      allocate(this%ut(0:msp, 0:nsp))
      allocate(this%vt(0:msp, 0:nsp))
      allocate(this%wt(0:msp, 0:nsp))
      allocate(this%Tt(0:msp, 0:nsp))

      allocate(this%uxx(0:msp, 0:nsp))
      allocate(this%vxx(0:msp, 0:nsp))
      allocate(this%wxx(0:msp, 0:nsp))
      allocate(this%Txx(0:msp, 0:nsp))

      allocate(this%uyy(0:msp, 0:nsp))
      allocate(this%vyy(0:msp, 0:nsp))
      allocate(this%wyy(0:msp, 0:nsp))
      allocate(this%Tyy(0:msp, 0:nsp))

      allocate(this%uzz(0:msp, 0:nsp))
      allocate(this%vzz(0:msp, 0:nsp))
      allocate(this%wzz(0:msp, 0:nsp))
      allocate(this%Tzz(0:msp, 0:nsp))

      allocate(this%uxy(0:msp, 0:nsp))
      allocate(this%vxy(0:msp, 0:nsp))
      allocate(this%uxz(0:msp, 0:nsp))
      allocate(this%wxz(0:msp, 0:nsp))
      allocate(this%vyz(0:msp, 0:nsp))
      allocate(this%wyz(0:msp, 0:nsp))

      allocate(this%FTerm(5, 0:msp, 0:nsp))

    end subroutine Create_dis_phy

    subroutine Finalize_dis_spt(this)

      implicit none
      class(dis_spt_type), intent(inout) :: this

      deallocate(this%rho)
      deallocate(this%u)
      deallocate(this%v)
      deallocate(this%w)
      deallocate(this%T)

      deallocate(this%rhox)
      deallocate(this%ux)
      deallocate(this%vx)
      deallocate(this%wx)
      deallocate(this%Tx)

      deallocate(this%rhoy)
      deallocate(this%uy)
      deallocate(this%vy)
      deallocate(this%wy)
      deallocate(this%Ty)

      deallocate(this%rhoz)
      deallocate(this%uz)
      deallocate(this%vz)
      deallocate(this%wz)
      deallocate(this%Tz)

      deallocate(this%rhot)
      deallocate(this%ut)
      deallocate(this%vt)
      deallocate(this%wt)
      deallocate(this%Tt)

      deallocate(this%uxx)
      deallocate(this%vxx)
      deallocate(this%wxx)
      deallocate(this%Txx)

      deallocate(this%uyy)
      deallocate(this%vyy)
      deallocate(this%wyy)
      deallocate(this%Tyy)

      deallocate(this%uzz)
      deallocate(this%vzz)
      deallocate(this%wzz)
      deallocate(this%Tzz)

      deallocate(this%uxy)
      deallocate(this%vxy)
      deallocate(this%uxz)
      deallocate(this%wxz)
      deallocate(this%vyz)
      deallocate(this%wyz)

    end subroutine Finalize_dis_spt

    subroutine Finalize_dis_phy(this)

      implicit none
      class(dis_phy_type), intent(inout) :: this

      deallocate(this%rho)
      deallocate(this%u)
      deallocate(this%v)
      deallocate(this%w)
      deallocate(this%T)

      deallocate(this%rhox)
      deallocate(this%ux)
      deallocate(this%vx)
      deallocate(this%wx)
      deallocate(this%Tx)

      deallocate(this%rhoy)
      deallocate(this%uy)
      deallocate(this%vy)
      deallocate(this%wy)
      deallocate(this%Ty)

      deallocate(this%rhoz)
      deallocate(this%uz)
      deallocate(this%vz)
      deallocate(this%wz)
      deallocate(this%Tz)

      deallocate(this%rhot)
      deallocate(this%ut)
      deallocate(this%vt)
      deallocate(this%wt)
      deallocate(this%Tt)

      deallocate(this%uxx)
      deallocate(this%vxx)
      deallocate(this%wxx)
      deallocate(this%Txx)

      deallocate(this%uyy)
      deallocate(this%vyy)
      deallocate(this%wyy)
      deallocate(this%Tyy)

      deallocate(this%uzz)
      deallocate(this%vzz)
      deallocate(this%wzz)
      deallocate(this%Tzz)

      deallocate(this%uxy)
      deallocate(this%vxy)
      deallocate(this%uxz)
      deallocate(this%wxz)
      deallocate(this%vyz)
      deallocate(this%wyz)

      deallocate(this%FTerm)

    end subroutine Finalize_dis_phy

    !>创建非线性项类型
    subroutine Create(this, jn, m_dim, n_dim, nonlinearTerm)

      implicit none
      class(nonlinear_type), intent(inout) :: this
      integer, intent(in) :: jn !<法向离散点数
      integer, intent(in) :: m_dim !<时间离散点
      integer, intent(in) :: n_dim !<展向空间离散点
      complex(R_P), intent(in), target :: nonlinearTerm(5, jn, 0:m_dim, -n_dim:n_dim)

      this%jn=jn
      this%mdim=m_dim
      this%ndim=n_dim

#ifdef DEBUG
        print*, 'nonlinear term parameters:'
        print*, jn, m_dim, n_dim
#endif
      !if(allocated(this%nonlinearTerm)) deallocate(this%nonlinearTerm)
      !if(.not. allocated(this%nonlinearTerm)) &
      !& allocate(this%nonlinearTerm(5, jn, 0:m_dim, -n_dim:n_dim))

      this%nonlinearTerm=>nonlinearTerm
      !
      ! block
      !   write(*, *)jn, m_dim, n_dim
      !   write(*, *) 'check bound'
      !   write(*, *)LBOUND(nonlinearTerm, dim=3)
      !   write(*, *)UBOUND(nonlinearTerm, dim=3)
      !   write(*, *)LBOUND(nonlinearTerm, dim=4)
      !   write(*, *)UBOUND(nonlinearTerm, dim=4)
      !   write(*, *)LBOUND(this%nonlinearTerm, dim=3)
      !   write(*, *)UBOUND(this%nonlinearTerm, dim=3)
      !   write(*, *)LBOUND(this%nonlinearTerm, dim=4)
      !   write(*, *)UBOUND(this%nonlinearTerm, dim=4)
      !   pause
      !
      !   endblock


      if(.not. this%SPT%isCreate) call this%SPT%Create(omega0, beta0, m_dim, n_dim)
      call this%dis_spt%Create(m_dim, n_dim, jn)
#ifdef DEBUG
      ! print*, size(this%dis_spt%rho)
      ! print*, size(this%dis_spt%rho, dim=1)
      ! print*, size(this%dis_spt%rho, dim=2)
      ! print*, size(this%dis_spt%rho, dim=3)
      ! pause
#endif

      call this%dis_phy%Create(m_dim, n_dim)
      allocate(this%IntergalAlf(0:this%mdim, -this%ndim:this%ndim))
      allocate(this%Alf(0:this%mdim, -this%ndim:this%ndim))

    end subroutine Create

    !> 求解非线性项
    subroutine SolveNonlinearTerm(this, DisNorm, DisNormFront, &
      & BFNorm, NormCoord, NormCoef, IntergalAlf, Alf)

      implicit none
      class(nonlinear_type), intent(inout) :: this
      type(lpse_dis_normal_type), target, intent(in) :: DisNorm(0:, -this%ndim:)
      type(lpse_dis_normal_type), target, intent(in) :: DisNormFront(0:, -this%ndim:)
      type(BF_normal_type), target, intent(in) :: BFNorm
      type(local_normal_coordinate_type), target, intent(in) :: NormCoord
      type(Norm_Diff_Coef_type), target, intent(in) :: NormCoef
      complex(R_P), intent(in) :: IntergalAlf(0:this%mdim, -this%ndim:this%ndim)
      complex(R_P), intent(in) :: Alf(0:this%mdim, -this%ndim:this%ndim)

      this%DisNorm => DisNorm
      this%DisNormFront => DisNormFront
      this%BFNorm => BFNorm
      this%NormCoord => NormCoord
      this%NormCoef => NormCoef

      this%IntergalAlf=IntergalAlf
      this%Alf=Alf

      ! print*, 'Alpha....'
      ! print*, DisNorm(1, 0)%GetAlpha()
      ! print*, DisNormFront(1, 0)%getAlpha()

      call this%Solve()

    end subroutine SolveNonlinearTerm

    !> 求解非线性项
    subroutine Solve(this)

      use mod_parameter, only: CPLI
      use mod_local_coordinate, only: local_coordinate_type

      implicit none
      class(nonlinear_type), intent(inout) :: this
      complex(R_P), allocatable, save :: DisFlux(:, :, :, :)
      complex(R_P), allocatable, save :: DisFluxFront(:, :, :, :)
      integer :: mCount, nCount, j, l
      real(R_P) :: coef(2)
      complex(R_P) :: Amp(0:this%mdim, -this%ndim:this%ndim)
      type(local_coordinate_type) :: Coord
      type(bf_point_type) :: BFPoint
      complex(R_P) :: AlfFront(0:this%mdim, -this%ndim:this%ndim)
      complex(R_P) :: dAlfdx(0:this%mdim, -this%ndim:this%ndim)

#ifdef DEBUG
      print*, 'subroutine check'
      block



      endblock

      pause
#endif

      if(.not. allocated(disflux)) then
        allocate(DisFlux(5, this%jn, 0:this%mdim, -this%ndim:this%ndim))
        allocate(DisFluxFront(5, this%jn, 0:this%mdim, -this%ndim:this%ndim))
      else
        DisFlux=0.0d0
        DisFluxFront=0.0d0
      endif

#ifdef DEBUG

	! print*, coef
#endif


      coef=this%NormCoef%Coef_di

#ifdef DEBUG
      ! pause
      ! print*, 'nonlinear start'
      ! pause
      ! !print*, size(this%dis_spt%rho)
      ! !print*, size(this%dis_spt%rho, dim=1)
      ! !print*, size(this%dis_spt%rho, dim=2)
      ! !print*, size(this%dis_spt%rho, dim=3)
      ! pause
#endif

      associate(dis_spt=>this%dis_spt)

      do nCount=-this%ndim, this%ndim
        do mCount=0, this%mdim
#ifdef DEBUG
     !       print*, '1'
     !       print*, this%DisNorm(mCount, nCount)%jn
	   ! pause 'this is done'
#endif
          call this%DisNorm(mCount, nCount)%Get(DisFlux(:, :, mCount, nCount))
          call this%DisNormFront(mCount, nCount)%Get(DisFluxFront(:, :, mCount, nCount))

          dis_spt%rho(mCount, nCount, :)=Disflux(1, :, mCount, nCount)
          dis_spt%u  (mCount, nCount, :)=Disflux(2, :, mCount, nCount)
          dis_spt%v  (mCount, nCount, :)=Disflux(3, :, mCount, nCount)
          dis_spt%w  (mCount, nCount, :)=Disflux(4, :, mCount, nCount)
          dis_spt%T  (mCount, nCount, :)=Disflux(5, :, mCount, nCount)



#ifdef DEBUG
           ! print*, '2'
           ! print*, size(dis_spt%rho)
           ! print*, size(dis_spt%rho, dim=1)
           ! print*, size(dis_spt%rho, dim=2)
           ! print*, size(dis_spt%rho, dim=3)
           ! print*, size(this%dis_spt%rho)
           ! print*, size(this%dis_spt%rho, dim=1)
           ! print*, size(this%dis_spt%rho, dim=2)
           ! print*, size(this%dis_spt%rho, dim=3)
           ! print*, dis_spt%rho(mCount, nCount, :)
           ! pause
#endif

           DisFlux(:, :, mCount, nCount)=Coef(1)*DisFluxFront(:, :, mCount, nCount) &
                        &              + Coef(2)*DisFlux(:, :, mCount, nCount)
           dis_spt%rhox(mCount, nCount, :)=DisFlux(1, :, mCount, nCount)
           dis_spt%ux  (mCount, nCount, :)=DisFlux(2, :, mCount, nCount)
           dis_spt%vx  (mCount, nCount, :)=DisFlux(3, :, mCount, nCount)
           dis_spt%wx  (mCount, nCount, :)=DisFlux(4, :, mCount, nCount)
           dis_spt%Tx  (mCount, nCount, :)=DisFlux(5, :, mCount, nCount)

           AlfFront(mCount, nCount)=this%DisNormFront(mCount, nCount)%getAlpha()
           dAlfdx(mCount, nCount)=Coef(1)*AlfFront(mCount, nCount) &
                        &        + Coef(2)*this%Alf(mCount, nCount)

           dis_spt%uxx(mCount, nCount, :)=CPLI*2.0d0*this%Alf(mCount, nCount)*dis_spt%ux(mCount, nCount, :) + &
                          &               CPLI*dAlfdx(mCount, nCount)*dis_spt%u(mCount, nCount, :)
           dis_spt%vxx(mCount, nCount, :)=CPLI*2.0d0*this%Alf(mCount, nCount)*dis_spt%vx(mCount, nCount, :) + &
                          &               CPLI*dAlfdx(mCount, nCount)*dis_spt%v(mCount, nCount, :)
           dis_spt%wxx(mCount, nCount, :)=CPLI*2.0d0*this%Alf(mCount, nCount)*dis_spt%wx(mCount, nCount, :) + &
                          &               CPLI*dAlfdx(mCount, nCount)*dis_spt%w(mCount, nCount, :)
           dis_spt%Txx(mCount, nCount, :)=CPLI*2.0d0*this%Alf(mCount, nCount)*dis_spt%Tx(mCount, nCount, :) + &
                          &               CPLI*dAlfdx(mCount, nCount)*dis_spt%T(mCount, nCount, :)


          dis_spt%rhox(mCount, nCount, :)=dis_spt%rhox(mCount, nCount, :)+&
                        & CPLI*this%Alf(mCount, nCount)*dis_spt%rho(mCount, nCount, :)
          dis_spt%ux  (mCount, nCount, :)=dis_spt%ux  (mCount, nCount, :)+&
                        & CPLI*this%Alf(mCount, nCount)*dis_spt%u  (mCount, nCount, :)
          dis_spt%vx  (mCount, nCount, :)=dis_spt%vx  (mCount, nCount, :)+&
                        & CPLI*this%Alf(mCount, nCount)*dis_spt%v  (mCount, nCount, :)
          dis_spt%wx  (mCount, nCount, :)=dis_spt%wx  (mCount, nCount, :)+&
                        & CPLI*this%Alf(mCount, nCount)*dis_spt%w  (mCount, nCount, :)
          dis_spt%Tx  (mCount, nCount, :)=dis_spt%Tx  (mCount, nCount, :)+&
                        & CPLI*this%Alf(mCount, nCount)*dis_spt%T  (mCount, nCount, :)


          ! block
          !
          !   print*, this%Alf(:, 0)
          !   pause
          !
          ! endblock

          dis_spt%rhoy(mCount, nCount, :)= &
                        & this%NormCoef%dy(dis_spt%rho(mCount, nCount, :))
          dis_spt%uy  (mCount, nCount, :)= &
                        & this%NormCoef%dy(dis_spt%u  (mCount, nCount, :))
          dis_spt%vy  (mCount, nCount, :)= &
                        & this%NormCoef%dy(dis_spt%v  (mCount, nCount, :))
          dis_spt%wy  (mCount, nCount, :)= &
                        & this%NormCoef%dy(dis_spt%w  (mCount, nCount, :))
          dis_spt%Ty  (mCount, nCount, :)= &
                        & this%NormCoef%dy(dis_spt%T  (mCount, nCount, :))

          dis_spt%uyy  (mCount, nCount, :)= &
                        & this%NormCoef%dyy(dis_spt%u  (mCount, nCount, :))
          dis_spt%vyy  (mCount, nCount, :)= &
                        & this%NormCoef%dyy(dis_spt%v  (mCount, nCount, :))
          dis_spt%wyy  (mCount, nCount, :)= &
                        & this%NormCoef%dyy(dis_spt%w  (mCount, nCount, :))
          dis_spt%Tyy  (mCount, nCount, :)= &
                        & this%NormCoef%dyy(dis_spt%T  (mCount, nCount, :))

!! @TODO_check correct!!!!
          dis_spt%uxx(mCount, nCount, :)= dis_spt%uxx(mCount, nCount, :) &
                & -this%Alf(mCount, nCount)**2*dis_spt%u(mCount, nCount, :)
          dis_spt%vxx(mCount, nCount, :)= dis_spt%vxx(mCount, nCount, :) &
                & -this%Alf(mCount, nCount)**2*dis_spt%v(mCount, nCount, :)
          dis_spt%wxx(mCount, nCount, :)= dis_spt%wxx(mCount, nCount, :) &
                & -this%Alf(mCount, nCount)**2*dis_spt%w(mCount, nCount, :)
          dis_spt%Txx(mCount, nCount, :)= dis_spt%Txx(mCount, nCount, :) &
                & -this%Alf(mCount, nCount)**2*dis_spt%T(mCount, nCount, :)

          dis_spt%uxy(mCount, nCount, :)=&
                & CPLI*this%Alf(mCount, nCount)*dis_spt%uy(mCount, nCount, :)
          dis_spt%vxy(mCount, nCount, :)=&
                & CPLI*this%Alf(mCount, nCount)*dis_spt%vy(mCount, nCount, :)
        enddo
      enddo

      do j=1, this%jn
        do mCount=0, this%mdim
          dis_spt%rhoz(mCount, :, j)=this%SPT%dz(dis_spt%rho(mCount, :, j))
          dis_spt%uz  (mCount, :, j)=this%SPT%dz(dis_spt%u  (mCount, :, j))
          dis_spt%vz  (mCount, :, j)=this%SPT%dz(dis_spt%v  (mCount, :, j))
          dis_spt%wz  (mCount, :, j)=this%SPT%dz(dis_spt%w  (mCount, :, j))
          dis_spt%Tz  (mCount, :, j)=this%SPT%dz(dis_spt%T  (mCount, :, j))

          dis_spt%uzz (mCount, :, j)=this%SPT%dzz(dis_spt%u  (mCount, :, j))
          dis_spt%vzz (mCount, :, j)=this%SPT%dzz(dis_spt%v  (mCount, :, j))
          dis_spt%wzz (mCount, :, j)=this%SPT%dzz(dis_spt%w  (mCount, :, j))
          dis_spt%Tzz (mCount, :, j)=this%SPT%dzz(dis_spt%T  (mCount, :, j))

          dis_spt%uxz(mCount, :, j)=this%SPT%dz(dis_spt%ux(mCount, :, j))
          dis_spt%wxz(mCount, :, j)=this%SPT%dz(dis_spt%wx(mCount, :, j))
          dis_spt%vyz(mCount, :, j)=this%SPT%dz(dis_spt%vy(mCount, :, j))
          dis_spt%wyz(mCount, :, j)=this%SPT%dz(dis_spt%wy(mCount, :, j))
        enddo
      enddo

      do j=1, this%jn
        do nCount=-this%ndim, this%ndim
          dis_spt%rhot(:, nCount, j)=this%SPT%dt(dis_spt%rho(:, nCount, j))
          dis_spt%ut  (:, nCount, j)=this%SPT%dt(dis_spt%u  (:, nCount, j))
          dis_spt%vt  (:, nCount, j)=this%SPT%dt(dis_spt%v  (:, nCount, j))
          dis_spt%wt  (:, nCount, j)=this%SPT%dt(dis_spt%w  (:, nCount, j))
          dis_spt%Tt  (:, nCount, j)=this%SPT%dt(dis_spt%T  (:, nCount, j))
        enddo
      enddo

      Amp=exp(CPLI*this%IntergalAlf)
#ifdef DEBUG
      print*, 'this%intergalAlf'
      do mCount=0, this%mdim
        do nCount=-this%ndim, this%ndim
            print*, this%IntergalAlf(mCount, nCount)
        enddo
      enddo
#endif
      do j=1, this%jn
        dis_spt%rho(:, :, j)=dis_spt%rho(:, :, j)*Amp
        dis_spt%u  (:, :, j)=dis_spt%u  (:, :, j)*Amp
        dis_spt%v  (:, :, j)=dis_spt%v  (:, :, j)*Amp
        dis_spt%w  (:, :, j)=dis_spt%w  (:, :, j)*Amp
        dis_spt%T  (:, :, j)=dis_spt%T  (:, :, j)*Amp


        dis_spt%rhox(:, :, j)=dis_spt%rhox(:, :, j)*Amp
        dis_spt%ux  (:, :, j)=dis_spt%ux  (:, :, j)*Amp
        dis_spt%vx  (:, :, j)=dis_spt%vx  (:, :, j)*Amp
        dis_spt%wx  (:, :, j)=dis_spt%wx  (:, :, j)*Amp
        dis_spt%Tx  (:, :, j)=dis_spt%Tx  (:, :, j)*Amp

        dis_spt%rhoy(:, :, j)=dis_spt%rhoy(:, :, j)*Amp
        dis_spt%uy  (:, :, j)=dis_spt%uy  (:, :, j)*Amp
        dis_spt%vy  (:, :, j)=dis_spt%vy  (:, :, j)*Amp
        dis_spt%wy  (:, :, j)=dis_spt%wy  (:, :, j)*Amp
        dis_spt%Ty  (:, :, j)=dis_spt%Ty  (:, :, j)*Amp

        dis_spt%rhoz(:, :, j)=dis_spt%rhoz(:, :, j)*Amp
        dis_spt%uz  (:, :, j)=dis_spt%uz  (:, :, j)*Amp
        dis_spt%vz  (:, :, j)=dis_spt%vz  (:, :, j)*Amp
        dis_spt%wz  (:, :, j)=dis_spt%wz  (:, :, j)*Amp
        dis_spt%Tz  (:, :, j)=dis_spt%Tz  (:, :, j)*Amp

        dis_spt%rhot(:, :, j)=dis_spt%rhot(:, :, j)*Amp
        dis_spt%ut  (:, :, j)=dis_spt%ut  (:, :, j)*Amp
        dis_spt%vt  (:, :, j)=dis_spt%vt  (:, :, j)*Amp
        dis_spt%wt  (:, :, j)=dis_spt%wt  (:, :, j)*Amp
        dis_spt%Tt  (:, :, j)=dis_spt%Tt  (:, :, j)*Amp

        dis_spt%uxx  (:, :, j)=dis_spt%uxx  (:, :, j)*Amp
        dis_spt%vxx  (:, :, j)=dis_spt%vxx  (:, :, j)*Amp
        dis_spt%wxx  (:, :, j)=dis_spt%wxx  (:, :, j)*Amp
        dis_spt%Txx  (:, :, j)=dis_spt%Txx  (:, :, j)*Amp

        dis_spt%uyy  (:, :, j)=dis_spt%uyy  (:, :, j)*Amp
        dis_spt%vyy  (:, :, j)=dis_spt%vyy  (:, :, j)*Amp
        dis_spt%wyy  (:, :, j)=dis_spt%wyy  (:, :, j)*Amp
        dis_spt%Tyy  (:, :, j)=dis_spt%Tyy  (:, :, j)*Amp

        dis_spt%uzz  (:, :, j)=dis_spt%uzz  (:, :, j)*Amp
        dis_spt%vzz  (:, :, j)=dis_spt%vzz  (:, :, j)*Amp
        dis_spt%wzz  (:, :, j)=dis_spt%wzz  (:, :, j)*Amp
        dis_spt%Tzz  (:, :, j)=dis_spt%Tzz  (:, :, j)*Amp

        dis_spt%uxy (:, :, j)= dis_spt%uxy (:, :, j)*Amp
        dis_spt%vxy (:, :, j)= dis_spt%vxy (:, :, j)*Amp
        dis_spt%uxz (:, :, j)= dis_spt%uxz (:, :, j)*Amp
        dis_spt%wxz (:, :, j)= dis_spt%wxz (:, :, j)*Amp
        dis_spt%vyz (:, :, j)= dis_spt%vyz (:, :, j)*Amp
        dis_spt%wyz (:, :, j)= dis_spt%wyz (:, :, j)*Amp

      enddo
#ifdef DEBUG
      print*, 'dis_spt_u'
      do mCount=0, this%mdim
        do nCount=-this%ndim, this%ndim
            print*, mCount, nCount, abs(amp(mCount, nCount)), maxval(abs(dis_spt%u(mCount, ncount, :)))
        ENDDO
      ENDDO

      pause
#endif

       ! block
    ! print*, 'r', maxval(abs(dis_spt%rho(1, 0, :))), maxloc(abs(dis_spt%rho(1, 0, :)))
		! print*, 'u', maxval(abs(dis_spt%u  (1, 0, :))), maxloc(abs(dis_spt%u  (1, 0, :)))
		! print*, 'v', maxval(abs(dis_spt%v  (1, 0, :))), maxloc(abs(dis_spt%v  (1, 0, :)))
		! print*, 'T', maxval(abs(dis_spt%t  (1, 0, :))), maxloc(abs(dis_spt%t  (1, 0, :)))
    !
		! print*, 'rx', maxval(abs(dis_spt%rhox(1, 0, :))), maxloc(abs(dis_spt%rhox(1, 0, :)))
		! print*, 'ux', maxval(abs(dis_spt%ux  (1, 0, :))), maxloc(abs(dis_spt%ux  (1, 0, :)))
		! print*, 'vx', maxval(abs(dis_spt%vx  (1, 0, :))), maxloc(abs(dis_spt%vx  (1, 0, :)))
		! print*, 'Tx', maxval(abs(dis_spt%tx  (1, 0, :))), maxloc(abs(dis_spt%tx  (1, 0, :)))
    !
		! print*, 'ry', maxval(abs(dis_spt%rhoy(1, 0, :))), maxloc(abs(dis_spt%rhoy(1, 0, :)))
		! print*, 'uy', maxval(abs(dis_spt%uy  (1, 0, :))), maxloc(abs(dis_spt%uy  (1, 0, :)))
		! print*, 'vy', maxval(abs(dis_spt%vy  (1, 0, :))), maxloc(abs(dis_spt%vy  (1, 0, :)))
		! print*, 'Ty', maxval(abs(dis_spt%ty  (1, 0, :))), maxloc(abs(dis_spt%ty  (1, 0, :)))
    !
		! print*, 'rz', maxval(abs(dis_spt%rhoz(1, 0, :))), maxloc(abs(dis_spt%rhoz(1, 0, :)))
		! print*, 'uz', maxval(abs(dis_spt%uz  (1, 0, :))), maxloc(abs(dis_spt%uz  (1, 0, :)))
		! print*, 'vz', maxval(abs(dis_spt%vz  (1, 0, :))), maxloc(abs(dis_spt%vz  (1, 0, :)))
		! print*, 'Tz', maxval(abs(dis_spt%tz  (1, 0, :))), maxloc(abs(dis_spt%tz  (1, 0, :)))
    !
		! print*, 'rt', maxval(abs(dis_spt%rhot(1, 0, :))), maxloc(abs(dis_spt%rhot(1, 0, :)))
		! print*, 'ut', maxval(abs(dis_spt%ut  (1, 0, :))), maxloc(abs(dis_spt%ut  (1, 0, :)))
		! print*, 'vt', maxval(abs(dis_spt%vt  (1, 0, :))), maxloc(abs(dis_spt%vt  (1, 0, :)))
		! print*, 'Tt', maxval(abs(dis_spt%tt  (1, 0, :))), maxloc(abs(dis_spt%tt  (1, 0, :)))
    !
		! print*, 'uxx', maxval(abs(dis_spt%uxx(1, 0, :))), maxloc(abs(dis_spt%uxx(1, 0, :)))
		! print*, 'vxx', maxval(abs(dis_spt%vxx(1, 0, :))), maxloc(abs(dis_spt%vxx(1, 0, :)))
		! print*, 'Txx', maxval(abs(dis_spt%txx(1, 0, :))), maxloc(abs(dis_spt%txx(1, 0, :)))
    !
		! print*, 'uyy', maxval(abs(dis_spt%uyy(1, 0, :))), maxloc(abs(dis_spt%uyy(1, 0, :)))
		! print*, 'vyy', maxval(abs(dis_spt%vyy(1, 0, :))), maxloc(abs(dis_spt%vyy(1, 0, :)))
		! print*, 'Tyy', maxval(abs(dis_spt%tyy(1, 0, :))), maxloc(abs(dis_spt%tyy(1, 0, :)))
    !
		! print*, 'uzz', maxval(abs(dis_spt%uzz(1, 0, :))), maxloc(abs(dis_spt%uzz(1, 0, :)))
		! print*, 'vzz', maxval(abs(dis_spt%vzz(1, 0, :))), maxloc(abs(dis_spt%vzz(1, 0, :)))
		! print*, 'Tzz', maxval(abs(dis_spt%tzz(1, 0, :))), maxloc(abs(dis_spt%tzz(1, 0, :)))
    !
		! print*, 'uxy', maxval(abs(dis_spt%uxy(1, 0, :))), maxloc(abs(dis_spt%uxy(1, 0, :)))
		! print*, 'uxz', maxval(abs(dis_spt%uxz(1, 0, :))), maxloc(abs(dis_spt%uxz(1, 0, :)))
		! print*, 'vxy', maxval(abs(dis_spt%vxy(1, 0, :))), maxloc(abs(dis_spt%vxy(1, 0, :)))
		! print*, 'vyz', maxval(abs(dis_spt%vyz(1, 0, :))), maxloc(abs(dis_spt%vyz(1, 0, :)))
		! print*, 'wxz', maxval(abs(dis_spt%wxz(1, 0, :))), maxloc(abs(dis_spt%wxz(1, 0, :)))
		! print*, 'wyz', maxval(abs(dis_spt%wyz(1, 0, :))), maxloc(abs(dis_spt%wyz(1, 0, :)))
    ! print*, '121 location'
    !
    ! print*, 'r', (dis_spt%rho(0, 0, 200)), abs(dis_spt%rho(0, 0, 200))
		! print*, 'u', (dis_spt%u  (0, 0, 200)), abs(dis_spt%u  (0, 0, 200))
		! print*, 'v', (dis_spt%v  (0, 0, 200)), abs(dis_spt%v  (0, 0, 200))
		! print*, 'T', (dis_spt%t  (0, 0, 200)), abs(dis_spt%t  (0, 0, 200))
    !
		! print*, 'rx', (dis_spt%rhox(0, 0, 200)), (abs(dis_spt%rhox(0, 0, 200)))
		! print*, 'ux', (dis_spt%ux  (0, 0, 200)), (abs(dis_spt%ux  (0, 0, 200)))
		! print*, 'vx', (dis_spt%vx  (0, 0, 200)), (abs(dis_spt%vx  (0, 0, 200)))
		! print*, 'Tx', (dis_spt%tx  (0, 0, 200)), (abs(dis_spt%tx  (0, 0, 200)))
    !
		! print*, 'ry', ((dis_spt%rhoy(0, 0, 200))), (abs(dis_spt%rhoy(0, 0, 200)))
		! print*, 'uy', ((dis_spt%uy  (0, 0, 200))), (abs(dis_spt%uy  (0, 0, 200)))
		! print*, 'vy', ((dis_spt%vy  (0, 0, 200))), (abs(dis_spt%vy  (0, 0, 200)))
		! print*, 'Ty', ((dis_spt%ty  (0, 0, 200))), (abs(dis_spt%ty  (0, 0, 200)))
    !
		! print*, 'rz', ((dis_spt%rhoz(0, 0, 200))), (abs(dis_spt%rhoz(0, 0, 200)))
		! print*, 'uz', ((dis_spt%uz  (0, 0, 200))), (abs(dis_spt%uz  (0, 0, 200)))
		! print*, 'vz', ((dis_spt%vz  (0, 0, 200))), (abs(dis_spt%vz  (0, 0, 200)))
		! print*, 'Tz', ((dis_spt%tz  (0, 0, 200))), (abs(dis_spt%tz  (0, 0, 200)))
    !
		! print*, 'rt', ((dis_spt%rhot(0, 0, 200))), (abs(dis_spt%rhot(0, 0, 200)))
		! print*, 'ut', ((dis_spt%ut  (0, 0, 200))), (abs(dis_spt%ut  (0, 0, 200)))
		! print*, 'vt', ((dis_spt%vt  (0, 0, 200))), (abs(dis_spt%vt  (0, 0, 200)))
		! print*, 'Tt', ((dis_spt%tt  (0, 0, 200))), (abs(dis_spt%tt  (0, 0, 200)))
    !
		! print*, 'uxx', ((dis_spt%uxx(0, 0, 200))), (abs(dis_spt%uxx(0, 0, 200)))
		! print*, 'vxx', ((dis_spt%vxx(0, 0, 200))), (abs(dis_spt%vxx(0, 0, 200)))
		! print*, 'Txx', ((dis_spt%txx(0, 0, 200))), (abs(dis_spt%txx(0, 0, 200)))
    !
		! print*, 'uyy', ((dis_spt%uyy(0, 0, 200))), (abs(dis_spt%uyy(0, 0, 200)))
		! print*, 'vyy', ((dis_spt%vyy(0, 0, 200))), (abs(dis_spt%vyy(0, 0, 200)))
		! print*, 'Tyy', ((dis_spt%tyy(0, 0, 200))), (abs(dis_spt%tyy(0, 0, 200)))
    !
		! print*, 'uzz', ((dis_spt%uzz(0, 0, 200))), (abs(dis_spt%uzz(0, 0, 200)))
		! print*, 'vzz', ((dis_spt%vzz(0, 0, 200))), (abs(dis_spt%vzz(0, 0, 200)))
		! print*, 'Tzz', ((dis_spt%tzz(0, 0, 200))), (abs(dis_spt%tzz(0, 0, 200)))
    !
		! print*, 'uxy', ((dis_spt%uxy(1, 0, 200))), (abs(dis_spt%uxy(0, 0, 200)))
		! print*, 'uxz', ((dis_spt%uxz(0, 0, 200))), (abs(dis_spt%uxz(0, 0, 200)))
		! print*, 'vxy', ((dis_spt%vxy(0, 0, 200))), (abs(dis_spt%vxy(0, 0, 200)))
		! print*, 'vyz', ((dis_spt%vyz(0, 0, 200))), (abs(dis_spt%vyz(0, 0, 200)))
		! print*, 'wxz', ((dis_spt%wxz(0, 0, 200))), (abs(dis_spt%wxz(0, 0, 200)))
		! print*, 'wyz', ((dis_spt%wyz(0, 0, 200))), (abs(dis_spt%wyz(0, 0, 200)))
    !
    !
    !
    !   endblock

      end associate

      ! block

        ! use mod_debug, only: jloc

      do j=1, this%jn

        ! jloc=j

        call this%spectral_to_physical(this%dis_spt, j, this%dis_phy)
        Coord=this%NormCoord%GetPoint(j)
        BFPoint=this%BFNorm%GetPoint(j)

        call this%ComputeF(Coord, BFPoint, j)
        call this%physical_to_spectral(this%dis_phy%FTerm, this%nonlinearTerm(:, j, :, :))

      enddo
    ! endblock
      !
      ! write(*, *) '............'
      ! write(*, *) 'Check nonlinearTerm1'
      ! write(*, *)this%nonlinearTerm(1, 200, 0, 0), abs(this%nonlinearTerm(1, 200, 0, 0))
      ! write(*, *)this%nonlinearTerm(2, 200, 0, 0), abs(this%nonlinearTerm(2, 200, 0, 0))
      ! write(*, *)this%nonlinearTerm(3, 200, 0, 0), abs(this%nonlinearTerm(3, 200, 0, 0))
      ! write(*, *)this%nonlinearTerm(4, 200, 0, 0), abs(this%nonlinearTerm(4, 200, 0, 0))
      ! write(*, *)this%nonlinearTerm(5, 200, 0, 0), abs(this%nonlinearTerm(5, 200, 0, 0))
      ! pause
      ! block
      !   integer :: m, n, l, j
      !   print*, this%mdim, this%ndim
      !   do m=0, this%mdim, 1
      !     do n=0, 0
      !       do l=1, 5
      !       !-this%ndim, this%ndim
      !       write(*, *)m, l, maxval(abs(this%nonlinearTerm(l, :, m, n))),maxloc(abs(this%nonlinearTerm(l, :, m, n)))
      !       enddo
      !     enddo
      !   enddo

      ! do j=1, this%jn
      !
      !   write(*, *)j, abs(this%nonlinearTerm(5, j, 0, 0))
      !
      ! enddo


      ! pause
      ! endblock

      do j=1, this%jn
        do l=1, 5
          this%nonlinearTerm(l, j, :, :)=this%nonlinearTerm(l, j, :, :)/Amp(:, :)
        enddo
      enddo
      !endblock
#ifdef DEBUG
      write(*, *) '............'
      write(*, *) 'Check nonlinearTerm2'
      write(*, *)this%nonlinearTerm(1, 200, 0, 0), abs(this%nonlinearTerm(1, 200, 0, 0))
      write(*, *)this%nonlinearTerm(2, 200, 0, 0), abs(this%nonlinearTerm(2, 200, 0, 0))
      write(*, *)this%nonlinearTerm(3, 200, 0, 0), abs(this%nonlinearTerm(3, 200, 0, 0))
      write(*, *)this%nonlinearTerm(4, 200, 0, 0), abs(this%nonlinearTerm(4, 200, 0, 0))
      write(*, *)this%nonlinearTerm(5, 200, 0, 0), abs(this%nonlinearTerm(5, 200, 0, 0))
      ! pause
#endif
    end subroutine solve

    !>非线性项F计算
    subroutine ComputeF1(this, Coord, BF, j)

      use mod_baseflow_org, only: BF_flux_org_ij_type
      use mod_local_coordinate
      use mod_lame
      use mod_basis
      use mod_gas, only: GAMMA, Pr, miu0=>miu, miu0T=>miuT, miu0TT=>miuTT
      use mod_parameter

      implicit none
      class(nonlinear_type), intent(inout) :: This
      type(local_coordinate_type), intent(in) :: Coord
      type(bf_point_type), intent(in) :: BF
      real(R_P), dimension(0:(this%mdim+1)*8-1, 0:(this%ndim+1)*8-1) :: &
        & divu, divux, divuy,divuz, &
        & s11, s22, s33, s12, s13, s23, &
        & s11x, s12y, s13z, &
        & s12x, s22y, s23z, &
        & s13x, s23y, s33z, &
        & miu, miux, miuy, miuz
      real(R_P) :: s011, s022, s033, s012, s013, s023
      real(R_P) :: rho0, u0, v0, w0, T0
      real(R_P) :: rho0x, u0x, v0x, w0x, T0x
      real(R_P) :: rho0y, u0y, v0y, w0y, T0y
      real(R_P) :: Pe0, c0, Vig
      real(R_P) :: GF, GF1, GF2M
!      real(R_P) :: u0, v0, w0

      real(R_P) :: hx, hy, hz
      real(R_P) :: d12, d32, d31
      type(lame_type) :: lame
      type(lame_grad_type) :: LameGrad
      type(basis_type) :: basis
      type(bf_flux_org_ij_type) :: flux, dxflux, dyflux, ddyflux

!# debug
      integer :: j

      call Coord%Get(basis, lame)
      call lame%Get(hx, hy, hz)
      call coord%GetGrad(LameGrad)
      call LameGrad%Get(d12, d32, d31)

      call BF%Get(flux, dxflux, dyflux, ddyflux)
      call Flux%get(rho0, u0, v0, w0, T0)
      call DxFlux%get(rho0x, u0x, v0x, w0x, T0x)
      call DyFlux%get(rho0y, u0y, v0y, w0y, T0y)

      ! call DDyFlux%get(rhoyy, uyy, vyy, wyy, Tyy)

      s011=u0x/hx+v0*d12
      s022=v0y
      s033=u0*d31+v0*d32
      s012=0.5d0*(v0x/hx-d12*u0+u0y)
      s013=0.5d0*(w0x/hx-d31*w0)
      s023=0.5d0*(w0y-w0*d32)

      Pe0=1.0d0/(GAMMA*Ma**2)
      GF=-1.0d0/GAMMA
      GF1=(GAMMA-1.0d0)/GAMMA
      GF2M=(GAMMA-1.d0)*Ma*Ma

      c0=U0/sqrt(T0)*Ma

      Vig=min(1.0d0, c0 **2)

      associate(dis_phy=>this%dis_phy)
        associate(rhox=>dis_phy%rhox, &
          &       ux=>dis_phy%ux, &
          &       vx=>dis_phy%vx, &
          &       wx=>dis_phy%wx, &
          &       Tx=>dis_phy%Tx, &
          &       rhoy=>dis_phy%rhoy, &
          &       uy=>dis_phy%uy, &
          &       vy=>dis_phy%vy, &
          &       wy=>dis_phy%wy, &
          &       Ty=>dis_phy%Ty, &
          &       rhoz=>dis_phy%rhoz, &
          &       uz=>dis_phy%uz, &
          &       vz=>dis_phy%vz, &
          &       wz=>dis_phy%wz, &
          &       Tz=>dis_phy%Tz, &
          &       rhot=>dis_phy%rhot, &
          &       ut=>dis_phy%ut, &
          &       vt=>dis_phy%vt, &
          &       wt=>dis_phy%wt, &
          &       Tt=>dis_phy%Tt, &
          &       rho=>dis_phy%rho, &
          &       u =>dis_phy%u , &
          &       v =>dis_phy%v , &
          &       w =>dis_phy%w , &
          &       T =>dis_phy%T , &
          &      uxx=>dis_phy%uxx, &
          &      uyy=>dis_phy%uyy, &
          &      uzz=>dis_phy%uzz, &
          &      vxx=>dis_phy%vxx, &
          &      vyy=>dis_phy%vyy, &
          &      vzz=>dis_phy%vzz, &
          &      wxx=>dis_phy%wxx, &
          &      wyy=>dis_phy%wyy, &
          &      wzz=>dis_phy%wzz, &
          &      Txx=>dis_phy%Txx, &
          &      Tyy=>dis_phy%Tyy, &
          &      Tzz=>dis_phy%Tzz, &
          &      uxy=>dis_phy%uxy, &
          &      uxz=>dis_phy%uxz, &
          &      vxy=>dis_phy%vxy, &
          &      vyz=>dis_phy%vyz, &
          &      wxz=>dis_phy%wxz, &
          &      wyz=>dis_phy%wyz, &
          &      FTerm=>dis_phy%FTerm)

          divu=ux/hx+vy+wz/hz+u*d31+v*(d12+d32)

        	s11=ux/hx+v*d12
        	s22=vy
        	s33=wz/hz+u*d31+v*d32
        	s12=0.5d0*(vx/hx-u*d12+uy)
        	s13=0.5d0*(wx/hx+uz/hz-w*d31)
        	s23=0.5d0*(wy+vz/hz-w*d32)

          divux=uxx/hx/hx+vxy/hx+wxz/hx/hz-d31*wz/hz+ &
        	&     d31*ux/hx+(d12+d32)*vx/hx
        	divuy=uxy/hx-d12*ux/hx+vyy+wyz/hz-d32*wz/hz+ &
          &     d31*uy+(d12+d32)*vy
        	divuz=uxz/hx/hz+vyz/hz+wzz/hz/hz+d31*uz/hz+ &
        	&     (d12+d32)*vz/hz

        	s11x=uxx/hx/hx+d12*vx/hx
        	s12y=0.5d0*(vxy/hx-d12*vx/hx+uyy-d12*uy)
        	s13z=0.5d0*(wxz/hx/hz+uzz/hz/hz-d31*wz/hz)

        	s12x=0.5d0*(vxx/hx/hx+uxy/hx-d12*ux/hx)
        	s22y=vyy
        	s23z=0.5d0*(wyz/hz+vzz/hz/hz-d32*wz/hz)

        	s13x=0.5d0*(wxx/hx/hx+uxz/hx/hz-d31*uz/hz-d31*wx/hx)
        	s23y=0.5d0*(wyy+vyz/hz-d32*vz/hz-d32*wy)
        	s33z=wzz/hz/hz+d31*uz/hz+d32*vz/hz

          miu=miu0T(T0)*T/Re
        	miux=(miu0TT(T0)*T0x/hx*T+miu0T(T0)*Tx/hx)/Re
        	miuy=(miu0TT(T0)*T0y*T+miu0T(T0)*Ty)/Re
        	miuz=(miu0T(T0)*Tz/hz)/Re

          FTerm(1, :, :) = -(rho*divu+u*rhox/hx+v*rhoy+w*rhoz/hz)

          FTerm(2, :, :) = -rho*(ut+ &
             &                   u0 *ux /hx+v0*uy+w0*uz/hz+u0*v*d12-w0*w*d31+ &
       	     &                   u  *u0x/hx+v*u0y         +d12*v0*u-d31*w0*w+ &
             &                   u  *ux /hx+v*uy +w *uz/hz+u *v*d12-w*w*d31)- &
       	     &              rho0*(u*ux/hx+v*uy+w*uz/hz+u*v*d12-w*w*d31)- &
       	     &              Pe0*(T*rhox/hx+rho*Tx/hx)*Vig+ &
      	     &  2.d0*(miux*s11+miu*s11x+miuy*s12+miu*s12y+miuz*s13+miu*s13z+  &
      	     &        miu*(2.d0*d12+d32)*s12+miu*s11*d31-miu*s33*d31)- &
      	     &  2.d0/3.d0*(miux*divu+miu*divux)

   ! FTerm(2, :, :) = (ut+ &
   !    &                   u0 *ux /hx+v0*uy+w0*uz/hz+u0*v*d12-w0*w*d31+ &
	 !     &                   u  *u0x/hx+v*u0y         +d12*v0*u-d31*w0*w)+ &
	 !     &              Pe0*(Tx/hx)*Vig
! block
!   use mod_debug
! if(jloc==121)then
!   print*, 'check baseflow'
!   print*, 'u0', u0
!   print*, 'v0', v0
!   print*, 'w0', w0
!   print*, 'u0x', u0x
!   print*, 'u0y', u0y
!   print*, 'Vig', vig
!   pause
!
! endif
! endblock
       	 FTerm(3, :, :) = -rho*(vt+ &
             &                  u0*vx /hx+v0*vy +w0*vz/hz-u0*u*d12-w0*w*d32+ &
       	     &                  u *v0x/hx+v *v0y         -u0*u*d12-w0*w*d32+ &
       	     &                  u *vx /hx+v *vy +w *vz/hz-u *u*d12-w *w*d32)- &
       	     &             rho0*(u*vx/hx+v*vy+w*vz/hz-u*u*d12-w*w*d32)-&
       	     &             Pe0 *(rho*Ty+T*rhoy)+        &
       	     &  2.d0*(miux*s12+miu*s12x+miuy*s22+miu*s22y+miuz*s23+miu*s23z+ &
       	     &        miu*s22*d12+miu*s12*d31+miu*s22*d32- &
             &        miu*s11*d12-miu*s33*d32)- &
       	     &  2.d0/3.d0*(miuy*divu+miu*divuy)

       	 FTerm(4, :, :)=-rho*(wt+ &
             &                u0*wx/hx+v0*wy+w0*wz/hz+w0*u*d31+w0*v*d32+ &
       	     &                u*w0x/hx+v*w0y         +d31*u0*w+d32*v0*w+ &
       	     &                u*wx /hx+v*wy +w *wz/hz+w *u*d31+w *v*d32)-&
       	     &          rho0*(u*wx/hx+v*wy+w*wz/hz+w*u*d31+w*v*d32)- &
       	     &          pe0*(rho*Tz/hz+T*rhoz/hz)+ &
       	     &  2.d0*(miux*s13+miu*s13x+miuy*s23+miu*s23y+miuz*s33+miu*s33z+ &
       	     &        miu*s23*d12+miu*s13*d31+miu*s23*d32+ &
             &        miu*s13*d31+miu*s23*d32)- &
       	     &  2.d0/3.d0*(miuz*divu+miu*divuz)


       	 FTerm(5, :, :)=rho*GF*(Tt+ &
             &                       (u *T0x/hx+v *T0y)+ &
             &                       (u0*Tx /hx+v0*Ty +w0*Tz/hz)+  &
       	     &                       (u *Tx /hx+v *Ty +w *Tz/hz))+ &
             &           rho0*GF*(u*Tx/hx+v*Ty+w*Tz/hz)+ &
       	     &          GF1*T*(rhot+ &
             &                      (u *rho0x/hx+v *rho0y) +         &
             &                      (u0*rhox /hx+v0*rhoy+w0*rhoz/hz)+&
       	     &                      (u *rhox /hx+v *rhoy+w *rhoz/hz))+ &
             &          GF1*T0*(     u *rhox /hx+v *rhoy+w *rhoz/hz)+ &
       	     &          ( miux*Tx/hx+miuy*Ty+miuz*Tz/hz+ &
             &            miu *(Txx/hx/hx+Tyy+Tzz/hz/hz)+ &
       	     &            miu *(d12*Ty+d32*Ty+d31*Tx/hx))/Pr+ &
       	     &GF2M*(4.d0*miu*(s011*s11+s022*s22+s033*s33+ &
             &                2.d0*(s012*s12+s013*s13+s023*s23)) &
             &      -4.d0/3.d0*miu*(u0x/hx+v0y+u0*d31+v0*(d12+d32))*divu)+&
       	     &     GF2M*(miu+miu0(T0)/Re)*(2.d0*(s11*s11+s22*s22+s33*s33+ &
             &                               2.d0*(s12*s12+s13*s13+s23*s23))- &
             &                         2.d0/3.d0*divu*divu)


             ! print*, (this%mdim+1)*4-1, (this%ndim+1)*4-1
             ! print*, j
             ! print*, 0, 15
             ! print*, Fterm(1, 0, 15)
             ! print*, Fterm(2, 0, 15)
             ! print*, Fterm(3, 0, 15)
             ! print*, Fterm(4, 0, 15)
             ! print*, Fterm(5, 0, 15)
             ! print*, 0, 16
             ! print*, Fterm(1, 0, 16)
             ! print*, Fterm(2, 0, 16)
             ! print*, Fterm(3, 0, 16)
             ! print*, Fterm(4, 0, 16)
             ! print*, Fterm(5, 0, 16)
             !
             ! pause



             !   block
             !     integer :: m, n, l
             !
             !     print*, size(Fterm(:, :, :), dim=2), size(Fterm(:, :, :), dim=3)
             !  do l=1, 5
             !   do m=0, size(Fterm(:, :, :), dim=2)-1
             !     do n=0, size(Fterm(:, :, :), dim=3)-1
             !       if(isnan(FTerm(5, m, n))) then
             !         print*, m, n, l, Fterm(l, m, n)
             !         ! print*, m, n, u(m, n)
             !         ! print*, m, n, v(m, n)
             !         ! print*, m, n, w(m, n)
             !         ! print*, m, n, T(m, n)
             !         ! print*, m, n, rhox(m, n)
             !         ! print*, m, n, ux(m, n)
             !         ! print*, m, n, vx(m, n)
             !         ! print*, m, n, wx(m, n)
             !         ! print*, m, n, Tx(m, n)
             !         ! print*, m, n, rhoy(m, n)
             !         ! print*, m, n, uy(m, n)
             !         ! print*, m, n, vy(m, n)
             !         ! print*, m, n, wy(m, n)
             !         ! print*, m, n, Ty(m, n)
             !         ! print*, m, n, rhoz(m, n)
             !         ! print*, m, n, uz(m, n)
             !         ! print*, m, n, vz(m, n)
             !         ! print*, m, n, wz(m, n)
             !         ! print*, m, n, Tz(m, n)
             !         ! print*, m, n, rhot(m, n)
             !         ! print*, m, n, ut(m, n)
             !         ! print*, m, n, vt(m, n)
             !         ! print*, m, n, wt(m, n)
             !         ! print*, m, n, Tt(m, n)
             !         pause
             !       endif
             !      enddo
             !   enddo
             ! enddo
             !
             ! end block

        end associate
      end associate

    end subroutine ComputeF1


    !>非线性项F计算
    subroutine ComputeF2(this, Coord, BF, j)

      use mod_baseflow_org, only: BF_flux_org_ij_type
      use mod_local_coordinate
      use mod_lame
      use mod_basis
      use mod_gas, only: GAMMA, Pr, miu0=>miu, miu0T=>miuT, miu0TT=>miuTT
      use mod_parameter

      implicit none
      class(nonlinear_type), intent(inout) :: This
      type(local_coordinate_type), intent(in) :: Coord
      type(bf_point_type), intent(in) :: BF
      real(R_P), dimension(0:(this%mdim+1)*8-1, 0:(this%ndim+1)*8-1) :: &
        & bMua,bMuTa,bMuTTa,&
        & bMut,bMuTt,bMuTTt,&
        & bMux,bMuTx,bMuTTx,&
        & bMuy,bMuTy,bMuTTy,&
    	& bra,bua,bva,bwa,bta,&
    	& brt,but,bvt,bwt,btt,brz,buz,bvz,bwz,btz,&
    	& brx,bux,bvx,bwx,btx,bry,buy,bvy,bwy,bty,&
        & buxx,bvxx,bwxx,btxx,buyy,bvyy,bwyy,btyy,&
    	& buzz,bvzz,bwzz,btzz,buxy,bvxy,buxz,bwxz,bvyz,bwyz,&
        & bfL1,bfL2,bfL3,bfL4,bfL5,afL1,afL2,afL3,afL4,afL5,&
        & miup,miup1,miup2,miup3,&
        & divup,s11p,s12p,s13p,s22p,s23p,s33p,&
        & divup1,divup2,divup3,s11p1,s12p2,s13p3,&
        & s12p1,s22p2,s23p3,s13p1,s23p2,s33p3,&
        & s11,s12,s13,s22,s23,s33
      real(R_P) :: s011, s022, s033, s012, s013, s023
      real(R_P) :: rho0, u0, v0, w0, T0
      real(R_P) :: rho0x, u0x, v0x, w0x, T0x
      real(R_P) :: rho0y, u0y, v0y, w0y, T0y
      real(R_P) :: Pe0, c0, Vig
      real(R_P) :: GF, GF1, GF2M
!      real(R_P) :: u0, v0, w0

      real(R_P) :: hx, hy, hz
      real(R_P) :: d12, d32, d31
      type(lame_type) :: lame
      type(lame_grad_type) :: LameGrad
      type(basis_type) :: basis
      type(bf_flux_org_ij_type) :: flux, dxflux, dyflux, ddyflux

      real(R_P) :: h10, h30
      real(R_P) :: d310, d320, d120
      real(R_P) :: Ta0, Ta0x, Ta0y
      real(R_P) :: r0, r0x, r0y
      real(R_P) :: pe, gmf, gmf1, gm1M2, paraR
      real(R_P) :: omv
      real(R_P) :: umu0

!# debug
      integer :: j

      call Coord%Get(basis, lame)
      call lame%Get(hx, hy, hz)
      call coord%GetGrad(LameGrad)
      call LameGrad%Get(d12, d32, d31)

      call BF%Get(flux, dxflux, dyflux, ddyflux)
      call Flux%get(rho0, u0, v0, w0, T0)
      call DxFlux%get(rho0x, u0x, v0x, w0x, T0x)
      call DyFlux%get(rho0y, u0y, v0y, w0y, T0y)

      c0=U0/sqrt(T0)*Ma

      Vig=min(1.0d0, c0 **2)
      omv=Vig


      gmf=-1.d0/Gamma
      gmf1=gmf+1.d0
      gm1M2=(Gamma-1.d0)*Ma*Ma/Re
      pe=1.d0/(Gamma*Ma*Ma)
      ParaR= gm1M2*0.4d0/3.0d0

      h10=hx
      h30=hz
      d310=d31
      d320=d32
      d120=d12

      ta0=miu0T(T0)
      ta0x=miu0TT(T0)*T0x
      ta0y=miu0TT(T0)*T0y
      umu0=miu0(T0)
      r0=rho0
      r0x=rho0x
      r0y=rho0y

      bra=this%dis_phy%rho
      bua=this%dis_phy%u
      bva=this%dis_phy%v
      bwa=this%dis_phy%w
      bta=this%dis_phy%T

      brt=this%dis_phy%rhot
      but=this%dis_phy%ut
      bvt=this%dis_phy%vt
      bwt=this%dis_phy%wt
      btt=this%dis_phy%tt

      brz=this%dis_phy%rhoz
      buz=this%dis_phy%uz
      bvz=this%dis_phy%vz
      bwz=this%dis_phy%wz
      btz=this%dis_phy%tz

      brx=this%dis_phy%rhox
      bux=this%dis_phy%ux
      bvx=this%dis_phy%vx
      bwx=this%dis_phy%wx
      btx=this%dis_phy%tx

      bry=this%dis_phy%rhoy
      buy=this%dis_phy%uy
      bvy=this%dis_phy%vy
      bwy=this%dis_phy%wy
      bty=this%dis_phy%ty

      buxx=this%dis_phy%uxx
      bvxx=this%dis_phy%Vxx
      bwxx=this%dis_phy%wxx
      btxx=this%dis_phy%txx

      buyy=this%dis_phy%uyy
      bvyy=this%dis_phy%Vyy
      bwyy=this%dis_phy%wyy
      btyy=this%dis_phy%tyy

      buzz=this%dis_phy%uzz
      bvzz=this%dis_phy%Vzz
      bwzz=this%dis_phy%wzz
      btzz=this%dis_phy%tzz

      buxy=this%dis_phy%uxy
      bvxy=this%dis_phy%vxy
      buxz=this%dis_phy%uxz
      bwxz=this%dis_phy%wxz
      bvyz=this%dis_phy%vyz
      bwyz=this%dis_phy%wyz

      divup=bux/H10 +bvy+bwz/H30 +bua*d310 +bva*(d120 +d320 )
    !	print*,divup
    !	stop
    	s11p=bux/H10 +bva*d120
    	s22p=bvy
    	s33p=bwz/H30 +bua*d310 +bva*d320
    	s12p=0.5d0*(bvx/H10 -bua*d120 +buy)
    	s13p=0.5d0*(bwx/H10 +buz/H30 -bwa*d310 )
    	s23p=0.5d0*(bwy+bvz/H30 -bwa*d320 )

    	divup1=buxx/H10 /H10 +bvxy/H10 +bwxz/H10 /H30 -d310 *bwz/H30 +&
    	&d310 *bux/H10 +(d120 +d320 )*bvx/H10
    	divup2=buxy/H10 -d120 *bux/H10 +bvyy+bwyz/H30 -d320 *bwz/H30 +&
    	&d310 *buy+(d120 +d320 )*bvy
    	divup3=buxz/H10 /H30 +bvyz/H30 +bwzz/H30 /H30 +d310 *buz/H30 +&
    	&(d120 +d320 )*bvz/H30

    	s11p1=buxx/H10 /H10 +d120 *bvx/H10
    	s12p2=0.5d0*(bvxy/H10 -d120 *bvx/H10 +buyy-d120 *buy)
    	s13p3=0.5d0*(bwxz/H10 /H30 +buzz/H30 /H30 -d310 *bwz/H30 )

    	s12p1=0.5d0*(bvxx/H10 /H10 +buxy/H10 -d120 *bux/H10 )
    	s22p2=bvyy
    	s23p3=0.5d0*(bwyz/H30 +bvzz/H30 /H30 -d320 *bwz/H30 )

    	s13p1=0.5d0*(bwxx/H10 /H10 +buxz/H10 /H30 -d310 *buz/H30 -d310 *bwx/H10 )
    	s23p2=0.5d0*(bwyy+bvyz/H30 -d320 *bvz/H30 -d320 *bwy)
    	s33p3=bwzz/H30 /H30 +d310 *buz/H30 +d320 *bvz/H30

    	miup=ta0 *bta
    	miup1=ta0x /H10 *bta+ta0 *btx/H10
    	miup2=ta0y *bta+ta0 *bty
    	miup3=ta0 *btz/H30

    	s11=u0x /H10 +v0 *d120
    	s22=v0y
    	s33=u0 *d310 +v0 *d320
    	s12=0.5d0*(v0x /H10 -d120 *u0 +u0y )
    	s13=0.5d0*(w0x /H10 -d310 *w0 )
    	s23=0.5d0*(w0y -w0 *d320 )


            !-----af0(1)----------------------

            bfL1=-(bra*bux/H10 +bua*brx/H10 +bra*bvy+bva*bry+bra*bwz/H30 +bwa*brz/H30 +&
    	     &bra*bva*(d120 +d320 )+bra*bua*d310 )

           this%dis_phy%FTerm(1, :, :)=bfL1





    	bfL2=-bra*(but+u0 *bux/H10 +v0 *buy+w0 *buz/H30 +u0 *bva*d120 -w0 *bwa*d310 +&
    	     &u0x /H10 *bua+u0y *bva+d120 *v0 *bua-d310 *w0 *bwa+&
                 &bua*bux/H10 +bva*buy+bwa*buz/H30 +bua*bva*d120 -bwa*bwa*d310 )-&
    	     &r0 *(bua*bux/H10 +bva*buy+bwa*buz/H30 +bua*bva*d120 -bwa*bwa*d310 )-&
    	     &pe*(bta*brx/H10 +bra*btx/H10 )*omv +&
    	     &2.d0*(miup1*s11p+miup*s11p1+miup2*s12p+miup*s12p2+miup3*s13p+miup*s13p3+&
    	     &miup*(2.d0*d120 +d320 )*s12p+miup*s11p*d310 -miup*s33p*d310 )/re-&
    	     &2.d0/3.d0*(miup1*divup+miup*divup1)/re

           this%dis_phy%Fterm(2, :, :)=bfL2

    	bfL3=-bra*(bvt+u0 *bvx/H10 +v0 *bvy+w0 *bvz/H30 -u0 *bua*d120 -w0 *bwa*d320 +&
    	     &bua*v0x /H10 +bva*v0y -d120 *u0 *bua-d320 *w0 *bwa+&
    	     &bua*bvx/H10 +bva*bvy+bwa*bvz/H30 -bua*bua*d120 -bwa*bwa*d320 )-&
    	     &r0 *(bua*bvx/H10 +bva*bvy+bwa*bvz/H30 -bua*bua*d120 -bwa*bwa*d320 )-&
    	     &pe*(bra*bty+bta*bry)+&
    	     &2.d0*(miup1*s12p+miup*s12p1+miup2*s22p+miup*s22p2+miup3*s23p+miup*s23p3+&
    	     &miup*s22p*d120 +miup*s12p*d310 +miup*s22p*d320 -miup*s11p*d120 -miup*s33p*d320 )/re-&
    	     &2.d0/3.d0*(miup2*divup+miup*divup2)/re
           this%dis_phy%Fterm(3, :, :)=bfL3

    	bfL4=-bra*(bwt+u0 *bwx/H10 +v0 *bwy+w0 *bwz/H30 +w0 *bua*d310 +w0 *bva*d320 +&
    	     &bua*w0x /H10 +bva*w0y +d310 *u0 *bwa+d320 *v0 *bwa+&
    	     &bua*bwx/H10 +bva*bwy+bwa*bwz/H30 +bwa*bua*d310 +bwa*bva*d320 )-&
    	     &r0 *(bua*bwx/H10 +bva*bwy+bwa*bwz/H30 +bwa*bua*d310 +bwa*bva*d320 )-&
    	     &pe*(bra*btz/H30 +bta*brz/H30 )+&
    	     &2.d0*(miup1*s13p+miup*s13p1+miup2*s23p+miup*s23p2+miup3*s33p+miup*s33p3+&
    	     &miup*s23p*d120 +miup*s13p*d310 +miup*s23p*d320 +miup*s13p*d310 +miup*s23p*d320 )/re-&
    	     &2.d0/3.d0*(miup3*divup+miup*divup3)/re
           this%dis_phy%Fterm(4, :, :)=bfL4

    	bfL5=bra*gmf*(btt+(bua*t0x /H10 +bva*t0y )+(u0 *btx/H10 +v0 *bty+w0 *btz/H30 )+&
    	     &(bua*btx/H10 +bva*bty+bwa*btz/H30 ))+ &
           & r0 *gmf*(bua*btx/H10 +bva*bty+bwa*btz/H30 )+&
    	     &gmf1*bta*(brt+(bua*r0x /H10 +bva*r0y )+(u0 *brx/H10 +v0 *bry+w0 *brz/H30 )+&
    	     &(bua*brx/H10 +bva*bry+bwa*brz/H30 )) &
           & +gmf1*t0 *(bua*brx/H10 +bva*bry+bwa*brz/H30 ) + &
    	     &(miup1*btx/H10 +miup2*bty+miup3*btz/H30 +miup*(btxx/H10 /H10 +btyy+btzz/H30 /H30 )+&
    	     &miup*(d120 *bty+d320 *bty+d310 *btx/H10 ))/re/pr+&
    	     &gm1M2*(4.d0*miup*(s11*s11p+s22*s22p+s33*s33p+2.d0*(s12*s12p+s13*s13p+&
    	     &s23*s23p))-4.d0/3.d0*miup*(u0x /H10 +v0y +u0 *d310 +v0 *(d120 +d320 ))*divup)+&
    	     &gm1M2*(miup+umu0 )*(2.d0*(s11p*s11p+s22p*s22p+s33p*s33p+2.d0*&
    	     &(s12p*s12p+s13p*s13p+s23p*s23p))-2.d0/3.d0*divup*divup)
           this%dis_phy%Fterm(5, :, :)=bfL5


    end subroutine ComputeF2

    subroutine physical_to_spectral(this, phy, spec)

      implicit none
      class(nonlinear_type), intent(inout) ::this
      real(R_P), intent(inout) :: phy(5, 0:(this%mdim+1)*8-1, 0:(this%ndim+1)*8-1)
      complex(R_P), intent(inout) :: spec(5, 0:this%mdim, -this%ndim:this%ndim)
      integer :: l

      ! print*, phy(1, :, :)
      ! pause

      do l=1, 5

        block

          use mod_debug, only:lloc

          lloc=l

        endblock

        call this%SPT%ForTrans(phy(l, :, :),spec(l, :, :))
      enddo

    end subroutine physical_to_spectral

    subroutine spectral_to_physical(this, spec, jloc, phy)

      implicit none
      class(nonlinear_type), intent(inout) :: this
      type(dis_spt_type), intent(inout) :: spec
      integer, intent(in) :: jloc
      type(dis_phy_type), intent(inout) :: phy


        call this%SPT%BackTrans(spec%rho(:, :, jloc), phy%rho(:, :))

        call this%SPT%BackTrans(spec%u  (:, :, jloc), phy%u  (:, :))
        call this%SPT%BackTrans(spec%v  (:, :, jloc), phy%v  (:, :))
        call this%SPT%BackTrans(spec%w  (:, :, jloc), phy%w  (:, :))
        call this%SPT%BackTrans(spec%T  (:, :, jloc), phy%T  (:, :))

        call this%SPT%BackTrans(spec%rhox(:, :, jloc), phy%rhox(:, :))
        call this%SPT%BackTrans(spec%ux  (:, :, jloc), phy%ux  (:, :))
        call this%SPT%BackTrans(spec%vx  (:, :, jloc), phy%vx  (:, :))
        call this%SPT%BackTrans(spec%wx  (:, :, jloc), phy%wx  (:, :))
        call this%SPT%BackTrans(spec%Tx  (:, :, jloc), phy%Tx  (:, :))

        call this%SPT%BackTrans(spec%rhoy(:, :, jloc), phy%rhoy(:, :))
        call this%SPT%BackTrans(spec%uy  (:, :, jloc), phy%uy  (:, :))
        call this%SPT%BackTrans(spec%vy  (:, :, jloc), phy%vy  (:, :))
        call this%SPT%BackTrans(spec%wy  (:, :, jloc), phy%wy  (:, :))
        call this%SPT%BackTrans(spec%Ty  (:, :, jloc), phy%Ty  (:, :))

        call this%SPT%BackTrans(spec%rhoz(:, :, jloc), phy%rhoz(:, :))
        call this%SPT%BackTrans(spec%uz  (:, :, jloc), phy%uz  (:, :))
        call this%SPT%BackTrans(spec%vz  (:, :, jloc), phy%vz  (:, :))
        call this%SPT%BackTrans(spec%wz  (:, :, jloc), phy%wz  (:, :))
        call this%SPT%BackTrans(spec%Tz  (:, :, jloc), phy%Tz  (:, :))

        call this%SPT%BackTrans(spec%rhot(:, :, jloc), phy%rhot(:, :))
        call this%SPT%BackTrans(spec%ut  (:, :, jloc), phy%ut  (:, :))
        call this%SPT%BackTrans(spec%vt  (:, :, jloc), phy%vt  (:, :))
        call this%SPT%BackTrans(spec%wt  (:, :, jloc), phy%wt  (:, :))
        call this%SPT%BackTrans(spec%Tt  (:, :, jloc), phy%Tt  (:, :))

        call this%SPT%BackTrans(spec%uxx(:, :, jloc), phy%uxx(:, :))
        call this%SPT%BackTrans(spec%vxx(:, :, jloc), phy%vxx(:, :))
        call this%SPT%BackTrans(spec%wxx(:, :, jloc), phy%wxx(:, :))
        call this%SPT%BackTrans(spec%Txx(:, :, jloc), phy%Txx(:, :))

        call this%SPT%BackTrans(spec%uyy(:, :, jloc), phy%uyy(:, :))
        call this%SPT%BackTrans(spec%vyy(:, :, jloc), phy%vyy(:, :))
        call this%SPT%BackTrans(spec%wyy(:, :, jloc), phy%wyy(:, :))
        call this%SPT%BackTrans(spec%Tyy(:, :, jloc), phy%Tyy(:, :))

        call this%SPT%BackTrans(spec%uzz(:, :, jloc), phy%uzz(:, :))
        call this%SPT%BackTrans(spec%vzz(:, :, jloc), phy%vzz(:, :))
        call this%SPT%BackTrans(spec%wzz(:, :, jloc), phy%wzz(:, :))
        call this%SPT%BackTrans(spec%Tzz(:, :, jloc), phy%Tzz(:, :))

        call this%SPT%BackTrans(spec%uxy(:, :, jloc), phy%uxy(:, :))
        call this%SPT%BackTrans(spec%vxy(:, :, jloc), phy%vxy(:, :))
        call this%SPT%BackTrans(spec%uxz(:, :, jloc), phy%uxz(:, :))
        call this%SPT%BackTrans(spec%wxz(:, :, jloc), phy%wxz(:, :))
        call this%SPT%BackTrans(spec%vyz(:, :, jloc), phy%vyz(:, :))
        call this%SPT%BackTrans(spec%wyz(:, :, jloc), phy%wyz(:, :))

    end subroutine spectral_to_physical

    !> 判断是否判断
    logical function IsCreate(this)
      implicit none
      class(nonlinear_type), intent(in) :: this

	      IsCreate=associated(this%nonlinearterm)

    end function IsCreate

    !>析构函数
    subroutine Finalize(this)

      implicit none
      class(nonlinear_type), intent(inout) :: this

      this%jn=0
      this%mdim=0
      this%ndim=0
      deallocate(this%nonlinearTerm)
      call this%dis_spt%Finalize
      call this%dis_phy%finalize
      deallocate(this%IntergalAlf)
      deallocate(this%alf)
      call this%SPT%Finalize()

    end subroutine Finalize


end module mod_nonlinear
