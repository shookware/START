!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_APSE_LPSE.f90
!> @file
!> @breif 先求解APSE，再求解LPSE.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: mod_APSE_LPSE
!> @breif 先求解APSE后求解LPSE.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2018-06-26 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2018-06-26
!------------------------------------------------------------------------------

module mod_APSE_LPSE

  use mod_solver
  use mod_lpse_eqn
  use mod_apse_eqn
  use mod_lpse_dis_normal
  use mod_lst_eqn_IR

  implicit none

  private

  type, extends(solver_type), public :: APSE_LPSE_type

    private
    type(lpse_eqn_type) :: lpse  !< PSE求解器
    type(apse_eqn_type) :: apse  !< PSE求解器
    type(dis_type), pointer :: Dis_LPSE !< 二维扰动场
    type(dis_type), pointer :: Dis_APSE !< 二维扰动场的伴随

    contains

    procedure :: SetDisInlet=>SetDisInlet_IR_APSE !<设置入口扰动
    procedure :: SetWholeDis !< 设置全局求解对应的扰动
    procedure :: Print !< 输出PSE的计算结果
    procedure :: Solve=>Solve_APSE_LPSE !< 推进求解PSE

    procedure, private :: SetDisInlet_IR_PSE !< 设置入口扰动
    procedure, private :: SetDisInlet_IR_APSE !<设置入口扰动伴随
    procedure, private :: PrintSigma => print_sigma !< 输出物理流向复波数

  end type APSE_LPSE_type

contains

  !> 设置求解对应的扰动
  subroutine SetWholeDis(this, Dis_LPSE, Dis_APSE)

      implicit none
      class(APSE_LPSE_type), intent(inout) :: This
      type(dis_type), target, intent(in) :: Dis_LPSE !< 二维扰动场
      type(dis_type), target, intent(in) :: Dis_APSE !< 二维扰动场的伴随

      if(.not. (associated(this%Dis_LPSE, Dis_LPSE))) this%Dis_LPSE => Dis_LPSE
      if(.not. (associated(this%Dis_APSE, Dis_APSE))) this%Dis_APSE => Dis_APSE

  end subroutine SetWholeDis

  !> LPSE推进求解
  subroutine Solve_APSE_LPSE(this)

      use mod_BF_normal
      use mod_dis_normal
      use mod_local_normal_coordinate
      use mod_difference

      implicit none
      class(APSE_LPSE_type), intent(inout) :: this
      type(lpse_dis_normal_type) :: DisNorm, DisNormFront, DisNormBack
      type(BF_normal_type), save :: BFNorm
      type(local_normal_coordinate_type), save :: NormCoord
      type(Norm_Diff_Coef_type), save :: NormCoef
      real(R_P), allocatable, save :: Eta(:)
      integer :: iloc, iend
      logical, parameter :: ISWITHALPHA=.false.


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

      !> APSE process

      call this%apse%SetSolver(method=1001)
      call this%apse%Create(this%Grid%GetJnSize())

#ifdef DEBUG
      print*, this%iloc
      pause
#endif
      do iloc=This%Iloc-1, this%IStart, -1
        call this%PrintILoc(iloc)
        call this%Dis_APSE%Get(iloc, DisNorm)

        call this%Dis_APSE%Get(iloc+1, DisNormBack )

        call BFNorm%Set(iloc, this%BaseFlow, this%Diff)
        call NormCoord%SetLameFromGrid(this%grid, iloc)
        call NormCoord%SetCurvature(this%Curvature%GetPoint(iloc))
        call this%Diff%GetLocalNormDiffCoef(iloc, NormCoef)
        call this%grid%Get_iy(iloc, Eta)

        DisNorm=DisNormBack

        call this%apse%SolvePSE(DisNorm, iloc, DisNormBack, &
        &   BFNorm, NormCoord, NormCoef, Eta, this%bctype)

        call this%Dis_APSE%Set(iloc, DisNorm)
      end do

      !> LPSE process

      call this%SetDisInlet_IR_PSE(DisNorm%GetWaveNum(), this%IStart)

      call this%lpse%SetSolver(method=1001)
      call this%lpse%Create(this%Grid%GetJnSize())

      do iloc=This%Istart+1, this%ILoc-1
        call this%PrintILoc(iloc)
        call this%Dis_APSE%Get(iloc, DisNorm)

        call this%Dis_LPSE%Get(iloc-1, DisNormFront )

        call BFNorm%Set(iloc, this%BaseFlow, this%Diff)
        call NormCoord%SetLameFromGrid(this%grid, iloc)
        call NormCoord%SetCurvature(this%Curvature%GetPoint(iloc))
        call this%Diff%GetLocalNormDiffCoef(iloc, NormCoef)
        call this%grid%Get_iy(iloc, Eta)

        !DisNorm=DisNormFront

        call this%lpse%SolvePSE(DisNorm, iloc, DisNormFront, &
        &   BFNorm, NormCoord, NormCoef, Eta, this%bctype, &
        &   isWithAlpha=ISWITHALPHA)

        call this%Dis_LPSE%Set(iloc, DisNorm)
      end do

      do iloc=This%iloc, this%Iend
        call this%PrintILoc(iloc)
        call this%Dis_LPSE%Get(iloc, DisNorm)

        call this%Dis_LPSE%Get(iloc-1, DisNormFront )

        call BFNorm%Set(iloc, this%BaseFlow, this%Diff)
        call NormCoord%SetLameFromGrid(this%grid, iloc)
        call NormCoord%SetCurvature(this%Curvature%GetPoint(iloc))
        call this%Diff%GetLocalNormDiffCoef(iloc, NormCoef)
        call this%grid%Get_iy(iloc, Eta)

        call this%lpse%SolvePSE(DisNorm, iloc, DisNormFront, &
        &   BFNorm, NormCoord, NormCoef, Eta, this%bctype)

        call this%Dis_LPSE%Set(iloc, DisNorm)
      end do

  end subroutine Solve_APSE_LPSE

  !> 输出PSE的计算结果
  subroutine Print(this, fn_surf)

      use stringifor
      implicit none
      class(APSE_LPSE_type), intent(in) :: this
      type(string), intent(in) :: fn_surf !< 文件名前缀
      type(string) :: fn_lpse
      integer, parameter :: LPSE=0, ALPSE=1, LST=2, ALPSE_LPSE=3, LPSE_ALPSE=4

      fn_lpse=fn_surf//'_LPSE'

      call this%Grid%Print(fn_lpse)
      call this%Dis_LPSE%Print(fn_lpse)
      call this%PrintSigma(fn_lpse)

      fn_lpse=fn_surf//'_APSE'

      call this%Grid%Print(fn_lpse)
      call this%Dis_APSE%Print(fn_lpse)

  end subroutine Print

  !> 给入口边界条件(反幂法求LST).
  !!
  !! 先求解入口处的LST问题,作为PSE入口边界条件
  subroutine SetDisInlet_IR_PSE(this, wavenum, iloc_in)

      use mod_dis_wavenum, only: dis_wavenum_type, dis_wavenum_lpse_type
      use mod_BF_normal
      use mod_dis_normal
      use mod_local_normal_coordinate
      use mod_difference
      use mod_lst_eqn_IR

      implicit none
      class(APSE_LPSE_type), intent(inout) :: this
      class(dis_wavenum_type), intent(in) :: wavenum !< 初始扰动色散特征
      integer, optional, intent(in) :: iloc_in
      type(dis_wavenum_lpse_type) :: wavenum_pse
      type(lpse_dis_normal_type) :: DisNorm
      type(BF_normal_type), save :: BFNorm
      type(local_normal_coordinate_type), save :: NormCoord
      type(Norm_Diff_Coef_type), save :: NormCoef
      real(R_P), allocatable, save :: Eta(:)
      integer :: iloc
      complex(R_P) :: sigma(7)

      type(lst_eqn_IR_type) :: LST

      if(.not. present(iloc_in)) then
        iloc=this%iloc
      else
        iloc=iloc_in
      endif

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

      call this%Dis_LPSE%Get(iloc, DisNorm)

      call wavenum_pse%set(wavenum)
      call DisNorm%SetWavenum(wavenum_pse)

      call BFNorm%Set(iloc, this%BaseFlow, this%Diff)
      call NormCoord%SetLameFromGrid(this%grid, iloc)
      call NormCoord%SetCurvature(this%Curvature%GetPoint(iloc))
      call this%Diff%GetLocalNormDiffCoef(iloc, NormCoef)
      call this%grid%Get_iy(iloc, Eta)

      call LST%SolveLST(DisNorm, iloc, &
      &   BFNorm, NormCoord, Eta, NormCoef)
      Sigma=DisNorm%GetAlpha()
      call DisNorm%SetSigma(Sigma)

      call this%Dis_LPSE%Set(iloc, DisNorm)

      call this%Dis_LPSE%Print(iloc, this%Grid)

  end subroutine SetDisInlet_IR_PSE

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
      class(APSE_LPSE_type), intent(inout) :: this
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
      if(iloc_in==this%Istart) this%iloc=this%Iend
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

      call this%Dis_APSE%Get(iloc, DisNorm)

      call wavenum_pse%set(wavenum)
      call DisNorm%SetWavenum(wavenum_pse)

      call BFNorm%Set(iloc, this%BaseFlow, this%Diff)
      call NormCoord%SetLameFromGrid(this%grid, iloc)
      call NormCoord%SetCurvature(this%Curvature%GetPoint(iloc))
      call this%Diff%GetLocalNormDiffCoef(iloc, NormCoef)
      call this%grid%Get_iy(iloc, Eta)

      call LST%SolveLST(DisNorm, iloc, &
      &   BFNorm, NormCoord, Eta, NormCoef, ISADJOINT)

      call this%Dis_APSE%Set(iloc, DisNorm)

      call this%Dis_APSE%Print(this%iloc, this%Grid)

  end subroutine SetDisInlet_IR_APSE

  !> 输出物理流向复波数
  subroutine print_sigma(this, fn_surf)

      use stringifor
      use mod_dis_flux
      use mod_baseflow_org

      implicit none
      class(APSE_LPSE_type), intent(in) :: this
      type(string), intent(in) :: fn_surf !< 文件名前缀
      integer :: i, l, j
      complex(R_P) :: sigma(7, This%Istart:This%Iend)
      real(R_P) :: GrowthRate(7, This%Istart:This%Iend)
      real(R_P) :: Nfact(7, This%Istart:This%Iend)
      real(R_P) :: Amp(7, This%Istart:This%Iend)
      real(R_P) :: Amp_pse(7, This%Istart:This%Iend)
      real(R_P) :: Amp_local(7)
      type(lpse_dis_normal_type) :: DisNorm
      complex(R_P) :: rho, u, v, w, T
      type(BF_flux_org_ij_type), allocatable :: BFFlux(:, :)
      type(dis_flux_ij_type), allocatable :: Flux(:)
      real(R_P) :: Rho0, U0, V0, W0, T0
      integer :: in, Jn
      real(R_P), allocatable :: xx(:)

      In=this%Grid%GetInSize()
      Jn=this%Grid%GetJnSize()

      allocate(BFFlux(1, Jn))
      allocate(Flux(jn))

      if(.not. allocated(xx)) allocate(xx(in))
      call this%Grid%Get_jx(1, xx)
      do i=This%istart, this%Iend
        call this%Dis_LPSE%Get(i, DisNorm)
        call DisNorm%GetSigma(sigma(:, i))
      end do
      GrowthRate=-aimag(sigma)

      Amp(:, This%Istart)=1.0d0
      Nfact(:, This%Istart)=0.0d0
      amp_pse(:, This%Istart)=1.0d0
      do i=This%Istart+1, This%Iend
        do l=1, 7
          Nfact(l, i)=Nfact(l, i-1)+(GrowthRate(l, i)+GrowthRate(l, i-1))*0.5d0 &
                                    *(xx(i)-xx(i-1))
          Amp(l, i)=exp(Nfact(l, i))
        enddo
        call this%Dis_LPSE%Get(i, DisNorm)
        DisNorm=Amp(6, i) .mx. DisNorm
        call this%BaseFlow%GetPart(i, i, 1, jn, BFFlux(:, :))
        call DisNorm%GetFlux(Flux)
        Amp_pse(:, i)=0.0d0
        do j=1, jn
          call BFFlux(1, j)%get(rho0, u0, v0, w0, T0)
          call Flux(j)%get(rho, u, v, w, T)
          Amp_local(1)=abs(rho*u0+rho0+u)
          Amp_local(2)=abs(rho)
          Amp_local(3)=abs(u)
          Amp_local(4)=abs(v)
          Amp_local(5)=abs(T)
          Amp_local(6)=0.0d0
          Amp_local(7)=abs(w)
          do l=1, 7
            if(amp_local(l)>=Amp_pse(l, i)) Amp_pse(l, i)=amp_local(l)
          enddo
        enddo
        amp_pse(6, i)=amp(6, i)
      end do

      open(997, file=fn_surf//'_Alf.plt', form='formatted')
      write(997, *)"variables=x, N_u, a_r, s_u, A_u"
      do i=This%Istart, this%Iend
        write(997, '(5ES20.8)')xx(i), Nfact(3, i), real(sigma(6, i)), -aimag(sigma(3, i)), Amp(3, i)
      enddo
      close(997)

      open(996, file=fn_surf//'_GrowthRate.plt', form='formatted')
      write(996, *)"variables=x, s_rhou, s_rho, s_u, s_v, s_T, s_E, s_w"
      do i=This%Istart, this%Iend
        write(996, '(8ES20.8)')xx(i), (GrowthRate(l, i), l=1, 7)
      enddo
      close(996)

      open(995, file=fn_surf//'_Amp.plt', form='formatted')
      write(995, *)"variables=x, A_rhou, A_rho, A_u, A_v, A_T, A_E, A_w"
      do i=This%Istart, this%Iend
        write(995, '(8ES20.8)')xx(i), (Amp(l, i), l=1, 7)
      enddo
      close(995)

      open(994, file=fn_surf//'_Nfact.plt', form='formatted')
      write(994, *)"variables=x, N_rhou, N_rho, N_u, N_v, N_T, N_E, N_w"
      do i=This%Istart, this%Iend
        write(994, '(8ES20.8)')xx(i), (Nfact(l, i), l=1, 7)
      enddo
      close(994)

      open(993, file=fn_surf//'_Amp_pse.plt', form='formatted')
      write(993, *)"variables=x, A_rhou, A_rho, A_u, A_v, A_T, A_E, A_w"
      do i=This%Istart, this%Iend
        write(993, '(8ES20.8)')xx(i), (Amp_pse(l, i), l=1, 7)
      enddo
      close(993)

      deallocate(BFFlux)
      deallocate(Flux)

  end subroutine print_sigma


end module mod_APSE_LPSE
