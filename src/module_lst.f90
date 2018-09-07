!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_lst.f90
!> @file
!> @breif LST计算全场的特征值.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: mod_lst
!> @breif 使用LST计算全场的特征值和特征函数.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2018-06-26 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2018-06-26
!------------------------------------------------------------------------------
module mod_lst

  use mod_solver
  use mod_lst_eqn_IR
  use mod_dis

  implicit none

  private

  type, extends(solver_type), public :: lst_type
    private
    type(lst_eqn_IR_type) :: lst  !< PSE求解器
    type(dis_type), pointer :: Dis !< 二维扰动场

    contains

    procedure :: SetWholeDis !< 设置全局求解对应的扰动
    procedure :: Print !< 输出LST的计算结果
    procedure :: Solve=>SolveLST !< 求解LST

    procedure, private :: PrintWave !< 输出物理流向复波数

  end type lst_type

contains

  !> 设置求解对应的扰动
  subroutine SetWholeDis(this, Dis)

      implicit none
      class(lst_type), intent(inout) :: This
      type(dis_type), target, intent(in) :: Dis !< 二维扰动场

      if(.not. (associated(this%Dis, Dis))) this%Dis => Dis

  end subroutine SetWholeDis

  !> 输出PSE的计算结果
  subroutine Print(this, fn_surf)

      use stringifor
      implicit none
      class(lst_type), intent(in) :: this
      type(string), intent(in) :: fn_surf !< 文件名前缀
      type(string) :: fn_lst
      integer, parameter :: LPSE=0, ALPSE=1, LST=2, ALPSE_LPSE=3, LPSE_ALPSE=4

      fn_lst=fn_surf//'_LST'

      call this%Grid%Print(fn_lst)
      call this%Dis%Print(fn_lst)
      call this%PrintWave(fn_lst)

  end subroutine Print

  !> 求解过程(整个计算域）
  subroutine SolveLST(this, iloc, wavenum)

    use mod_dis_wavenum, only: dis_wavenum_type, dis_wavenum_lpse_type
    use mod_BF_normal
    use mod_dis_normal
    use mod_local_normal_coordinate
    use mod_difference
    use mod_lpse_dis_normal

    implicit none
    class(lst_type), intent(inout) :: this
    class(dis_wavenum_type), intent(in) :: wavenum !< 初始扰动色散特征
    integer, intent(in) :: iloc
    type(dis_wavenum_lpse_type) :: wavenum_
    type(dis_wavenum_lpse_type) :: wavenum_pse
    type(lpse_dis_normal_type) :: DisNorm
    type(BF_normal_type), save :: BFNorm
    type(local_normal_coordinate_type), save :: NormCoord
    type(Norm_Diff_Coef_type), save :: NormCoef
    real(R_P), allocatable, save :: Eta(:)
    integer :: iend, iCount

    this%iloc=iloc


    if(this%iend==-1)then
      iend=this%Grid%GetInSize()
      this%iend=iend
    end if

    if(this%istart==-1) this%istart=1

    if(.not. BFNorm%HasCreated()) then
      call BFNorm%Create(this%Grid%GetJnSize())
    endif

    if(.not. NormCoord%HasCreated()) then
      call NormCoord%Create(this%Grid%GetJnSize())
    endif

    if(.not. allocated(eta)) then
      allocate(Eta(this%Grid%GetJnSize()))
    endif

    !if(.not. this%LSTHasCreated()) then
      call this%LST%Create(this%Grid%GetJnSize())
    !endif

    call wavenum_pse%set(wavenum)
    wavenum_=wavenum_pse
    do iCount=this%iloc, this%iend

      write(*, *)'iLoc=', iCount
      call this%Dis%Get(iCount, DisNorm)

      call DisNorm%SetWavenum(wavenum_pse)

      call BFNorm%Set(iCount, this%BaseFlow, this%Diff)
      call NormCoord%SetLameFromGrid(this%grid, iCount)
      call NormCoord%SetCurvature(this%Curvature%GetPoint(iCount))
      call this%Diff%GetLocalNormDiffCoef(iCount, NormCoef)
      call this%grid%Get_iy(iCount, Eta)

      call this%LST%SolveLST(DisNorm, iCount, &
      &   BFNorm, NormCoord, Eta, NormCoef)

      call this%Dis%Set(iCount, DisNorm)
      wavenum_pse=DisNorm%GetWaveNum()
      write(*, *)wavenum_pse%getAlpha()

      !call this%Dis%Print(this%iloc, this%Grid)
    enddo
    wavenum_pse=wavenum_
    do iCount=this%iloc-1, this%istart, -1

      write(*, *)'iLoc=', iCount
      call this%Dis%Get(iCount, DisNorm)

      call DisNorm%SetWavenum(wavenum_pse)

      call BFNorm%Set(iCount, this%BaseFlow, this%Diff)
      call NormCoord%SetLameFromGrid(this%grid, iCount)
      call NormCoord%SetCurvature(this%Curvature%GetPoint(iCount))
      call this%Diff%GetLocalNormDiffCoef(iCount, NormCoef)
      call this%grid%Get_iy(iCount, Eta)

      call this%LST%SolveLST(DisNorm, iCount, &
      &   BFNorm, NormCoord, Eta, NormCoef)

      call this%Dis%Set(iCount, DisNorm)
      wavenum_pse=DisNorm%GetWaveNum()
      write(*, *)wavenum_pse%getAlpha()

      !call this%Dis%Print(this%iloc, this%Grid)
    enddo

  end subroutine SolveLST

  ! !> 使用Muller法求全场的线性稳定性
  ! subroutine SolveLST_muller(this, iloc, wavenum)
  !
  !   use mod_dis_wavenum, only: dis_wavenum_type, dis_wavenum_lpse_type
  !   use mod_BF_normal
  !   use mod_dis_normal
  !   use mod_local_normal_coordinate
  !   use mod_difference
  !   use mod_lpse_dis_normal
  !   use mod_lst_eqn_muller
  !
  !   implicit none
  !   class(lst_type), intent(inout) :: this
  !   class(dis_wavenum_type), intent(in) :: wavenum !< 初始扰动色散特征
  !   integer, intent(in) :: iloc
  !   type(dis_wavenum_lpse_type) :: wavenum_
  !   type(dis_wavenum_lpse_type) :: wavenum_pse
  !   type(lpse_dis_normal_type) :: DisNorm
  !   type(BF_normal_type), save :: BFNorm
  !   type(local_normal_coordinate_type), save :: NormCoord
  !   type(Norm_Diff_Coef_type), save :: NormCoef
  !   real(R_P), allocatable, save :: Eta(:)
  !   integer :: iend, iCount
  !   type(lst_eqn_muller_type) :: LST
  !   integer, parameter :: SPATIAL_LST=1, TEMPORAL_LST=2
  !
  !   this%iloc=iloc
  !
  !
  !   if(this%iend==-1)then
  !     iend=this%Grid%GetInSize()
  !     this%iend=iend
  !   end if
  !
  !   if(this%istart==-1) this%istart=1
  !
  !   if(.not. BFNorm%HasCreated()) then
  !     call BFNorm%Create(this%Grid%GetJnSize())
  !   endif
  !
  !   if(.not. NormCoord%HasCreated()) then
  !     call NormCoord%Create(this%Grid%GetJnSize())
  !   endif
  !
  !   if(.not. allocated(eta)) then
  !     allocate(Eta(this%Grid%GetJnSize()))
  !   endif
  !
  !   !if(.not. this%LSTHasCreated()) then
  !     !call LST%Initialize()
  !     !call LST%SetMode(SPATIAL_LST)
  !   !endif
  !   call LST%Create(this%grid%GetJnSize())
  !
  !   call wavenum_pse%set(wavenum)
  !   wavenum_=wavenum_pse
  !   do iCount=this%iloc, this%iend
  !
  !     write(*, *)'iLoc=', iCount
  !     call this%Dis%Get(iCount, DisNorm)
  !
  !     call DisNorm%SetWavenum(wavenum_pse)
  !
  !     call BFNorm%Set(iCount, this%BaseFlow, this%Diff)
  !     call this%grid%Get_iy(iCount, Eta)
  !
  !     call NormCoord%SetCurvature(this%Curvature%GetPoint(iCount))
  !
  !     !call LST%SetCurvature(this%Curvature%GetPoint(iCount))
  !     !call LST%SetEta(this%Grid%GetJnSize(), Eta)
  !     !call LST%SetBaseflow(BFNorm)
  !     call DisNorm%SetWavenum(wavenum_pse)
  !     !call LST%Solve(DisNorm)
  !     call LST%SolveLST(DisNorm, iCount, &
  !     &   BFNorm, NormCoord, Eta)
  !
  !     call DisNorm%SetILoc(iCount)
  !
  !     call this%Dis%Set(iCount, DisNorm)
  !     wavenum_pse=DisNorm%GetWaveNum()
  !     write(*, *)wavenum_pse%getAlpha()
  !
  !     !call this%Dis%Print(this%iloc, this%Grid)
  !   enddo
  !   wavenum_pse=wavenum_
  !   do iCount=this%iloc-1, this%istart, -1
  !
  !     write(*, *)'iLoc=', iCount
  !     call this%Dis%Get(iCount, DisNorm)
  !
  !     call DisNorm%SetWavenum(wavenum_pse)
  !
  !     call BFNorm%Set(iCount, this%BaseFlow, this%Diff)
  !     call this%grid%Get_iy(iCount, Eta)
  !
  !     !call LST%SetCurvature(this%Curvature%GetPoint(iCount))
  !     !call LST%SetEta(this%Grid%GetJnSize(), Eta)
  !     call NormCoord%SetCurvature(this%Curvature%GetPoint(iCount))
  !     !call LST%SetBaseflow(BFNorm)
  !     call DisNorm%SetWavenum(wavenum_pse)
  !     !call LST%Solve(DisNorm)
  !     call LST%SolveLST(DisNorm, iCount, &
  !     &   BFNorm, NormCoord, Eta)
  !     call DisNorm%SetILoc(iCount)
  !
  !     call this%Dis%Set(iCount, DisNorm)
  !     wavenum_pse=DisNorm%GetWaveNum()
  !     write(*, *)wavenum_pse%getAlpha()
  !
  !     !call this%Dis%Print(this%iloc, this%Grid)
  !   enddo
  !
  ! end subroutine SolveLST_muller


  !> 输出物理流向复波数
  subroutine PrintWave(this, fn_surf)

      use stringifor
      !use mod_dis_flux
      use mod_lpse_dis_normal
      use mod_dis_wavenum

      implicit none
      class(lst_type), intent(in) :: this
      type(string), intent(in) :: fn_surf !< 文件名前缀
      integer :: i, l, j
      complex(R_P) :: Alpha(This%Istart:This%Iend)
      complex(R_P) :: Beta(This%Istart:This%Iend)
      complex(R_P) :: Omega(This%Istart:This%Iend)
      real(R_P) :: GrowthRate(This%Istart:This%Iend)
      real(R_P) :: Nfact(This%Istart:This%Iend)
      real(R_P) :: Amp(This%Istart:This%Iend)
      type(lpse_dis_normal_type) :: DisNorm
      type(dis_wavenum_lpse_type) :: wavenum
      integer :: in, Jn
      real(R_P), allocatable :: xx(:)

      In=this%Grid%GetInSize()
      Jn=this%Grid%GetJnSize()

      if(.not. allocated(xx)) allocate(xx(in))
      call this%Grid%Get_jx(1, xx)
      do i=This%istart, this%Iend
        call this%Dis%Get(i, DisNorm)
        wavenum=DisNorm%GetWaveNum()
        call wavenum%get(alpha(i), beta(i), omega(i))
      end do
      GrowthRate=-aimag(alpha)

      Amp(This%Istart)=1.0d0
      Nfact(This%Istart)=0.0d0
      do i=This%Istart+1, This%Iend
          Nfact(i)=Nfact(i-1)+(GrowthRate(i)+GrowthRate(i-1))*0.5d0 &
                                    *(xx(i)-xx(i-1))
          Amp(i)=exp(Nfact(i))
      end do

      open(997, file=fn_surf//'_WaveNum.plt', form='formatted')
      write(997, *)"variables=x, w, a_r, a_i, c, b"
      do i=This%Istart, this%Iend
        write(997, '(6ES20.8)')xx(i), real(omega(i)), real(alpha(i)), &
        &         aimag(alpha(i)), real(omega(i))/real(alpha(i)), real(beta(i))
      enddo
      close(997)

      open(996, file=fn_surf//'_Growthrate.plt', form='formatted')
      write(996, *)"variables=x, s"
      do i=This%Istart, this%Iend
        write(996, '(8ES20.8)')xx(i), GrowthRate(i)
      enddo
      close(996)

      open(995, file=fn_surf//'_Amp.plt', form='formatted')
      write(995, *)"variables=x, Amp, N"
      do i=This%Istart, this%Iend
        write(995, '(8ES20.8)')xx(i), Amp(i), Nfact(i)
      enddo
      close(995)

  end subroutine PrintWave


end module mod_lst
