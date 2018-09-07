!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_npse.f90
!> @file
!> @breif NPSE求解器文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: mod_npse
!> @breif NPSE求解器模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2018-06-22 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2018-06-22
!------------------------------------------------------------------------------
module mod_npse

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

    use mod_cfgio_adapter, only: fn_surf=>prefix

    use mod_dis_wavenum
    use mod_npse_eqn

    implicit none

    private

    !> NPSE求解器类
    type, extends(solver_type), public :: npse_type

        private
        type(npse_eqn_type), private :: npse  !< NPSE求解器
        integer, private :: m_dim !<时间方向离散点 0:m_dim
        integer, private :: n_dim !<展向离散点 -n_dim:n_dim
        type(lpse_dis_normal_type), allocatable, private :: DisNorm(:, :) !<当前站位扰动
        type(lpse_dis_normal_type), allocatable, private :: DisNormFront(:, :)!<前一站位扰动
        logical, allocatable, private :: isWithAlpha(:, :)

        type(dis_type), pointer :: Dis !< 二维扰动场(仅输出用，非做应用使用)

        contains

        generic :: SetDisInlet_NPSE=> SetDisInlet_IR_NPSE_single, & !< 设置入口扰动
                                      &   SetDisInlet_IR_NPSE_multi

        procedure :: SetDim  !< 设置时间展向空间离散点数
        procedure :: SetWholeDis !< 设置全局求解对应的扰动
        procedure :: Print !< 输出PSE的计算结果
        procedure :: Solve=>Solve_NPSE !< 推进求解PSE

        procedure, private :: PrintSigma => print_sigma !< 输出物理流向复波数
        procedure, private :: SetDisInlet_IR_NPSE_single !< 设置入口扰动（单个波）
        procedure, private :: SetDisInlet_IR_NPSE_multi  !< 设置入口扰动（多个波）
        procedure, private :: CopyToDis !<将DisNorm写入到扰动变量
        procedure, private :: Write !< 将当前扰动写入到文件中

    end type npse_type

    contains

    !> 时间展向空间离散点数设置.
    subroutine SetDim(this, m_dim, n_dim)

        use mod_parameter, only: EPS_REL
        implicit none
        class(npse_type), intent(inout) :: this
        integer, intent(in) :: m_dim !< 时间法向离散点数
        integer, intent(in) :: n_dim !< 空间方向离散点数
        integer :: m, n

        this%m_dim=m_dim-1  !< 时间方向8个离散点则对应0:7
        this%n_dim=n_dim-1  !< 空间法向8个离散点则对应-7：7

        write(*, *)'The dimension in time and spanwise direction is:'
        write(*, *)'0:', this%m_dim, ',   -',this%n_dim,':', this%n_dim

#IFDEF DEBUG
        print*, this%m_dim, this%n_dim, 1
#ENDIF
        if(allocated(this%DisNorm)) deallocate(this%DisNorm)
        if(allocated(this%DisNormFront)) deallocate(this%DisNormFront)
        if(allocated(this%isWithAlpha)) deallocate(this%isWithAlpha)
        if(.not. allocated(this%DisNorm)) &
        & ALLOCATE(this%DisNorm(0:this%m_dim, -this%n_dim:this%n_dim))
#IFDEF DEBUG
        print*, this%m_dim, this%n_dim, 2
#ENDIF
        if(.not. allocated(this%DisNormFront)) &
        & ALLOCATE(this%DisNormFront(0:this%m_dim, -this%n_dim:this%n_dim))
        if(.not. allocated(this%isWithAlpha)) &
        & allocate(this%isWithalpha(0:this%m_dim, -this%n_dim:this%n_dim))
        this%isWithAlpha=.False.
#IFDEF DEBUG
        print*, this%m_dim, this%n_dim, 3
#ENDIF
        do n=-this%n_dim, this%n_dim
          do m=0, this%m_dim
            call this%disNorm(m, n)%create(this%grid%GetJnSize())
            call this%disNormFront(m, n)%create(This%grid%GetJnSize())
          enddo
        enddo

        !EPS_REL=1d-6

    end subroutine SetDim

    !> 设置求解对应的扰动
    subroutine SetWholeDis(this, Dis)

        implicit none
        class(npse_type), intent(inout) :: This
        type(dis_type), target, intent(in) :: Dis !< 二维扰动场

        if(.not. (associated(this%Dis, Dis))) this%Dis => Dis

    end subroutine SetWholeDis

    !> LPSE推进求解
    subroutine Solve_NPSE(this)

        use mod_BF_normal
        use mod_dis_normal
        use mod_local_normal_coordinate
        use mod_difference

        implicit none
        class(npse_type), intent(inout) :: this
        type(BF_normal_type), save :: BFNorm
        type(local_normal_coordinate_type), save :: NormCoord
        type(Norm_Diff_Coef_type), save :: NormCoef
        real(R_P), allocatable, save :: Eta(:)
        integer :: iloc, iend

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

        !call this%npse%SetSolver(method=1001)
        call this%npse%Create(this%m_dim, this%n_dim, this%Grid%GetJnSize())

        do iloc=This%Istart+1, this%Iend

          call this%PrintILoc(iloc)

          this%DisNorm=this%DisNormFront
          call BFNorm%Set(iloc, this%BaseFlow, this%Diff)
          call NormCoord%SetLameFromGrid(this%grid, iloc)
          call NormCoord%SetCurvature(this%Curvature%GetPoint(iloc))
          call this%Diff%GetLocalNormDiffCoef(iloc, NormCoef)
          call this%grid%Get_iy(iloc, Eta)

          call this%npse%SolveNPSE(this%DisNorm, iloc, this%DisNormFront, &
          &               BFNorm, NormCoord, NormCoef, Eta, this%BCtype, &
          &               this%isWithAlpha)

          call this%write(iloc)
          this%DisNormFront=this%DisNorm

        end do

    end subroutine Solve_NPSE

    !> 输出PSE的计算结果
    subroutine Print(this, fn_surf)

        use stringifor
        use mod_cfgio_adapter
        implicit none
        class(npse_type), intent(inout) :: this
        type(string), intent(in) :: fn_surf !< 文件名前缀
        type(string) :: fn_npse, fn_npse_sub
        integer :: m, n, l


        fn_npse=fn_surf//'_NPSE'

        call this%Grid%Print(fn_npse)

        do n=-this%n_dim, this%n_dim
         do m=0, this%m_dim
           print*, 'Writing....  m=', m, 'n=', n
            call this%CopyToDis(m, n, fn_npse)
            fn_npse_sub=GenFileName(fn_npse, m=m, n=n)
            do l=1, size(m_index)
              if((m==m_index(l) .and. n==n_index(l)) .or. (m==0 .and. n==0))then
                call this%Dis%Print(fn_npse_sub)
              endif
            enddo
            call this%PrintSigma(fn_npse_sub)
          enddo
        enddo

    end subroutine Print

    function GenFileName(fn, m, n, iloc) result(fn_out)

      implicit none
      type(string), intent(in) :: fn
      integer, intent(in), optional :: m
      integer, intent(in), optional :: n
      integer, intent(in), optional :: iloc
      type(string) :: fn_out
      type(string) :: strm, strn, striloc

      fn_out=fn
      if(present(m)) then
        strm=m; strm='_'//strm
      endif
      if(present(n)) then
        strn=n; strn='_'//strn
      endif
      if(present(iloc)) then
        striloc=iloc; striloc='_'//striloc
      endif

      fn_out=fn_out//strm//strn//striloc

    end function GenFileName

    !> 将对应的结果写入到Dis变量
    subroutine CopyToDis(this, m, n, fn)

      implicit none
      class(npse_type), intent(inout) :: this
      integer, intent(in) :: m, n
      type(string), intent(in) :: fn
      type(string) :: fn_dis
      integer :: i
      type(string) :: command

      do i=this%istart, this%iend
        fn_dis=GenFileName(fn, m=m, n=n, iloc=i)
        call this%DisNorm(m, n)%Read(fn_dis)
        command='rm -f '//fn_dis//'_DisNorm.dat'
        call execute_command_line(command//'', wait=.false.)
        call this%Dis%Set(i, this%DisNorm(m, n))
      enddo

    end subroutine CopyToDis

    !> 给入口边界条件(反幂法求LST)--单个波.
    !!
    !! 先求解入口处的LST问题,作为PSE入口边界条件
    subroutine SetDisInlet_IR_NPSE_single(this, wavenum, m_index, n_index, &
                                          & amp, iloc_in)

        use mod_dis_wavenum, only: dis_wavenum_type, dis_wavenum_lpse_type
        use mod_BF_normal
        use mod_dis_normal
        use mod_local_normal_coordinate
        use mod_difference

        implicit none
        class(npse_type), intent(inout) :: this
        class(dis_wavenum_type), intent(in) :: wavenum !< 初始扰动色散特征
        integer, optional, intent(in) :: iloc_in !<流向坐标
        integer, intent(in) :: m_index !<时间序号
        integer, intent(in) :: n_index !<展向空间序号
        real(R_P), intent(in) :: amp !< 扰动初始幅值
        type(dis_wavenum_lpse_type) :: wavenum_pse
        type(lpse_dis_normal_type) :: DisNorm
        type(BF_normal_type), save :: BFNorm
        type(local_normal_coordinate_type), save :: NormCoord
        type(Norm_Diff_Coef_type), save :: NormCoef
        real(R_P), allocatable, save :: Eta(:)
        integer :: iloc

        type(lst_eqn_ir_type) :: LST
        complex(R_P) :: sigma(7)

        if(present(iloc_in)) then
          this%iloc=iloc_in
        else
          this%iloc=this%Istart
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

        DisNorm=this%DisNorm(m_index, n_index)

        call wavenum_pse%set(wavenum)
        call DisNorm%SetWavenum(wavenum_pse)

        call BFNorm%Set(iloc, this%BaseFlow, this%Diff)
        call NormCoord%SetLameFromGrid(this%grid, iloc)
        call NormCoord%SetCurvature(this%Curvature%GetPoint(iloc))
        call this%Diff%GetLocalNormDiffCoef(iloc, NormCoef)
        call this%grid%Get_iy(iloc, Eta)

        call LST%SolveLST(DisNorm, iloc, &
        &   BFNorm, NormCoord, Eta, NormCoef)

        wavenum_pse=DisNorm%GetWaveNum()
        DisNorm=amp .mx. DisNorm
        call DisNorm%SetWavenum(wavenum_pse)

        this%DisNorm(m_index, n_index)=DisNorm
        this%DisNormFront(m_index, n_index)=DisNorm
        Sigma=DisNorm%GetAlpha()
        call DisNorm%SetSigma(Sigma)

        wavenum_pse=DisNorm%GetWaveNum()

        write(*, *)'The Disturbance added in the flow is:'
        write(*, *)m_index, n_index, real(wavenum_pse%getOmega()), &
                  real(wavenum_pse%getBeta()), wavenum_pse%getAlpha()

        call this%write(iloc)

        this%isWithAlpha(m_index, n_index)=.True.

    end subroutine SetDisInlet_IR_NPSE_single

    !> 给入口边界条件(反幂法求LST)--多个波.
    !!
    !! 先求解入口处的LST问题,作为PSE入口边界条件
    subroutine SetDisInlet_IR_NPSE_multi(this, wavenum, m_index, n_index, &
                                         amp, iloc_in)

        use mod_dis_wavenum, only: dis_wavenum_type, dis_wavenum_lpse_type
        use mod_BF_normal
        use mod_dis_normal
        use mod_local_normal_coordinate
        use mod_difference
        use mod_cfgio_adapter, only: omega0, beta0

        implicit none
        class(npse_type), intent(inout) :: this
        class(dis_wavenum_type), intent(inout) :: wavenum(:) !< 初始扰动色散特征
        integer, intent(in), optional :: iloc_in !<流向坐标
        integer, intent(in) :: m_index(:) !<时间序号
        integer, intent(in) :: n_index(:) !<展向空间序号
        real(R_P), intent(in) :: amp(:)
        type(dis_wavenum_lpse_type) :: wavenum_pse
        type(lpse_dis_normal_type) :: DisNorm
        type(BF_normal_type), save :: BFNorm
        type(local_normal_coordinate_type), save :: NormCoord
        type(Norm_Diff_Coef_type), save :: NormCoef
        real(R_P), allocatable, save :: Eta(:)
        integer :: i, iloc
        type(lst_eqn_ir_type) :: LST
        type(string) :: fn
        complex(R_P) :: omega, beta
        complex(R_P) :: sigma(7)

        if(present(iloc_in)) then
          this%iloc=iloc_in
        else
          this%iloc=this%Istart
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

        write(*, *)'The Disturbance added in the flow is:'
        write(*, '(A, 2X, A, 10X, A, 14X, 3(A, 14X))') &
        &     'm_index', 'n_index', 'Amp', 'omega', 'beta', 'alpha'

        do i=1, size(wavenum, dim=1)

          DisNorm=this%DisNorm(m_index(i), n_index(i))
          beta=dble(n_index(i))*beta0
          omega=dble(m_index(i))*omega0

          call wavenum(i)%SetBeta(beta)
          call wavenum(i)%SetOmega(omega)

          call wavenum_pse%set(wavenum(i))
          call DisNorm%SetWavenum(wavenum_pse)
          call BFNorm%Set(iloc, this%BaseFlow, this%Diff)
          call NormCoord%SetLameFromGrid(this%grid, iloc)
          call NormCoord%SetCurvature(this%Curvature%GetPoint(iloc))
          call this%Diff%GetLocalNormDiffCoef(iloc, NormCoef)
          call this%grid%Get_iy(iloc, Eta)
#IFDEF DEBUG
print*, 'lst solver start'
#ENDIF
          call LST%SolveLST(DisNorm, iloc, &
          &   BFNorm, NormCoord, Eta, NormCoef)
          wavenum_pse=DisNorm%GetWaveNum()
          DisNorm=amp(i) .mx. DisNorm
          call DisNorm%SetWavenum(wavenum_pse)
#IFDEF DEBUG
print*, 'lst solver end'
#ENDIF

          this%DisNorm(m_index(i), n_index(i))=DisNorm
          this%DisNormFront(m_index(i), n_index(i))=DisNorm
          Sigma=DisNorm%GetAlpha()
          call DisNorm%SetSigma(Sigma)

          write(*, '(I5, 4X, I5, 3(3X, E16.8), 3X,"(",E16.8,",",2X,E16.8,")")') &
          &     m_index(i), n_index(i), amp(i), real(wavenum_pse%getOmega()), &
          &     real(wavenum_pse%getBeta()), wavenum_pse%getAlpha()
          if(i==1) this%isWithAlpha(m_index(i), n_index(i)) =.True.

        enddo

        call this%Write(iloc)

    end subroutine SetDisInlet_IR_NPSE_multi

    subroutine Write(this, iloc)
      implicit none
      class(npse_type) :: this
      integer, intent(in) :: iloc
      type(string) :: fn
      integer :: nCount, mCount

      do nCount=-this%n_dim, this%n_dim
        do mCount=0, this%m_dim

          fn=fn_surf//'_NPSE'
          fn=GenFileName(fn,m=mCount,n=nCount,iloc=iloc)

          call this%DisNorm(mCount, nCount)%write(fn)
        enddo
      enddo

    end subroutine Write



    !> 输出物理流向复波数
    subroutine print_sigma(this, fn_surf)

        use stringifor
        use mod_dis_flux

        implicit none
        class(npse_type), intent(in) :: this
        type(string), intent(in) :: fn_surf !< 文件名前缀
        integer :: i, l, j
        complex(R_P) :: Alf(This%Istart:This%Iend)
        real(R_P) :: GrowthRate(This%Istart:This%Iend)
        real(R_P) :: Nfact(This%Istart:This%Iend)
        real(R_P) :: Amp(This%Istart:This%Iend)
        real(R_P) :: Amp_pse(7, This%Istart:This%Iend)
        real(R_P) :: Amp_local(7)
        type(lpse_dis_normal_type) :: DisNorm
        complex(R_P) :: rho, u, v, w, T
        type(BF_flux_org_ij_type), allocatable :: BFFlux(:, :)
        type(dis_flux_ij_type), allocatable :: Flux(:)
        real(R_P) :: Rho0, U0, V0, W0, T0
        integer :: in, Jn
        real(R_P), allocatable :: xx(:), yy(:)
        real(R_P) :: E_local0, E_local1

        In=this%Grid%GetInSize()
        Jn=this%Grid%GetJnSize()

        allocate(BFFlux(1, Jn))
        allocate(Flux(jn))

        if(.not. allocated(xx)) allocate(xx(in))
        if(.not. allocated(yy)) allocate(yy(jn))
        call this%Grid%Get_jx(1, xx)
        do i=This%istart, this%Iend
          call this%Dis%Get(i, DisNorm)
          Alf(i)=DisNorm%GetAlpha()
        end do
        GrowthRate=-aimag(Alf)

        Amp(This%Istart)=1.0d0
        Nfact(This%Istart)=0.0d0
        amp_pse(:, This%Istart)=1.0d0
        do i=This%Istart+1, This%Iend
          Nfact(i)=Nfact(i-1)+(GrowthRate(i)+GrowthRate(i-1))*0.5d0 &
                                    *(xx(i)-xx(i-1))
          Amp(i)=exp(Nfact(i))
        enddo
        do i=This%Istart, This%Iend
          call this%Grid%get_iy(i, yy)
          call this%Dis%Get(i, DisNorm)
          DisNorm=Amp(i) .mx. DisNorm
          call this%BaseFlow%GetPart(i, i, 1, jn, BFFlux(:, :))
          call DisNorm%GetFlux(Flux)

          Amp_local(7)=0.0d0
          call BFFlux(1, 1)%get(rho0, u0, v0, w0, T0)
          call Flux(1)%get(rho, u, v, w, T)
          Amp_local(1)=abs(rho*u0+rho0+u)
          Amp_local(2)=abs(rho)
          Amp_local(3)=abs(u)
          Amp_local(4)=abs(v)
          Amp_local(5)=abs(T)
          Amp_local(6)=0.0d0
          Amp_local(7)=abs(w)
          E_local0=rho0*(u*conjg(u)+v*conjg(v))+w*conjg(w)
          Amp_pse(:, i)=0.0d0

          do j=2, jn
            call BFFlux(1, j)%get(rho0, u0, v0, w0, T0)
            call Flux(j)%get(rho, u, v, w, T)
            Amp_local(1)=abs(rho*u0+rho0+u)
            Amp_local(2)=abs(rho)
            Amp_local(3)=abs(u)
            Amp_local(4)=abs(v)
            Amp_local(5)=abs(T)
            E_local1=rho0*(u*conjg(u)+v*conjg(v))+w*conjg(w)
            Amp_local(6)=Amp_local(6)+0.5d0*(E_local1+E_local0)*(yy(j)-yy(j-1))
            E_local0=E_local1
            Amp_local(7)=abs(w)
            do l=1, 7
              if(amp_local(l)>=Amp_pse(l, i)) Amp_pse(l, i)=amp_local(l)
            enddo
          enddo
          amp_pse(6, i)=amp_local(6)
        end do

        open(993, file=fn_surf//'_Alf.plt', form='formatted')
        write(993, *)"variables=x, <greek>a</greek><sub>r</sub>, <greek>a</greek><sub>i</sub>"
        do i=This%Istart, this%Iend
          write(993, '(8ES20.8)')xx(i), real(Alf(i)), AIMAG(Alf(i))
        enddo
        close(993)

        open(993, file=fn_surf//'_Amp_pse.plt', form='formatted')
        write(993, *)"variables=x, A_rhou, A_rho, A_u, A_v, A_T, A_E, A_w"
        do i=This%Istart, this%Iend
          write(993, '(8ES20.8)')xx(i), (Amp_pse(l, i), l=1, 7)
        enddo
        close(993)

        deallocate(BFFlux)
        deallocate(Flux)

    end subroutine print_sigma


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

end module mod_npse
