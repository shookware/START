!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_lpse_dis_normal.f90
!> @file
!> @breif pse某站位处的扰动文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: lpse_dis_normal
!> @breif 某站位pse扰动模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-04 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-04
!------------------------------------------------------------------------------
module mod_lpse_dis_normal

    use mod_dis_normal
    use mod_lpse_dis_OP_normal
    use mod_lpse_dis_OP_point
    use mod_dis_flux
    use penf, only: R_P
    use stringifor

    implicit none

    private
    public :: lpse_dis_normal_type

    !> PSE在一个站位的扰动类
    type, extends(dis_normal_type) :: lpse_dis_normal_type

        private
        complex(R_P), private :: sigma(7) !< PSE的物理复波数(流向波数和增长率),分别对应于动量,密度,流向速度,法向速度,温度,动能以及展向速度

        Contains

        generic :: operator(.add.) => addDis !< 扰动相加
        generic :: operator(.mx.) => mul_scalardble_type !< 扰动乘法     !mul_disop,
        generic :: SetSigma => set_sigma_multi, set_sigma_single !<设置物理复波数
        generic :: GetSigma => get_sigma_multi, get_sigma_single !<获得物理复波数

        procedure :: Write=>WriteDisNorm !<存储数据
        procedure :: Read=>ReadDisNorm !<读取数据

!        procedure, pass(this), private :: mul_disop !<扰动算子与其相乘
        procedure, pass(this), private :: mul_scalardble_type !<标量与其相乘
        procedure, pass(this), private :: addDis !< 扰动加法
        procedure, pass(this), private :: set_sigma_single !< 设置单个的物理复波数
        procedure, pass(this), private :: set_sigma_multi !<设置全部的物理复波数
        procedure, pass(this), private :: get_sigma_single !<获得单个的物理复波数
        procedure, pass(this), private :: get_sigma_multi !< 获得全部的物理复波数

    end type lpse_dis_normal_type

    contains

     !>读取LPSE法向扰动分布
     subroutine ReadDisNorm(this, fn_surf)

       use mod_dis_shape
       use mod_dis_wavenum
       implicit none
       class(lpse_dis_normal_type), intent(inout) :: this
       type(string), intent(in) :: fn_surf
       integer :: ios
       type(dis_wavenum_lpse_type) :: wavenum
       type(dis_shape_type) :: DisShape

       open(unit=999, file=fn_surf//'_DisNorm.dat', iostat=ios, &
            form='unformatted', access='stream')
       if ( ios /= 0 ) then
         print*, "Error opening file ", fn_surf//'_DisNorm.dat'
         stop
       endif
       read(999)this%sigma
       call wavenum%Read(999)

       call DisShape%ReadDisShape(999)
       call this%SetDisShape(DisShape)
       call this%SetWavenum(Wavenum)

     end subroutine ReadDisNorm

     !>写入LPSE法向扰动分布
     subroutine WriteDisNorm(this, fn_surf)

       use mod_dis_shape
       use mod_dis_wavenum
       implicit none
       class(lpse_dis_normal_type), intent(in) :: this
       type(string), intent(in) :: fn_surf
       integer :: ios
       type(dis_wavenum_lpse_type) :: wavenum
       type(dis_shape_type) :: DisShape

       open(unit=999, file=fn_surf//'_DisNorm.dat', iostat=ios, &
            form='unformatted', access='stream')
       if ( ios /= 0 ) then
         print*, "Error opening file ", fn_surf//'_DisNorm.dat'
         stop
       endif

       write(999)this%sigma
       wavenum=this%GetWaveNum()
       DisShape=this%GetDisShape()

       call Wavenum%Write(999)
       call DisShape%WriteDisShape(999)
       close(999)

     end subroutine WriteDisNorm

     !> 获得单个的物理复波数
     subroutine get_sigma_single(this, char_sign, sigma)

        implicit none
        class(lpse_dis_normal_type), intent(in) :: this
        character, intent(in) :: char_sign !< 波数类型
        complex(R_P), intent(out) :: sigma !< 物理复波数

        select case (char_sign)
            case ("m")
              sigma=this%sigma(1)
            case ("r")
              sigma=this%sigma(2)
            case ("u")
              sigma=this%sigma(3)
            case ("v")
              sigma=this%sigma(4)
            case ("T")
              sigma=this%sigma(5)
            case ("E")
              sigma=this%sigma(6)
            case ("w")
              sigma=this%sigma(7)
            case default
              Write(*, *)"Please Input a correct flag for defining &
              &           the growth rate!"
        end select

    end subroutine get_sigma_single

    !> 获得全部的物理复波数
    subroutine get_sigma_multi(this, sigma)

        implicit none
        class(lpse_dis_normal_type), intent(in) :: this
        complex(R_P), intent(out) :: sigma(7) !< 物理复波数

        sigma=this%sigma

    end subroutine get_sigma_multi

    !> 设置单个的物理复波数
    subroutine set_sigma_single(this, char_sign, sigma)

        implicit none
        class(lpse_dis_normal_type), intent(inout) :: this
        character, intent(in) :: char_sign !< 波数类型
        complex(R_P), intent(in) :: sigma !<物理复波数

        select case (char_sign)
            case ("m")
              this%sigma(1)=sigma
            case ("r")
              this%sigma(2)=sigma
            case ("u")
              this%sigma(3)=sigma
            case ("v")
              this%sigma(4)=sigma
            case ("T")
              this%sigma(5)=sigma
            case ("E")
              this%sigma(6)=sigma
            case ("w")
              this%sigma(7)=sigma
            case default
              Write(*, *)"Please Input a correct flag for defining the growth rate!"
        end select

    end subroutine set_sigma_single

    !> 设置全部的复波数
    subroutine set_sigma_multi(this, sigma)

        implicit none
        class(lpse_dis_normal_type), intent(inout) :: this
        complex(R_P), intent(in) :: sigma(7) !< 物理复波数

        this%sigma=sigma

    end subroutine set_sigma_multi

!    !> PSE扰动的算子与扰动相乘\f$\varphi=L\phi\f$, 波数不做变化
!    function mul_disop(disop, this) result(DisNorm)
!
!        use mod_difference
!        use mod_dis_shape
!        implicit none
!        type(lpse_dis_op_normal_type), target, intent(in) :: DisOP !< PSE扰动算子\f$L\f$
!        class(lpse_dis_normal_type), intent(in) :: this
!        type(lpse_dis_normal_type) :: DisNorm !< 乘积\f$\varphi\f$
!        type(lpse_dis_op_point_type) :: DisOPPoint
!        type(dis_flux_ij_type) :: FluxPoint, DyFlux, DyyFlux
!        type(dis_shape_diff_type) :: DiffShape
!        type(dis_shape_type) :: DisShape
!        complex(R_P), dimension(5, 5) :: A, B, D, Vyy
!!        type(difference_2D_type) :: Diff
!        real(R_P), pointer :: Coef(:)
!        integer, parameter :: XI1=1
!
!        integer :: Jn
!        integer :: iloc, j
!
!        Jn= this%getJn()
!        DisNorm=this
!        Coef=> disop%Coef
!
!        iloc= this%GetILoc()
!        call this%GetDiffShape(DiffShape)
!        do j=1, Jn
!            DisOPPoint=DisOP%GetPoint(j)
!            call DisOPPoint%get(A, B, D, Vyy)
!            call this%GetFluxPoint(j, FluxPoint)
!            DyFlux=DiffShape%GetPartialEta(j)
!            DyyFlux=DiffShape%GetPartialEta2(j)
!            FluxPoint=((A*Coef(size(Coef))) .mx. FluxPoint)+(B .mx. DyFlux)+(D .mx. FluxPoint)-(Vyy .mx. DyyFlux)
!            call DisNorm%SetFluxPoint(J, FluxPoint)
!        enddo
!
!    end function mul_disop

    !> 标量乘:\f$\varphi=\alpha\phi\f$
    function mul_scalardble_type(scalar, this) result(mul)

        use mod_dis_wavenum
        use mod_dis_shape
        implicit none
        class(lpse_dis_normal_type), intent(in) :: this
        real(R_P), intent(in) :: scalar !< 标量\f$\alpha\f$
        type(lpse_dis_normal_type) :: mul !< 乘积\f$\varphi\f$

        type(dis_wavenum_lpse_type) :: wavenum
        type(dis_shape_type) :: DisShape

        call mul%Create(this%getJn())
!        mul=this
        wavenum=this%GetWaveNum()
        call mul%SetWavenum(scalar*wavenum)
       !! DisShape=scalar*this%GetDisShape()
        call mul%SetDisShape(scalar*this%GetDisShape())

    end function mul_scalardble_type

    !> 扰动相加
    function addDis(this, obj) result(add)

        use mod_dis_wavenum
        use mod_dis_shape
        implicit none
        class(lpse_dis_normal_type), intent(in) :: this
        type(lpse_dis_normal_type), intent(in) :: obj
        type(lpse_dis_normal_type) :: add

        type(dis_wavenum_lpse_type) :: wavenum1, wavenum2, wavenum
        type(dis_shape_type) :: DisShape

        wavenum1=this%GetWaveNum(); wavenum2= obj%GetWaveNum()

        DisShape=this%GetDisShape()+obj%GetDisShape()
        call add%create(this%getJn())
        call add%SetDisShape(DisShape)
        wavenum=wavenum1+wavenum2
        call add%SetWavenum(wavenum)

    end function addDis

!    Subroutine finalize_lpse_dis_normal(this)

!        Implicit None
!        Type(lpse_dis_normal_type), intent(Inout) :: this


!    End Subroutine finalize_lpse_dis_normal

end module mod_lpse_dis_normal
