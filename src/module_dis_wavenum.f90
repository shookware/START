!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_dis_wavenum.f90
!> @file
!> @breif 扰动色散参数相关文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: dis_wavenum
!> @breif 扰动色散参数相关模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-03 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-03
!------------------------------------------------------------------------------`
module mod_dis_wavenum

    use penf, only: R_P
    implicit none
    private

    !> 扰动色散关系类型
    type, public :: dis_wavenum_type

        complex(R_P), private :: alpha !< 流向波数\f$\alpha\f$
        complex(R_P), private :: beta !< 展向波数\f$\beta\f$
        complex(R_P), private :: omega !< 频率\f$\omega\f$

        Contains

          procedure :: SetAlpha !< 设置流向波数
          procedure :: SetOmega=> SetOmega_complex !< 设置频率
          procedure :: SetBeta=> SetBeta_complex !< 设置频率
          procedure :: getAlpha => get_alpha !< 获得流向波数
          procedure :: getBeta => get_beta !< 获得展向波数
          procedure :: getOmega => get_omega !< 获得频率
          procedure :: Print !< 输出波数信息
          procedure :: Write !< 写入波数信息
          procedure :: Read !< 读取波数信息

          generic :: Set => set_wavenum, set_alpha_beta_omega      !< 设置扰动色散信息
          generic :: Get => get_wavenum, get_alpha_beta_omega !< 获得扰动色散信息

          procedure, private :: set_alpha_beta_omega !< 设置波数和频率信息
          procedure, private :: get_alpha_beta_omega !< 获得波数和频率信息
          procedure, private :: set_wavenum !< 设置扰动色散信息
          procedure, private :: get_wavenum !<获得扰动色散信息

    end type dis_wavenum_type

    !> 适用于PSE的扰动色散信息
    !!
    !! 添加了一个字段表示流向波数随流向的导数
    type, extends(dis_wavenum_type), public :: dis_wavenum_lpse_type

        private
        complex(R_P), private :: Dx_alpha=0.0d0 !< 流向波数沿流向的导数

        contains
        procedure :: SetDxAlpha => set_Dx_alpha !< 设置流向波数沿流向的导数
        procedure :: GetDxAlpha => Get_Dx_Alpha !< 获得流向波数沿流向的导数

        generic :: operator(*) => mul_scalar_lpse_type !< 数乘(标量乘)
        generic :: operator(+) => add_wavenum_lpse   !< 加法操作
        procedure, pass(this), private :: mul_scalar_lpse_type !<标量乘类型
        procedure, pass(this), private :: add_wavenum_lpse !< 加法操作

    end type

    contains

    !> 设置流向波数
    !! @param[in] alpha 流向波数(复数)
    subroutine SetAlpha(this, alpha)

        implicit none
        class(dis_wavenum_type), intent(inout) :: this
        complex(R_P), intent(in) :: alpha

        this%alpha=alpha

    end subroutine SetAlpha

    !> 设置频率
    !! @param[in] omega 频率(复数)
    subroutine SetOmega_complex(this, omega)

        implicit none
        class(dis_wavenum_type), intent(inout) :: this
        complex(R_P), intent(in) :: omega

        this%omega=omega

    end subroutine SetOmega_complex

    !> 设置展向波数
    !! @param[in] beta 展向波数(复数)
    subroutine SetBeta_complex(this, beta)

        implicit none
        class(dis_wavenum_type), intent(inout) :: this
        complex(R_P), intent(in) :: beta

        this%beta=beta

    end subroutine SetBeta_complex

    !> 获得扰动色散关系
    !! @param[out] 扰动色散关系
    subroutine get_wavenum(this, wavenum)

        implicit none
        class(dis_wavenum_type), intent(in) :: this
        type(dis_wavenum_type), intent(out) :: wavenum

        wavenum%alpha = this%alpha
        wavenum%beta  = this%beta
        wavenum%omega = this%omega

    end subroutine get_wavenum

    !> 设置扰动的波数和频率信息
    !! @param[in] alpha 流向波数\f$\alpha\f$
    !! @param[in] beta 展向波数\f$\beta\f$
    !! @param[in] omega 频率\f$\omega\f$
    subroutine set_alpha_beta_omega(this, alpha, beta, omega)

        implicit none
        complex(R_P), intent(in) :: alpha, beta, omega
        class(dis_wavenum_type), intent(inout) :: this

        this%alpha=alpha; this%beta=beta; this%omega=omega

    end subroutine set_alpha_beta_omega

    !> 设置扰动色散信息
    !! @param[in] wavenum 扰动色散信息
    subroutine set_wavenum(this, wavenum)

        implicit none
        type(dis_wavenum_type), intent(in) :: wavenum
        class(dis_wavenum_type), intent(inout) :: this

        this%alpha=wavenum%alpha
        this%beta=wavenum%beta
        this%omega=wavenum%omega

    end subroutine set_wavenum

    !> 设置扰动流向波数沿流向的导数
    !! @param[in] Dx_alpha 流向波数沿流向的导数
    subroutine set_Dx_alpha(this, Dx_alpha)

        implicit none
        complex(R_P), intent(in) :: Dx_alpha
        class(dis_wavenum_lpse_type), intent(inout) :: this

        this%Dx_alpha=Dx_alpha

    end subroutine set_Dx_alpha

    !> 获得扰动流向波数沿流向的导数\f$ \partial{\alpha} / \partial{x}\f$
    !! @return 流向波数沿流向的导数
    function Get_Dx_Alpha(this) result(Dx_alpha)

        implicit none
        class(dis_wavenum_lpse_type), intent(in) :: this
        complex(R_P) :: Dx_alpha

        Dx_alpha=this%Dx_alpha

    end function Get_Dx_Alpha

    !> 获得扰动的波数和频率信息
    !! @param[out] alpha 流向波数\f$\alpha\f$
    !! @param[out] beta 展向波数\f$\beta\f$
    !! @param[out] omega 频率\f$\omega\f$
    subroutine get_alpha_beta_omega(this, alpha, beta, omega)

        implicit none
        class(dis_wavenum_type), intent(in) :: this
        complex(R_P), intent(inout) :: alpha, beta, omega

        alpha=this%alpha; beta=this%beta; omega=this%omega

    end subroutine get_alpha_beta_omega

    !> 获得流向波数\f$\alpha\f$
    !! @return 流向波数\f$\alpha\f$
    pure complex(R_P) function get_alpha(this)

        implicit none
        class(dis_wavenum_type), intent(in) :: this

        get_alpha=this%alpha

    end function get_alpha

    !> 获得展向波数\f$\beta\f$
    !! @return 展向波数\f$\beta\f$
    complex(R_P) function get_beta(this)

        implicit none
        class(dis_wavenum_type), intent(in) :: this

        get_beta=this%beta

    end function get_beta

    !> 获得频率\f$\omega\f$
    !! @return 频率\f$\omega\f$
    complex(R_P) function get_omega(this)

        implicit none
        class(dis_wavenum_type), intent(in) :: this

        get_omega=this%omega

    end function get_omega

    !> 输出波数信息
    subroutine Print(this)

        implicit none
        class(dis_wavenum_type),intent(inout) :: this

        write(*,"(1x, 'omega=', 2F20.8)") this%omega
        write(*,"(1x, 'alpha=', 2F20.8)") this%alpha
        write(*,"(1x, ' beta=', 2F20.8)") this%beta

    end subroutine Print

    !> 写入波数信息
    subroutine write(this, unit)
      implicit none
      class(dis_wavenum_type), intent(in) :: this
      integer, intent(in) :: unit

      write(unit)this%omega, this%alpha, this%beta

    end subroutine write

    !> 读取波数信息
    subroutine read(this, unit)
      implicit none
      class(dis_wavenum_type), intent(inout) :: this
      integer, intent(in) :: unit

      read(unit)this%omega, this%alpha, this%beta

    end subroutine read

    !> 标量乘
    function mul_scalar_lpse_type(scalar, this) result(mul)

        implicit none
        real(R_P), intent(in) :: scalar
        class(dis_wavenum_lpse_type), intent(in) :: this
        type(dis_wavenum_lpse_type) :: mul

        call mul%set(scalar*this%alpha, scalar*this%beta, scalar*this%omega)

    end function mul_scalar_lpse_type

    !> 类型相加
    function add_wavenum_lpse(this, obj) result(add)

        implicit none
        class(dis_wavenum_lpse_type), intent(in) :: this
        type(dis_wavenum_lpse_type), intent(in) :: obj
        type(dis_wavenum_lpse_type) :: add

        call add%set(this%alpha+obj%alpha, this%beta+obj%beta, this%omega+obj%omega)

    end function add_wavenum_lpse

end module mod_dis_wavenum
