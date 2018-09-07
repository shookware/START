!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_dis_normal.f90
!> @file
!> @breif 某流向展向站位处的扰动信息(形函数和波数)文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: dis_normal
!> @breif 某站位扰动信息(形函数和波数)模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-03 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-03
!------------------------------------------------------------------------------
module mod_dis_normal

    use mod_dis_shape
    use mod_dis_wavenum
    use mod_dis_flux
    use penf, only: R_P

    implicit none

    private

    !> 局部扰动类
    type, public, extends(dis_shape_type) :: dis_normal_type

        private
        type(dis_wavenum_lpse_type), private :: wavenum !< 扰动色散信息
!        type(dis_shape_diff_type), private :: DiffShape !< 扰动形函数的导数

        Contains

          Procedure :: SetWavenum => set_wavenum !< 设置扰动色散信息
          procedure :: GetWaveNum !< 获得扰动色散信息
          procedure :: SetAlpha => set_alpha !< 设置流向波数
          procedure :: SetOmega => set_omega !< 设置频率
          procedure :: SetBeta => set_beta !< 设置展向波数
          procedure :: GetAlpha => get_alpha !< 读取流向波数
          procedure :: GetOmega => get_omega !< 读取频率
          procedure :: GetBeta => get_beta !< 读取频率
          !procedure :: UpdateDiff !< 更新扰动形函数导数
          !procedure :: SetDiffIloc !< 设置扰动形函数导数的流向站位
          !procedure :: CreateDiff !< 创建扰动形函数导数, 分配内存
          !procedure :: GetDiffShape !<获得扰动形函数的导数信息
!          Final     ::

    end type dis_normal_type

    contains

    !> 设置扰动色散信息
    !! @param[in] wavenum 扰动色散信息
    subroutine set_wavenum(this, wavenum)

        implicit none
        type(dis_wavenum_lpse_type), intent(in) :: wavenum
        class(dis_normal_type), intent(inout) :: this

        this%wavenum=wavenum

    end subroutine set_wavenum


        !> 获得扰动色散信息
        !! @return 扰动色散信息
    function GetWaveNum(this) result(WaveNum)

        implicit none
        class(dis_normal_type), intent(in) :: this
        type(dis_wavenum_lpse_type) :: wavenum

        wavenum=this%wavenum

    end function GetWaveNum

    !> 获得扰动流向波数
  elemental function Get_alpha(this) result(alpha)

      implicit none
      class(dis_normal_type), intent(in) :: this
      complex(R_P) :: alpha
      type(dis_wavenum_lpse_type) :: wavenum

      wavenum=this%wavenum
      alpha=wavenum%getAlpha()

  end function Get_alpha

  !> 获得扰动频率
  function Get_omega(this) result(omega)

      implicit none
      class(dis_normal_type), intent(in) :: this
      complex(R_P) :: omega
      type(dis_wavenum_lpse_type) :: wavenum

      wavenum=this%wavenum
      omega=wavenum%getOmega()

  end function Get_omega

  !> 获得扰动展向波数
  function Get_beta(this) result(beta)

      implicit none
      class(dis_normal_type), intent(in) :: this
      complex(R_P) :: beta
      type(dis_wavenum_lpse_type) :: wavenum

      wavenum=this%wavenum
      beta=wavenum%getBeta()

  end function Get_beta

    !> 设置流向波数.
    subroutine set_alpha(this, alpha)

        use penf, only: R_P
        implicit none
        class(dis_normal_type), intent(inout) :: this
        complex(R_P), intent(in) :: alpha !< 流向波数

        call this%wavenum%setalpha(alpha)

    end subroutine set_alpha

    !> 设置频率.
    subroutine set_omega(this, omega)

        implicit none
        class(dis_normal_type), intent(inout) :: this
        complex(R_P), intent(in) :: omega !< 频率

        call this%wavenum%SetOmega(omega)

    end subroutine set_omega

    !> 设置展向波数.
    subroutine set_beta(this, beta)

        implicit none
        class(dis_normal_type), intent(inout) :: this
        complex(R_P), intent(in) :: beta !< 频率

        call this%wavenum%SetBeta(beta)

    end subroutine set_beta

    !!> 创建扰动形函数导数, 分配内存
    !!! @param[in] jn 法向点数
    !elemental subroutine CreateDiff(this, jn)
    !
    !    implicit none
    !    class(dis_normal_type), intent(inout) :: this
    !    integer, intent(in) :: jn
    !
    !    call this%DiffShape%Create(jn)
    !
    !end subroutine CreateDiff

    !!> 设置扰动形函数导数的流向站位
    !!! @param[in] iloc 流向站位
    !subroutine SetDiffIloc(this, iloc)
    !
    !    implicit none
    !    integer, intent(in) :: iloc
    !    class(dis_normal_type), intent(inout) :: this
    !
    !    call this%DiffShape%SetILoc(iloc)
    !
    !end subroutine SetDiffIloc

    !!> 更新扰动形函数导数
    !!! @param[in] Diff 扰动的场差分基函数
    !subroutine UpdateDiff(this, Diff)
    !
    !    use mod_difference
    !    implicit none
    !    class(dis_normal_type), intent(inout) :: this
    !    type(difference_2D_type), intent(in) :: Diff
    !
    !    call this%DiffShape%Set(this, Diff)
    !
    !end subroutine UpdateDiff

    !!> 获得扰动形函数的导数信息
    !!! @param[out] DiffShape 扰动形函数的导数
    !subroutine GetDiffShape(this, DiffShape)
    !
    !    implicit none
    !    class(dis_normal_type), intent(in) :: this
    !    type(dis_shape_diff_type), intent(inout) :: DiffShape
    !
    !    DiffShape=this%DiffShape
    !
    !end subroutine GetDiffShape

end module mod_dis_normal
