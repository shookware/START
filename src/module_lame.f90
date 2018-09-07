!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: /media/psf/Home/Fortran/START/src/module_lame.f90 
!> @file
!> @breif Lame系数相关文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: Lame
!> @breif Lame系数模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-01 - Initial Version
!> @author
!> Liu Jianxin
!> @date 2017-08-01
!------------------------------------------------------------------------------
module mod_lame

    use penf, only: R_P
    use mod_vector
    implicit none

    private

    !> Lame系数类
    type, public :: lame_type

        private
        type(vector_type) :: lame    !< Lame系数(\f$ h_1, h_2, h_3 \f$), 对应于流向法向展向三个方向

        Contains

          Procedure :: Set => set_lame  !< 设置Lame系数
          procedure :: Get => get_lame  !< 获得Lame系数
          procedure :: finalize => finalize_lame   !< 析构函数

    end type lame_type

    !> Lame系数随坐标梯度
    type, public :: lame_grad_type

        private
        real(R_P), private :: d12   !< \f$ \frac 1 {h_1} \frac {\partial {h_1}} {\partial{h_2}} \f$
        real(R_P), private :: d32   !< \f$ \frac 1 {h_3} \frac {\partial {h_3}} {\partial{h_2}} \f$
        real(R_P), private :: d31   !< \f$ \frac 1 {h_{1}h_{3}}\frac {\partial {h_3}} {\partial{h_1}} \f$

        Contains

          Procedure :: Set => Set_lame_grad   !< 设置Lame系数梯度
          procedure :: Get => Get_lame_grad   !< 获得Lame系数梯度

    end type lame_grad_type

    contains

    !> 设置Lame系数.
    !!
    !! @param[in] hx \f$ h_1 \f$
    !! @param[in] hy \f$ h_2 \f$
    !! @param[in] hz \f$ h_3 \f$
    subroutine set_lame(this, hx, hy, hz)

        implicit none
        real(R_P)          , intent(in) :: hx, hy, hz
        class(lame_type), intent(inout) :: this

        this%lame=Vector(hx, hy, hz)

    end subroutine set_lame

    !> 设置Lame系数梯度.
    !! @param[in] d12 \f$ \frac 1 {h_1} \frac {\partial {h_1}} {\partial{h_2}} \f$
    !! @param[in] d32 \f$ \frac 1 {h_3} \frac {\partial {h_3}} {\partial{h_2}} \f$
    !! @param[in] d31 \f$ \frac 1 {h_{1}h_{3}}\frac {\partial {h_3}} {\partial{h_1}} \f$
    subroutine set_lame_grad(this, d12, d32, d31)

        implicit none
        real(R_P)          , intent(in) :: d12, d32, d31
        class(lame_grad_type), intent(inout) :: this

        this%d12=d12; this%d32=d32; this%d31=d31

    end subroutine set_lame_grad

    !> 获得Lame系数.
    !!
    !! @param[out] hx \f$ h_1 \f$
    !! @param[out] hy \f$ h_2 \f$
    !! @param[out] hz \f$ h_3 \f$
    subroutine get_lame(this, hx, hy, hz)

        implicit none
        class(lame_type), intent(in) :: this
        real(R_P)          , intent(inout) :: hx, hy, hz

        call this%lame%Get(hx, hy, hz)

    end subroutine get_lame

    !> 获得Lame系数梯度.
    !! @param[out] d12 \f$ \frac 1 {h_1} \frac {\partial {h_1}} {\partial{h_2}} \f$
    !! @param[out] d32 \f$ \frac 1 {h_3} \frac {\partial {h_3}} {\partial{h_2}} \f$
    !! @param[out] d31 \f$ \frac 1 {h_{1}h_{3}}\frac {\partial {h_3}} {\partial{h_1}} \f$
    subroutine get_lame_grad(this, d12, d32, d31)

        implicit none
        class(lame_grad_type), intent(in) :: this
        real(R_P)          , intent(inout) :: d12, d32, d31

        d12=this%d12; d32=this%d32; d31=this%d31

    end subroutine get_lame_grad

    !> 析构函数
    Subroutine finalize_lame(this)

        Implicit None
        class(lame_type), intent(Inout) :: this

        call this%lame%Set(1.0d0, 1.0d0, 1.0d0)

    End Subroutine finalize_lame

end module mod_lame
