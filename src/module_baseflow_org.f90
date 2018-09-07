!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: /media/psf/Home/Fortran/START/src/module_baseflow_org.f90
!> @file
!> @breif 基本流点通量文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: baseflow_org
!> @breif 基本流点通量模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-01 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-08-01
!------------------------------------------------------------------------------

module mod_baseflow_org !!BF=BaseFlow

    use mod_vector
    use penf, only: R_P
    implicit none

    private
    public :: BF_FLUX_NULL, BF_FLUXIJ_NULL !< 空常量

    !> 基本流通量\f$ \left( \rho, U, V, W, T\right) \f$类(不带序号).
    type, public :: bf_flux_org_type

        private
        real(R_P), private :: rho   !< 密度
        type(vector_type), private :: vel  !< 速度
        real(R_P), private :: T        !< 温度

        Contains

          Procedure :: Set => set_BF_flux_org !< 设置基本流通量
          procedure :: Zeros !< 置零函数
          generic :: Get => get_BF_flux_org, get_BF_flux_org_vel !< 获得基本流通量
!          generic  :: operator(*) => multiply_double_type, multiply_type_double
          generic  :: operator(+) => add_type                !< 通量加法

          procedure, pass(this), private :: get_BF_flux_org        !<获得基本流通量
          procedure, pass(this), private :: get_BF_flux_org_vel    !<获得基本流通量
          procedure, pass(this), private :: multiply_double_type   !<放缩基本流通量
          procedure, pass(this), private :: multiply_type_double   !<放缩基本流通量
          procedure, pass(obj1), private :: add_type              !<基本流通量相加

    end type bf_flux_org_type

    !> 基本流通量(带序号).
    type, public, extends(bf_flux_org_type) :: bf_flux_org_ij_type

        private
        integer, private :: i !<i方向序号
        integer, private :: j !<j方向序号

        Contains

          Procedure :: SetIJ => set_BF_flux_org_ij !<设置通量I和J
          procedure :: GetIJ => get_BF_flux_org_ij  !<获得通量I和J
          procedure :: Zeros => Zeros_ij !< 置零函数
          generic  :: operator(*) => multiply_double_typeij, &
                                  &   multiply_typeij_double   !< 通量乘法
          generic  :: operator(+) => add_typeij               !< 通量加法

          procedure, pass(this), private :: multiply_double_typeij   !< 通量放缩
          procedure, pass(this), private :: multiply_typeij_double    !<通量放缩
          procedure, pass(obj1), private :: add_typeij                !<通量相加

    end type bf_flux_org_ij_type

    type(bf_flux_org_type), parameter :: BF_FLUX_NULL= &
                                bf_flux_org_type(0.0d0, VEC_NULL, 0.0d0)
    type(bf_flux_org_ij_type), parameter :: BF_FLUXIJ_NULL= &
        bf_flux_org_ij_type(0.0d0, VEC_NULL, 0.0d0, 0, 0)

    contains

    !> 设置基本流通量.
    !! @param[in] rho 密度
    !! @param[in] vel 速度矢量
    !! @param[in] T   温度
    elemental subroutine set_BF_flux_org(this, rho, Vel, T)

        implicit none
        real(R_P), intent(in) :: rho, T
        type(vector_type), intent(in) :: Vel
        class(bf_flux_org_type), intent(inout) :: this

        this%rho=rho; this%Vel=Vel; this%T=T

    end subroutine set_BF_flux_org

    !> 获得基本流通量.
    !! @param[out] rho 密度
    !! @param[out] vel 速度矢量
    !! @param[out] T   温度
    elemental subroutine get_BF_flux_org(this, rho, Vel, T)

        implicit none
        class(bf_flux_org_type), intent(in) :: this
        real(R_P), intent(inout) :: rho, T
        type(vector_type), intent(inout) :: Vel

        rho=this%rho; Vel=this%Vel; T=this%T

    end subroutine get_BF_flux_org

    !> 获得基本流通量.
    !! @param[out] rho 密度
    !! @param[out] U 流向速度
    !! @param[out] V 法向速度
    !! @param[out] W 展向速度
    !! @param[out] T   温度
    subroutine get_BF_flux_org_vel(this, rho, U, V, W, T)

        implicit none
        class(bf_flux_org_type), intent(in) :: this
        real(R_P), intent(inout) :: rho, T, U, V, W

        rho=this%rho; T=this%T
        call this%Vel%Get(U, V, W)

    end subroutine get_BF_flux_org_vel

    !< 基本流通量乘标量(基本流通量放缩).
    !! @param[in] scalar 放缩系数
    !! @return 缩放之后的基本流通量
    type(bf_flux_org_type) function multiply_double_type(this, scalar)

        implicit none
        class(bf_flux_org_type), intent(in) :: this
        real(R_P), intent(in) :: scalar

        multiply_double_type%rho=scalar*this%rho
        multiply_double_type%Vel=scalar*this%Vel
        multiply_double_type%T  =scalar*this%T

    end function multiply_double_type

    !< 基本流通量乘标量(基本流通量放缩).
    !! @param[in] scalar 放缩系数
    !! @return 缩放之后的基本流通量
    type(bf_flux_org_type) function multiply_type_double(scalar, this)

        implicit none
        class(bf_flux_org_type), intent(in) :: this
        real(R_P), intent(in) :: scalar

        multiply_type_double%rho=scalar*this%rho
        multiply_type_double%Vel=scalar*this%Vel
        multiply_type_double%T  =scalar*this%T

    end function multiply_type_double

    !< 基本流通量乘标量(基本流通量放缩).
    !! @param[in] scalar 放缩系数
    !! @return 缩放之后的基本流通量
    type(bf_flux_org_ij_type) function multiply_double_typeij(this, scalar)

        implicit none
        class(bf_flux_org_ij_type), intent(in) :: this
        real(R_P), intent(in) :: scalar

        multiply_double_typeij%rho=scalar*this%rho
        multiply_double_typeij%Vel=scalar*this%Vel
        multiply_double_typeij%T  =scalar*this%T
        multiply_double_typeij%i  =this%i
        multiply_double_typeij%j  =this%j

    end function multiply_double_typeij

    !< 基本流通量乘标量(基本流通量放缩).
    !! @param[in] scalar 放缩系数
    !! @return 缩放之后的基本流通量
    type(bf_flux_org_ij_type) function multiply_typeij_double(scalar, this)

        implicit none
        class(bf_flux_org_ij_type), intent(in) :: this
        real(R_P), intent(in) :: scalar

        multiply_typeij_double%rho=scalar*this%rho
        multiply_typeij_double%Vel=scalar*this%Vel
        multiply_typeij_double%T  =scalar*this%T
        multiply_typeij_double%i  =this%i
        multiply_typeij_double%j  =this%j

    end function multiply_typeij_double

    !> 基本流通量相加 \f$U_3=U_1+U_2\f$.
    !! @param[in] obj1 通量\f$U_1\f$
    !! @param[in] obj2 通量\f$U_2\f$
    !! @return 通量\f$U_3\f$
    type(bf_flux_org_type) function add_type(obj1, obj2)

        implicit none
        class(bf_flux_org_type), intent(in) :: obj1
        type(bf_flux_org_type), intent(in) :: obj2

        add_type%rho=obj1%rho+obj2%rho
        add_type%Vel=obj1%Vel+obj2%Vel
        add_type%T  =obj1%T+obj2%T

    end function add_type

    !> 基本流通量相加 \f$U_3=U_1+U_2\f$.
    !! @param[in] obj1 通量\f$U_1\f$
    !! @param[in] obj2 通量\f$U_2\f$
    !! @return 通量\f$U_3\f$
    type(bf_flux_org_ij_type) function add_typeij(obj1, obj2)

        implicit none
        class(bf_flux_org_ij_type), intent(in) :: obj1
        type(bf_flux_org_ij_type), intent(in) :: obj2

        add_typeij%rho=obj1%rho+obj2%rho
        add_typeij%Vel=obj1%Vel+obj2%Vel
        add_typeij%T  =obj1%T+obj2%T
        add_typeij%i=obj1%i
        add_typeij%j=obj1%j

    end function add_typeij

    !> 设置基本流通量的序号
    !! @param[in] i i方向序号
    !! @param[in] j j方向序号
    subroutine set_BF_flux_org_ij(this, i, j)

        implicit none
        integer, intent(in) :: i, j
        class(bf_flux_org_ij_type), intent(inout) :: this

        this%i=i; this%j=j

    end subroutine set_BF_flux_org_ij

    !> 获得基本流通量的序号
    !! @param[out] i i方向序号
    !! @param[out] j j方向序号
    subroutine get_BF_flux_org_ij(this, i, j)

        implicit none
        integer, intent(inout) :: i, j
        class(bf_flux_org_ij_type), intent(in) :: this

        i=this%i; j=this%j

    end subroutine get_BF_flux_org_ij

    !> 置零函数.
    elemental subroutine Zeros(this)

        implicit none
        class(bf_flux_org_type), intent(inout) :: this

        this%rho=0.0d0
        call this%vel%zeros()
        this%T=0.0d0

    end subroutine Zeros

    !> 置零函数.
    elemental subroutine Zeros_ij(this)

        implicit none
        class(bf_flux_org_ij_type), intent(inout) :: this

        this%rho=0.0d0
        call this%vel%Zeros()
        this%T=0.0d0
        this%i=0; this%j=0

    end subroutine Zeros_ij


end module mod_baseflow_org
