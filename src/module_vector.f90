!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: /media/psf/Home/Fortran/START/src/module_vector.f90
!> @file
!> @breif 实型矢量类文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: vector
!> @breif 矢量模块(实型).
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-07-31 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-07-31
!------------------------------------------------------------------------------

module mod_vector

    use penf, only: R_P
    implicit none
    public  :: Vector, VEC_NULL, VECTOR100, VECTOR010, VECTOR001
    private

    real(R_P), parameter :: ONE=1.0d0, ZERO=0.0d0

    !> 矢量类(实型)
    type, public :: vector_type

        real(R_P), private :: a(3)                                              !<矢量分量

        Contains

          Procedure :: Set => set_vector !< 设置矢量
          Procedure :: Get => get_vector   !< 获得矢量
          procedure :: Zeros
          generic  :: operator(*) => multiply_scalar_type, multiply_type_scalar !<标量乘矢量
          generic  :: operator(+) => add, add_array, array_add                  !<矢量加法

          procedure, pass(this), private :: multiply_scalar_type !<标量乘矢量
          procedure, pass(this), private :: multiply_type_scalar !<矢量乘标量
          procedure, pass(obj1), private :: add                  !<矢量相加
          procedure, pass(obj1), private :: add_array            !<矢量相加
          procedure, pass(obj2), private :: array_add            !<矢量相加

          procedure :: finalize => finalize_Vector                 !<析构函数

    end type vector_type

    !>空矢量
    type(vector_type), parameter :: VEC_NULL=vector_type((/zero, zero, zero/))
    !>(1, 0, 0)
    type(vector_type), parameter :: VECTOR100=vector_type((/one, zero, zero/))
    !>(0, 1, 0)
    type(vector_type), parameter :: VECTOR010=vector_type((/zero, one, zero/))
    !>(0, 0, 1)
    type(vector_type), parameter :: VECTOR001=vector_type((/zero, zero, one/))

    contains

    !> 构造矢量.
    !!
    !! @param[in] a 第一分量
    !! @param[in] b 第二分量
    !! @param[in] c 第三分量
    !! @return 向量
    elemental type(vector_type) function Vector(a, b, c)

        implicit none
        real(R_P), intent(in) :: a, b, c
        Vector%a=(/a, b, c/)

    end function Vector

    !> 设置向量.
    !!
    !! @param[in] a 第一分量
    !! @param[in] b 第二分量
    !! @param[in] c 第三分量
    !! @return 向量(a, b, c)
    subroutine set_vector(this, a, b, c)

        implicit none
        real(R_P), intent(in) :: a, b, c
        class(vector_type), intent(inout) :: this

        this%a=(/a, b, c/)

    end subroutine set_vector

    !> 获得向量.
    !!
    !! @param[out] a 第一分量
    !! @param[out] b 第二分量
    !! @param[out] c 第三分量
    subroutine Get_vector(this, a, b, c)

        implicit none
        real(R_P), intent(inout) :: a, b, c
        class(vector_type), intent(in) :: this

        a=this%a(1)
        b=this%a(2)
        c=this%a(3)

    end subroutine get_vector


    !> 矢量加法\f$\mathbf c=\mathbf a+\mathbf b\f$.
    !!
    !! @param[in] obj1 矢量\f$ \mathbf a\f$
    !! @param[in] obj2 矢量\f$ \mathbf b\f$
    !! @retval add 和矢量\f$ \mathbf c\f$
    !! @return 矢量的和
    type(vector_type) function add(obj1, obj2)

        implicit none
        class(vector_type), intent(in) :: obj1, obj2

        add%a=obj1%a+obj2%a

    end function add

    !> 矢量加法\f$\mathbf c=\mathbf a+\mathbf b\f$.
    !!
    !! @param[in] obj1 矢量\f$ \mathbf a\f$
    !! @param[in] obj2 矢量\f$ \mathbf b\f$
    !! @retval add 和矢量\f$ \mathbf c\f$
    !! @return 矢量的和
    type(vector_type) function add_array(obj1, obj2)

        implicit none
        class(vector_type), intent(in) :: obj1
        real(R_P), intent(in) ::obj2(3)

        add_array%a=obj1%a+obj2

    end function add_array

    !> 矢量加法\f$\mathbf c=\mathbf a+\mathbf b\f$.
    !!
    !! @param[in] obj1 矢量\f$ \mathbf a\f$
    !! @param[in] obj2 矢量\f$ \mathbf b\f$
    !! @retval add 和矢量\f$ \mathbf c\f$
    !! @return 矢量的和
    type(vector_type) function array_add(obj1, obj2)

        implicit none
        real(R_P), intent(in) ::obj1(3)
        class(vector_type), intent(in) :: obj2

        array_add%a=obj2%a+obj1

    end function array_add

    !> 标量乘矢量: \f$ \mathbf c=\alpha \mathbf b\f$.
    !!
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return \f$ \mathbf c\f$
    type(vector_type) function multiply_scalar_type(this, scalar)

        implicit none
        class(vector_type), intent(in) :: this
        real(R_P), intent(in) :: scalar

        multiply_scalar_type%a=scalar*this%a

    end function multiply_scalar_type

    !> 标量乘矢量: \f$ \mathbf c=\alpha \mathbf b\f$.
    !!
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return \f$ \mathbf c\f$
    type(vector_type) function multiply_type_scalar(scalar, this)

        implicit none
        class(vector_type), intent(in) :: this
        real(R_P), intent(in) :: scalar

        multiply_type_scalar%a=scalar*this%a

    end function multiply_type_scalar

    !> 置零函数.
    elemental subroutine Zeros(this)

        implicit none
        class(vector_type), intent(inout) :: this

        this%a=0.0d0

    end subroutine Zeros

    !> 析构函数
    Subroutine finalize_Vector(this)

        Implicit None
        class(vector_type), intent(Inout) :: this

        this%a=0.0d0

    End Subroutine finalize_Vector

end module mod_vector
