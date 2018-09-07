!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: vector
!> @breif 矢量模块(复型).
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
module mod_vector_cmplx

    use penf, only: R_P
    implicit none
    public :: VECTOR_CMPLX, VEC_CMPLX_NULL
    private

    complex(R_P), parameter :: ZERO=(0.0d0, 0.0d0)

    !> 矢量类(复型).
    type, public :: vector_cmplx_type

        complex(R_P), private :: a(3)     !<矢量分量

    Contains

        Procedure :: Set            => set_vector_cmplx !< 设置矢量
        Procedure :: Get            => get_vector_cmplx  !< 获得矢量
        procedure :: SelfCmplxDot !<自内点乘
        procedure :: Cojg         !<共轭
        procedure :: Zeros    !< 置零函数

        generic  :: operator(*)     => multiply_double_type, &
            &  multiply_type_double, &
            &  multiply_cmplx_type, &
            &  multiply_type_cmplx, &
            &  multiply_cmplxarray_type, &
            &  multiply_dblearray_type                     !<点乘法
        generic  :: operator(+)     => add, add_array, array_add   !<矢量加法
        generic  :: operator(-)     => minus !<矢量减法
        generic  :: operator(.dot.) => dotmul  !<点乘

        procedure, pass(this), private :: multiply_double_type!<标量乘矢量
        procedure, pass(this), private :: multiply_type_double!<标量乘矢量
        procedure, pass(this), private :: multiply_cmplx_type !<标量乘矢量
        procedure, pass(this), private :: multiply_type_cmplx!<标量乘矢量
        procedure, pass(this), private :: multiply_cmplxarray_type !<矢量点乘
        procedure, pass(this), private :: multiply_dblearray_type  !<矢量点乘

        procedure, pass(obj1), private :: add                  !<矢量相加
        procedure, pass(obj1), private :: add_array            !<矢量相加
        procedure, pass(obj2), private :: array_add            !<矢量相加

        procedure, pass(obj1), private :: minus
        procedure, pass(this), private :: dotmul

        procedure    :: finalize    => finalize_vector_cmplx

    end type vector_cmplx_type

    !>空矢量
    type(vector_cmplx_type), parameter :: VEC_CMPLX_NULL=vector_cmplx_type((/ZERO, ZERO, ZERO/))

contains

    !> 矢量点乘\f$ c = \mathbf{a} \cdot \mathbf{b}\f$.
    !!
    !! @param[in] obj 被乘矢量 \f$\mathbf b\f$
    !! @retval Dot 矢量\f$\mathbf a\f$与被乘矢量\f$\mathbf b\f$的点积
    !! @return 点乘值 *c*
    function dotmul(this, obj) result(Dot)

        implicit none
        class(vector_cmplx_type), intent(in) :: this
        type(vector_cmplx_type), intent(in) :: obj
        complex(R_P) :: Dot

        Dot=dot_product(this%a, obj%a)

    end function dotmul

    !> 矢量求共轭.
    !!
    !! @return 共轭矢量

    function cojg(this) result(ReValue)

        implicit none
        class(vector_cmplx_type), intent(in) :: this
        type(vector_cmplx_type) :: Revalue

        Revalue%a=conjg(this%a)

    end function cojg

    !> 矢量自身内积.
    !!
    !! @return 矢量自内积
    function SelfCmplxDot(this) result(Norm)

        implicit none
        class(vector_cmplx_type), intent(in) :: this
        complex(R_P) :: Norm

        Norm=sum(this%a*conjg(this%a))

    end function SelfCmplxDot

    !> 构造矢量.
    !!
    !! @param[in] a 第一分量
    !! @param[in] b 第二分量
    !! @param[in] c 第三分量
    !! @return 向量
    type(vector_cmplx_type) function vector_cmplx(a, b, c)

        implicit none
        real(R_P), intent(in) :: a, b, c
        Vector_cmplx%a=(/a, b, c/)

    end function vector_cmplx

    !> 设置向量.
    !!
    !! @param[in] a 第一分量
    !! @param[in] b 第二分量
    !! @param[in] c 第三分量
    !! @return 向量(a, b, c)
    subroutine set_vector_cmplx(this, a, b, c)

        implicit none
        complex(R_P), intent(in) :: a, b, c
        class(vector_cmplx_type), intent(inout) :: this

        this%a=(/a, b, c/)

    end subroutine set_vector_cmplx

    !> 获得向量.
    !!
    !! @param[out] a 第一分量
    !! @param[out] b 第二分量
    !! @param[out] c 第三分量
    pure  subroutine Get_vector_cmplx(this, a, b, c)

        implicit none
        complex(R_P), intent(inout) :: a, b, c
        class(vector_cmplx_type), intent(in) :: this

        a=this%a(1); b=this%a(2); c=this%a(3)

    end subroutine get_vector_cmplx

    !> 矢量加法\f$\mathbf c=\mathbf a+\mathbf b\f$.
    !!
    !! @param[in] obj1 矢量\f$ \mathbf a\f$
    !! @param[in] obj2 矢量\f$ \mathbf b\f$
    !! @retval add 和矢量\f$ \mathbf c\f$
    !! @return 矢量的和
    type(vector_cmplx_type) function add(obj1, obj2)

        implicit none
        class(vector_cmplx_type), intent(in) :: obj1, obj2

        add%a=obj1%a+obj2%a

    end function add

    !> 矢量加法\f$\mathbf c=\mathbf a+\mathbf b\f$.
    !!
    !! @param[in] obj1 矢量\f$ \mathbf a\f$
    !! @param[in] obj2 矢量\f$ \mathbf b\f$
    !! @retval add 和矢量\f$ \mathbf c\f$
    !! @return 矢量的和
    type(vector_cmplx_type) function add_array(obj1, obj2)

        implicit none
        class(vector_cmplx_type), intent(in) :: obj1
        complex(R_P), intent(in) ::obj2(3)

        add_array%a=obj1%a+obj2

    end function add_array

    !> 矢量加法\f$\mathbf c=\mathbf a+\mathbf b\f$.
    !!
    !! @param[in] obj1 矢量\f$ \mathbf a\f$
    !! @param[in] obj2 矢量\f$ \mathbf b\f$
    !! @retval add 和矢量\f$ \mathbf c\f$
    !! @return 矢量的和
    type(vector_cmplx_type) function array_add(obj1, obj2)

        implicit none
        complex(R_P), intent(in) ::obj1(3)
        class(vector_cmplx_type), intent(in) :: obj2

        array_add%a=obj2%a+obj1

    end function array_add


    !> 矢量减法\f$\mathbf c=\mathbf a-\mathbf b\f$.
    !!
    !! @param[in] obj1 矢量\f$ \mathbf a\f$
    !! @param[in] obj2 矢量\f$ \mathbf b\f$
    !! @retval add 差矢量\f$ \mathbf c\f$
    !! @return 矢量的差
    pure function minus(obj1, obj2) result(dis)

        implicit none
        class(vector_cmplx_type), intent(in) :: obj1, obj2
        type(vector_cmplx_type) :: dis

        dis%a=obj1%a-obj2%a

    end function minus

    !> 标量乘矢量: \f$ \mathbf c=\alpha \mathbf b\f$.
    !!
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return \f$ \mathbf c\f$
    pure function multiply_double_type(this, scalar) result(mul)

        implicit none
        class(vector_cmplx_type), intent(in) :: this
        real(R_P), intent(in) :: scalar
        type(vector_cmplx_type) :: mul

        mul%a=scalar*this%a

    end function multiply_double_type

    !> 标量乘矢量: \f$ \mathbf c=\alpha \mathbf b\f$.
    !!
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return \f$ \mathbf c\f$
    pure function multiply_type_double(scalar, this) result(mul)

        implicit none
        class(vector_cmplx_type), intent(in) :: this
        real(R_P), intent(in) :: scalar
        type(vector_cmplx_type) :: mul

        mul%a=scalar*this%a

    end function multiply_type_double

    !> 标量乘矢量: \f$ \mathbf c=\alpha \mathbf b\f$.
    !!
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return \f$ \mathbf c\f$
    pure function multiply_cmplx_type(this, scalar) result(mul)

        implicit none
        class(vector_cmplx_type), intent(in) :: this
        complex(R_P), intent(in) :: scalar
        type(vector_cmplx_type) :: mul

        mul%a=scalar*this%a

    end function multiply_cmplx_type

    !> 标量乘矢量: \f$ \mathbf c=\alpha \mathbf b\f$.
    !!
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return \f$ \mathbf c\f$
    pure function multiply_type_cmplx(scalar, this) result(mul)

        implicit none
        class(vector_cmplx_type), intent(in) :: this
        complex(R_P), intent(in) :: scalar
        type(vector_cmplx_type) :: mul

        mul%a=scalar*this%a

    end function multiply_type_cmplx

    !> 矢量点乘: \f$ c=\mathbf a \cdot \mathbf b\f$.
    !!
    !! @param[in] scalar 标量\f$\mathbf a\f$
    !! @return \f$ c\f$
    pure function multiply_dblearray_type(array, this) result(dotmult)

        implicit none
        class(vector_cmplx_type), intent(in) :: this
        real(R_P), dimension(3), intent(in) :: array
        complex(R_P) :: dotmult

        dotmult=sum(array*this%a)

    end function multiply_dblearray_type

    !> 矢量点乘: \f$ c=\mathbf a \cdot \mathbf b\f$.
    !!
    !! @param[in] scalar 标量\f$\mathbf a\f$
    !! @return \f$ c\f$
    pure function multiply_cmplxarray_type(array, this) result(dotmult)

        implicit none
        class(vector_cmplx_type), intent(in) :: this
        complex(R_P), dimension(3), intent(in) :: array
        complex(R_P) :: dotmult

        dotmult=sum(array*this%a)

    end function multiply_cmplxarray_type

    !> 置零函数.
    elemental subroutine Zeros(this)

        implicit none
        class(vector_cmplx_type), intent(inout) :: this

        this%a=0.0d0

    end subroutine Zeros

    !> 析构函数
    Subroutine finalize_vector_cmplx(this)

        Implicit None
        class(vector_cmplx_type), intent(Inout) :: this

        this%a=zero

    End Subroutine finalize_vector_cmplx

end module mod_vector_cmplx
