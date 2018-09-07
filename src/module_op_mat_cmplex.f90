!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: mod_op_mat_cmplx.f90
!> @file
!> @breif 矩阵操作模块文件_复型.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: op_mat
!> @breif 矩阵操作模块(复型).
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-03 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-03
!------------------------------------------------------------------------------
module mod_op_mat_cmplx

    use penf, only: R_P
    implicit none
    private

    !> 矩阵操作类
    type, public :: op_mat_cmplx_type

        private

        complex(R_P), private :: Mat(5, 5) !< 复型的5*5的矩阵


        Contains

          Procedure, pass(this) :: Set => Set !< 设置矩阵的值
          procedure, pass(this) :: Get => get !< 获得矩阵的值
          procedure, pass(this) :: trans => trans !< 实型矩阵变换成复型
          procedure, pass(this) :: Print => PrintMat !< 输出矩阵
          generic :: assignment(=) => set, get, Trans !<定义赋值操作
          generic :: operator(+) => add !< 加法
          generic :: operator(-) => Subtract, minus !< 减法
          generic :: operator(*) => multCmplx1, multCmplx2, MultReal1, MultReal2 !<矩阵标量乘
          generic :: operator(.t.) => transpose_op_mat !<转置

          procedure, pass(this), private :: MultCmplx1 => mult_type_cmplx !<矩阵标量乘
          procedure, pass(this), private :: MultCmplx2 => mult_cmplx_type !<矩阵标量乘
          procedure, pass(this), private :: MultReal1 => mult_real_type !<矩阵标量乘
          procedure, pass(this), private :: MultReal2 => mult_type_real !<矩阵标量乘
          procedure, pass(obj1), private :: Subtract => subtract !<矩阵减法
          procedure, pass(obj1), private :: Add => Add !<矩阵加法
          procedure, pass(this), private :: transpose_op_mat !<矩阵转置
          procedure, pass(this), private :: minus !<求相反数

    end type op_mat_cmplx_type

    contains

    !< 设置矩阵的值
    !! @param[in] Mat 5*5的矩阵
    subroutine Set(this, Mat)

        implicit none
        complex(R_P), intent(in) :: Mat(5, 5)
        class(op_mat_cmplx_type), intent(inout) :: this

        this%Mat=Mat

    end subroutine Set

    !< 获得矩阵的值
    !! @param[out] Mat 5*5的矩阵
    subroutine Get( Mat, this)

        implicit none
        complex(R_P), intent(inout) :: Mat(5, 5)
        class(op_mat_cmplx_type), intent(in) :: this

        Mat=this%Mat

    end subroutine Get

    !> 实型矩阵转换乘复型矩阵
    !! @param[in] obj 实型矩阵
    subroutine Trans(this, obj)

        use mod_op_mat
        implicit none
        type(op_mat_type), intent(in) :: obj
        class(op_mat_cmplx_type), intent(inout) :: this
        real(R_P) :: tmp(5, 5)

        call obj%Get(tmp)
        this%Mat=tmp

    end subroutine

    !> 矩阵标量乘
    !! @param[in] scalar 标量
    !! @return 标量乘的结果
    type(op_mat_cmplx_type) function mult_type_cmplx(this, scalar)

        implicit none
        class(op_mat_cmplx_type), intent(in) :: this
        complex(R_P) , intent(in) :: scalar

        mult_type_cmplx%mat=this%Mat*scalar

    end function mult_type_cmplx

    !> 矩阵标量乘
    !! @param[in] scalar 标量
    !! @return 标量乘的结果
    type(op_mat_cmplx_type) function mult_cmplx_type(scalar, this)

        implicit none
        class(op_mat_cmplx_type), intent(in) :: this
        complex(R_P) , intent(in) :: scalar

        mult_cmplx_type%mat=this%mat*scalar

    end function mult_cmplx_type

    !> 矩阵标量乘
    !! @param[in] scalar 标量
    !! @return 标量乘的结果
    type(op_mat_cmplx_type) function mult_type_real(this, scalar)

        implicit none
        class(op_mat_cmplx_type), intent(in) :: this
        real(R_P) , intent(in) :: scalar

        mult_type_real%mat=this%mat*scalar

    end function mult_type_real

    !> 矩阵标量乘
    !! @param[in] scalar 标量
    !! @return 标量乘的结果
    type(op_mat_cmplx_type) function mult_real_type(scalar, this)

        implicit none
        class(op_mat_cmplx_type), intent(in) :: this
        real(R_P) , intent(in) :: scalar

        mult_real_type%mat=this%mat*scalar

    end function mult_real_type

    !> 矩阵加法
    !! @return 矩阵的和
    type(op_mat_cmplx_type) function add(obj1, obj2)

        implicit none
        class(op_mat_cmplx_type), intent(in) :: obj1, obj2

        add%mat=obj1%mat+obj2%mat

    end function add

    !> 矩阵减法
    !! @return 矩阵的差
    type(op_mat_cmplx_type) function subtract(obj1, obj2)

        implicit none
        class(op_mat_cmplx_type), intent(in) :: obj1, obj2

        subtract%mat=obj1%mat-obj2%mat

    end function subtract

    !> 矩阵转置
    !! @return 矩阵的转置矩阵
    function transpose_op_mat(this) result(obj)

        implicit none
        class(op_mat_cmplx_type), intent(in) :: this
        type(op_mat_cmplx_type) :: obj

        obj%Mat=transpose(this%Mat)

    end function transpose_op_mat

    !> 矩阵求相反数.
    !!
    !! 矩阵前加负号
    function minus(this) result(obj)

        implicit none
        class(op_mat_cmplx_type), intent(in) :: this
        type(op_mat_cmplx_type) :: obj

        obj%Mat=-this%Mat

    end function minus

    !> 输出矩阵
    subroutine PrintMat(this)

        use mod_parameter, only: eno => ENO_DIS
        implicit none
        class(op_mat_cmplx_type), intent(in) :: this
        integer :: i, j

        do j=1, 5
          do i=1, 5
            write(eno, *)i, j, this%Mat(i, j)
          end do
        end do

    end subroutine PrintMat

end module mod_op_mat_cmplx
