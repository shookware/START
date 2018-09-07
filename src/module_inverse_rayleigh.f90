!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_inverse_rayleigh.f90
!> @file
!> @breif 反幂法求特征值文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: inverse_rayleigh
!> @breif 反幂法求特征值模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-03 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-03
!------------------------------------------------------------------------------
module mod_inverse_rayleigh

    use mod_sparse_matrix
    use mod_pardiso
    use penf, only: R_P
    implicit none

    private

    !> 反幂法求特征值解算器类
    type, public :: inverse_rayleigh_type

        private

        integer :: ndim !< 矩阵维数
        complex(R_P) :: lambda=0.0d0 !< 特征值
        complex(R_P) :: delta=0.0d0 !< 位移量
        type(mat_csr_type) :: A !< 广义问题左矩阵
        type(mat_csr_type) :: B !< 广义问题右矩阵
        type(mat_csr_type) :: matA !< 广义问题左矩阵(谱平移变换前）
        complex(R_P), allocatable :: phi(:)  !< eigenfunction
        complex(R_P), allocatable :: psai(:) !< adjoint eigenfunction
        integer :: inner_loop=1 !< 内循环次数
        integer :: outer_loop=0 !< 外循环次数
        type(pardiso_type) :: solver !< pardiso解算器

        Contains

            generic :: SetMatrix=> SetMatrix_coo !< 设置矩阵
            procedure :: SetEigenValueGuess !< 设置特征值初值
            procedure :: solve !< 求解过程
            procedure :: GetEigenFunc !< 获得特征函数
            procedure :: GetEigenFuncAdjoint !< 获得特征函数的伴随
            procedure :: GetEigenValue !< 获得特征值

            procedure, private :: UpdateEigenfun !< 更新特征函数
            procedure, private :: UpdateLambda !< 更新特征值
            procedure, private :: check !< 检查特征值方程两边是否相等
!            procedure, private :: SetMatrix_bsr !< 设置矩阵(BSR)
            procedure, private :: SetMatrix_coo !< 设置矩阵(COO)

    end type inverse_rayleigh_type

    real(R_P), parameter :: EPS=1.0d-9

    contains

    !> 设置矩阵
    !! @param[in] ndim 矩阵维数
    !! @param[in] A 左矩阵
    !! @param[in] B 右矩阵
    subroutine SetMatrix_coo(this, ndim, A, B)

        implicit none
        class(inverse_rayleigh_type), intent(inout) :: this
        integer, intent(in) :: ndim
        type(mat_coo_type), intent(inout) :: A
        type(mat_coo_type), intent(inout), optional :: B
        integer :: i, l
        complex(R_P), parameter :: one=(1.0d0, 0.0d0)
        real(R_P), parameter :: pi=3.141592653590d0
        type(mat_coo_type) :: B_identify
        real(R_P) :: random(ndim)

        this%ndim=ndim
        call this%A%set(A)

        if( present(B)) then
            call this%B%set(B)
        else
            call B_identify%create(ndim, ndim)
            do i=1, ndim
                call B_identify%Set(i, i, one)
            end do
            call this%B%Set(B_identify)
        endif
        if(.not. allocated(this%phi)) then
            allocate(this%phi(ndim), this%psai(ndim))
        endif

        this%phi=0.0d0
        this%psai=0.0d0

        do i=1, this%ndim
            call random_number(random)
            this%phi=random
            call random_number(random)
            this%psai=random
        end do

        this%outer_loop=0
!        call this%solver%SetNDim(ndim)
        call this%solver%Initialize(13, 0)

    end subroutine SetMatrix_coo

!    !> 设置矩阵
!    !! @param[in] ndim 矩阵维数
!    !! @param[in] A 左矩阵
!    !! @param[in] B 右矩阵
!    subroutine SetMatrix_bsr(this, ndim, blocksize, A, B)
!
!        implicit none
!        class(inverse_rayleigh_type), intent(inout) :: this
!        integer, intent(in) :: ndim
!        integer, intent(in) :: blocksize
!        type(mat_bsr_type), intent(inout) :: A
!        type(mat_bsr_type), intent(inout), optional :: B
!        integer :: i, l
!        complex(R_P), parameter :: one=(1.0d0, 0.0d0)
!        real(R_P), parameter :: pi=3.141592653590d0
!        type(mat_coo_type) :: B_identify
!        real(R_P) :: random(ndim)
!
!        this%ndim=ndim
!        call this%A%set(A)
!
!        if( present(B)) then
!            call this%B%set(B)
!        else
!            call B_identify%create(ndim, ndim)
!            do i=1, ndim
!                call B_identify%Set(i, i, one)
!            end do
!            call this%B%Set(B_identify)
!        endif
!        if(.not. allocated(this%phi)) then
!            allocate(this%phi(ndim), this%psai(ndim))
!        endif
!
!        this%phi=0.0d0
!        this%psai=0.0d0
!
!        do i=1, this%ndim
!            call random_number(random)
!            this%phi=random
!            call random_number(random)
!            this%psai=random
!        end do
!
!        this%outer_loop=0
!!        call this%solver%SetNDim(ndim)
!        call this%solver%Initialize(13, 0)
!
!    end subroutine SetMatrix_bsr



    !> 设置特征值初值
    !! @param[in] lambda 特征值初值
    subroutine SetEigenValueGuess(this, lambda)

        implicit none
        class(inverse_rayleigh_type), intent(inout) :: this
        complex(R_P), intent(in) :: lambda
        type(mat_csr_type) :: matA

        this%lambda=lambda
        this%delta=lambda

        matA=this%A-lambda*this%B
        this%A=matA
        this%lambda=0.0d0

    end subroutine SetEigenValueGuess

    !> 反幂法求解过程
    subroutine solve(this)

        implicit none
        class(inverse_rayleigh_type), intent(inout) :: this
        complex(R_P) :: dlambda, lambda_old, lambda_new

        dlambda=1.0d9

        do while (abs(dlambda)>=EPS)
            lambda_old=this%lambda
            call this%UpdateEigenfun()
            call this%UpdateLambda()
            lambda_new=this%lambda
            dlambda=lambda_new-lambda_old
        !    print*, dlambda
        end do

        !call this%check()
        call this%A%clear
        call this%B%clear

    end subroutine solve

    !> 检查特征值问题等号两边是否相等
    subroutine check(this)

        implicit none
        class(inverse_rayleigh_type), intent(inout) :: this
        complex(R_P) :: outvalue(this%ndim)

        outvalue =(this%A .mx. this%phi) - (this%lambda* (this%B.mx. this%phi))
        print*, this%lambda+this%delta

    end subroutine check

    !> 更新特征值
    subroutine UpdateLambda(this)

        implicit none
        class(inverse_rayleigh_type), intent(inout) :: this

        complex(R_P) :: dot1, dot2

        dot1=dot_product(conjg(this%psai), this%A .mx. this%phi)
        dot2=dot_product(conjg(this%psai), this%B .mx. this%phi)
        !dot1=dot_product(this%A .mx. this%phi, this%B .mx. this%phi)
        !dot2=dot_product(this%B .mx. this%phi, this%B .mx. this%phi)

        this%lambda=dot1/dot2

        this%outer_loop=this%outer_loop+1

    end subroutine UpdateLambda

    !> 更新特征函数和它的伴随
    subroutine UpdateEigenfun(this)

        implicit none
        class(inverse_rayleigh_type), intent(inout) :: this
        integer :: i
        complex(R_P) :: rhs(this%ndim)
        type(mat_csr_type) :: matA
        complex(R_P) :: tmp(this%ndim)

        rhs=0.0d0

        matA=this%A-this%lambda*this%B
        call this%solver%SetMatA(matA)

        do i=1, this%inner_loop
            rhs=this%B .mx. this%phi
            call this%solver%SetRHS(rhs)
            call this%solver%Solve(0)
            call this%solver%GetX(this%phi)

            rhs=this%B .tmx. this%psai
            call this%solver%SetRHS(rhs)
            call this%solver%Solve(2)
            call this%solver%GetX(this%psai)
        end do
        call normalize(this%ndim, this%phi)
        call normalize(this%ndim, this%psai)

    end subroutine UpdateEigenfun

    !> 获得特征函数
    !! @param[out] phi 特征函数
    subroutine GetEigenFunc(this, phi)

        implicit none
        class(inverse_rayleigh_type), intent(inout) :: this
        complex(R_P), intent(out) :: phi(this%ndim)

        phi=this%phi

    end subroutine GetEigenFunc

    !> 获得特征函数的伴随
    !! @param[out] psai 特征函数的伴随
    subroutine GetEigenFuncAdjoint(this, psai)

        implicit none
        class(inverse_rayleigh_type), intent(inout) :: this
        complex(R_P), intent(out) :: psai(this%ndim)

        psai=this%psai

    end subroutine GetEigenFuncAdjoint

    !> 获得特征值
    !! @param[out] lambda 特征值
    subroutine GetEigenValue(this, lambda)

        implicit none
        class(inverse_rayleigh_type), intent(inout) :: this
        complex(R_P) :: lambda

        lambda=this%lambda+this%delta

    end subroutine GetEigenValue

    !> 对向量归一化
    !! @param[in,out] 待归一向量
    subroutine normalize(n, array)

        implicit none
        integer, intent(in) :: n
        complex(R_P), intent(inout) :: array(n)
        real(R_P) :: absarray(n)
        integer :: imax

        absarray=abs(array)
        imax=maxloc(absarray, 1)

        array=array/array(imax)

    end subroutine normalize

end module mod_inverse_rayleigh
