!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_pardiso_adapter.f90
!> @file
!> @breif paridiso适配器文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: pardiso_adapter
!> @breif pardiso适配器模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-02 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-08-02
!------------------------------------------------------------------------------
module mod_pardiso_adapter

    use mod_pardiso
    use penf, only: R_P

    private

    !> pardiso适配器类型
    type, public :: pardiso_adapter_type

        private
        integer :: ndim !< 方程维数 @rivate
        integer :: nnz  !< 系数矩阵非0元素个数   @private
        integer, pointer :: ia(:)=>null()  !< csr格式系数矩阵行信息@private
        integer, pointer :: ja(:)=>null()  !< csr格式系数矩阵列信息@private
        complex(R_P), pointer :: ma(:)=>null()  !< csr格式系数矩阵非0元素@private
        complex(R_P), pointer :: rhs(:)=>null() !< 线性方程右端项@private
        complex(R_P), allocatable :: X(:)   !< 线性方程组的解@private
        type(pardiso_type) :: solver    !< 解算器@private

        Contains

          Procedure :: Initialize => initialize_ap !< 初始化
          procedure :: Analysis => analysis_ap     !< 方程分析
          generic :: Solve => solve_ap, solve_ap_rhs_x  !< 求解方程
          procedure :: Finalize => finalized_ap         !< 析构函数

          generic :: SetMatA => set_matA_coo, set_matA_bsr !< 设置系数矩阵
          Procedure :: SetNDim => set_ndim              !< 设置方程维数
          procedure :: SetNNZ => set_nnz                !< 设置系数矩阵非零元素个数
          procedure :: SetRHS => set_rhs                !< 设置方程右端项
          procedure :: SetX => set_X                    !< 设置方程的解
          generic :: GetX => get_X                      !< 获得方程的解

          procedure, private :: solve_ap               !< 解方程
          procedure, private :: solve_ap_rhs_x         !< 解方程给rhs,获得x
          procedure, private :: get_X                  !< 获得方程的解X

          procedure, private :: Set_matA_coo  !< 设置系数矩阵(COO)
          procedure, private :: Set_matA_bsr  !< 设置系数矩阵(BSR)

    end type pardiso_adapter_type

    contains

    !> 初始化
    !! @param[in] method 方程类型:(正常问题0,转置问题1,共轭转置问题2)
    subroutine initialize_ap(this, method)

        implicit none
        class(pardiso_adapter_type), intent(inout) :: this
        integer, intent(in) :: method
        integer, parameter :: mtype=13

        call this%solver%Initialize(mtype, method)

    end subroutine initialize_ap

    !> 分析方程
    subroutine analysis_ap(this)

        implicit none
        class(pardiso_adapter_type), intent(inout) :: this

        call this%solver%Analysis()

    end subroutine analysis_ap

    !> 求解方程
    subroutine solve_ap(this)

        implicit none
        class(pardiso_adapter_type), intent(inout) :: this

        call this%solver%solve()

    end subroutine solve_ap

    !> 给RHS,求解方程,同时返回X
    !! @param[in] RHS 方程的右端项b
    !! @param[out] RHS 方程的解X
    subroutine solve_ap_rhs_x(this, RHS)

        implicit none
        class(pardiso_adapter_type), intent(inout) :: this
        complex(R_P), intent(inout) :: RHS(:)

        call this%SetRHS(RHS)
        call this%Solve()
        call this%GetX(RHS)

    end subroutine solve_ap_rhs_x

    !> 设置方程的维数
    !! @param[in] ndim 方程维数
    subroutine set_ndim(this, ndim, blocksize)

        implicit none
        integer, intent(in) :: ndim
        integer, intent(in) :: blocksize
        class(pardiso_adapter_type), intent(inout) :: this
        

        this%ndim=ndim*blocksize
!        call this%solver%SetNDim(ndim*blocksize)
        if(.not. allocated(this%X)) allocate(this%X(ndim*blocksize))
!        if(.not. allocated(this%ia)) then
!          allocate(this%ia(ndim+1))
!        endif

    end subroutine set_ndim

    !> 设置系数矩阵非零元素个数
    !! @param[in] nnz 系数矩阵非零元素个数
    subroutine set_nnz(this, nnz, blocksize)

        implicit none
        integer, intent(in) :: nnz
        integer, intent(in) :: blocksize
        class(pardiso_adapter_type), intent(inout) :: this

!        if(.not. allocated(this%ja)) then
          this%nnz=nnz*blocksize*blocksize
!          allocate(this%ja(nnz), this%ma(nnz))
!        endif

    end subroutine set_nnz

    !>设置系数矩阵(COO格式)
    !! @param[in] MatAcoo 系数矩阵A(COO格式)
    subroutine set_matA_coo(this, MatAcoo)

        use mod_sparse_matrix
        implicit none
        class(pardiso_adapter_type), intent(inout) :: this
        type(mat_coo_type) :: MatAcoo
        type(mat_csr_type), save :: MatA

        if(.not. MatA%IsCreate()) then
          call MatA%Create(this%ndim, this%nnz)
        else
          call MatA%zeros()
        endif

        call MatA%set(MatAcoo)
        call MatA%get(this%ia, this%ja, this%ma)
        call this%solver%SetMatA(this%ma, this%ia, this%ja)

    end subroutine set_matA_coo

    !>设置系数矩阵(BSR格式)
    !! @param[in] MatAcoo 系数矩阵A(BSR格式)
    subroutine set_matA_bsr(this, MatAbsr)

        use mod_sparse_matrix
        implicit none
        class(pardiso_adapter_type), intent(inout) :: this
        type(mat_bsr_type) :: MatAbsr

        call MatAbsr%get(this%ia, this%ja, this%ma)
        call this%solver%SetMatA(this%ma, this%ia, this%ja)

    end subroutine set_matA_bsr

    !> 设置方程右端项
    !! @param[in] RHS 方程右端项b
    subroutine set_rhs(this, RHS)


        implicit none
        complex(R_P), intent(in), target :: RHS(:)
        class(pardiso_adapter_type), intent(inout) :: this

        if((associated(this%rhs))) this%rhs=>null() 
        this%rhs => rhs

        call this%solver%SetRHS(this%rhs)

    end subroutine set_rhs

    !> 设置方程的解的初值
    !! @param[in] X 方程解初值X
    subroutine set_X(this, X)

    implicit none
    class(pardiso_adapter_type), intent(inout) :: this
    complex(R_P), intent(in) :: X(:)

        if(allocated(this%X)) allocate(this%X(size(X)))
        this%X = X

        call this%solver%SetX(this%X)

    end subroutine set_X

    !> 获得方程的解
    !! @param[out] X 方程的解X
    subroutine get_X(this, X)

        implicit none
        complex(R_P), intent(inout) :: X(:)
        class(pardiso_adapter_type), intent(inout) :: this

        call this%solver%GetX(this%X)
        X=this%X

    end subroutine get_X

    !> 析构函数
    subroutine finalized_ap(this)

        implicit none
        class(pardiso_adapter_type), intent(inout) :: this

        call this%solver%Finalize()
        nullify(this%ja, this%ia, this%ma, this%rhs)
        deallocate(this%x)

    end subroutine finalized_ap


end module mod_pardiso_adapter
