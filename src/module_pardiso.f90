!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: /media/psf/Home/Fortran/START/src/module_pardiso.f90
!> @file
!> @breif pardiso模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: pardiso
!> @breif pardiso模块.
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
module mod_pardiso

    use penf, only: R_P
    implicit none

    private

    !> pardiso类型
    !!
    !! 使用pardiso库计算线性方程组AX=b.
    !!  方程类型|method参数
    !!  --------|:-----:
    !!  \f$ AX=b \f$|0
    !! \f$ A^HX=b \f$|1
    !! \f$ A^TX=b \f$|2
    type, public :: pardiso_type

        private
        logical, private :: Flg_analysis=.False. !< 是否已经执行过分析
        logical, private :: Flg_create=.False.   !< 是否计算容器创建
        integer(8), private :: pt(64)             !< pardiso计算器指针
        integer, private :: mtype  !< 系数矩阵类型
        integer, private :: method !< 求解类型:一般问题,转置问题,共轭转置问题
        integer, private :: ierr   !< 出错信号
        integer, private :: msglvl !< 信息显示
        integer, private :: idum(1)     !<站位临时数组
        integer, private :: iparm(64)   !<pardiso运行参数
        integer, private :: ndim        !<方程维数
        integer, private :: blocksize !< 块大小
        complex(R_P), private, pointer :: ma(:)  !<数组非0元素
        integer, private, pointer :: ia(:)       !<数组元素行个数
        integer, private, pointer :: ja(:)       !<数组元素列序号
        complex(R_P), private, pointer :: rhs(:) !方程右向量b
        complex(R_P), private, allocatable :: X(:) !方程的解X

        Contains

          Procedure :: Initialize => initialize_pardiso       !<初始化
          procedure :: Analysis => analysis_pardiso           !<分析问题
          procedure :: Solve => solve_pardiso                 !<求解方程
          procedure :: SetRHS => set_rhs                      !<设置右向量b
          procedure :: SetX => set_X                          !<设置方程的解X的猜测值
          procedure :: GetX => get_X                          !<获得方程的解X
          procedure :: Finalize => finalize_pardiso           !<析构函数

          generic :: SetMatA => set_matA, set_matA_csr, set_matA_bsr         !<设置方程的系数矩阵A

          procedure, private :: set_matA_csr                  !< 直接给csr矩阵设置A
          procedure, private :: set_matA_bsr                  !< 直接给bsr矩阵设置A
          procedure, private :: set_matA                      !< 分量形式设置A
          procedure, private :: SetNDim => set_ndim_pardiso            !<设置方程维数

    end type pardiso_type

    contains

    !> 初始化
    !! @param[in] mtype 系数矩阵类型
    !! @param[in] method 线性方程类型
    subroutine initialize_pardiso(this, mtype, method)

        implicit none
        class(pardiso_type), intent(inout) :: this
        integer, intent(in) :: mtype, method

        if(.not.this%Flg_create) then
            this%mtype=mtype
            this%method=method
            call pardisoinit(this%pt, this%mtype, this%iparm)
        endif

    end subroutine initialize_pardiso

    !> 设置方程维数
    !! @param[in] ndim 方程维数
    subroutine set_ndim_pardiso(this, ndim)

        implicit none
        integer, intent(in) :: ndim
        class(pardiso_type), intent(inout) :: this

        this%ndim=ndim
        if(.not. this%Flg_create)then
            allocate(this%X(ndim))
            this%Flg_create=.True.
        endif

    end subroutine set_ndim_pardiso

    !> 分量形式给系数矩阵A
    !! @param[in] ma csr形式矩阵非0元素
    !! @param[in] ia csr形式矩阵每行非0元素个数
    !! @param[in] ja csr形式矩阵非0元素列位置
    subroutine set_matA(this, ma, ia, ja)

        implicit none
        integer, intent(in), target :: ia(:), ja(:)
        complex(R_P), intent(in), target :: ma(:)
        class(pardiso_type), intent(inout) :: this
        integer :: ndim, nnz, blocksize, nblock

        this%ma => ma; this%ia => ia; this%ja => ja
        nnz=size(ma)
        ndim=size(this%ia)-1
        nblock=this%ia(ndim+1)-1
        blocksize=nint(sqrt(real(nnz/nblock, R_P)))
        this%blocksize=blocksize
        call this%SetNDim(ndim*blocksize)

    end subroutine set_matA

    !> csr形式直接设置系数矩阵A
    !! @param[in] matA csr形式系数矩阵A
    !! @param[in] method 线性问题类型(可选)
    subroutine set_matA_csr(this, matA, method)

        use mod_sparse_matrix
        implicit none
        class(pardiso_type), intent(inout) :: this
        type(mat_csr_type), intent(in) :: matA
        integer, intent(in), optional :: method
        integer, pointer :: ia(:), ja(:)
        complex(R_P), pointer :: ma(:)
        integer :: ndim, nnz, blocksize, nblock

        call matA%get(ia, ja, ma)

        this%ma => ma; this%ia => ia; this%ja => ja
        nnz=size(ma)
        ndim=size(this%ia)-1
        nblock=this%ia(ndim)-1
        blocksize=nint(sqrt(real(nnz/nblock, R_P)))
        this%blocksize=blocksize
        call this%SetNDim(ndim*blocksize)

        !this%iparm(2)=3
        !this%iparm(4)=91 !61
        !this%iparm(12)=this%method
        !if(present(method)) &
        !    & this%iparm(12)=method
        !this%msglvl=0
        !
        !call pardiso(this%pt, 1, 1, this%mtype, 12, this%ndim, this%ma, this%ia, this%ja, this%idum, 1 &
        !            , this%iparm, this%msglvl, this%rhs, this%X, this%ierr)

    end subroutine set_matA_csr

    !> bsr形式直接设置系数矩阵A
    !! @param[in] matA bsr形式系数矩阵A
    !! @param[in] method 线性问题类型(可选)
    subroutine set_matA_bsr(this, matA, method)

        use mod_sparse_matrix
        implicit none
        class(pardiso_type), intent(inout) :: this
        type(mat_bsr_type), intent(in) :: matA
        integer, intent(in), optional :: method
        integer, pointer :: ia(:), ja(:)
        complex(R_P), pointer :: ma(:)
        integer :: ndim, nnz, blocksize, nblock

        call matA%get(ia, ja, ma)

        this%ma => ma; this%ia => ia; this%ja => ja
        nnz=size(ma)
        ndim=size(this%ia)-1
        nblock=this%ia(ndim)-1
        blocksize=nint(sqrt(real(nnz/nblock, R_P)))
        this%blocksize=blocksize
        call this%SetNDim(ndim*blocksize)


    end subroutine set_matA_bsr

    !> 设置方程右向量b
    !! @param[in] rhs 方程右向量b
    subroutine set_RHS(this, rhs)

        implicit none
        complex(R_P), intent(in), target :: rhs(:)
        class(pardiso_type), intent(inout) :: this

        this%rhs => rhs

    end subroutine set_RHS

    !> 设置方程解初值X
    !! @param[in] X 方程的解X
    subroutine set_X(this, X)

        implicit none
        complex(R_P), intent(in) :: X(:)
        class(pardiso_type), intent(inout) :: this

        this%X=X

    end subroutine set_X

    !> 获得方程解初值X
    !! @param[out] X 方程的解X
    subroutine get_X(this, X)

        implicit none
        class(pardiso_type), intent(inout) :: this
        complex(R_P), intent(out) :: X(this%ndim)

        X=this%X

    end subroutine get_X

    !> 分析线性方程
    !! @param[in] method 方程类型(可选)
    subroutine analysis_pardiso(this, method)

        implicit none
        class(pardiso_type), intent(inout) :: this
        integer , intent(in), optional :: method
        !complex(R_P) :: ma((2498-1)*25)
        !integer(4) :: ia(this%ndim/5+1)
        !integer(4) :: ja((2498-1))
        !integer :: job(6)
        !integer :: info

        this%Flg_analysis=.True.

        !job=0
        !job(1)=0
        !job(2)=1
        !job(3)=1
        !job(6)=1
        !
        !call mkl_zcsrbsr(job, this%ndim, 5, 5*5, this%ma, this%ja, this%ia, ma, ja, ia, info)

        this%iparm(2)=3
        this%iparm(4)=0 !61
        this%iparm(12)=this%method
        if(this%blocksize==1) then
            this%iparm(37)=0
        else
            this%iparm(37)=this%blocksize
        endif

        this%iparm(11)=0
        this%iparm(13)=0
        if(present(method)) &
            & this%iparm(12)=method
        this%msglvl=0
        call pardiso(this%pt, 1, 1, this%mtype, 11, this%ndim/this%blocksize, this%ma, this%ia, this%ja, this%idum, 1 &
                    , this%iparm, this%msglvl, this%idum, this%idum, this%ierr)
!        call pardiso(this%pt, 1, 1, this%mtype, 11, this%ndim/5, ma, ia, ja, this%idum, 1 &
!                    , this%iparm, this%msglvl, this%idum, this%idum, this%ierr)

    end subroutine analysis_pardiso

    !> 求解线性方程组
    !! @param method 线性问题类型(可选)
    subroutine solve_pardiso(this, method)

        implicit none
        class(pardiso_type), intent(inout) :: this
        integer, intent(in), optional :: method
        integer i
        integer :: jstep
        !complex(R_P) :: ma(this%ia(this%ndim+1)-1)
        !integer(4) :: ia(this%ndim/5+1)
        !integer(4) :: ja((this%ia(this%ndim+1)-1)/25)
        !integer :: job(6)
        !integer :: info
        !

        !print*, this%ia(1), this%ia(this%ndim+1)
        !print*, this%ndim/5+1
        !pause
        if(.not. this%Flg_analysis) call this%Analysis()


        !job=0
        !job(1)=0
        !job(2)=1
        !job(3)=1
        !job(6)=1
        !
        !call mkl_zcsrbsr(job, this%ndim, 5, 5*5, this%ma, this%ja, this%ia, ma, ja, ia, info)
        !

        this%iparm(2)=3
        this%iparm(4)=0 !61
        if(this%blocksize==1) then
            this%iparm(37)=0
        else
            this%iparm(37)=this%blocksize
        endif
        this%iparm(11)=0
        this%iparm(13)=0
        this%iparm(12)=this%method
        if(present(method)) &
            & this%iparm(12)=method
        this%msglvl=0

        call pardiso(this%pt, 1, 1, this%mtype, 23, this%ndim/this%blocksize, this%ma, this%ia, this%ja, this%idum, 1 &
                    , this%iparm, this%msglvl, this%rhs, this%X, this%ierr)

        !this%iparm(2)=3
        !this%iparm(4)=0 !61
        !this%iparm(37)=5
        !this%iparm(11)=0
        !this%iparm(13)=0
        !this%iparm(12)=this%method
        !if(present(method)) &
        !    & this%iparm(12)=method
        !this%msglvl=1
        !
        !
        !call pardiso(this%pt, 1, 1, this%mtype, 13, this%ndim/5, ma, ia, ja, this%idum, 1 &
        !            , this%iparm, this%msglvl, this%rhs, this%X, this%ierr)
        !
        !stop
    end subroutine solve_pardiso


    !> 析构函数
    subroutine finalize_pardiso(this)

        implicit none
        class(pardiso_type), intent(inout) :: this

        call pardiso(this%pt, 1, 1, this%mtype, -1, this%ndim, this%ma, this%ia, this%ja, this%idum, 1 &
                , this%iparm, this%msglvl, this%idum, this%idum, this%ierr)
        deallocate(this%X)
        this%Flg_create=.False.

    end subroutine finalize_pardiso

end module mod_pardiso
