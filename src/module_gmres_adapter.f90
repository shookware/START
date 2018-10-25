                                                                                                           !------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_gmres_adapter.f90
!> @file
!> @breif paridiso�������ļ�.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: gmres_adapter
!> @breif gmres������ģ��.
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
module mod_gmres_adapter

    use mod_gmres
    use mod_solver
    use penf, only: R_P

    private

    !> gmres����������
    type, public :: gmres_adapter_type

        private
 !       complex(R_P), pointer :: rhs(:) !< ���Է����Ҷ���@private
 !       complex(R_P), pointer :: X(:)   !< ���Է������Ľ�@private
        type(mkl_gmres) :: solver    !< ������@private

        Contains

!          Procedure :: Initialize => initialize_ap !< ��ʼ��
!          procedure :: Analysis => analysis_ap     !< ���̷���
!          generic :: Solve => solve_ap, solve_ap_rhs_x  !< ���ⷽ��
!          procedure :: Finalize => finalized_ap         !< ��������
!
!          Procedure :: SetNDim => set_ndim              !< ���÷���ά��
!          procedure :: SetNNZ => set_nnz                !< ����ϵ����������Ԫ�ظ���
          procedure :: initial=> initialize_ap
          procedure :: Set
          procedure :: SetMatA => set_matA              !< ����ϵ������
          procedure :: SetOP
          procedure :: SetPC
!          procedure :: SetRHS => set_rhs                !< ���÷����Ҷ���
!          procedure :: SetX => set_X                    !< ���÷��̵Ľ�
          generic :: GetX => get_X                      !< ���÷��̵Ľ�
          generic :: solve => solve_ap

          procedure, private :: solve_ap               !< �ⷽ��
!          procedure, private :: solve_ap_rhs_x         !< �ⷽ�̸�rhs,����x
          procedure, private :: get_X                  !< ���÷��̵Ľ�X

    end type gmres_adapter_type

        integer, pointer :: ia(:)  !< csr��ʽϵ����������Ϣ@private
        integer, pointer :: ja(:)  !< csr��ʽϵ����������Ϣ@private
        complex(R_P), pointer :: ma(:)  !< csr��ʽϵ��������0Ԫ��@private
        integer :: ndim !< ����ά�� @rivate
        integer :: nnz  !< ϵ��������0Ԫ�ظ���   @private

    abstract interface

        subroutine proc_real(content, x, rhs) ! Lx=rhs
            import
            Implicit None
            class(solver_type), intent(in) :: content
            real(8), intent(in) :: x(:)
            real(8), intent(out) :: rhs(:)

        end subroutine proc_real

    End interface

    contains

    subroutine Set(this, SolutionGuess, RHS)

        implicit none
        class(gmres_adapter_type),intent(inout) :: this
        complex(R_P), intent(in) :: SolutionGuess(:)
        complex(R_P), intent(in) :: RHS(:)

        call this%solver%SetGuess(SolutionGuess)
        call this%solver%SetRHS(RHS)

    end subroutine Set

    !> ��ʼ��
    !! @param[in] method ��������:(��������0,ת������1,����ת������2)
    subroutine initialize_ap(this, n, SolutionVec, RHS, content)

        use mod_solver
        implicit none
        class(gmres_adapter_type), intent(inout) :: this
        integer, intent(in) :: n
        complex(R_P), intent(in) :: solutionVec(:)
        complex(R_P), intent(in) :: RHS(:)
        class(solver_type) :: content

        this%solver=mkl_gmres(n, solutionVec, rhs)
        call this%solver%SetContent(content)
        call this%solver%check()

    end subroutine initialize_ap

    !> ���ⷽ��
    subroutine solve_ap(this)

        implicit none
        class(gmres_adapter_type), intent(inout) :: this

        call this%solver%solve()
      !call mv(this%solver%content, this%solver%rhs, this%solver%solution)

    end subroutine solve_ap

    !>����ϵ������
    !! @param[in] MatAcoo ϵ������A(COO��ʽ)
    subroutine set_matA(this, MatAcoo)

        use mod_sparse_matrix
        use mod_sparsekit
        implicit none
        class(gmres_adapter_type), intent(inout) :: this
        type(mat_coo_type) :: MatAcoo
        type(mat_csr_type), save :: MatA
        integer :: nz_num

        if(.not. MatA%IsCreate()) then
          call matAcoo%GetNRow(ndim)
          call matAcoo%GetNnz(nnz)
          call MatA%Create(ndim, nnz)
        else
          call MatA%zeros()
        endif

        call MatA%set(MatAcoo)

        call MatA%get(ia, ja, ma)
!        call this%solver%setop(av)
!        call this%solver%setpc(mv)
        call zilu( ndim, nnz, ia, ja, ma )
        call this%solver%setpc(mv_ilu)

    end subroutine set_matA

    subroutine setop(this, op)

        implicit none
        class(gmres_adapter_type),intent(inout) :: this
        procedure(proc_real) :: op

        call this%solver%setop(op)

    end subroutine setop

    subroutine setpc(this, op)

        implicit none
        class(gmres_adapter_type), intent(inout) :: this
        procedure(proc_real) :: op

        call this%solver%setpc(op)

    end subroutine setpc


!    subroutine av(content, x, b)
!    use mod_solver
!    implicit none
!    class(solver_type), intent(in) :: content
!    real(R_P), intent(in) ::x(:)
!    real(R_P), intent(out) :: b(:)
!    complex(R_P), allocatable :: ArrayX(:), Arrayb(:)
!    integer :: n
!
!        n=size(x)/2
!
!        allocate(ArrayX(n), Arrayb(n))
!        call array_complex(x, ArrayX)
!
!        !! OP*ArrayX=Arrayb
!!---------------------------------------------------------------------------
!!Initialize variables and the right hand side through matrix-vector product
!!---------------------------------------------------------------------------
!      call MKL_ZCSRGEMV('N', n, mA, IA, JA, arrayx, arrayb)  !AX=RHS
!!---------------------------------------------------------------------------
!
!        call array_real(Arrayb, b)
!        deallocate(ArrayX, Arrayb)
!
!    end subroutine av

    subroutine mv(content, b, x)

    use mod_solver
    use mod_pardiso
    implicit none
    class(solver_type), intent(in) :: content
    real(R_P), intent(in) :: b(:)
    real(R_P), intent(out) :: x(:)
    complex(R_P), allocatable :: ArrayX(:), Arrayb(:)
    !real(R_P), allocatable :: residual(:)
    !complex(R_P), allocatable :: zresidual(:)
    integer :: n
    type(pardiso_type), save :: psolver


        n=size(x)/2

        allocate(ArrayX(n), Arrayb(n))
        !allocate(zresidual(n))
        !allocate(residual(n*2))
        call array_complex(b, ArrayB)

        !! OP*ArrayX=Arrayb
        call psolver%Initialize(13, 0)
        call psolver%SetMatA(ma, ia, ja)
!        call psolver%SetNDim(ndim)
        call psolver%SetRHS(Arrayb)
        call psolver%SetX(arrayX)
        call psolver%Solve(13)
        call psolver%GetX(arrayX)

        call array_real(Arrayx, X)

        deallocate(ArrayX, Arrayb)

    end subroutine mv

    subroutine mv_ilu(content, b, x)

    use mod_solver
    use mod_pardiso
    implicit none
    class(solver_type), intent(in) :: content
    real(R_P), intent(in) :: b(:)
    real(R_P), intent(out) :: x(:)
    complex(R_P), allocatable :: ArrayX(:), Arrayb(:)
    complex(R_P), allocatable :: TRVEC(:)
    !real(R_P), allocatable :: residual(:)
    !complex(R_P), allocatable :: zresidual(:)
    integer :: n
    type(pardiso_type), save :: psolver


        n=size(x)/2

        allocate(ArrayX(n), Arrayb(n), trvec(n))
        !allocate(zresidual(n))
        !allocate(residual(n*2))
        call array_complex(b, ArrayB)

        !! OP*ArrayX=Arrayb
       CALL MKL_zCSRTRSV('L','N','U',N,ma,IA,JA,arrayb,TRVEC)
       CALL MKL_zCSRTRSV('U','N','N',N,ma,IA,JA,TRVEC,arrayx)

        call array_real(Arrayx, X)

        deallocate(ArrayX, Arrayb, trvec)

    end subroutine mv_ilu

    !> ���÷��̵Ľ�
    !! @param[out] X ���̵Ľ�X
    subroutine get_X(this, X)

        implicit none
        complex(R_P), intent(inout) :: X(:)
        class(gmres_adapter_type), intent(inout) :: this

        call this%solver%get(X)

    end subroutine get_X

end module mod_gmres_adapter
