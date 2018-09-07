!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  Module: mkl��gmres��װ��
!
!  DESCRIPTION:
!> @breif mkl��gmres�㷨.
!>
!!
!  REVISION HISTORY:
!  2017-09-21 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-09-20
!------------------------------------------------------------------------------
module mod_gmres

    use mod_solver, only: solver_type
    use penf, only: R_P
    implicit none

    private
    public :: mkl_gmres, array_real, array_complex

    integer, parameter, private :: SIZE_=128

    type :: mkl_gmres

        private
        integer :: ipar(size_)
        real*8 :: dpar(size_)
        integer :: rci_request
        integer :: itercount
        integer :: n !dimenision of the problem
        real*8, allocatable :: tmp(:)
        real*8, allocatable, public :: rhs(:)
        real*8, allocatable, public :: solution(:)
        procedure(proc_real), nopass, pointer, private :: proc
        procedure(proc_real), nopass, pointer, private :: pcproc
        class(solver_type), pointer, public :: content
        logical :: IsInitial=.False.
!         complex(R_P), private, pointer :: ma(:)  !<csr������0Ԫ��
!        integer, private, pointer :: ia(:)       !<csr����Ԫ���и���
!       integer, private, pointer :: ja(:)       !<csr����Ԫ��������

        Contains

          procedure, private :: initial
          procedure :: SetContent
          procedure :: setop
          procedure :: setpc
          procedure :: check
          procedure :: solve
          procedure :: SetRHS
          procedure :: SetGuess

          generic :: get => get_real, get_complex
          procedure, private :: get_real
          procedure, private :: get_complex

          !generic :: SetMatA => set_matA_csr, set_matA         !<���÷��̵�ϵ������A
          !
          !procedure, private :: set_matA_csr                  !< ֱ�Ӹ�csr��������A
          !procedure, private :: set_matA

    end type mkl_gmres

    abstract interface

        subroutine proc_real(content, x, rhs) ! Lx=rhs
            import
            Implicit None
            class(solver_type), intent(in) :: content
            real*8, intent(in) :: x(:)
            real*8, intent(out) :: rhs(:)

        end subroutine proc_real

        subroutine proc_complex(content, x, rhs) ! Lx=rhs
            import
            Implicit None
            class(solver_type), intent(in) :: content
            real*8, intent(in) :: x(:)
            real*8, intent(out) :: rhs(:)
        end subroutine proc_complex

    End interface

    interface mkl_gmres

        module procedure construct_real
        module procedure construct_complex

    End interface mkl_gmres

    contains

    subroutine SetGuess(this, Guess)

        implicit none
        class(mkl_gmres),intent(inout) :: this
        complex(R_P), intent(in) :: Guess(:)

        call array_real(guess, this%solution)

    end subroutine SetGuess

    subroutine SetRHS(this, RHS)

        implicit none
        class(mkl_gmres),intent(inout) :: this
        complex(R_P), intent(in) :: RHS(:)

        call array_real(rhs, this%rhs)

    end subroutine SetRHS


    subroutine SetContent(this, Content)

        implicit none
        class(mkl_gmres), intent(inout) :: this
        class(solver_type), target, intent(in) :: Content

        this%Content=> content

    end subroutine SetContent

    function construct_real(n, solution_guess, rhs) result(this)

        implicit none
        type(mkl_gmres) :: this
        integer, intent(in) :: n
        real*8, intent(in) :: solution_guess(n)
        real*8, intent(in) :: rhs(n)

        this%n=n
        allocate(this%tmp(this%n*(2*this%n+1)+(this%n*(this%n+9))/2+1))
        this%tmp=0.0d0
        allocate(this%rhs(this%n))
        this%rhs=rhs
        allocate(this%solution(this%n))
        this%solution=solution_guess
        call this%initial()

    end function construct_real

    function construct_complex(n, solution_guess, rhs) result(this)

        implicit none
        type(mkl_gmres) :: this
        integer, intent(in) :: n
        complex*16, intent(in) :: solution_guess(n)
        complex*16, intent(in) :: rhs(n)

        this%n=n*2
        if(.not. allocated(this%tmp)) &
        & allocate(this%tmp(this%n*(2*this%n+1)+(this%n*(this%n+9))/2+1))
        if(.not. allocated(this%rhs)) allocate(this%rhs(this%n))
        if(.not. allocated(this%solution)) allocate(this%solution(this%n))
        this%tmp=0.0d0
        call array_real(rhs, this%rhs)
        call array_real(solution_guess, this%solution)
        if(.not. This%IsInitial) call this%initial()

    end function construct_complex

    subroutine initial(this)

        implicit none
        class(mkl_gmres),intent(inout) :: this

        call DFGMRES_INIT(this%n, this%solution, this%rhs, this%rci_request, this%ipar, &
     &                   this%dpar, this%tmp)

        if (this%rci_request.ne.0) call throw_error(this%rci_request)

        this%ipar(15) = min( this%n/2, 150)     ! do the restart after 2 iterations
        this%ipar( 8) = 0            ! do not do the stopping test for the maximal number of iterations
        this%ipar( 9) = 0            ! do not do the stopping test for the maximal number of iterations
        this%ipar(11) = 0            ! precondition
        this%dpar( 1) = 0.0d0        ! tolerance
        this%IsInitial=.True.

    end subroutine initial

    subroutine array_real(array_in, array_out)

        implicit none
        complex*16 :: array_in(:)
        real*8 :: array_out(:)
        integer :: i

        do i=1, size(array_in)
            array_out((i-1)*2+1)=real(array_in(i))
            array_out((i-1)*2+2)=aimag(array_in(i))
        end do

    end subroutine array_real

    subroutine array_complex(array_in, array_out)

        implicit none
        real*8 :: array_in(:)
        complex*16 :: array_out(:)
        integer :: i, num

        num=size(array_in)
        do i=1, num, 2
            array_out(i/2+1)=cmplx(array_in(i), array_in(i+1), 8)
        end do

    end subroutine array_complex

    subroutine check(this)

        implicit none
        class(mkl_gmres),intent(inout) :: this

    !---------------------------------------------------------------------------
    ! Check the correctness and consistency of the newly set parameters
    !---------------------------------------------------------------------------
          call DFGMRES_CHECK(this%N, this%solution, this%rhs, this%rci_request, &
         & this%ipar, this%dpar, this%tmp)
          if (this%rci_request.ne.0) call throw_error(this%rci_request)
    !---------------------------------------------------------------------------
    ! print the info about the RCI FGMRES method
    !---------------------------------------------------------------------------
!          print *, ''
!          print *,'Some info about the current run of RCI FGMRES method:'
!          print *, ''
!          if (this%ipar(8).ne.0) then
!             write(*,'(A,I1,A,A)') 'As ipar(8)=',this%ipar(8),', the automatic', &
!         & ' test for the maximal number of iterations will be'
!             print *,'performed'
!          else
!            write(*,'(A,I1,A,A)') 'As ipar(8)=',this%ipar(8),', the automatic', &
!         & ' test for the maximal number of iterations will be'
!            print *,'skipped'
!          endif
!          print *,'+++'
!          if (this%ipar(9).ne.0) then
!            write(*,'(A,I1,A,A)') 'As ipar(9)=',this%ipar(9),', the automatic', &
!         & ' residual test will be performed'
!          else
!            write(*,'(A,I1,A,A)') 'As ipar(9)=',this%ipar(9),', the automatic', &
!         & ' residual test will be skipped'
!          endif
!          print *,'+++'
!          if (this%ipar(10).NE.0) then
!            write(*,'(A,I1,A,A)') 'As ipar(10)=',this%ipar(10),', the', &
!         & ' user-defined stopping test will be requested via'
!            print *,'rci_request=2'
!          else
!            write(*,'(A,I1,A,A)') 'As ipar(10)=',this%ipar(10),', the', &
!         & ' user-defined stopping test will not be requested, thus,'
!            print *,'rci_request will not take the value 2'
!          endif
!          print *,'+++'
!          if (this%ipar(11).NE.0) then
!            write(*,'(A,I1,A,A)') 'As ipar(11)=',this%ipar(11),', the', &
!         & ' Preconditioned FGMRES iterations will be performed, thus,'
!            write(*,'(A,A)') 'the preconditioner action will be requested', &
!         & ' via rci_request=3'
!          else
!            write(*,'(A,I1,A,A)') 'As ipar(11)=',this%ipar(11),', the', &
!         & ' Preconditioned FGMRES iterations will not be performed,'
!            write(*,'(A)') 'thus, rci_request will not take the value 3'
!          endif
!          print *,'+++'
!          if (this%ipar(12).NE.0) then
!            write(*,'(A,I1,A,A)') 'As ipar(12)=',this%ipar(12),', the automatic', &
!         & ' test for the norm of the next generated vector is'
!            write(*,'(A,A)') 'not equal to zero up to rounding and', &
!         & ' computational errors will be performed,'
!            print *,'thus, rci_request will not take the value 4'
!          else
!            write(*,'(A,I1,A,A)') 'As ipar(12)=',this%ipar(12),', the automatic', &
!         & ' test for the norm of the next generated vector is'
!            write(*,'(A,A)') 'not equal to zero up to rounding and', &
!         & ' computational errors will be skipped,'
!            write(*,'(A,A)') 'thus, the user-defined test will be requested', &
!         & ' via rci_request=4'
!          endif
!          print *,'+++'

    end subroutine check

    subroutine setop(this, proc)

        implicit none
        class(mkl_gmres),intent(inout) :: this
        procedure(proc_real) :: proc

        this%proc=>proc

    end subroutine setop

    subroutine setpc(this, proc)

        implicit none
        class(mkl_gmres),intent(inout) :: this
        procedure(proc_real) :: proc

        this%pcproc=>proc
        this%ipar(11) = 1     ! precondition

    end subroutine setpc

    subroutine solve(this)

        implicit none
        class(mkl_gmres),intent(inout) :: this
        integer :: i
        real*8 :: b(this%n), residual(this%n)
        real*8 :: dvar
        logical :: loop
                    real*8 :: B1(this%n) , B2(This%n), Bres(this%n)

        !---------------------------------------------------------------------------
! An external BLAS function is taken from MKL BLAS to use
! with the RCI (P)FGMRES solver
!---------------------------------------------------------------------------
        real*8 DNRM2
        external DNRM2

        loop=.True.

    do while (loop)
    !---------------------------------------------------------------------------
    ! Compute the solution by RCI (P)FGMRES solver with preconditioning
    ! Reverse Communication starts here
    !---------------------------------------------------------------------------
    call DFGMRES(this%n, this%solution, this%rhs, this%rci_request, this%ipar, &
         & this%dpar, this%tmp)

         select case (this%rci_request)
    !---------------------------------------------------------------------------
    ! If rci_request=0, then the solution was found with the required precision
    !---------------------------------------------------------------------------
             case (0)
                 loop=.False.
    !---------------------------------------------------------------------------
    ! If rci_request=1, then compute the vector A*tmp(ipar(22))
    ! and put the result in vector tmp(ipar(23))
    !---------------------------------------------------------------------------
             case (1)
                call this%proc(this%content, this%tmp(this%ipar(22): this%ipar(22)+this%n-1), &
                               this%tmp(this%ipar(23): this%ipar(23)+this%n-1))
                loop=.True.
    !---------------------------------------------------------------------------
    ! If rci_request=2, then do the user-defined stopping test
    ! The residual stopping test for the computed solution is performed here
    !---------------------------------------------------------------------------
    ! NOTE: from this point vector B(N) is no longer containing the right-hand
    ! side of the problem! It contains the current FGMRES approximation to the
    ! solution. If you need to keep the right-hand side, save it in some other
    ! vector before the call to DFGMRES routine. Here we saved it in vector
    ! RHS(N). The vector B is used instead of RHS to preserve the original
    ! right-hand side of the problem and guarantee the proper restart of FGMRES
    ! method. Vector B will be altered when computing the residual stopping
    ! criterion!
    !---------------------------------------------------------------------------
              case (2)
        ! Request to the DFGMRES_GET routine to put the solution into B(N) via ipar(13)
                this%ipar(13)=1
        ! Get the current FGMRES solution in the vector B(N)
                call DFGMRES_GET(this%n, this%solution, b, this%rci_request, this%ipar, &
             & this%dpar, this%tmp, this%itercount)
        ! Compute the current true residual via MKL (Sparse) BLAS routines
                call this%proc(this%content, B, residual)

                residual=this%rhs-residual
                dvar=dnrm2(this%n, residual, 1)
                IF (dvar.lt. 1.0E-9) then
                   loop=.False.
                else
                   loop=.True.
                endif
    !---------------------------------------------------------------------------
    ! If rci_request=3, then apply the preconditioner on the vector
    ! tmp(ipar(22)) and put the result in vector tmp(ipar(23))
    !---------------------------------------------------------------------------
              case (3)
                call this%pcproc(this%content, this%tmp(this%ipar(22): this%ipar(22)+this%n-1), &
                               this%tmp(this%ipar(23): this%ipar(23)+this%n-1))
                loop=.True.
    !---------------------------------------------------------------------------
    ! If rci_request=4, then check if the norm of the next generated vector is
    ! not zero up to rounding and computational errors. The norm is contained
    ! in DPAR(7) parameter
    !---------------------------------------------------------------------------
              case (4)
                if (this%dpar(7).lt.1.0D-16) then
                   loop=.False.
                else
                   loop=.True.
                endif
    !---------------------------------------------------------------------------
    ! If rci_request=anything else, then DFGMRES subroutine failed
    ! to compute the solution vector: COMPUTED_SOLUTION(N)
    !---------------------------------------------------------------------------
             case default
                call throw_error(this%rci_request)

         end select

        enddo

    end subroutine solve

    subroutine get_real(this, solution)

        implicit none
        class(mkl_gmres),intent(inout) :: this
        real*8, intent(out) :: solution(this%n)
        integer :: i

!---------------------------------------------------------------------------
! Reverse Communication ends here
! Get the current iteration number and the FGMRES solution. (DO NOT FORGET to
! call DFGMRES_GET routine as computed_solution is still containing
! the initial guess!). Request to DFGMRES_GET to put the solution into
! vector COMPUTED_SOLUTION(N) via ipar(13)
!---------------------------------------------------------------------------
      this%ipar(13)=0
      call DFGMRES_GET(this%N, this%SOLUTION, this%RHS, this%rci_request, this%ipar, &
     & this%DPAR, this%tmp, this%ITERCOUNT)
      solution=this%solution

    end subroutine get_real

    subroutine get_complex(this, solution)

        implicit none
        class(mkl_gmres),intent(inout) :: this
        complex*16, intent(out) :: solution(this%n/2)
        real*8 :: b(this%n)
        real*8 :: residual(this%n), dvar
        integer :: i
        real*8 DNRM2
        external DNRM2
!---------------------------------------------------------------------------
! Reverse Communication ends here
! Get the current iteration number and the FGMRES solution. (DO NOT FORGET to
! call DFGMRES_GET routine as computed_solution is still containing
! the initial guess!). Request to DFGMRES_GET to put the solution into
! vector COMPUTED_SOLUTION(N) via ipar(13)
!---------------------------------------------------------------------------
!      this%ipar(13)=0
!      call DFGMRES_GET(this%N, this%SOLUTION, this%RHS, this%rci_request, this%ipar, &
!     & this%DPAR, this%tmp, this%ITERCOUNT)

        ! Request to the DFGMRES_GET routine to put the solution into B(N) via ipar(13)
                this%ipar(13)=0
        ! Get the current FGMRES solution in the vector B(N)
                call DFGMRES_GET(this%n, this%solution, b, this%rci_request, this%ipar, &
             & this%dpar, this%tmp, this%itercount)
        ! Compute the current true residual via MKL (Sparse) BLAS routines
              !  call this%proc(this%content, this%solution, residual)
              !  residual=this%rhs-residual
              !  dvar=dnrm2(this%n, residual, 1)

      call array_complex(this%solution, solution)

    end subroutine get_complex

    subroutine throw_error(rci_request)

        implicit none
        integer :: rci_request

        print *,'The solver has returned the ERROR code ', rci_request
        stop

    end subroutine throw_error

    end module mod_gmres

!    module ops
!  use mod_gmres, only: array_real, array_complex
!  use mod_solver
!  contains
!    subroutine   av(content, x, b)
!        implicit none
!        type(solver_type), intent(in) :: content
!        real*8, intent(in) :: x(:)
!        real*8, intent(out) :: b(:)
!        integer :: n=5
!!---------------------------------------------------------------------------
!! Define arrays for the upper triangle of the coefficient matrix
!! Compressed sparse row storage is used for sparse representation
!!---------------------------------------------------------------------------
!      integer IA(6)
!      data IA /1,3,6,9,12,14/
!      integer JA(13)
!      data JA    /  1,        3,            &
!     &              1,   2,        4,       &
!     &                   2,   3,        5,  &
!     &                        3,   4,   5,  &
!     &                             4,   5  /
!      real*8 A(13)
!      DATA A     / 1.0,     -1.0,               &
!     &             -1.0, 1.0,     -1.0,         &
!     &                    1.0,-2.0,      1.0,   &
!     &                        -1.0, 2.0,-1.0,   &
!     &                             -1.0,-3.0 /
!
!
!!---------------------------------------------------------------------------
!!Initialize variables and the right hand side through matrix-vector product
!!---------------------------------------------------------------------------
!      call MKL_DCSRGEMV('N', n, A, IA, JA, x, b)  !AX=RHS
!!---------------------------------------------------------------------------
!
!    end subroutine
!
!    subroutine av_complex(content, x, b)
!        implicit none
!        type(solver_type) :: content
!        real*8, intent(in) :: x(:)
!        real*8, intent(out) :: b(:)
!
!        complex*16, allocatable :: ArrayX(:)
!        complex*16, allocatable :: Arrayb(:)
!        integer :: n
!        integer :: i
!
!        n=size(x)/2
!
!        allocate(ArrayX(n), Arrayb(n))
!        call array_complex(x, ArrayX)
!
!        !! OP*ArrayX=Arrayb
!
!        call array_real(Arrayb, b)
!        deallocate(ArrayX, Arrayb)
!
!    end subroutine av_complex
!
!    subroutine   mv(content, x, b)
!        implicit none
!        type(solver_type), intent(in) :: content
!        real*8, intent(in) :: x(:)
!        real*8, intent(out) :: b(:)
!        integer, parameter :: n=5
!        real*8 :: diag(n)
!        REAL*8 :: XTMP(N)
!        integer :: i
!!---------------------------------------------------------------------------
!! Define arrays for the upper triangle of the coefficient matrix
!! Compressed sparse row storage is used for sparse representation
!!---------------------------------------------------------------------------
!      integer IA(6)
!      data IA /1,3,6,9,12,14/
!      integer JA(13)
!      data JA    /  1,        3,            &
!     &              1,   2,        4,       &
!     &                   2,   3,        5,  &
!     &                        3,   4,   5,  &
!     &                             4,   5  /
!      real*8 A(13)
!      DATA A     / 1.0,     -1.0,               &
!     &             -1.0, 1.0,     -1.0,         &
!     &                    1.0,-2.0,      1.0,   &
!     &                        -1.0, 2.0,-1.0,   &
!     &                             -1.0,-3.0 /
!      integer IA_L(6)
!      data IA_L /1,2,4,6,8,10/
!      integer JA_L(9)
!      data JA_L  /  1,                      &
!     &              1,   2,                 &
!     &                   2,   3,            &
!     &                        3,   4,       &
!     &                             4,   5  /
!      real*8 A_L(9)
!      DATA A_L   / 1.0,                         &
!     &             -1.0, 1.0,                   &
!     &                    1.0,-2.0,             &
!     &                        -1.0, 2.0,        &
!     &                             -1.0,-3.0 /
!      integer IA_u(6)
!      data IA_u /1,3,5,7,9,10/
!      integer JA_u(9)
!      data JA_u    /  1,        3,            &
!     &                   2,        4,       &
!     &                        3,        5,  &
!     &                             4,   5,  &
!     &                                  5  /
!      real*8 A_U(9)
!      DATA A_U     / 1.0,     -1.0,               &
!     &                    1.0,      -1.0,         &
!     &                         -2.0,      1.0,    &
!     &                               2.0,-1.0,    &
!     &                                   -3.0 /
!
!!---------------------------------------------------------------------------
!!Initialize variables and the right hand side through matrix-vector product
!!---------------------------------------------------------------------------
!      call MKL_DCSRGEMV('N', n, A_U, IA_U, JA_U, x, b)  !(D+U)X=RHS
!      xtmp=1.0/[1.0, 1.0, -2.0, 2.0, -3.0]*b            !D^(-1)(D+U)
!      call MKL_DCSRGEMV('N', n, A_L, IA_L, JA_L, xtmp, b)  !AX=RHS  !(D+L)D^(-1)(D+U)
!      !b=[1.0, 1.0, -2.0, 2.0, -3.0]*x
!!      call MKL_DCSRGEMV('N', n, A, IA, JA, x, b)  !AX=RHS  !(D+L)D^(-1)(D+U)
!
!!---------------------------------------------------------------------------
!
!    end subroutine
!
!    subroutine mv_complex(x, b)
!        implicit none
!        real*8, intent(in) :: x(:)
!        real*8, intent(out) :: b(:)
!        complex*16, allocatable :: ArrayX(:)
!        complex*16, allocatable :: Arrayb(:)
!        integer :: n
!        integer :: i
!
!        n=size(x)/2
!
!        allocate(ArrayX(n), Arrayb(n))
!        call array_complex(x, ArrayX)
!
!        !! OP*ArrayX=Arrayb
!
!        call array_real(Arrayb, b)
!        deallocate(ArrayX, Arrayb)
!
!    end subroutine mv_complex
!
!!*******************************************************************************
!!                              INTEL CONFIDENTIAL
!!   Copyright(C) 2005-2008 Intel Corporation. All Rights Reserved.
!!   The source code contained  or  described herein and all documents related to
!!   the source code ("Material") are owned by Intel Corporation or its suppliers
!!   or licensors.  Title to the  Material remains with  Intel Corporation or its
!!   suppliers and licensors. The Material contains trade secrets and proprietary
!!   and  confidential  information of  Intel or its suppliers and licensors. The
!!   Material  is  protected  by  worldwide  copyright  and trade secret laws and
!!   treaty  provisions. No part of the Material may be used, copied, reproduced,
!!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!!   in any way without Intel's prior express written permission.
!!   No license  under any  patent, copyright, trade secret or other intellectual
!!   property right is granted to or conferred upon you by disclosure or delivery
!!   of the Materials,  either expressly, by implication, inducement, estoppel or
!!   otherwise.  Any  license  under  such  intellectual property  rights must be
!!   express and approved by Intel in writing.
!!*******************************************************************************
!!  Content:
!!  Intel MKL RCI (P)FGMRES ((Preconditioned) Flexible Generalized Minimal
!!                                                       RESidual method) example
!!*******************************************************************************
!
!!---------------------------------------------------------------------------
!!  Example program for solving non-symmetric indefinite system of equations
!!  Fully advanced case: full functionality of RCI FGMRES solver is exploited
!!---------------------------------------------------------------------------
!     subroutine  test()
!
!      use mod_gmres
!
!      implicit none
!
!      include "mkl_rci.fi"
!
!      integer n, m
!      parameter(n=5)
!      integer :: i
!!---------------------------------------------------------------------------
!! Allocate storage for the ?par parameters and the solution/rhs/residual vectors
!!---------------------------------------------------------------------------
!      real*8 expected_solution(n)
!      data expected_solution /-1.0,1.0,0.0,1.0,-1.0/
!      real*8 computed_solution(n), rhs(n)
!      type(solver_type) :: content
!!---------------------------------------------------------------------------
!! Some additional variables to use with the RCI (P)FGMRES solver
!!---------------------------------------------------------------------------
!
!      type(mkl_gmres) :: gmres
!
!      print *,'--------------------------------------------------------'
!      print *,'The FULLY ADVANCED example of usage of RCI FGMRES solver'
!      print *,'   to solve a non-symmetric indefinite non-degenerate'
!      print *,'          algebraic system of linear equations'
!      print *,'--------------------------------------------------------'
!
!
!      call av(content, expected_solution, rhs)
!! Save the right-hand side in vector B for future use
!!---------------------------------------------------------------------------
!      !call dcopy(n, rhs, 1, b, 1)  ! b=rhs
!!---------------------------------------------------------------------------
!! Initialize the initial guess
!!---------------------------------------------------------------------------
!      do i=1, n
!         computed_solution(i)=1.0
!      enddo
!
!!---------------------------------------------------------------------------
!! Initialize the solver
!!---------------------------------------------------------------------------
!      !call gmres%initial(n, computed_solution, rhs)
!      gmres=mkl_gmres(n, computed_solution, rhs)
!      call gmres%check()
!      call gmres%setop(av)
!      call gmres%setpc(mv)
!      call gmres%solve()
!      call gmres%get(computed_solution)
!
!!---------------------------------------------------------------------------
!! print solution vector: COMPUTED_SOLUTION(N) and
!! the number of iterations: ITERCOUNT
!!---------------------------------------------------------------------------
!      print *, ''
!      print *,' The system has been SUCCESSFULLY solved'
!      print *, ''
!      print *,' The following solution has been obtained:'
!      DO I=1,N
!         write(*,'(A18,I1,A2,E10.3)') 'COMPUTED_SOLUTION(',I,')=',  &
!     & computed_solution(I)
!      ENDDO
!      print *, ''
!!      print *,' The expected solution is:'
!!      DO I=1,this%N
!!         write(*,'(A18,I1,A2,E10.3)') 'EXPECTED_SOLUTION(',I,')=', &
!!     & EXPECTED_SOLUTION(I)
!!      ENDDO
!!      print *, ''
!!      print *,' Number of iterations: ', this%ITERCOUNT
!    end   subroutine test
!
!end module
!
