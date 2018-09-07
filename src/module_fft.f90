include 'mkl_dfti.f90'

module mod_fft

    use  MKL_DFTI
    !use penf, only: R_P
    implicit none

    private

    integer :: status
    integer, parameter :: r_p=8

    type, public :: fft_type

        private
        type(DFTI_DESCRIPTOR), pointer :: handle
        integer :: rstrides(3)
        integer :: cstrides(3)
        logical :: isInitial
        logical :: isRZ

        Contains

          !Procedure :: Initial
          procedure :: Finalize
          procedure :: HasInitial

          generic :: fft=> fft_rz, fft_zz, fft_zz_inplace, fft_rz_2d
          generic :: ifft=> ifft_zr, ifft_zz, ifft_zz_inplace, ifft_zr_2d
          generic :: Initial=> Initial_1d, initial_2d

          Procedure, private :: Initial_1d
          procedure, private :: initial_2d
          procedure, private :: fft_rz
          procedure, private :: fft_zz
          procedure, private :: fft_zz_inplace
          procedure, private :: fft_rz_2d
          procedure, private :: ifft_zr
          procedure, private :: ifft_zz
          procedure, private :: ifft_zr_2d
          procedure, private :: ifft_zz_inplace

    end type fft_type


    contains

    function HasInitial(this) result(Has)

        implicit none
        class(fft_type), intent(inout) :: this
        logical :: Has

        has=this%isInitial

    end function HasInitial

    subroutine initial_1d(this, n)

        implicit none
        class(fft_type),intent(inout) :: this
        integer, intent(in) :: n

        call this%finalize()

        status = DftiCreateDescriptor(this%handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, n)
        status = DftiSetValue(this%handle, DFTI_BACKWARD_SCALE, 1.0_R_P/real(n, R_P))

        status = DftiCommitDescriptor(this%handle)

        this%isInitial=.True.

    end subroutine initial_1d

    subroutine initial_2d(this, n1, n2)

        implicit none
        class(fft_type),intent(inout) :: this
        integer, intent(in) :: n1, n2

        call this%finalize()

        status = DftiCreateDescriptor(this%handle, DFTI_DOUBLE, DFTI_REAL, 2, [n1, n2])
        status = DftiSetValue(this%handle, DFTI_FORWARD_SCALE, 1.0_R_P/real(n1*n2, R_P))
        status = DftiSetValue(this%handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        status = DftiSetValue(this%handle, DFTI_CONJUGATE_EVEN_STORAGE,                   &
          &                    DFTI_COMPLEX_COMPLEX)

        status = DftiCommitDescriptor(this%handle)

        this%rstrides = [0, 1, n1]
        this%cstrides = [0, 1, INT(n1/2)+1]


        status = DftiSetValue(this%handle, DFTI_INPUT_STRIDES, this%rstrides)
        status = DftiSetValue(this%handle, DFTI_OUTPUT_STRIDES, this%cstrides)
        status = DftiCommitDescriptor(this%handle)
        this%isRZ=.True.
        this%isInitial=.True.

    end subroutine initial_2d

    subroutine fft_rz_2d(this, Ain, Aout)

      implicit none
      class(fft_type), intent(inout) :: this
      real(R_P), intent(in) :: Ain(:, :)
      complex(R_P), intent(out) :: Aout(:, :)

      if(.not. this%isRZ) then
        status = DftiSetValue(this%handle, DFTI_INPUT_STRIDES, this%rstrides)
        status = DftiSetValue(this%handle, DFTI_OUTPUT_STRIDES, this%cstrides)
        status = DftiCommitDescriptor(this%handle)
        this%isRZ=.True.
      endif

      if(size(Ain, dim=1) /= this%rstrides(3)) then
          write(*,*) 'the first rank in real Ain is wrong: n1=',  size(Ain, dim=1), &
      &   'not', this%rstrides(3)
          stop
      endif
      if(size(Aout, dim=1) /= this%cstrides(3)) then
          write(*,*) 'the first rank in complex Aout is wrong: n1=',  size(Aout, dim=1), &
      &   'not', this%cstrides(3)
          stop
      endif
      status = DftiComputeForward(this%handle, Ain(:,1), Aout(:,1))

    end subroutine fft_rz_2d

    subroutine ifft_zr_2d(this, Ain, Aout)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(R_P), intent(in) :: Ain(:, :)
      real(R_P), intent(out) :: Aout(:, :)

      if(this%isRZ) then
        status = DftiSetValue(this%handle, DFTI_INPUT_STRIDES, this%cstrides)
        status = DftiSetValue(this%handle, DFTI_OUTPUT_STRIDES, this%rstrides)
        status = DftiCommitDescriptor(this%handle)
        this%isRZ=.False.
      endif

      if(size(Ain, dim=1) /= this%cstrides(3)) then
          write(*,*) 'the first rank in complex Ain is wrong: n1=',  size(Ain, dim=1), &
      &   'not', this%cstrides(3)
          stop
      endif
      if(size(Aout, dim=1) /= this%rstrides(3)) then
          write(*,*) 'the first rank in real Aout is wrong: n1=',  size(Aout, dim=1), &
      &   'not', this%rstrides(3)
          stop
      endif
      status = DftiComputeBackward(this%handle, Ain(:,1), Aout(:,1))

    end subroutine ifft_zr_2d

    subroutine fft_zz(this, Ain, Aout)

        implicit none
        class(fft_type), intent(inout) :: this
        complex(R_P), intent(in) :: Ain(:)
        complex(R_P), intent(out) :: Aout(size(Ain))

        Aout=Ain
        status = DftiComputeBackward(this%handle, Aout)

    end subroutine fft_zz

    subroutine ifft_zz(this, Ain, Aout)

        implicit none
        class(fft_type), intent(inout) :: this
        complex(R_P), intent(in) :: Ain(:)
        complex(R_P), intent(out) :: Aout(size(Ain))

        Aout=Ain
        status = DftiComputeForward(this%handle, Aout)

    end subroutine ifft_zz

    subroutine fft_zz_inplace(this, Ain)

        implicit none
        class(fft_type), intent(inout) :: this
        complex(R_P), intent(inout) :: Ain(:)

        status = DftiComputeBackward(this%handle, Ain)

    end subroutine fft_zz_inplace

    subroutine ifft_zz_inplace(this, Ain)

        implicit none
        class(fft_type), intent(inout) :: this
        complex(R_P), intent(inout) :: Ain(:)

        status = DftiComputeForward(this%handle, Ain)

    end subroutine ifft_zz_inplace

    subroutine fft_rz(this, Ain, Aout)

        implicit none
        class(fft_type), intent(inout) :: this
        real(R_P), intent(in) :: Ain(:)
        complex(R_P), intent(out) :: Aout(size(Ain))

        Aout=Ain
        status = DftiComputeBackward(this%handle, Aout)

    end subroutine fft_rz

    subroutine ifft_zr(this, Ain, Aout)

        implicit none
        class(fft_type), intent(inout) :: this
        complex(R_P), intent(in) :: Ain(:)
        real(R_P), intent(out) :: Aout(size(Ain))
        complex(R_P) :: tmp(size(Ain))

        tmp=Ain
        status = DftiComputeForward(this%handle, tmp)
        Aout=real(tmp, R_P)

    end subroutine ifft_zr

    subroutine finalize(this)
        implicit none
        class(fft_type), intent(inout) :: this

        status = DftiFreeDescriptor(this%handle)
        this%isInitial=.False.

    end subroutine finalize

end module mod_fft
