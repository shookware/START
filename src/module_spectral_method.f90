!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_spectral_method.f90
!> @file
!> @breif 谱方法导数计算子程序.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: mod_spectral_method
!> @breif 谱方法求导模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2018-07-08 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2018-07-08
!------------------------------------------------------------------------------
module mod_spectral_method

  use penf, only: R_P
  use mod_parameter, only: CPLI
  use mod_fft
  implicit none

  public :: spectral_method_type
  private

  logical :: checkflag=.False.
  integer, save :: icheck=0

  type spectral_method_type

  !  real(R_P), private :: alpha0=0.0d0
    integer, private :: mdim
    integer, private :: ndim
    integer, private :: msp
    integer, private :: nsp
    integer, private :: nsp4
    integer, private :: msp4

    real(R_P), private :: omega0=0.0d0
    real(R_P), private :: beta0=0.0d0
    type(fft_type), private :: fft_solver
    logical, public :: isCreate=.False.

  contains

!    generic :: dx=> dx_real
    generic :: dz=> dz_real, dz_complex
    generic :: dt=> dt_real, dt_complex
    generic :: dzz=> dzz_real, dzz_complex
    procedure :: ForTrans=>ForTrans_foo
    procedure :: BackTrans=>BackTrans_foo

    procedure :: Create
    procedure :: Finalize

    !procedure, private :: dx_real
    procedure, private :: dz_real
    procedure, private :: dt_real
    procedure, private :: dz_complex
    procedure, private :: dt_complex
    procedure, private :: dzz_real
    procedure, private :: dzz_complex

  end type spectral_method_type

contains

  subroutine Create(this, omega, beta, m, n)

    implicit none
    class(spectral_method_type), intent(inout) :: this
    real(R_P), intent(in) :: beta
    real(R_P), intent(in) :: omega
    integer, intent(in) :: m, n

    this%beta0=beta
    this%omega0=omega
    this%mdim=m
    this%ndim=n
    if(this%mdim .ne. 0) then
      this%msp=(m+1)*2
      this%msp4=this%msp*4
    else
      this%msp=1
      this%msp4=1
    endif
    if(this%ndim .ne. 0) then
      this%nsp=(n+1)*2
      this%nsp4=this%nsp*4
    else
      this%nsp=1
      this%nsp4=1
    endif
    !
    ! print*, 'm=', m, "n=", n
    ! print*, 'this%msp4=', this%msp4, 'this%nsp4=', this%nsp4
    ! pause

    call this%fft_solver%initial(this%msp4, this%nsp4)

    this%isCreate=.True.

  end subroutine Create

  ! elemental function dx_real(this, array) result(dx)
  !
  !   implicit none
  !   class(spectral_method_type), intent(in) :: this
  !   real(R_P), intent(in) :: array
  !   real(R_P) :: dx
  !
  !     dx=CPLI*dble(1)*this%alpha0*array
  !
  ! end function  dx_real

  pure function dz_real(this, array) result(dz)

    implicit none
    class(spectral_method_type), intent(in) :: this
    real(R_P), intent(in) :: array(:)
    real(R_P), allocatable :: dz(:)

    integer :: dim
    integer :: i

    dim=(size(array)-1)/2
    allocate(dz(-dim:dim))
    dz=0.0d0

    do i=-dim, dim
      dz(i)=CPLI*dble(i)*this%beta0*array(i+dim+1)
    enddo

  end function  dz_real

  pure function dzz_real(this, array) result(dzz)

    implicit none
    class(spectral_method_type), intent(in) :: this
    real(R_P), intent(in) :: array(:)
    real(R_P), allocatable :: dzz(:)

    integer :: dim
    integer :: i

    dim=(size(array)-1)/2
    allocate(dzz(-dim:dim))
    dzz=0.0d0

    do i=-dim, dim
      dzz(i)=-(dble(i)*this%beta0)**2*array(i+dim+1)
    enddo

  end function  dzz_real

  pure function dz_complex(this, array) result(dz)

    implicit none
    class(spectral_method_type), intent(in) :: this
    complex(R_P), intent(in) :: array(:)
    complex(R_P), allocatable :: dz(:)

    integer :: dim
    integer :: i

    dim=(size(array)-1)/2
    allocate(dz(-dim:dim))
    dz=0.0d0

    do i=-dim, dim
      dz(i)=CPLI*dble(i)*this%beta0*array(i+dim+1)
    enddo

  end function  dz_complex

  pure function dzz_complex(this, array) result(dzz)

    implicit none
    class(spectral_method_type), intent(in) :: this
    complex(R_P), intent(in) :: array(:)
    complex(R_P), allocatable :: dzz(:)

    integer :: dim
    integer :: i

    dim=(size(array)-1)/2
    allocate(dzz(-dim:dim))
    dzz=0.0d0

    do i=-dim, dim
      dzz(i)=-(dble(i)*this%beta0)**2*array(i+dim+1)
    enddo

  end function  dzz_complex

  pure function dt_real(this, array) result(dt)

    implicit none
    class(spectral_method_type), intent(in) :: this
    real(R_P), intent(in) :: array(0:)
    real(R_P), allocatable :: dt(:)

    integer :: dim
    integer :: i

    dim=(size(array)-1)
    allocate(dt(0:dim))
    dt=0.0d0

    do i=0, dim
      dt(i)=-CPLI*dble(i)*this%omega0*array(i)
    enddo

  end function  dt_real

  pure function dt_complex(this, array) result(dt)

    implicit none
    class(spectral_method_type), intent(in) :: this
    complex(R_P), intent(in) :: array(0:)
    complex(R_P), allocatable :: dt(:)

    integer :: dim
    integer :: i

    dim=(size(array)-1)
    allocate(dt(0:dim))
    dt=0.0d0

    do i=0, dim
      dt(i)=-CPLI*dble(i)*this%omega0*array(i)
    enddo

  end function  dt_complex

  subroutine BackTrans_foo(this, x_in, BackTrans)

    implicit none
    class(spectral_method_type), intent(inout) :: this
    complex(R_P), intent(in) :: x_in(:, :)
    !real(R_P), allocatable :: BackTrans(:, :)
    real(R_P), intent(out) :: BackTrans(0:this%msp4-1, 0:this%nsp4-1)
    !integer :: mdim, ndim
    complex(R_P) :: x2(0:this%msp4/2, 0:this%nsp4-1)
    complex(R_P) :: xt(0:this%mdim+1, -this%ndim: this%ndim)
    integer :: m, n, dim2(2)
    !integer :: m, n

    !mdim=(this%mdim+1)*4
    !ndim=(this%ndim+1)*4
    x2=0.0d0
    xt=0.0d0
    !allocate(BackTrans(0:mdim, 0:ndim*2-1))
    BackTrans=0.0d0
    !
    ! x2(0, 0:this%ndim)=x_in(0, 0:this%ndim)
    ! x2(0, ndim*4-1:ndim*4-this%ndim:-1)=x_in(0, -1:-this%ndim:-1)
    ! x2(1:this%mdim,        1:          this%ndim)=conjg(x_in(1:this%mdim, -1:-this%ndim:-1))
    ! x2(1:this%mdim, ndim*2-1:ndim*2-this%ndim:-1)=conjg(x_in(1:this%mdim,  1:this%ndim))
    ! x2(1:this%mdim, 0)=conjg(x_in(1:this%mdim, 0))
    dim2=shape(x_in)
    m=dim2(1); n=dim2(2)
    !
    ! print*, shape(x_in)
    ! print*, m, n
    ! print*, this%mdim, this%ndim

    xt(0:0, :) = x_in(1:1, :)
    xt(1:this%mdim, 0)= conjg(x_in(2:this%mdim+1, this%ndim+1))
    xt(1:this%mdim, 1:this%ndim)= conjg(x_in(2:this%mdim+1, this%ndim:1:-1))
    xt(1:this%mdim, -1:-this%ndim:-1)=conjg(x_in(2:this%mdim+1, this%ndim+2:n))
    ! print*, this%mdim, this%ndim
    ! print*, this%ndim+2, n
    ! pause
    ! print*, 'numbers....'
    ! print*, 'xt...1', 1, this%mdim, 1, this%ndim
    ! print*, xt(1:this%mdim, 1:this%ndim)
    ! print*, 'x_in....1', this%mdim, this%ndim
    ! print*, x_in(2:this%mdim+1, this%ndim:1:-1)
    ! print*, 'xt...2', this%mdim, this%ndim
    ! print*, xt(1:this%mdim, -1:-this%ndim:-1)
    ! print*, 'x_in....2', this%mdim, n-(this%ndim+2)+1
    ! print*, x_in(2:this%mdim+1, this%ndim+2:n)
    !
    ! pause 'numbers....'
    !

    x2(0:this%mdim, 0: this%ndim)= xt(0:this%mdim, 0:this%ndim)
    x2(0:this%mdim, this%nsp4-this%ndim:this%nsp4-1)=xt(0:this%mdim, -this%ndim:-1)


    !!!!!!!!!!TEST!!!!!!!!!

    !x2=0
    !x2(1, 0)=(0.0d0, 0.5d0)

!    print*, 'x2....'
!    print*, maxval(abs(x2))


!    pause

    call this%fft_solver%ifft(x2, BackTrans)

!     block
!
!       integer :: m, n
!       if(icheck/1000*1000-icheck==0)then
!       if( checkflag==.true.)then
!     do m=0, 0
!       do n=0, this%nsp4-1
!     print*, m, n, Backtrans(m, n)
!   enddo
! enddo
! pause
! endif
! endif
! endblock

!    print*, maxval(abs(BackTrans))
!    pause

    !!!!!!!!TEST END!!!!!!!!!!
    !stop

  end subroutine BackTrans_foo


  subroutine ForTrans_foo(this, x_in, ForTrans)

    implicit none
    class(spectral_method_type), intent(inout) :: this
    real(R_P), intent(in) :: x_in(0:this%msp4-1, 0:this%nsp4-1)
    !complex(R_P), allocatable :: ForTrans(:, :)
    complex(R_P), intent(inout) :: Fortrans(0: this%mdim, -this%ndim:this%ndim)
    complex(R_P) :: xt(0:this%mdim+1, -this%ndim:this%ndim)
    !integer :: mdim, ndim
    complex(R_P) :: x2(0:this%msp4/2, 0:this%nsp4-1)

    !### debug
    !integer :: m, n

    !mdim=(this%mdim+1)*2
    !ndim=(this%ndim+1)*2
    x2=0.0d0
    !allocate(ForTrans(0:this%mdim, -this%ndim:this%ndim))
    ForTrans=0.0d0

    ! print*, x_in
    ! pause

    call this%fft_solver%fft(x_in, x2)

    !
    ! ForTrans(0, 0:this%ndim)=x2(0, 0:this%ndim)
    ! ForTrans(0, -1:-this%ndim:-1)=x2(0, ndim*2-1:ndim*2-this%ndim:-1)
    ! ForTrans(1:this%mdim, -1:-this%ndim:-1)=&
    !   & conjg(x2(1:this%mdim, 1:this%ndim))
    ! ForTrans(1:this%mdim, 1:this%ndim)= &
    !   & conjg(x2(1:this%mdim, ndim*2-1:ndim*2-this%ndim:-1))
    ! ForTrans(1:this%mdim, 0)=conjg(x2(1:this%mdim, 0))

    xt=0.0d0
    xt(0: this%mdim, 0: this%ndim)= x2(0:this%mdim, 0:this%ndim)
    xt(0: this%mdim, -this%ndim: -1)= &
    &   x2(0: this%mdim, this%nsp4-this%ndim:this%nsp4-1)

    ForTrans(0, :) = xt(0, :)
    ForTrans(1:this%mdim, 0)= conjg(xt(1:this%mdim, 0))
    ForTrans(1:this%mdim, -1:-this%ndim:-1)= conjg(xt(1:this%mdim, 1:this%ndim))
    ForTrans(1:this%mdim, 1:this%ndim)=conjg(xt(1:this%mdim, -1:-this%ndim:-1))


  !  block
  !  use mod_debug, only: m=>mloc, n=>nloc, file_log, jloc, lloc
  !  if(ii>(150*301+123)*5 .and. mod(ii, 301*5)>123*5 .and. mod(ii, 301*5)<=124*5)then

    !   if((jloc<=10 .or. jloc>290 .or. (jloc==123)) .and. lloc .ne. 4)then
    !
    !     write(file_log, *)'j=', jloc,  'l=', lloc
    !     write(file_log, *)'xin......'
    !     do m=0, this%msp4-1
    !       write(file_log, *)m, 0, x_in(m, 0)
    !     enddo
    !   ! write(file_log, *)'x2.........'
    !   ! do m=0, this%msp4/2
    !   !   write(file_log, *)m, 0, x2(m ,0)
    !   ! enddo
    !   ! write(file_log, *)'xt.....'
    !   ! do m=0, this%mdim
    !   !   write(file_log, *)m, 0, xt(m ,0)
    !   ! enddo
    !   write(file_log, *)'fortrans..'
    !   do m=0, this%mdim
    !     write(file_log, *)m, 0, fortrans(m, 0)
    !   enddo
    !   endif
    !
    ! end block



!     checkflag=.true.
!
!     block
!       integer :: m, n
!       icheck=icheck+1
!       if(icheck/1000*1000-icheck==0)then
!       print*, icheck
!       do n=0, this%nsp4-1
!         print*, 0, x_in(0, n)
!       end do
!       pause
!     endif
!
!     endblock
    ! print*, maxval(abs(x_in))
    if(maxval(abs(Fortrans))>=15e0) then
    print*, maxval(abs(Fortrans)), maxloc(abs(Fortrans))
    print*, 'xt...'
    print*, xt
    print*, 'x_in..'
    block
      integer :: m, n
    do m=0, this%mdim
      do n=0, this%nsp4-1
    print*, m, n, x_in(m, n)
  enddo
enddo

stop
endblock
  !  pause
!
    endif

  end subroutine ForTrans_foo

  subroutine Finalize(this)

    implicit none
    class(spectral_method_type), intent(inout) :: this

    !this%alpha0=0.0d0
    this%beta0=0.0d0
    this%omega0=0.0d0
    this%isCreate=.False.

  end subroutine Finalize
end module mod_spectral_method
