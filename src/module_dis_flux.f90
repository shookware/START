!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_dis_flux.f90
!> @file
!> @breif 扰动形函数通量点模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: dis_flux
!> @breif 扰动形函数通量点模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-02 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-02
!------------------------------------------------------------------------------
module mod_dis_flux

    use mod_vector_cmplx
    use penf, only: R_P
    implicit none

    private
    public :: DIS_FLUX_NULL, DIS_FLUXIJ_NULL, fft, ifft

    complex(R_P), parameter :: ZERO=(0.0d0, 0.0d0)

    !> 扰动形函数通量点\f$ {\widehat\phi}^\prime=\left(\widehat{\rho}^\prime, \widehat {u}^\prime,
    !! \widehat {v}^\prime, \widehat {w}^\prime, \widehat {T}^\prime\right)\f$类.
    type, public :: dis_flux_type

        private
        complex(R_P), private :: rho !< 密度扰动形函数\f$\widehat{\rho}^\prime\f$
        type(vector_cmplx_type), private :: Vel !<速度扰动形函数\f$\left(\widehat{u}^\prime,
                                        !! \widehat{v}^\prime, \widehat{w}^\prime\right)\f$
        complex(R_P), private :: T !< 温度扰动形函数\f$\widehat{T}^\prime\f$

        Contains

          procedure :: GetVel   !< 获得速度扰动形函数点
          procedure :: Print    !< 输出扰动
          generic :: Set => set_dis_flux, set_dis_flux_vel, set_dis_flux_array !< 设置扰动形函数通量点
          generic :: Get => get_dis_flux, get_dis_flux_vel, get_dis_flux_array !< 获得扰动形函数通量点
          generic :: operator(+) => add_type !< 扰动形函数相加

!          Final     ::
          procedure, private :: add_type     !< 形函数相加
          procedure, private :: get_dis_flux  !< 获得扰动形函数
          procedure, private :: set_dis_flux  !< 设置扰动形函数
          procedure, private :: set_dis_flux_vel !< 设置扰动形函数(速度以分量形式给出)
          procedure, private :: get_dis_flux_vel !< 获得扰动形函数(速度以分量形式给出)
          procedure, private :: set_dis_flux_array !< 设置扰动形函数(以数组形式给出)
          procedure, private :: get_dis_flux_array !< 设置扰动形函数(以数组形式给出)

    end type dis_flux_type

    !> 扰动形函数通量点\f$ {\widehat\phi}^\prime=\left(\widehat{\rho}^\prime, \widehat {u}^\prime,
    !! \widehat {v}^\prime, \widehat {w}^\prime, \widehat {T}^\prime\right)\f$类(带序号).
    type, public, extends(dis_flux_type) :: dis_flux_ij_type

        private
        integer, private :: i  !< 流向坐标序号
        integer, private :: j  !< 法向坐标序号

        Contains

          Procedure :: SetIJ => set_dis_flux_ij !< 设置序号
          procedure :: GetIJ => get_dis_flux_ij !< 获得序号
          procedure :: GetI => get_dis_flux_i  !< 获得流向坐标序号
          generic :: operator(+) => add_typeij, add_array_type, add_type_array   !< 形函数通量相加
          generic :: operator(-) => minus_typeij !< 形函数通量相减
          generic :: operator(*) => multiply_double_typeij, &
          &                         multiply_typeij_double, &
                                    multiply_cmplx_typeij,  &
                                    multiply_typeij_cmplx, &
                                    multiply_dblearray_typeij, &
                                    multiply_dblearray3_typeij, &
                                    multiply_cmplxarray_typeij, &
                                    multiply_cmplxarray3_typeij  !< 形函数通量乘法(标量乘和点乘矢量)
          generic :: operator(.mx.) => multiply_cmplxmatrix_typeij, &
          &                            multiply_dblematrix_typeij !< 形函数通量乘法(矩阵乘)
          generic :: operator(.div.) => div_cmplxmatrix_typeij   !< 形函数通量除法
          generic :: assignment(=) => assignment_typeij

          procedure, pass(obj1), private :: add_typeij  !< 形函数通量相加
          procedure, pass(obj1), private :: minus_typeij !< 形函数通量相减
          procedure, pass(this), private :: multiply_double_typeij !< 标量乘
          procedure, pass(this), private :: multiply_typeij_double !< 标量乘
          procedure, pass(this), private :: multiply_cmplx_typeij !< 标量乘
          procedure, pass(this), private :: multiply_typeij_cmplx !< 标量乘
          procedure, pass(this), private :: multiply_dblearray_typeij !<数组乘
          procedure, pass(this), private :: multiply_dblearray3_typeij !< 数组乘
          procedure, pass(this), private :: multiply_cmplxarray_typeij !< 数组乘
          procedure, pass(this), private :: multiply_cmplxarray3_typeij !< 数组乘
          procedure, pass(this), private :: multiply_dblematrix_typeij !<矩阵乘
          procedure, pass(this), private :: multiply_cmplxmatrix_typeij !< 矩阵乘
          procedure, pass(this), private :: div_cmplxmatrix_typeij !<矩阵逆乘形函数
!          procedure, pass(This), private :: multiply_disop_typeij
          procedure, pass(this), private :: assignment_typeij
          procedure, pass(obj1), private :: add_type_array !< 形函数相加
          procedure, pass(obj1), private :: add_array_type !< 形函数相加

    end type dis_flux_ij_type

    type(dis_flux_type), parameter :: DIS_FLUX_NULL= &
                        &   dis_flux_type(ZERO, VEC_CMPLX_NULL, ZERO) !< 空扰动形函数
    type(dis_flux_ij_type), parameter :: DIS_FLUXIJ_NULL= &
                        &   dis_flux_ij_type(ZERO, VEC_CMPLX_NULL, ZERO, 0, 0) !< 空扰动形函数

    interface fft
       module procedure fft_zz
       module procedure FFT_rz
    end interface

    interface ifft
       module procedure ifft_zz
       module procedure ifft_zr
    end interface

    contains

    !> 获得扰动速度形函数
    !! @return 速度扰动形函数
    function GetVel(this) Result(Vel)

        implicit none
        class(dis_flux_type), intent(in) :: this
        type(vector_cmplx_type) :: Vel

        Vel=this%Vel

    end function GetVel

    !> 设置扰动形函数
    !! @param[in] rho 密度扰动形函数
    !! @param[in] Vel 速度扰动形函数
    !! @param[in] T 温度扰动形函数
    elemental subroutine set_dis_flux(this, rho, Vel, T)

        implicit none
        complex(R_P), intent(in) :: rho, T
        type(vector_cmplx_type), intent(in) :: Vel
        class(dis_flux_type), intent(inout) :: this

        this%rho=rho; this%Vel=Vel; this%T=T

    end subroutine set_dis_flux

    !> 设置扰动形函数
    !! @param[in] rho 密度扰动形函数
    !! @param[in] u 流向速度扰动形函数
    !! @param[in] v 法向速度扰动形函数
    !! @param[in] w 展向速度扰动形函数
    !! @param[in] T 温度扰动形函数
    subroutine set_dis_flux_vel(this, rho, u, v, w, T)

        implicit none
        complex(R_P), intent(in) :: rho, T, u, v, w
        type(vector_cmplx_type) :: Vel
        class(dis_flux_type), intent(inout) :: this

        call Vel%set(u, v, w)
        this%rho=rho; this%Vel=Vel; this%T=T

    end subroutine set_dis_flux_vel

    !> 设置扰动形函数
    !! @param[in] flux 扰动形函数
    subroutine set_dis_flux_array(this, flux)

        implicit none
        class(dis_flux_type), intent(inout) :: this
        complex(R_P), intent(in) :: flux(5)
        type(vector_cmplx_type) :: Vel

        call Vel%set(flux(2), flux(3), flux(4))
        this%rho=flux(1); this%Vel=Vel; this%T=flux(5)

    end subroutine set_dis_flux_array

    !> 获得扰动形函数
    !! @param[out] rho 密度扰动形函数
    !! @param[out] Vel 速度扰动形函数
    !! @param[out] T 温度扰动形函数
    elemental subroutine get_dis_flux(this, rho, Vel, T)

        implicit none
        class(dis_flux_type), intent(in) :: this
        complex(R_P), intent(inout) :: rho, T
        type(vector_cmplx_type), intent(inout) :: Vel

        rho=this%rho; Vel=this%Vel; T=this%T

    end subroutine get_dis_flux

    !> 获得速度扰动形函数
    !! @param[out] rho 密度扰动形函数
    !! @param[out] u 流向速度扰动形函数
    !! @param[out] v 法向速度扰动形函数
    !! @param[out] w 展向速度扰动形函数
    !! @param[out] T 温度扰动形函数
    pure subroutine get_dis_flux_vel(this, rho, u, v, w, T)

        implicit none
        class(dis_flux_type), intent(in) :: this
        complex(R_P), intent(inout) :: rho, T, u, v, w
        type(vector_cmplx_type) :: Vel

        rho=this%rho; Vel=this%Vel; T=this%T
        call Vel%get(U, v, w)

    end subroutine get_dis_flux_vel

    !> 获得速度扰动形函数
    !! @param[out] flux 扰动形函数
    pure subroutine get_dis_flux_array(this, flux)

        implicit none
        class(dis_flux_type), intent(in) :: this
        complex(R_P), intent(inout) :: flux(5)
        type(vector_cmplx_type) :: Vel

        flux(1)=this%rho; Vel=this%Vel; flux(5)=this%T
        call Vel%get(flux(2), flux(3), flux(4))

    end subroutine get_dis_flux_array

    !> 设置扰动形函数通量的序号
    !! @param[in] i 流向序号
    !! @param[in] j 法向序号
    subroutine set_dis_flux_ij(this, i, j)

        implicit none
        integer, intent(in) :: i, j
        class(dis_flux_ij_type), intent(inout) :: this

        this%i=i; this%j=j

    end subroutine set_dis_flux_ij

    !> 获得扰动形函数通量的序号
    !! @param[out] i 流向序号
    !! @param[out] j 法向序号
    subroutine get_dis_flux_ij(this, i, j)

        implicit none
        integer, intent(inout) :: i, j
        class(dis_flux_ij_type), intent(in) :: this

        i=this%i; j=this%j

    end subroutine get_dis_flux_ij

    !> 获得扰动形函数通量的流向序号
    !! @param[out] i 流向序号
    subroutine get_dis_flux_i(this, i)

        implicit none
        integer, intent(inout) :: i
        class(dis_flux_ij_type), intent(in) :: this

        i=this%i

    end subroutine get_dis_flux_i

    !> 标量乘 \f$\widehat\varphi^\prime=\alpha\widehat\phi^\prime\f$
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return 乘的结果\f$\widehat\varphi^\prime\f$
    function multiply_double_type(scalar, this) result(mul)

        implicit none
        class(dis_flux_type), intent(in) :: this
        real(R_P), intent(in) :: scalar
        type(dis_flux_type) :: mul

        mul%rho=scalar*this%rho
        mul%Vel=scalar*this%Vel
        mul%T  =scalar*this%T

    end function multiply_double_type

    !> 标量乘 \f$\widehat\varphi^\prime=\widehat\phi^\prime\alpha\f$
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return 乘的结果\f$\widehat\varphi^\prime\f$
    type(dis_flux_type) function multiply_type_double(this, scalar)

        implicit none
        class(dis_flux_type), intent(in) :: this
        real(R_P), intent(in) :: scalar

        multiply_type_double%rho=scalar*this%rho
        multiply_type_double%Vel=scalar*this%Vel
        multiply_type_double%T  =scalar*this%T

    end function multiply_type_double

    !> 标量乘 \f$\widehat\varphi^\prime=\alpha\widehat\phi^\prime\f$
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return 乘的结果\f$\widehat\varphi^\prime\f$
    type(dis_flux_ij_type) function multiply_double_typeij(scalar, this)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        real(R_P), intent(in) :: scalar

        multiply_double_typeij%rho=scalar*this%rho
        multiply_double_typeij%Vel=scalar*this%Vel
        multiply_double_typeij%T  =scalar*this%T
        multiply_double_typeij%i  =this%i
        multiply_double_typeij%j  =this%j

    end function multiply_double_typeij

    !> 标量乘 \f$\widehat\varphi^\prime=\alpha\widehat\phi^\prime\f$
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return 乘的结果\f$\widehat\varphi^\prime\f$
    type(dis_flux_ij_type) function multiply_cmplx_typeij(scalar, this)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        complex(R_P), intent(in) :: scalar

        multiply_cmplx_typeij%rho=scalar*this%rho
        multiply_cmplx_typeij%Vel=scalar*this%Vel
        multiply_cmplx_typeij%T  =scalar*this%T
        multiply_cmplx_typeij%i  =this%i
        multiply_cmplx_typeij%j  =this%j

    end function multiply_cmplx_typeij

    !> 标量乘 \f$\widehat\varphi^\prime=\widehat\phi^\prime\alpha\f$
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return 乘的结果\f$\widehat\varphi^\prime\f$
    type(dis_flux_ij_type) function multiply_typeij_double(this, scalar)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        real(R_P), intent(in) :: scalar

        multiply_typeij_double%rho=scalar*this%rho
        multiply_typeij_double%Vel=scalar*this%Vel
        multiply_typeij_double%T  =scalar*this%T
        multiply_typeij_double%i  =this%i
        multiply_typeij_double%j  =this%j

    end function multiply_typeij_double

    !> 标量乘 \f$\widehat\varphi^\prime=\widehat\phi^\prime\alpha\f$
    !! @param[in] scalar 标量\f$\alpha\f$
    !! @return 乘的结果\f$\widehat\varphi^\prime\f$
    type(dis_flux_ij_type) function multiply_typeij_cmplx(this, scalar)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        complex(R_P), intent(in) :: scalar

        multiply_typeij_cmplx%rho = scalar*this%rho
        multiply_typeij_cmplx%Vel = scalar*this%Vel
        multiply_typeij_cmplx%T   = scalar*this%T
        multiply_typeij_cmplx%i   = this%i
        multiply_typeij_cmplx%j   = this%j

    end function multiply_typeij_cmplx

    !> 向量乘(内积) \f$c^\prime=\overrightarrow a \cdot \widehat\phi^\prime\f$
    !! @param[in] scalar 向量\f$\overrightarrow a\f$
    !! @return 乘的结果\f$c^\prime\f$
    complex(R_P) function multiply_dblearray_typeij(array, this)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        real(R_P), dimension(5), intent(in) :: array

        multiply_dblearray_typeij=array(1)*this%rho + &
                                  array(2:4)*this%Vel +&
                                  array(5)*this%T


    end function multiply_dblearray_typeij

    !> 同时做向量乘(内积) \f$c_i^\prime=\overrightarrow a_i \cdot \widehat\phi^\prime, i=1,2,3\f$
    !! @param[in] scalar 向量\f$\overrightarrow a\f$
    !! @return 乘的结果\f$c_i^\prime\f$
    function multiply_dblearray3_typeij(array, this) result(dotmult)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        real(R_P), dimension(3, 5), intent(in) :: array
        type(vector_cmplx_type) :: dotmult
        complex(R_P) :: a(3)
        integer :: i

        do i=1, 3
          a(i)=array(i, :)*this
        end do

        call dotmult%Set(a(1), a(2), a(3))

    end function multiply_dblearray3_typeij

    !> 向量乘(内积) \f$c^\prime=\overrightarrow a \cdot \widehat\phi^\prime\f$
    !! @param[in] scalar 向量\f$\overrightarrow a\f$
    !! @return 乘的结果\f$c^\prime\f$
    complex(R_P) function multiply_cmplxarray_typeij(array, this)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        complex(R_P), dimension(5), intent(in) :: array

        multiply_cmplxarray_typeij=array(1)*this%rho + &
                                  array(2:4)*this%Vel +&
                                  array(5)*this%T

    end function multiply_cmplxarray_typeij

    !> 同时做向量乘(内积) \f$c_i^\prime=\overrightarrow a_i \cdot \widehat\phi^\prime, i=1,2,3\f$
    !! @param[in] scalar 向量\f$\overrightarrow a\f$
    !! @return 乘的结果\f$c_i^\prime\f$
    function multiply_cmplxarray3_typeij(array, this) result(dotmult)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        complex(R_P), dimension(3, 5), intent(in) :: array
        type(vector_cmplx_type) :: dotmult
        complex(R_P) :: a(3)
        integer :: i

        do i=1, 3
          a(i)=array(i, :)*this
        end do

        call dotmult%Set(a(1), a(2), a(3))

    end function multiply_cmplxarray3_typeij

    !> 矩阵乘 \f$\widehat\varphi^\prime=A\widehat\phi^\prime\f$
    !! @param[in] matrix 矩阵\f$A\f$
    !! @return 乘的结果\f$\widehat\varphi^\prime\f$
    type(dis_flux_ij_type) function multiply_dblematrix_typeij(matrix, this)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        real(R_P), dimension(5, 5), intent(in) :: matrix

        multiply_dblematrix_typeij%rho = matrix(1  , :)*this
        multiply_dblematrix_typeij%Vel = matrix(2:4, :)*this
        multiply_dblematrix_typeij%T   = matrix(5  , :)*this
        multiply_dblematrix_typeij%i   = this%i
        multiply_dblematrix_typeij%j   = this%j

    end function multiply_dblematrix_typeij

    !> 矩阵乘 \f$\widehat\varphi^\prime=A\widehat\phi^\prime\f$
    !! @param[in] matrix 矩阵\f$A\f$
    !! @return 乘的结果\f$\widehat\varphi^\prime\f$
    type(dis_flux_ij_type) function multiply_cmplxmatrix_typeij(matrix, this)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        complex(R_P), dimension(5, 5), intent(in) :: matrix

        multiply_cmplxmatrix_typeij%rho = matrix(1  , :)*this
        multiply_cmplxmatrix_typeij%Vel = matrix(2:4, :)*this
        multiply_cmplxmatrix_typeij%T   = matrix(5  , :)*this
        multiply_cmplxmatrix_typeij%i   = this%i
        multiply_cmplxmatrix_typeij%j   = this%j

    end function multiply_cmplxmatrix_typeij

    !> 矩阵除 \f$\widehat\varphi^\prime=A^{-1}\widehat\phi^\prime\f$
    !! @param[in] matrix 矩阵\f$A\f$
    !! @return 除的结果\f$\widehat\varphi^\prime\f$
    type(dis_flux_ij_type) function div_cmplxmatrix_typeij(matrix, this)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        complex(R_P), dimension(5, 5), intent(in) :: matrix
        complex(R_P), dimension(5) :: Flux(5)
        complex(R_P), dimension(5,5) :: mf
        integer :: ipiv(5), info
        complex(R_P), dimension(25) :: ma
        integer :: ia(6), ja(25)
        integer :: ii, ij
        integer(8) :: pt(64)
        integer :: iparm(64)
        complex(R_P), dimension(5) :: rhs, x
        integer :: mtype, idum(1), msglvl, ierr

        call this%get(flux)

        mf=matrix

        ! do ii=1, 5
        !   do ij=1, 5
        !     ma(ij+(ii-1)*5)=matrix(ii, ij)
        !     ja(ij+(ii-1)*5)=ij
        !   enddo
        !   ia(ii)=1+(ii-1)*5
        ! enddo
        ! ia(6)=26
        !
        ! mtype=13
        ! msglvl=0
        ! call pardisoinit(pt, mtype, iparm)
        ! rhs=flux
        ! call pardiso(pt, 1, 1, mtype, 13, 5, ma, ia, ja, idum, 1, iparm, msglvl, rhs, x, ierr)
        ! flux=x


        ! call zgetrf(5, 5, mf, 5, ipiv, info)
        ! call zgetrs('N', 5, 1, mf, 5, ipiv, Flux, 5, info)

        call solve_matrix(5, matrix, flux, x)
        flux=x

        div_cmplxmatrix_typeij%rho = Flux(1)
        call div_cmplxmatrix_typeij%Vel%set(Flux(2), Flux(3), Flux(4))
        div_cmplxmatrix_typeij%T   = Flux(5)
        div_cmplxmatrix_typeij%i   = this%i
        div_cmplxmatrix_typeij%j   = this%j

    end function div_cmplxmatrix_typeij

    subroutine solve_matrix(n, a, rhs, x)

      implicit none
      integer, intent(in) :: n
      complex(R_P), intent(in) :: a(n, n), rhs(n)
      complex(R_P), intent(out) :: x(n)
      complex(R_P) :: mf(5,5)
      integer :: ipiv(5), info

      mf=a(5:1:-1, 5:1:-1)
      call zgetrf(n, n, mf, 5, ipiv, info)
      x=rhs(5:1:-1)
      call zgetrs('N', 5, 1, mf, 5, ipiv, x, 5, info)
      x=x(5:1:-1)

    end subroutine solve_matrix


    !> 通量相加
    type(dis_flux_type) function add_type(obj1, obj2)

        implicit none
        class(dis_flux_type), intent(in) :: obj1
        type(dis_flux_type), intent(in) :: obj2

        add_type%rho = obj1%rho+obj2%rho
        add_type%Vel = obj1%Vel+obj2%Vel
        add_type%T   = obj1%T+obj2%T

    end function add_type

    !> 通量赋值
    subroutine assignment_typeij(obj, this)

        implicit none
        class(dis_flux_ij_type), intent(in) :: this
        complex(R_P), intent(inout) :: obj(5)

        obj(1)=this%rho; obj(2:4)=0.0d0; obj(5)=this%T

    end subroutine assignment_typeij


    !> 通量相加
    type(dis_flux_ij_type) function add_typeij(obj1, obj2)

        implicit none
        class(dis_flux_ij_type), intent(in) :: obj1
        type(dis_flux_ij_type), intent(in) :: obj2

        add_typeij%rho = obj1%rho+obj2%rho
        add_typeij%Vel = obj1%Vel+obj2%Vel
        add_typeij%T   = obj1%T+obj2%T
        add_typeij%i   = obj1%i
        add_typeij%j   = obj1%j

    end function add_typeij

    !> 通量相加
    type(dis_flux_ij_type) function add_type_array(obj1, obj2)

        implicit none
        class(dis_flux_ij_type), intent(in) :: obj1
        complex(R_P), intent(in) :: obj2(5)

        add_type_array%rho = obj1%rho+obj2(1)
        add_type_array%Vel = obj1%Vel+obj2(2:4)
        add_type_array%T   = obj1%T+obj2(5)

    end function add_type_array

    !> 通量相加
    type(dis_flux_ij_type) function add_array_type(obj2, obj1)

        implicit none
        class(dis_flux_ij_type), intent(in) :: obj1
        complex(R_P), intent(in) :: obj2(5)

        add_array_type%rho = obj1%rho+obj2(1)
        add_array_type%Vel = obj1%Vel+obj2(2:4)
        add_array_type%T   = obj1%T+obj2(5)

    end function add_array_type


    !> 通量相减
    type(dis_flux_ij_type) function minus_typeij(obj1, obj2)

        implicit none
        class(dis_flux_ij_type), intent(in) :: obj1
        type(dis_flux_ij_type), intent(in) :: obj2

        minus_typeij%rho = obj1%rho-obj2%rho
        minus_typeij%Vel = obj1%Vel-obj2%Vel
        minus_typeij%T   = obj1%T-obj2%T
        minus_typeij%i   = obj1%i
        minus_typeij%j   = obj1%j

    end function minus_typeij

    !> 输出扰动.
        subroutine print(this, iounit)

            implicit none
            class(dis_flux_type), intent(in) :: this
            integer, intent(in) :: iounit !< 设备编号
            complex(R_P) :: dis(5)

            call this%get(dis)
            write(iounit, *)'rho', dis(1)
            write(iounit, *)'u', dis(2)
            write(iounit, *)'v', dis(3)
            write(iounit, *)'w', dis(4)
            write(iounit, *)'T', dis(5)

    end subroutine print

!    type(dis_flux_ij_type) function multiply_disop_typeij(Disop, this)
!
!        use mod_lpse_dis_OP_point
!        use mod_vector_cmplx
!        implicit none
!        class(dis_flux_ij_type), intent(in) :: this
!        type(lpse_dis_op_point_type), intent(in) :: disop
!
!        call disop%get()
!
!    end function multiply_disop_type
!
!    !> X=A^{-1}b
!    function div_typeij(A, this) result(x)
!
!        implicit none
!        complex*16, dimension(5,5), intent(in) :: A
!        class(dis_flux_ij_type), intent(in) :: this
!        type(dis_flux_ij_type) :: x
!        complex*16 :: b(5)
!
!        x%i=this%i; x%j=this%j
!        b=this
!        call x%Set(zsolve(5, A, b))
!
!    end function div_typeij
!
!
!    function zsolve(n, a, b)  result(x)
!
!        implicit none
!        integer, intent(in) :: n
!        complex(R_P), intent(in) :: a(n,n),b(n)
!        integer :: lda,ldb,ipiv(n),info
!        complex(R_P) ::ra(n,n), x(n)
!
!          ra=a
!          lda=n
!          ldb=n
!          x=b
!          call zgetrf(n, n, ra, lda, ipiv, info)
!          if(info/=0) write(0,*) 'Error occured in zgetrf!'
!          call zgetrs('N', n, 1, ra, lda, ipiv, x, ldb, info)
!          if(info/=0) write(0,*) 'Error occured in zgetri!'
!
!    endfunction zsolve

    !> FFT变换 物理空间到谱空间
    subroutine FFT_zz(ZI, ZO)

        use mod_fft, only: fft_type
        implicit none
        class(dis_flux_type), intent(in) :: ZI(:)
        class(dis_flux_type), intent(out) :: ZO(:)
        complex(R_P) :: zvalue(5, size(ZI))

        integer :: n
        integer :: i, l
        type(fft_type), save :: handle

        n=size(ZI)
        if(.not. handle%HasInitial()) call handle%Initial(n)

        do i=1, n
            call ZI(i)%get(zvalue(:, i))
        enddo

        do l=1, 5
            call handle%fft(zvalue(l, :))
        end do

        do i=1, n
            call ZO(i)%set(zvalue(:, i))
        end do



    end subroutine FFT_zz

    !> FFT变换 物理空间到谱空间
    subroutine iFFT_zz(ZI, ZO)

        use mod_fft, only: fft_type
        implicit none
        class(dis_flux_type), intent(in) :: ZI(:)
        class(dis_flux_type), intent(out) :: ZO(:)
        complex(R_P) :: zvalue(5, size(ZI))

        integer :: n
        integer :: i, l
        type(fft_type), save :: handle

        n=size(ZI)
        if(.not. handle%HasInitial()) call handle%Initial(n)

        do i=1, n
            call ZI(i)%get(zvalue(:, i))
        enddo

        do l=1, 5
            call handle%ifft(zvalue(l, :))
        end do

        do i=1, n
            call ZO(i)%set(zvalue(:, i))
        end do

    end subroutine iFFT_zz

    !> FFT变换 物理空间到谱空间
    subroutine FFT_rz(ZI, ZO)

        use mod_fft, only: fft_type
        use mod_baseflow_org, only :bf_flux_org_type
        implicit none
        class(bf_flux_org_type), intent(in) :: ZI(:)
        class(dis_flux_type), intent(out) :: ZO(:)
        real(R_P) :: rvalue(5, size(ZI))
        complex(R_P) :: zvalue(5, size(ZI))

        integer :: n
        integer :: i, l
        type(fft_type), save :: handle

        n=size(ZI)
        if(.not. handle%HasInitial()) call handle%Initial(n)

        do i=1, n
            call ZI(i)%Get(rvalue(1, i), rvalue(2, i), rvalue(3, i), rvalue(4, i), rvalue(5, i))
        enddo

        do l=1, 5
            call handle%fft(rvalue(l, :), zvalue(l, :))
        end do

        do i=1, n
            call ZO(i)%set(zvalue(:, i))
        end do

    end subroutine FFT_rz

    !> FFT变换 谱空间到物理空间
    subroutine iFFT_zr(ZI, ZO)

        use mod_fft, only: fft_type
        use mod_baseflow_org, only : bf_flux_org_type
        use mod_vector
        implicit none
        class(dis_flux_type), intent(in) :: ZI(:)
        class(bf_flux_org_type), intent(out) :: ZO(:)
        real(R_P) :: rvalue(5, size(ZI))
        complex(R_P) :: zvalue(5, size(ZI))

        integer :: n
        integer :: i, l
        type(fft_type), save :: handle

        n=size(ZI)
        if(.not. handle%HasInitial()) call handle%Initial(n)

        do i=1, n
            call ZI(i)%get(zvalue(:, i))
        enddo

        do l=1, 5
            call handle%ifft(zvalue(l, :), rvalue(l, :))
        end do

        do i=1, n
            call ZO(i)%set(rvalue(1, i), Vector(rvalue(2, i), rvalue(3, i), rvalue(4, i)), rvalue(5, i))
        end do

    end subroutine iFFT_zr

end module mod_dis_flux
