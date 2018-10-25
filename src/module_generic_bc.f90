!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_generic_bc.f90
!> @file
!> @breif 通用边界条件模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  Module: generic_bc
!
!  DESCRIPTION:
!> @breif 通用边界条件模块.
!>
!!
!  REVISION HISTORY:
!  2017-08-04 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-08-04
!------------------------------------------------------------------------------
module mod_generic_bc

    use mod_lns_OP_point
    use mod_local_coordinate
    use mod_dis_wavenum
    use penf, only: R_P

    implicit none
    complex(R_P), parameter :: zi=(0.0d0, 1.0d0)

    private

    !> 通用边界条件模块
    type, public :: bctype_type

        private
        integer, private :: bctype !< 边界条件类型
        complex(R_P) :: A(5)=0.0d0 !< 流向导数系数\f$A_i\f$
        complex(R_P) :: B(5)=0.0d0 !< 法向导数系数\f$B_i\f$
        complex(R_P) :: D(5)=0.0d0 !< 0阶项前系数\f$C_i\f$
        complex(R_P) :: Vyy(5)=0.0d0 !< 粘性项前系数(右端项)\f$V_yy_i\f$
        complex(R_P) :: rhs !<右端项的值
        type(local_coordinate_type), private :: Coord !< 局部坐标系
        type(dis_wavenum_lpse_type) :: wavenum !< 波数色散关系
        type(lpse_BF_OP_point_type), pointer :: BFOP=>null() !< 基本流场算子

        Contains

        generic :: Set => set_bctype_lpse !< 设置边界条件类型
        generic :: Get => get_bctype, get_bctype_rhs !< 获得边界上方程系数、右端项

        procedure, private :: set_from_bctype !< 从边界条件类型设置边界条件
        procedure, private :: set1000 !< 密度扰动\f$\rho^\prime=0\f$边界条件
        procedure, private :: set1004 !< 连续性方程边界条件
!        procedure, private :: set1094
        procedure, private :: set2000 !< 流向速度扰动无滑移边界条件
        procedure, private :: set3000 !< 法向速度扰动无滑移边界条件
        procedure, private :: set3001 !< 法向速度扰动法向梯度边界条件
        procedure, private :: set4000 !< 展向速度扰动无滑移边界条件
        procedure, private :: set5000 !< 温度扰动等温无滑移边界条件
        procedure, private :: set5001 !< 温度扰动绝热无滑移边界条件
        procedure, private :: set_bctype_lpse !< 设置边界条件类型, 给定基本信息
        procedure, private :: get_bctype !< 获得边界上方程系数矩阵
        procedure, private :: get_bctype_rhs !<获得边界上的右端项
!        procedure, private :: set_bctype_alpse

    end type bctype_type

    contains

    !> 设置边界条件类型
    !!@param[in] bctype  边界条件类型
    !!@param[in] BFOP 基本流算子
    !!@param[in] wavenum 扰动色散关系
    !!@param[in] Coord 当地坐标系
    subroutine set_bctype_lpse(this, bctype, BFOP, wavenum, Coord)

        implicit none
        integer, intent(in) :: bctype
        type(local_coordinate_type), intent(in) :: Coord
        type(lpse_BF_OP_point_type), target, intent(in) :: bfop
        type(dis_wavenum_lpse_type) :: wavenum
        class(bctype_type), intent(inout) :: this

        this%bctype=bctype
        this%Coord = Coord
        this%wavenum= wavenum
        if((associated(this%BFOP))) this%BFOP=>null()
        this%BFOP => bfop

        call this%set_from_bctype

    end subroutine set_bctype_lpse

    !> 获得边界上方程系数
    !!@param[out] A 扰动形函数流向一阶导数前系数
    !!@param[out] B 扰动形函数法向一阶导数前系数
    !!@param[out] D 扰动形函数前系数
    !!@param[out] Vyy 扰动形函数法向二阶导数前系数
    subroutine get_bctype(this, A, B, D, Vyy)

        implicit none
        class(bctype_type), intent(in) :: this
        complex(R_P), dimension(5), intent(inout) :: A, B, D, Vyy

        A=this%A; B=this%B; D=this%D; Vyy=this%Vyy

    end subroutine get_bctype

    !> 获得边界条件的右端项
    subroutine get_bctype_rhs(this, rhs)

      implicit none
      class(bctype_type), intent(in) :: this
      complex(R_P), intent(out) :: rhs !<右端项

      rhs=this%rhs
    end subroutine get_bctype_rhs

    !> 从边界条件类型设置边界条件
    subroutine set_from_bctype(this)

        implicit none
        !, intent(in) ::
        class(bctype_type), intent(inout) :: this

        select case ( this%bctype)
            case (1000)  !! 密度边界条件为0
              call this%set1000
            case (1004)  !! 密度边界条件 连续性方程
              call this%set1004
            !case (1094)
            !  call this%set1094
            case (2000)  !! u 边界条件 为0
              call this%set2000
            case (3000)  !! v 边界条件 为0
              call this%set3000
            case (3001)
              call this%set3001
            case (4000)  !! w 边界条件 为0
              call this%set4000
            case (5000)  !! T 边界条件 为0
              call this%set5000
            case (5001)  !! T 边界条件 导数为0
              call this%set5001
            case default
                stop "Please input a correct bctype"
        end select

    end subroutine set_from_bctype

    !> 密度扰动\f$\rho^\prime=0\f$边界条件
    subroutine set1000(this)

        implicit none
        class(bctype_type), intent(inout) :: this

        this%A=0.0d0
        this%B=0.0d0
        this%D=0.0d0
        this%D(1)=1.0d0
        this%Vyy=0.0d0
        this%rhs=0.0d0

    end subroutine set1000

    !> 连续性方程边界条件
    subroutine set1004(this)

        implicit none
        class(bctype_type), intent(inout) :: this
        real(R_P), dimension(5, 5) :: G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
        complex(R_P) :: alpha, beta, omega, DxAlpha
        type(basis_type) :: basis
        type(lame_type) :: lame
        real(R_P) :: hx, hy, hz
        real(R_P) :: hxhx, hxhz, hzhz
        real(R_P) :: InvHx, InvHz, InvHxHx, InvHxHz, InvHzHz

        call this%BFOP%Get(G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz)
        call this%wavenum%Get(alpha, beta, omega)

        This%A=A(1, 1:5)
        This%B=B(1, 1:5)
        This%D=D(1, 1:5)-zi*omega*G(1, 1:5)+zi*alpha*A(1, 1:5)+zi*C(1, 1:5)*beta
        This%Vyy=0.0d0
        this%rhs=1.0d17

    end subroutine set1004

    !subroutine set1094(this)
    !
    !    implicit none
    !    class(bctype_type), intent(inout) :: this
    !    real(R_P), dimension(5, 5) :: G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
    !    real(R_P), dimension(5, 5) :: DxA, DyB, DyyVyy, DyVyy, DyVxy, DyVyz
    !    complex(R_P) :: alpha, beta, omega, DxAlpha
    !    type(type_basis) :: basis
    !    type(type_lame) :: lame
    !    real(R_P) :: hx, hy, hz
    !    real(R_P) :: hxhx, hxhz, hzhz
    !    real(R_P) :: InvHx, InvHz, InvHxHx, InvHxHz, InvHzHz
    !
    !    call this%BFOP%Get(G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz)
    !    select type(BFOP => this%BFOP)
    !    type is (type_alpse_BF_OP_point)
    !        call BFOP%GetGrad(DxA, DyB, DyyVyy, DyVyy, DyVxy, DyVyz)
    !    end select
    !    call this%wavenum%Get(alpha, beta, omega)
    !
    !    This%A=A(1, 1:5)
    !    This%B=-B(1, 1:5)
    !    This%D=D(1, 1:5)-zi*omega*G(1, 1:5)+zi*alpha*A(1, 1:5)+zi*C(1, 1:5)*beta-DyB(1, 1:5)
    !    This%Vyy=0.0d0
    !
    !end subroutine set1094

    !> 流向速度扰动无滑移边界条件
    subroutine set2000(this)

        implicit none
        class(bctype_type), intent(inout) :: this

        this%A=0.0d0
        this%B=0.0d0
        this%D=0.0d0
        this%D(2)=1.0d0
        this%Vyy=0.0d0
        this%rhs=0.0d0

    end subroutine set2000

    !> 法向速度扰动无滑移边界条件
    subroutine set3000(this)

        implicit none
        class(bctype_type), intent(inout) :: this

        this%A=0.0d0
        this%B=0.0d0
        this%D=0.0d0
        this%D(3)=1.0d0
        this%Vyy=0.0d0
        this%rhs=0.0d0

    end subroutine set3000

    !> 法向速度扰动梯度为零边界条件
    subroutine set3001(this)

        implicit none
        class(bctype_type), intent(inout) :: this

        this%A=0.0d0
        this%B=0.0d0
        this%B(3)=1.0d0
        this%D=0.0d0
        this%Vyy=0.0d0
        this%rhs=0.0d0

    end subroutine set3001

    !> 展向速度扰动无滑移边界条件
    subroutine set4000(this)

        implicit none
        class(bctype_type), intent(inout) :: this

        this%A=0.0d0
        this%B=0.0d0
        this%D=0.0d0
        this%D(4)=1.0d0
        this%Vyy=0.0d0
        this%rhs=0.0d0

    end subroutine set4000

    !> 温度扰动等温无滑移边界条件
    subroutine set5000(this)

        implicit none
        class(bctype_type), intent(inout) :: this

        this%A=0.0d0
        this%B=0.0d0
        this%D=0.0d0
        this%D(5)=1.0d0
        this%Vyy=0.0d0
        this%rhs=0.0d0

    end subroutine set5000

    !> 温度扰动绝热无滑移边界条件
    subroutine set5001(this)

        implicit none
        class(bctype_type), intent(inout) :: this

        this%A=0.0d0
        this%B=0.0d0
        this%B(5)=1.0d0
        this%D=0.0d0
        this%Vyy=0.0d0
        this%rhs=0.0d0

    end subroutine set5001

end module mod_generic_bc
