!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_lpse_op_point.f90
!> @file
!> @breif PSE在法向某点的扰动线性算子.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: lpse_op_point
!> @breif PSE在法向某点处的扰动线性系数算子模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-04 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-04
!------------------------------------------------------------------------------
module mod_lpse_dis_OP_point

    use mod_op_mat_cmplx
    use mod_local_coordinate
!    use mod_lpse_BF_OP_point
    use mod_lns_OP_point
    use mod_dis_wavenum
    use mod_op_mat
    use penf, only: R_P
    implicit none

    private

    complex(R_P), parameter :: zi=(0.0d0, 1.0d0) !< 复数i

    !> PSE在法向某点处的扰动线性系数算子类
    type, public :: lpse_dis_op_point_type

        private

        integer, private :: jloc !< 法向位置
        type(op_mat_cmplx_type), private :: A !< 流向一阶导数系数算子
        type(op_mat_cmplx_type), private :: B !< 法向一阶导数系数算子
        type(op_mat_cmplx_type), private :: D !< 零阶项系数算子
        type(op_mat_cmplx_type), private :: Vyy !< 法向二阶导数系数算子
        type(op_mat_cmplx_type), private :: DxA !<  流向一阶导数系数算子在流向的一阶导数
        type(local_coordinate_type), private :: Coord !<当地坐标系
        type(dis_wavenum_lpse_type) :: wavenum !< 扰动色散特征
        type(lpse_BF_OP_point_type) :: BFOP !< LNS方程算子(基本流部分算子)

        Contains

          Procedure :: Set !< 设置基本流、色散特征以及当地坐标系信息
          procedure :: SetBC => set_bc !< 设置边界上的算子
          procedure :: SetCoord => set_coord !< 设置坐标系
          procedure :: GetSpt !<获得谱半径
          generic :: get => getall, getOP!<获得PSE扰动系数算子

          procedure, private :: GetAll !< 获得PSE扰动全部系数算子
          procedure, private :: GetOP !< 获得指定的PSE扰动系数算子
          procedure, private :: compute !< 计算PSE扰动系数算子

    end type lpse_dis_op_point_type

    contains

    !> 设置基本流、色散特征以及当地坐标系信息,计算系数算子
    subroutine Set(this, BFOP, wavenum, Coord)

        implicit none
        type(dis_wavenum_lpse_type), intent(in) :: wavenum !< 扰动色散特征
        type(lpse_BF_OP_point_type), intent(in) :: BFOP !< LNS方程算子(基本流部分算子)
        type(local_coordinate_type), intent(in) :: Coord !< 当地坐标系特征
        class(lpse_dis_op_point_type), intent(inout) :: this

        this%BFOP=BFOP
        this%wavenum=wavenum
        this%Coord = Coord
        call this%Compute

    end subroutine Set

    !> 设定系数算子在边界上的值
    subroutine set_bc(this, LpseBC)

        use mod_lpse_bc
        implicit none
        class(lpse_dis_op_point_type), intent(inout) :: this
        type(lpse_bc_type), intent(in) :: LpseBC !< 边界条件处的算子系数特征
        complex(R_P), dimension(5, 5) :: A, B, D, Vyy

        call LpseBC%get(A, B, D, Vyy)
        this%A=A
        this%B=B
        this%D=D
        this%Vyy=Vyy

    end subroutine set_bc

    !> 计算PSE扰动系数算子
    subroutine compute(this)

        implicit none
        class(lpse_dis_op_point_type), intent(inout) :: this
        type(op_mat_cmplx_type) :: G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz, DxA

        complex(R_P) :: alpha, beta, omega, DxAlpha
        type(basis_type) :: basis
        type(lame_type) :: lame
        real(R_P) :: hx, hy, hz
        real(R_P) :: hxhx, hxhz, hzhz
        real(R_P) :: InvHx, InvHz, InvHxHx, InvHxHz, InvHzHz
        real(R_P) :: VigneronA(2)
        complex(R_P) :: matA(5, 5)

        call this%BFOP%get(G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz)
        call this%BFOP%getDxA(DxA)
        call this%wavenum%Get(alpha, beta, omega)
        DxAlpha=this%wavenum%GetDxAlpha()
        call this%Coord%Get(basis, lame)
        call lame%Get(hx, hy, hz)
        hxhx=hx**2; hxhz=hx*hz; hzhz=hz**2
        InvHx=1.0d0/hx; InvHz=1.0d0/hz
        InvHxHx=1.0d0/Hxhx; InvHxHz=1.0d0/HxHz; InvHzHz=1.0d0/HzHz

        This%B=B-zi*alpha*Vxy-Vyz*zi*beta
        This%D=D-zi*omega*G+zi*alpha*A+zi*C*beta+ &
               (alpha**2)*Vxx+Vzz*beta**2+ & !-zi*Dxalpha
               Vxz*alpha*beta
        This%Vyy=Vyy
        call this%BFOP%GetVigneronA(VigneronA)
        call A%Get(matA)
        matA(2, 1)=matA(2, 1)-VigneronA(1)
        matA(2, 5)=matA(2, 5)-VigneronA(2)
        call A%Set(matA)
        This%A=A-zi*2.0d0*alpha*Vxx-Vxz*zi*beta
        this%DxA=DxA


    end subroutine compute

    !> 获得PSE扰动全部系数算子
    subroutine GetAll(this, A, B, D, Vyy)

        implicit none
        complex(R_P), dimension(5, 5), intent(inout) :: A, B, D, Vyy
        class(lpse_dis_op_point_type), intent(in) :: this

        A=this%A
        B=this%B
        D=this%D
        Vyy=this%Vyy

    end subroutine GetAll

    !>获得指定的PSE扰动系数算子
    subroutine GetOP(this, A, charsign)

        implicit none
        complex(R_P), dimension(5, 5), intent(inout) :: A
        character :: charsign
        class(lpse_dis_op_point_type), intent(in) :: this

        select case (trim(charsign))
            case ('A')
                A=this%A
            case ('B')
                A=this%B
            case ('D')
                A=this%D
            case ('V')
                A=this%Vyy
            case ('X')
                A=this%DxA
            case default
                write(*, *)'the charsign must be the one of A, B, D, V!'
        end select

   end subroutine GetOP

   subroutine GetSpt(this, rhox, rhoy, rhoz)

       implicit none
       class(lpse_dis_op_point_type),intent(inout) :: this
       real(R_P), intent(out) :: rhox, rhoy, rhoz

       call this%BFOP%GetSpt(rhox, rhoy, rhoz)

   end subroutine GetSpt


    !> 设置当地坐标系
    subroutine set_coord(this, Coord)

        implicit none
        type(local_coordinate_type), intent(in) :: Coord
        class(lpse_dis_op_point_type), intent(inout) :: this

        this%coord=Coord

    end subroutine set_coord

end module mod_lpse_dis_OP_point
