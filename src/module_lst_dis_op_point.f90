!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_lst_dis_op_point.f90
!> @file
!> @breif LST方程法向某点扰动算子系数文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: lst_dis_op_point
!> @breif 线性稳定性理论(LST)法向某点处的扰动系数算子模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-05 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-05
!------------------------------------------------------------------------------
module mod_lst_dis_OP_point

    use mod_op_mat_cmplx
    use mod_local_coordinate
    use mod_lns_OP_point
    use mod_dis_wavenum
    use mod_op_mat
    use penf, only: R_P

    implicit none

    private

    complex(R_P), parameter :: zi=(0.0d0, 1.0d0)

    !> 性稳定性理论(LST)法向某点处的扰动算子类
    !!
    !! @todo 需要完善方程部分说明文档
    type, public :: lst_dis_op_point_type

        private

        type(op_mat_cmplx_type), private :: BL !<左侧扰动法向一阶导数算子系数\f$B_L\f$
        type(op_mat_cmplx_type), private :: BR !<右侧扰动法向一阶导数算子系数\f$B_R\f$
        type(op_mat_cmplx_type), private :: DL !<左侧零阶扰动算子系数\f$D_L\f$
        type(op_mat_cmplx_type), private :: DR1 !<右侧零阶扰动算子系数\f$D_{R1}\f$
        type(op_mat_cmplx_type), private :: DR2 !<右侧零阶扰动算子系数\f$D_{R2}\f$
        type(op_mat_cmplx_type), private :: Vyy !<右侧零阶扰动算子系数\f$V_{yy}\f$
        type(local_coordinate_type), private :: Coord !< 当地坐标系
        type(dis_wavenum_type) :: wavenum !< 扰动色散特征
        type(lst_bf_op_point_type) :: BFOP !< LNS方程算子系数

        Contains

          Procedure :: Set !< 设置基本流、色散特征以及当地坐标系信息
          procedure :: SetCoord => set_coord !< 设置当地坐标系
          procedure :: SetBCFarField !< 设置远场边界系数算子
          procedure :: SetBCWall !<设置壁面边界系数算子
          generic :: Get => getall !< 获得PSE扰动全部系数算子

          procedure, private :: GetAll !<获得PSE扰动全部系数算子
          procedure, private :: compute !< 计算PSE扰动系数算子

    end type lst_dis_op_point_type

    contains

    !> 设置基本流、色散特征以及当地坐标系信息,计算系数算子
    subroutine Set(this, BFOP, wavenum, Coord) !

        implicit none
        class(lst_dis_op_point_type), intent(inout) :: this
        type(dis_wavenum_type), intent(in) :: wavenum !< 扰动色散信息
        type(lst_bf_op_point_type), intent(in) :: BFOP !< LNS方程扰动系数算子
        type(local_coordinate_type), intent(in) :: Coord !< 当地坐标系

        this%BFOP=BFOP
        this%wavenum=wavenum
        this%Coord = Coord
        call this%Compute

    end subroutine Set

    !> 设置壁面边界系数算子
    subroutine SetBCWall(this, BFOP, wavenum, bctype)

        implicit none
        class(lst_dis_op_point_type), intent(inout) :: this
        type(lst_bf_op_point_type), intent(in) :: BFOP !< LNS方程扰动系数算子
        type(dis_wavenum_type), intent(in) :: wavenum
        integer :: bctype(5)
        real(R_P), dimension(5, 5) :: G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
        complex(R_P), dimension(5, 5) :: matB, matDL, matDR1
        complex(R_P) :: alpha, beta, omega
        complex(R_P), parameter :: zero(5, 5)=0.0d0
        complex(R_P), parameter :: one(5, 5)=reshape([1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                                    &  0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, &
                                    &  0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, &
                                    &  0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, &
                                    &  0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0], [5,5])
        call BFOP%Get(G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz)

        call wavenum%Get(alpha, beta, omega)

        matB=zero
        matB(1, :)=B(1, :)
        matDL=one
        matDL(1, :)=D(1, :)-zi*omega*G(1, :)+zi*C(1, :)*beta
        matDR1=zero
        matDR1(1, :)=-zi*A(1, :)

        call This%BL%Set(matB)
        call this%BR%Set(zero)
        call This%DL%set(matDL)
        call this%DR1%Set(matDR1)
        call this%DR2%set(zero)
        call this%Vyy%set(zero)

    end subroutine SetBCWall

    !> 设置原场边界系数算子
    subroutine SetBCFarField(this, BFOP, wavenum, bctype)

        implicit none
        class(lst_dis_op_point_type), intent(inout) :: this
        type(lst_bf_op_point_type), intent(in) :: BFOP !< LNS方程扰动系数算子
        type(dis_wavenum_type), intent(in) :: wavenum
        integer, intent(in) :: bctype(:)
        complex(R_P), parameter :: zero(5, 5)=0.0d0
        complex(R_P), parameter :: one(5, 5)=reshape([1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                                        &  0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, &
                                        &  0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, &
                                        &  0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, &
                                        &  0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0], [5,5])

        call This%BL%set(zero)
        call This%BR%set(zero)
        call This%DL%Set(one)
        call this%DR1%set(zero)
        call this%DR2%set(zero)
        call This%Vyy%set(zero)

    end subroutine SetBCFarField

    !> 计算PSE扰动系数算子
    subroutine compute(this)

        implicit none
        class(lst_dis_op_point_type), intent(inout) :: this
        type(op_mat_cmplx_type) :: G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz

        complex(R_P) :: alpha, beta, omega
        type(basis_type) :: basis
        type(lame_type) :: lame

        call this%BFOP%get(G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz)
        call this%wavenum%Get(alpha, beta, omega)

        This%BL=B-Vyz*zi*beta
        this%BR=zi*Vxy
        This%DL=D-zi*omega*G+zi*C*beta+Vzz*beta**2
        this%DR1=-zi*A-Vxz*beta
        this%DR2=-Vxx
        This%Vyy=Vyy

    end subroutine compute


    subroutine GetAll(this, BL, BR, DL, DR1, DR2, Vyy)

        implicit none
        complex(R_P), dimension(5, 5), intent(inout) :: BL, BR, DL, DR1, DR2, Vyy
        class(lst_dis_op_point_type), intent(in) :: this

        BL=this%BL
        BR=this%BR
        DL=this%DL
        DR1=this%DR1
        DR2=this%DR2
        Vyy=this%Vyy

    end subroutine GetAll

    !> 设置当地坐标系
    subroutine set_coord(this, Coord)

        implicit none
        type(local_coordinate_type), intent(in) :: Coord !< 当地坐标系
        class(lst_dis_op_point_type), intent(inout) :: this

        this%coord=Coord

    end subroutine set_coord

end module mod_lst_dis_OP_point
