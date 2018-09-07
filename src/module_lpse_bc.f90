!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_lpse_bc.f90
!> @file
!> @breif PSE边界条件文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  Module: lpse_bc
!
!  DESCRIPTION:
!> @breif PSE边界条件类型模块.
!>
!!
!  REVISION HISTORY:
!  2017-08-04 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-08-04
!------------------------------------------------------------------------------
 module mod_lpse_bc

     use mod_generic_bc
     use mod_local_coordinate
     use mod_dis_wavenum
     use mod_lns_OP_point
     use penf, only: R_P

     implicit none

     private

     !> PSE边界条件矩阵类型.
     !!
     !! 包含PSE方程在边界条件处的方程系数矩阵信息
     type, public :: lpse_bc_type

        private
        type(bctype_type) :: BC(5) !< 法向某点处的边界条件

        Contains

          generic :: set => set_bc_lpse!, set_bc_alpse !< 设置边界条件系数矩阵
          generic :: get => get_bc, get_bc_rhs !< 获得边界条件系数矩阵和值

          procedure, private :: set_bc_lpse !< 设置边界条件系数矩阵
!          procedure, private :: set_bc_alpse
          procedure, private :: get_bc !< 获得边界条件系数矩阵
          procedure, private :: get_bc_rhs !<获得边界条件右端项的值

     end type lpse_bc_type

     contains

      !> 设置边界条件系数矩阵
      !!@param[in] bctype 边界条件类型
      !!@param[in] BFOP 基本流算子
      !!@param[in] wavenum 扰动波色散信息
      !!@param[in] Coord 当地坐标系信息
      subroutine set_bc_lpse(this, bctype, BFOP, wavenum, Coord)

          implicit none
          integer, intent(in) :: bctype(5)
          type(local_coordinate_type), intent(in) :: Coord
          type(lpse_BF_OP_point_type), intent(in) :: bfop
          type(dis_wavenum_lpse_type) :: wavenum
          class(lpse_bc_type), intent(inout) :: this
          integer :: icount

          do icount=1, 5
#IFDEF    DEBUG
!            print*, 'iCount=', icount
#ENDIF
            call this%bc(icount)%Set(bctype(icount), BFOP, wavenum, Coord)
          end do

     end subroutine set_bc_lpse

     !subroutine set_bc_alpse(this, bctype, BFOP, wavenum, Coord)
     !
     !     implicit none
     !     integer, intent(in) :: bctype(5)
     !     type(local_coordinate_type), intent(in) :: Coord
     !     type(alpse_BF_OP_point_type), intent(in) :: bfop
     !     type(dis_wavenum_lpse_type) :: wavenum
     !     class(lpse_bc_type), intent(inout) :: this
     !     integer :: icount
     !
     !     do icount=1, 5
     !       call this%bc(icount)%Set(bctype(icount), BFOP, wavenum, Coord)
     !     end do
     !
     !end subroutine set_bc_alpse

    !> 获得边界条件系数矩阵
    !!@param[out] A 扰动形函数流向一阶导数前系数
    !!@param[out] B 扰动形函数法向一阶导数前系数
    !!@param[out] D 扰动形函数前系数
    !!@param[out] Vyy 扰动形函数法向二阶导数前系数
     subroutine get_bc(this, A, B, D, Vyy)

         implicit none
         class(lpse_bc_type), intent(in) :: this
         complex(R_P), dimension(5, 5), intent(inout) :: A, B, D, Vyy
         integer :: icount

         do icount=1, 5
            call this%BC(icount)%Get(A(icount, :), B(icount, :), D(icount, :), Vyy(icount, :))
         end do

     end subroutine get_bc

     !>获得边界条件右端项的值
     subroutine get_bc_rhs(this, rhs)

       implicit none
       class(lpse_bc_type), intent(in) :: this
       complex(R_P), intent(out) :: rhs(5)
       integer :: iCount

       do iCount=1, 5
         call this%bc(iCount)%get(rhs(iCount))
       enddo

     end subroutine get_bc_rhs

 end module mod_lpse_bc
