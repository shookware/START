!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_baseflow_normal.f90
!> @file
!> @breif 基本流法向分布模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: baseflow_normal
!> @breif 法向基本流类型模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-03 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-03
!------------------------------------------------------------------------------
module mod_BF_normal

    use mod_baseflow_org
    use mod_grid
    use mod_baseflow
    use mod_difference
    use penf, only: R_P

    implicit none

    private

    !> 基本流通量点(含导数信息)类型
    type, public :: bf_point_type

        private
        type(bf_flux_org_ij_type) :: Flux                                       !< 基本流通量\private
        type(bf_flux_org_ij_type) :: Dx_Flux                                    !< 基本流通量流向一阶导数\private
        type(bf_flux_org_ij_type) :: Dy_Flux                                    !< 基本流通量法向一阶导数\private
        type(bf_flux_org_ij_type) :: DDy_Flux                                   !< 基本流通量法向二阶导数\private
        type(bf_flux_org_ij_type) :: Dz_Flux                                    !< 基本流通量展向一阶导数\private
        type(bf_flux_org_ij_type) :: Dzz_Flux                                   !< 基本流通量展向二阶导数\private

        Contains

        procedure :: Get => get_Flux

    end type bf_point_type

    !> 法向基本流通量类型.
    !!
    !! 表示了在某流向展向位置沿法向一列基本流的信息
    type, public :: BF_normal_type

        private
        integer, private :: iloc    !< 流向站位
        integer, private :: jn      !< 法向点数
        type(bf_point_type), allocatable, private :: BF_normal(:) !< 基本流通量法向分布

        Contains

          generic :: Set=> Set_diff    !< 设置法向基本流相关信息
          Procedure :: Create             !< 创建一个法向基本流分布类,分配内存
          procedure :: GetPoint !< 获得法向基本流点通量的基本流通量以及各导数通量
          procedure :: GetJn   !< 获得法向点数
          procedure :: HasCreated !< 是否分配空间
          procedure :: finalize => finalize_BF_nromal
          procedure :: GetBFPoint !< 获得法向点基本流点通量

!          procedure :: Set_dbf
          procedure :: Set_diff

    end type BF_normal_type

    contains

    !> 判断是否分配内存
    function HasCreated(this) result(Flag)
        class(BF_normal_type), intent(inout) :: this
        logical :: Flag

        Flag=ALLOCATED(this%BF_normal)

    end function HasCreated

    !> 获得法向点数
    !! @return 法向点数
    integer function GetJn(this) result(jn)

        implicit none
        class(BF_normal_type), intent(in) :: this

        jn=this%jn

    end function GetJn

    !> 创建一个法向基本流分布类,分配内存
    !! @param[in] jn 法向点数
    subroutine create(this, jn)

        implicit none
        integer, intent(in) :: jn
        class(BF_normal_type), intent(inout) :: this

        this%jn=jn
        if(.not. allocated(this%BF_normal)) &
        & allocate(this%BF_normal(jn))

    end subroutine create

    !!> 设置法向基本流相关信息
    !!! @param[in] iloc 流向站位
    !!! @param[in] Flow 二维基本流场信息
    !!! @param[in] DiffFlow 二维基本流场导数信息
    !subroutine set_dbf(this, iloc, Flow, DiffFlow)
    !
    !    implicit none
    !    integer, intent(in) :: iloc
    !    type(Baseflow2D_type), intent(in) :: Flow
    !    type(Baseflow2D_Diff_type), intent(in) :: DiffFlow
    !    class(BF_normal_type), intent(inout) :: this
    !    integer :: j
    !
    !    this%iloc = iloc
    !    Do j=1, this%jn
    !      this%BF_normal(j)%Flux=Flow%GetPoint(iloc, j)
    !      this%BF_normal(j)%Dx_Flux=DiffFlow%GetPartial_Xi(iloc, j)
    !      this%BF_normal(j)%Dy_Flux=DiffFlow%GetPartial_Eta(iloc, j)
    !      this%BF_normal(j)%DDy_Flux=DiffFlow%GetPartial_Eta2(iloc, j)
    !    enddo
    !
    !end subroutine set_dbf

    !> 设置法向基本流相关信息
    !! @param[in] iloc 流向站位
    !! @param[in] Flow 二维基本流场信息
    !! @param[in] Diff 二维差分算子信息
    subroutine set_diff(this, iloc, Flow, Diff)

        use mod_baseflow_org
        implicit none
        integer, intent(in) :: iloc
        type(Baseflow2D_type), target, intent(in) :: Flow
        type(difference_2D_type), intent(in) :: Diff
        class(BF_normal_type), intent(inout) :: this
        type(bf_flux_org_ij_type) :: Dy_Flux(this%jn), Dyy_Flux(this%jn)
        type(bf_flux_org_ij_type) :: Dx_Flux(this%jn)
        integer :: j

        this%iloc = iloc
        call Diff%Partial_Eta_BFFluxIJ_iloc(iloc, Flow%Flux_org(:,:), Dy_Flux(:))
        call Diff%Partial_Eta2_BFFluxIJ_iloc(iloc, Flow%Flux_org(:,:), Dyy_Flux(:))
        call Diff%Partial_XI_BFFluxIJ_iloc(1, this%jn, iloc, Flow%Flux_org(:,:), Dx_Flux( :))

        Do j=1, this%jn
          this%BF_normal(j)%Flux=Flow%GetPoint(iloc, j)
          this%BF_normal(j)%Dx_Flux=Dx_flux(j)
          this%BF_normal(j)%Dy_Flux=Dy_flux(j)
          this%BF_normal(j)%DDy_Flux=Dyy_flux(j)
        enddo


    end subroutine set_diff

    !> 获得法向基本流中一个法向位置点处的信息
    !! @param[in] j 法向位置
    !! @return j点的基本流信息以及其导数信息
    function getpoint(this, j) result(bfpoint)

        implicit none
        class(BF_normal_type), intent(in) :: this
        integer, intent(in) :: j
        type(bf_point_type) :: bfpoint

        bfpoint=this%BF_normal(j)

    end function getpoint

    !> 获得法向基本流中一个法向位置点处的信息
    !! @param[in] j 法向位置
    !! @return j点的基本流信息以及其导数信息
    function getBFpoint(this, j) result(BFFLUX)

        implicit none
        class(BF_normal_type), intent(in) :: this
        integer, intent(in) :: j
        type(BF_flux_org_ij_type) :: BFFlux
        type(bf_point_type) :: bfpoint

        bfpoint=this%BF_normal(j)
        BFFlux=bfpoint%Flux

    end function getBFpoint

    !> 获得法向基本流点通量的基本流通量以及各导数通量
    !! @param[out] flux 基本流通量
    !! @param[out] dx_flux 基本流通量沿流向的一阶导数
    !! @param[out] dy_flux 基本流通量沿法向的一阶导数
    !! @param[out] ddy_flux 基本流通量沿法向的二阶导数
    subroutine get_flux(this, flux, dx_flux, dy_flux, ddy_flux, dz_flux, dzz_flux)

        implicit none
        class(bf_point_type), intent(in) :: this
        type(bf_flux_org_ij_type), intent(inout) :: flux, dx_flux, dy_flux, ddy_flux
        type(bf_flux_org_ij_type), intent(inout), optional :: dz_flux, dzz_flux

        flux    =this%flux
        dx_flux =this%dx_flux
        dy_flux =this%dy_flux
        ddy_flux=this%ddy_flux

        if(present(dz_flux))then
            dz_flux=this%Dz_Flux
            dzz_flux=this%Dzz_Flux
        end if

    end subroutine get_flux

    !> 析构函数
    Subroutine finalize_BF_nromal(this)

        Implicit None
        class(BF_normal_type), intent(Inout) :: this

        this%iloc=0
        this%jn=0
        deallocate(this%BF_normal)

    End Subroutine finalize_BF_nromal

end module mod_BF_normal
