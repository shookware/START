!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: /media/psf/Home/Fortran/START/src/module_coordinate.f90
!> @file
!> @breif 当地坐标系相关文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: local_coordinate
!> @breif 当地点坐标系模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-01 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-08-01
!------------------------------------------------------------------------------
module mod_local_coordinate

    use mod_basis
    use mod_lame
    use mod_grid
    implicit none

      !> 当地点坐标系类型
      type, public :: local_coordinate_type

          type(basis_type), private :: basis  !< 坐标系基矢量
          type(lame_type), private :: lame    !< 坐标系Lame系数
          type(lame_grad_type), private :: lameGrad !< 坐标系Lame系数梯度

          Contains

            Procedure :: Set => set_local_coordinate   !< 设置当地坐标系信息
            procedure :: Get => get_local_coordinate   !< 获得当地坐标系信息
            procedure :: GetGrad => get_lame_grad      !< 设置当地坐标系Lame系数梯度
            procedure :: SetGrad => set_lame_grad      !< 获得当地坐标系Lame系数梯度

      end type local_coordinate_type

    contains

    !> 设置当地坐标系信息
    !! @param[in] basis 当地坐标系基矢量
    !! @param[in] lame 当地坐标系Lame系数
    subroutine set_local_coordinate(this, basis, lame)

        implicit none
        type(basis_type), intent(in) :: basis
        type(lame_type), intent(in) :: lame
        class(local_coordinate_type), intent(inout) :: this

        this%basis=basis
        this%lame=lame

    end subroutine set_local_coordinate

    !> 获得当地坐标系信息
    !! @param[out] basis 当地坐标系基矢量
    !! @param[out] lame 当地坐标系Lame系数
    subroutine get_local_coordinate(this, basis, lame)

        implicit none
        class(local_coordinate_type), intent(in) :: this
        type(basis_type), intent(inout) :: basis
        type(lame_type), intent(inout) :: lame

        basis=this%basis
        lame=this%lame

    end subroutine get_local_coordinate

    !> 设置当地坐标系Lame系数梯度
    !! @param[in] Lame 当地坐标系Lame系数梯度
    subroutine set_lame_grad(this, LameGrad)

        implicit none
        class(local_coordinate_type), intent(inout) :: this
        type(lame_grad_type), intent(in) :: LameGrad

        this%lameGrad=LameGrad

    end subroutine set_lame_grad

    !> 获得当地坐标系Lame系数梯度
    !! @param[out] Lame 当地坐标系Lame系数梯度
    subroutine get_lame_grad(this, LameGrad)

        implicit none
        class(local_coordinate_type), intent(in) :: this
        type(lame_grad_type), intent(inout) :: LameGrad

        LameGrad=this%lameGrad

    end subroutine get_lame_grad

end module mod_local_coordinate

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: Curvature
!> @breif 点曲率相关模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-01 - Initial Version
!> @author
!> Liu Jianxin
!> @date 2017-08-01
!------------------------------------------------------------------------------
module mod_curvature

    use penf, only: R_P
    implicit none

    private

    !> 点曲率类.
    type, public :: curvature_type

          private

          real(R_P), private :: curvature_i=0.0d0  !< 流向曲率
          real(R_P), private :: curvature_k=0.0d0  !< 展向曲率

          contains

          procedure :: Set => set_K   !< 设置曲率
          procedure :: GetKi => get_K_i   !< 获得流向曲率
          procedure :: GetKk => get_K_k   !< 获得展向曲率

    end type curvature_type

    contains

    !> 获得流向曲率
    !! @return 流向曲率
    real(R_P) function get_K_i(this)

      implicit  none
      class(curvature_type), intent(in) :: this

      get_k_i=this%curvature_i

    end function get_K_i

    !> 获得展向曲率
    !! @return 展向曲率
    real(R_P) function get_K_k(this)

      implicit  none
      class(curvature_type), intent(in) :: this

      get_k_k=this%curvature_k

    end function get_K_k

    !> 设置流向和展向曲率
    !! @param[in] K_i 流向曲率
    !! @param[in] K_k 展向曲率
    elemental subroutine set_K(this, K_i, K_k)

        implicit none
        real(R_P), intent(in) :: K_i, K_k
        class(curvature_type), intent(inout) :: this

        this%curvature_i=K_i; this%curvature_k=K_k

    end subroutine set_K

end module mod_curvature

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: Curvature2D
!> @breif 二维网格的壁面曲率模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-01 - Initial Version
!  TODO_ 2017-08-01 - TODO_由网格生成当地曲率算法
!> @author
!> Liu Jianxin
!> @date 2017-08-01
!------------------------------------------------------------------------------
module mod_curvature_2d

    use mod_curvature
    use penf, only: R_P

    implicit none

    !> 二维壁面曲率类.
    type, public :: curvature_2d_type

        integer :: in  !< 二维流场流向点数
        type(curvature_type), allocatable :: curvature_2d(:) !< 二维场壁面曲率

        contains

        procedure :: CreateFromPlot3D    !> 从Plot3D文件直接创建二维曲率
        procedure :: CreateFromASCII     !> 从文本文件直接创建二维曲率
        procedure :: CreateWithNoCurvature  !> 创建曲率为0的情况
        procedure :: GetPoint               !> 获得其中i点位置的曲率

        procedure :: finalize => finalize_curvature_2d      !> 析构函数

    end type

    contains

    !> 从Plot3D文件创建二维曲率
    !! @param[in] filename Plot3d文件名
    subroutine CreateFromPlot3D(this, filename)

        use stringifor
        implicit none

        type(string) :: filename
        class(curvature_2d_type), intent(inout) :: this
        integer :: in
        real(R_P), allocatable :: drxz(:, :)
        integer :: i

        open(11, file=filename//'', form='unformatted')
        read(11)
        read(11)in
        if(.not. allocated(drxz)) allocate(drxz(in, 2))
        read(11)drxz
        close(11)

        this%in=in

        if(.not. allocated(this%curvature_2d)) &
        &   allocate(this%curvature_2d(in))
        do i=1, in
            call this%curvature_2d%Set(drxz(i, 1), drxz(i, 2))
        enddo
        deallocate(drxz)

    end subroutine CreateFromPlot3d

    !> 创建壁面无曲率的情况(平板)
    !! @param[in] in 流向网格点数
    subroutine CreateWithNoCurvature(this, in)

        implicit none

        class(curvature_2d_type), intent(inout) :: this
        integer :: in
        real(R_P), allocatable :: drxz(:, :)
        integer :: i

        if(.not. allocated(drxz)) allocate(drxz(in, 2))
        drxz=0.0d0

        this%in=in

        if(.not. allocated(this%curvature_2d)) &
          &   allocate(this%curvature_2d(in))
        do i=1, in
            call this%curvature_2d%Set(drxz(i, 1), drxz(i, 2))
        enddo
        deallocate(drxz)

    end subroutine CreateWithNoCurvature

    !> 从文本文件创建二维曲率
    !! @param[in] filename 文本文件文件名
    subroutine CreateFromASCII(this, filename)

        use stringifor
        implicit none

        type(string) :: filename
        class(curvature_2d_type), intent(inout) :: this
        integer :: in
        real(R_P), allocatable :: drxz(:, :)
        real(R_P) :: tmpx
        integer :: i

        open(11, file=filename//'', form='formatted')
        read(11, *)in
        if(.not. allocated(drxz)) allocate(drxz(in, 2))
        do i=1, in
            read(11, *)tmpx, drxz(i, 1), drxz(i, 2)
        enddo
        close(11)

        this%in=in

        if(.not. allocated(this%curvature_2d)) allocate(this%curvature_2d(in))
        do i=1, in
            call this%curvature_2d(i)%Set(drxz(i, 1), drxz(i, 2))
        enddo
        deallocate(drxz)

    end subroutine CreateFromASCII

    !> 获得流向某点处的壁面曲率.
    !! @param[in] iloc 流向位置序号
    !! @return 该位置曲率
    function GetPoint(this, iloc) result(point)

        implicit none
        class(curvature_2d_type), intent(in) :: this
        integer, intent(in) :: iloc
        type(curvature_type) :: point

        point=this%curvature_2d(iloc)

    end function GetPoint

    !> 析构函数
    subroutine finalize_curvature_2d(this)

        implicit none
        class(curvature_2d_type), intent(inout) :: this

        deallocate(this%curvature_2d)
        this%in=0

    end subroutine finalize_curvature_2d


end module

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: local_normal_coordinate.
!> @breif 壁面某点法向坐标系模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-01 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-08-01
!------------------------------------------------------------------------------
module mod_local_normal_coordinate

    use mod_grid
    use mod_curvature
    use mod_local_coordinate
    use penf, only: R_P
    implicit none

    private

      !> 壁面某点法向坐标系类.
      type, public :: local_normal_coordinate_type

          private
          integer, private :: jn !< 法向点数
          type(local_coordinate_type), allocatable, private :: localCoord(:) !< 法向上每点的坐标系
          type(curvature_type), private :: curvature    !< 当地壁面曲率

          Contains

            Procedure :: Create => Create_local_normal_coordinate_from_grid  !< 创建类型,分配内存
            procedure :: SetLameFromGrid => &                                !< 从网格设置Lame系数的值
                                &   Set_local_normal_coordinate_lame_from_grid
            procedure :: SetCurvature => set_curvature         !< 设置壁面曲率
            procedure :: GetCurvature => get_curvature !<获得壁面曲率
            procedure :: GetPoint => get_local_coord            !< 获得当地法向某位置上的当地点坐标系
            procedure :: finalize => finalize_local_normal_coordinate   !< 析构函数
            procedure :: HasCreated !< 判断类型是否分配

      end type local_normal_coordinate_type

    contains


    !>判断类型是否分配
    !! @return 类型是否分配
    function HasCreated(this) result(Flag)
        class(local_normal_coordinate_type), intent(inout) :: this
        logical :: Flag

        Flag=allocated(this%localCoord)

    end function HasCreated


    !> 获得当地某法向位置上的当地点坐标系
    !! @param[in] Jindex 法向点位置序号
    !! @return 对应Jindex点的当地点坐标系
    function get_local_coord(this, Jindex) result(Coord)

        implicit none
        class(local_normal_coordinate_type), intent(in) :: this
        integer, intent(in) :: Jindex
        type(local_coordinate_type) :: Coord

        Coord=this%localCoord(Jindex)

    end function get_local_coord

    !> 从网格创建当地法向坐标系类型,分配内存
    !! @param[in] jn 法向点数
    subroutine Create_local_normal_coordinate_from_grid(this, jn)

        implicit none
        integer, intent(in) :: jn
        class(local_normal_coordinate_type), intent(inout) :: this
        integer :: in

        this%jn=jn
        if(.not. allocated(this%localCoord)) allocate(this%localCoord(jn))

    end subroutine Create_local_normal_coordinate_from_grid

    !> 设置壁面曲率
    !! @param[in] curvature 壁面曲率
    subroutine set_curvature(this, curvature)

        implicit none
        type(curvature_type), intent(in) :: curvature
        class(local_normal_coordinate_type), intent(inout) :: this

        this%curvature=curvature

    end subroutine set_curvature

    !> 获得壁面曲率
    !! @param[in] curvature 壁面曲率
    subroutine get_curvature(this, curvature)

        implicit none
        type(curvature_type), intent(out) :: curvature
        class(local_normal_coordinate_type), intent(in) :: this

        curvature=this%curvature

    end subroutine get_curvature

    !> 从网格设置Lame系数的值.
    !! @param[in] grid 二维网格
    !! @param[in] i_loc 当前位置的流向站位序号
    subroutine Set_local_normal_coordinate_lame_from_grid(this, grid, i_loc)

        implicit none
        integer, intent(in) :: i_loc
        type(grid2d_type), intent(in) :: grid !! must be a grid bodyfixxed
        class(local_normal_coordinate_type), intent(inout) :: this
        type(lame_type) :: lame
        type(lame_grad_type) :: LameGrad
        real(R_P) :: hx, hy, hz
        real(R_P) :: yy(this%jn)
        integer :: j

        call grid%Get_iy(i_loc, yy)

        associate(curvature_i => this%curvature%GetKi(), curvature_k => this%curvature%GetKk())
        do j=1, this%jn
          hx=1.0d0+curvature_i*yy(j); hy=1.0d0; hz=1.0d0+curvature_k*yy(j)
          call lame%Set(hx, hy, hz)
          call LameGrad%Set(curvature_i/hx, curvature_k/hz, 0.0d0) ! d31 value must be calculated by the difference
          call this%localcoord(j)%Set(BASIS_NULL, lame)
          call this%localcoord(j)%SetGrad(LameGrad)
        enddo
        end associate

    end subroutine Set_local_normal_coordinate_lame_from_grid

    !> 析构函数
    Subroutine finalize_local_normal_coordinate(this)

        Implicit None
        class(local_normal_coordinate_type), intent(Inout) :: this

        this%jn=0
        deallocate(this%localCoord)

    End Subroutine finalize_local_normal_coordinate

end module mod_local_normal_coordinate
