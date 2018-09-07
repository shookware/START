!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: /media/psf/Home/Fortran/START/src/module_grid.f90
!> @file
!> @breif 网格模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: Grid
!> @breif 网格模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-01 - Initial Version
!  TODO 2017-08-01 - 3D网格类 - 编写三维网格类
!> @author
!> Liu Jianxin
!> @date 2017-08-01
!------------------------------------------------------------------------------
Module mod_grid

    Use mod_point
    use penf, only: R_P

    Implicit None

    private

    type, public :: grid_object_type
        integer, private :: In=0 !< 流向网格点数
        integer, private :: jn=0 !< 法向网格点数
        logical, private :: isBodyFixed =.False.   !< 是否贴体网格
    contains
        procedure :: GetJnSize => get_jn_size               !<获得网格法向网格点数
        procedure :: GetInSize => get_in_size               !<获得网格流向网格点数
        procedure :: HasBodyFixed                !< 检查是否贴体网格

    end type

    !> 二维网格类型
    Type, extends(grid_object_type), public :: grid2d_type

!        Integer, private :: In=0 !< 流向网格点数
!        Integer, private :: jn=0  !<法向网格点数
        Type(point_type), allocatable, private :: points(:, :) !网格点坐标
!        logical, private :: isBodyFixed =.False.   !< 是否贴体网格

    Contains

        Procedure :: Create => create_grid2d     !< 分配网格空间
        procedure :: HasCreate=> HasCreate_2d               !< 检查网格空间是否分配
        procedure :: CreateFromPLOT3D => createfromPLOT3D_2d   !<从Plot3D格式构造网格
        Procedure :: Set=> Set_2d                                    !<设置网格点信息
        procedure :: TransBodyFix => trans_body_fix_2d         !<将网格转换为随体网格
        procedure :: Get_iy => Get_Grid2d_iy                !<获得i站位y信息
        procedure :: Get_jx => Get_Grid2d_jx                !<获得j站位x信息
        procedure :: GetSize => get_grid2d_size             !<获得网格流向和法向网格点数
        procedure :: Print=> print_2d                                !<输出网格
        procedure :: Zeros=> zeros_2d !< 置零函数
        procedure :: Finalize => finalizeGrid2d             !<析构函数

    End Type grid2d_type

    !> 三维网格类型
    type, extends(grid_object_type), public :: Grid3d_type

        private
        integer :: kn !< 展向点数
        Type(point_type), allocatable, private :: points(:, :, :) !网格点坐标


        Contains

          Procedure :: Create => Create_Grid3d
          procedure :: HasCreate=> HasCreate_3d
          procedure :: CreateFromPLOT3D => createfromPLOT3D_3d   !<从Plot3D格式构造网格
          procedure :: Set=> Set_3d !<设置网格点信息
          procedure :: TransBodyFix=> trans_body_fix_3d !<将网格转换为随体网格
          procedure :: Get_ijz=> Get_Grid3d_ijz !< 获得3D网格ij点站位沿k方向的z坐标点.
          procedure :: Get_iky=> Get_Grid3d_iky !< 获得3D网格ik点站位沿j方向的y坐标点.
          procedure :: Get_jkx=> Get_Grid3d_jkx !< 获得3D网格jk点站位沿k方向的x坐标点.
          procedure :: GetSize=> get_grid3d_size !< 获得网格流向法向和展向的网格点数.
          procedure :: Print=> print_3d  !<输出网格
          procedure :: Zeros=> zeros_3d !< 置零函数

    end type Grid3d_type

Contains

    !> 判断网格空间是否分配.
    !! @return 网格空间是否分配
    !! @retval TRUE 已分配
    !! @retval FALSE 未分配
    function HasCreate_2d(this) result(is)

        implicit none
        class(grid2d_type), intent(in) :: this
        logical :: is

        is = allocated(this%points)

    end function HasCreate_2d

    !> 判断网格空间是否分配.
    !! @return 网格空间是否分配
    !! @retval TRUE 已分配
    !! @retval FALSE 未分配
    function HasCreate_3d(this) result(is)

        implicit none
        class(grid3d_type), intent(in) :: this
        logical :: is

        is = allocated(this%points)

    end function HasCreate_3d

    !> 判断网格是否贴体网格.
    !! @return 网格是否贴体
    !! @retval TRUE 是贴体网格
    !! @retval FALSE 不是贴体网格
    function HasBodyFixed(this) result(is)

        implicit none
        class(grid_object_type), intent(inout) :: this
        logical :: is

        is = this%isBodyFixed

    end function HasBodyFixed

    !>   为2D网格分配内存.
    !!   @param[in] in 流向网格点数
    !!   @param[in] jn 法向网格点数
    Subroutine  create_grid2d(this, in, jn)

        Implicit None
        Integer in, jn
        Class(grid2d_type) :: this

        this%In=in
        this%jn=jn

        if(.not. allocated(this%points)) Allocate(this%points(in, jn))

    End Subroutine create_grid2d

    !>   为3D网格分配内存.
    !!   @param[in] in 流向网格点数
    !!   @param[in] jn 法向网格点数
    !!   @param[in] kn 展向网格点数
    Subroutine  create_grid3d(this, in, jn, kn)

        Implicit None
        Integer, intent(in) :: in, jn, kn
        Class(grid3d_type) :: this

        this%In=in
        this%jn=jn
        this%kn=kn

        if(.not. allocated(this%points)) Allocate(this%points(in, jn, kn))

    End Subroutine create_grid3d

    !> 从Plot3D格式文件创建网格并设置网格点信息
    !! @param[in] filename Plot3D网格文件名

    subroutine createfromPLOT3D_2d(this, filename)

        use stringifor
        implicit none
        type(string), intent(in) :: filename
        class(grid2d_type), intent(inout) :: this
        integer :: in , jn, kn
        real(R_P), allocatable :: xx(:, :, :), yy(:, :, :), zz(:, :, :)
        integer :: l

        open(11, file=filename//'', form='unformatted')
        read(11)
        read(11)in, jn!, kn
        kn=1
        print*, in, jn, kn

        if(kn .ne. 1) stop "kn must be 1 for 2D problem"
        call this%Create(in, jn)
        if(.not. allocated(xx)) &
        &   allocate(xx(in, jn, 1), yy(in, jn, 1), zz(in, jn, 1) )
        read(11)xx(:, :, :), yy(:, :, :)!, zz(:, :, :)
        zz=0.0d0
        close(11)
        call this%Set(xx(:, :, 1), yy(:, :, 1), zz(:, :, 1))
        deallocate(xx, yy, zz)

    end subroutine createfromPLOT3D_2d

    !> 从Plot3D格式文件创建网格并设置网格点信息
    !! @param[in] filename Plot3D网格文件名

    subroutine createfromPLOT3D_3d(this, filename)

        use stringifor
        implicit none
        type(string), intent(in) :: filename
        class(grid3d_type), intent(inout) :: this
        integer :: in , jn, kn
        real(R_P), allocatable :: xx(:, :, :), yy(:, :, :), zz(:, :, :)

        open(11, file=filename//'', form='unformatted')
        read(11)
        read(11)in, jn!, kn
        kn=1
        print*, in, jn, kn
        if(kn == 1) stop "kn must be > 1 for 3D problem"
        call this%Create(in, jn, kn)
        if(.not. allocated(xx)) &
        &   allocate(xx(in, jn, 1), yy(in, jn, kn), zz(in, jn, kn) )
        read(11)xx(:, :, :), yy(:, :, :)!, zz(:, :, :)
        zz=0.0d0
        close(11)
        call this%Set(xx(:, :, :), yy(:, :, :), zz(:, :, :))
        deallocate(xx, yy, zz)

    end subroutine createfromPLOT3D_3d

    !>   设置2D网格网格点信息.
    !!
    !!   @param[in] xx 网格x点坐标
    !!   @param[in] yy 网格y点坐标
    !!   @param[in] zz 网格z点坐标
    Subroutine set_2d(this, xx, yy, zz)

        Implicit None
        Class(grid2d_type) :: this
        real(R_P), Dimension(:, :) :: xx, yy, zz

        Call this%points%Set(xx, yy, zz)

    End Subroutine set_2d

    !>   设置3D网格网格点信息.
    !!
    !!   @param[in] xx 网格x点坐标
    !!   @param[in] yy 网格y点坐标
    !!   @param[in] zz 网格z点坐标
    Subroutine set_3d(this, xx, yy, zz)

        Implicit None
        Class(grid3d_type) :: this
        real(R_P), Dimension(:, :, :) :: xx, yy, zz

        Call this%points%Set(xx, yy, zz)

    End Subroutine set_3d

    !>   获得2D网格网格点信息.
    !!
    !!   @param[out] xx 网格x点坐标
    !!   @param[out] yy 网格y点坐标
    !!   @param[out] zz 网格z点坐标
    Subroutine get_grid2d(this, xx, yy, zz)
        Implicit None
        Class(grid2d_type) :: this
        real(R_P) :: xx(:, :), yy(:, :), zz(:, :)

        Call this%Points%Get(xx, yy, zz)

    End Subroutine get_grid2d

    !> 获得2D网格i点站位沿j方向的y坐标点.
    !!
    !! @param[in] i i站位序号
    !! @param[out] yy 沿j方向的y坐标点
    subroutine Get_Grid2d_iy(this, i, yy)

        implicit none
        integer, intent(in) :: i
        class(grid2d_type), intent(in) :: this
        real(R_P) :: yy(:), xx, zz
        integer :: J

        do j=1, size(yy)
            call this%Points(i, j)%Get(xx, yy(j), zz)
        enddo

    end subroutine Get_Grid2d_iy

    !> 获得2D网格j点站位沿i方向的x坐标点.
    !!
    !! @param[in] j j站位序号
    !! @param[out] xx 沿i方向的x坐标点
    subroutine Get_Grid2d_jx(this, j, xx)

        implicit none
        integer, intent(in) :: j
        class(grid2d_type), intent(in) :: this
        integer :: i
        real(R_P) :: xx(:), yy, zz

        do i=1, size(xx)
            call this%Points(i, j)%Get(xx(i), yy, zz)
        enddo

    end subroutine Get_Grid2d_jx

    !> 获得3D网格ik点站位沿j方向的y坐标点.
    !!
    !! @param[in] i i站位序号
    !! @param[in] k k站位序号
    !! @param[out] yy 沿j方向的y坐标点
    subroutine Get_Grid3D_iky(this, i, k, yy)

        implicit none
        integer, intent(in) :: i, k
        class(grid3d_type), intent(in) :: this
        real(R_P) :: yy(:), xx, zz
        integer :: J

        do j=1, size(yy)
            call this%Points(i, j, k)%Get(xx, yy(j), zz)
        enddo

    end subroutine Get_Grid3D_iky

    !> 获得3D网格jk点站位沿i方向的x坐标点.
    !!
    !! @param[in] j j站位序号
    !! @param[in] k k站位序号
    !! @param[out] xx 沿i方向的x坐标点
    subroutine Get_Grid3d_jkx(this, j, k, xx)

        implicit none
        integer, intent(in) :: j, k
        class(grid3d_type), intent(in) :: this
        integer :: i
        real(R_P) :: xx(:), yy, zz

        do i=1, size(xx)
            call this%Points(i, j, k)%Get(xx(i), yy, zz)
        enddo

    end subroutine Get_Grid3d_jkx

    !> 获得3D网格ij点站位沿k方向的z坐标点.
    !!
    !! @param[in] i i站位序号
    !! @param[in] j j站位序号
    !! @param[out] zz 沿k方向的z坐标点
    subroutine Get_Grid3d_ijz(this, i, j, zz)

        implicit none
        integer, intent(in) :: i, j
        class(grid3d_type), intent(in) :: this
        integer :: k
        real(R_P) :: xx, yy, zz(:)

        do k=1, size(zz)
            call this%Points(i, j, k)%Get(xx, yy, zz(k))
        enddo

    end subroutine Get_Grid3d_ijz

    !> 获得2D网格流向和法向点数.
    !>
    !! @param[out] in i方向点数
    !! @param[out] jn j方向点数
    subroutine get_grid2d_size(this, in, jn)

        implicit none
        class(grid2d_type), intent(in) :: this
        integer , intent(inout) :: in, jn

        in=this%In
        jn=this%jn

    end subroutine get_grid2d_size

    !> 获得3D网格流向和法向点数.
    !>
    !! @param[out] in i方向点数
    !! @param[out] jn j方向点数
    !! @param[out] kn k方向点数
    subroutine get_grid3d_size(this, in, jn, kn)

        implicit none
        class(grid3d_type), intent(in) :: this
        integer , intent(inout) :: in, jn, kn

        in=this%In
        jn=this%jn
        kn=this%kn

    end subroutine get_grid3d_size

    !> 获得2D网格法向点数.
    !>
    !! @retval jn j方向点数
    function get_jn_size(this) result(jn)

        implicit none
        class(grid_object_type), intent(in) :: this
        integer :: jn

        jn=this%jn

    end function get_jn_size

    !> 获得2D网格流向点数.
    !>
    !! @retval in i方向点数
    function get_in_size(this) result(in)

        implicit none
        class(grid_object_type), intent(in) :: this
        integer :: in

        in=this%in

    end function get_in_size

    !> 将网格变为随体网格表达
    !!
    !! @return 随体网格
    function trans_body_fix_2d(this) result(grid)

        implicit none
        class(grid2d_type), intent(in) :: this
        type(grid2d_type), pointer :: grid
        integer :: i, j
        real(R_P) :: hypot_i, hypot_j
        type(point_type) :: PointHypot

        i=0; j=0

        if(associated(grid)) grid=>null()
        allocate(grid)
        grid%in=this%in; grid%jn=this%jn
        allocate(grid%points(grid%in, grid%jn))

        call grid%points(1, :)%Set(0.0d0, 0.0d0, 0.0d0)
        call grid%points(:, 1)%Set(0.0d0, 0.0d0, 0.0d0)

        do j=1, grid%jn
            do i=2, grid%in
                hypot_i=.abs.(this%points(i, j) - this%points(i-1, j))
                call PointHypot%Set(hypot_i, 0.0d0, 0.0d0)
                grid%points(i, j)=grid%points(i-1, j) + PointHypot
            enddo
        end do
        do j=2, grid%jn
            do i=1, grid%in
                hypot_j=.abs.(this%points(i, j) - this%points(i, j-1))
                call PointHypot%Set(0.0d0, hypot_j, 0.0d0)
                grid%points(i, j)=grid%points(i, j-1) + PointHypot
            enddo
        end do

        grid%isBodyfixed=.True.

    end function trans_body_fix_2d

    !> 将网格变为随体网格表达
    !!
    !! @return 随体网格
    function trans_body_fix_3d(this) result(grid)

        implicit none
        class(grid3d_type), intent(in) :: this
        type(grid3d_type), pointer :: grid
        integer :: i, j, k
        real(R_P) :: hypot_i, hypot_j, hypot_k
        type(point_type) :: PointHypot

        i=0; j=0; k=0

        if(associated(grid)) grid=>null()
        allocate(grid)
        grid%in=this%in; grid%jn=this%jn; grid%kn=this%kn
        allocate(grid%points(grid%in, grid%jn, grid%kn))

        call grid%points(1, :, 1)%Set(0.0d0, 0.0d0, 0.0d0)
        call grid%points(:, 1, 1)%Set(0.0d0, 0.0d0, 0.0d0)
        call grid%points(1, 1, :)%Set(0.0d0, 0.0d0, 0.0d0)

        do k=1, grid%kn
            do j=1, grid%jn
                do i=2, grid%in
                    hypot_i=.abs.(this%points(i, j,  k) - this%points(i-1, j, k))
                    call PointHypot%Set(hypot_i, 0.0d0, 0.0d0)
                    grid%points(i, j, k)=grid%points(i-1, j, k) + PointHypot
                enddo
            end do
            do j=2, grid%jn
                do i=1, grid%in
                    hypot_j=.abs.(this%points(i, j, k) - this%points(i, j-1, k))
                    call PointHypot%Set(0.0d0, hypot_j, 0.0d0)
                    grid%points(i, j, k)=grid%points(i, j-1, k) + PointHypot
                enddo
            end do
        enddo
        do j=1, grid%jn
            do i=1, grid%in
                do k=2, grid%kn
                    hypot_k=.abs.(this%points(i, j, k) - this%points(i, j, k-1))
                    call PointHypot%Set(0.0d0, 0.0d0, hypot_j)
                    grid%points(i, j, k)=grid%points(i, j, k-1) + PointHypot
                end do
            end do
        end do

        grid%isBodyfixed=.True.

    end function trans_body_fix_3d

    !> 输出网格到文件.
    !! @param[in] fn_surf 文件前缀
    !!
    subroutine print_2d(this, fn_surf)

        use stringifor
        implicit none
        class(grid2d_type), intent(in) :: this
        type(string), intent(in) :: fn_surf
        type(string) :: fn
        real(R_P), dimension(this%in, this%jn) :: xx, yy, zz
        integer :: i, j

        call this%points%Get(xx, yy, zz)
        fn=fn_surf//'_Grid.dat'
        open(998, file=fn%chars(), form='unformatted')
        write(998)1
        write(998)this%in, this%jn
        write(998)((xx(i, j), i=1, this%in), j=1, this%jn), &
            ((yy(i, j), i=1, this%in), j=1, this%jn)
        close(998)

    end subroutine print_2d

    !> 输出网格到文件.
    !! @param[in] fn_surf 文件前缀
    !!
    subroutine print_3d(this, fn_surf)

        use stringifor
        implicit none
        class(grid3d_type), intent(in) :: this
        type(string), intent(in) :: fn_surf
        type(string) :: fn
        real(R_P), dimension(this%in, this%jn, this%kn) :: xx, yy, zz
        integer :: i, j, k

        call this%points%Get(xx, yy, zz)
        fn=fn_surf//'_Grid.dat'
        open(998, file=fn%chars(), form='unformatted')
        write(998)1
        write(998)this%in, this%jn, this%kn
        write(998)(((xx(i, j, k), i=1, this%in), j=1, this%jn), k=1, this%kn), &
                  (((yy(i, j, k), i=1, this%in), j=1, this%jn), k=1, this%kn), &
                  (((zz(i, j, k), i=1, this%in), j=1, this%jn), k=1, this%kn)
        close(998)

    end subroutine print_3d

    !> 置零函数.
    elemental subroutine Zeros_2d(this)

        implicit none
        class(grid2d_type), intent(inout) :: this

        call this%points%Zeros()

    end subroutine Zeros_2d

    !> 置零函数.
    elemental subroutine Zeros_3d(this)

        implicit none
        class(grid3d_type), intent(inout) :: this

        call this%points%Zeros()

    end subroutine Zeros_3d

    !> 析构函数
    Subroutine finalizeGrid2d(this)

        Implicit None
        class(grid2d_type) :: this
        this%In=0; this%jn=0; this%isBodyFixed=.False.
        call this%points%Finalize()

    End Subroutine finalizeGrid2d

End Module mod_grid
