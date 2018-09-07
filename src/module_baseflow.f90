!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_baseflow.f90
!> @file
!> @breif 基本流类文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  Module: baseflow_org
!
!  DESCRIPTION:
!> @breif 二维基本流模块.
!>
!!
!  REVISION HISTORY:
!  2017-07-25 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-07-25
!------------------------------------------------------------------------------
module mod_baseflow

    use mod_baseflow_org
    !    use mod_baseflow_grad
    use mod_grid
    use mod_vector
    use penf, only: R_P
    implicit none

    private

    !> 基本Baseflow类型
    type :: baseflow_object_type

        private
        integer                               :: in                             !< 流向点数\private
        integer                               :: jn                             !< 法向点数\private

    end type baseflow_object_type


    !>  二维基本流.
    !!
    type, extends(baseflow_object_type), public :: Baseflow2D_type

        private
        type(bf_flux_org_ij_type), public, allocatable :: Flux_org(:, :)                  !< 基本流通量\private

    Contains

        Procedure :: CreateFromGrid  => create_baseflow2d_From_Grid           !< 从网格分配流场内存
        procedure :: Create          => create_baseflow2d                     !<从in, jn分配流场内存
        procedure :: CreateFromPLOT3D => create_from_PLOT3D                    !<从Plot3D文件读取流场
        procedure :: IsCreate        => IsCreate_BF                           !<检查是否创建流场
        procedure :: Get             => get_baseflow2d                        !< 读取全部二维基本流
        procedure :: GetPoint        => get_PointFlux                         !< 读取某一个点的基本流
        procedure :: GetPart         => get_baseflow2dPartFlux                !<读取某一部分的基本流
!        procedure :: GetIlocFlux     => get_iloc                              !<读取流向站位序号的流场
        procedure :: GetSize         => get_size                              !<读取in和jn
        procedure :: finalize        => finalize_baseflow2d                   !<析构类型

        generic  :: Set              => set_baseflow2d, set_baseflow2d_dble    !<设置流场通量值

        procedure, private :: set_baseflow2d                                   !<设置流场通量
        procedure, private :: set_baseflow2d_dble                              !<设置流场通量

    end type Baseflow2D_type

    !>  二维基本流导数.
    !!
    type, private :: Baseflow2D_Diff_type

        private
        type(bf_flux_org_ij_type), allocatable :: DxFlux_org(:, :)                  !<流向一阶导数\private
        type(bf_flux_org_ij_type), allocatable :: DyFlux_org(:, :)                  !<法向一阶导数\private
        type(bf_flux_org_ij_type), allocatable :: DyyFlux_org(:, :)                 !<法向二阶导数\private

    contains

        procedure :: CreateDiff       => Create_Diff                              !<创建基本流导数
        procedure :: IsCreate         => IsCreate_BFDiff                          !<检查是否创建流场导数
        procedure :: CreateDiffFromFile                                            !<从plot3d文件读取基本流导数
        procedure :: GetPartial_Xi    => Partial_Xi_2d                               !<读取基本流流向一阶导数
        procedure :: GetPartial_Eta   => Partial_Eta_2d                              !<读取基本流法向一阶导数
        procedure :: GetPartial_Eta2  => Partial_Eta2_2d                             !<读取基本流法向二阶导数
        procedure :: finalize         => finalize_baseflow2d_diff                 !<析构类型

    end type Baseflow2D_Diff_type

    !>  三维基本流.
    !!
    type, extends(baseflow_object_type), private :: Baseflow3D_type

        private
        integer :: kn !< 展向点数
        type(bf_flux_org_ij_type), public, allocatable :: Flux_org(:, :, :)                  !< 基本流通量\private

    Contains

        Procedure :: CreateFromGrid  => create_baseflow3d_From_Grid           !< 从网格分配流场内存
        procedure :: Create          => create_baseflow3d                     !<从in, jn, kn分配流场内存
        procedure :: CreateFromPLOT3D => create_from_PLOT3D_3d                    !<从Plot3D文件读取流场
        procedure :: IsCreate        => IsCreate_BF_3d                           !<检查是否创建流场
        procedure :: Get             => get_baseflow3d                        !< 读取全部三维基本流
        procedure :: GetPoint        => get_PointFlux_3d                         !< 读取某一个点的基本流
        procedure :: GetPart         => get_baseflow3dPartFlux                !<读取某一部分的基本流
!        procedure :: GetIloc         => get_iloc_3d                              !<读取流向站位序号
        procedure :: GetSize         => get_size_3d                              !<读取in, jn和kn
        procedure :: finalize        => finalize_baseflow3d                   !<析构类型

        generic  :: Set              => set_baseflow3d, set_baseflow3d_dble    !<设置流场通量值

        procedure, private :: set_baseflow3d                                   !<设置流场通量
        procedure, private :: set_baseflow3d_dble                              !<设置流场通量

    end type Baseflow3D_type

    !>  三维基本流导数.
    !!
    type, public :: Baseflow3D_Diff_type

        private
        type(bf_flux_org_ij_type), allocatable :: DxFlux_org(:, :, :)                  !<流向一阶导数\private
        type(bf_flux_org_ij_type), allocatable :: DyFlux_org(:, :, :)                  !<法向一阶导数\private
        type(bf_flux_org_ij_type), allocatable :: DzFlux_org(:, :, :)                  !<展向一阶导数\private
        type(bf_flux_org_ij_type), allocatable :: DyyFlux_org(:, :, :)                 !<法向二阶导数\private
        type(bf_flux_org_ij_type), allocatable :: DyzFlux_org(:, :, :)                 !<法向展向二阶导数\private
        type(bf_flux_org_ij_type), allocatable :: DzzFlux_org(:, :, :)                 !<展向二阶导数\private

    contains

        procedure :: CreateDiff       => Create_Diff_3D                              !<创建基本流导数
        procedure :: IsCreate         => IsCreate_BFDiff_3D                          !<检查是否创建流场导数
        procedure :: CreateDiffFromFile_3D                                           !<从plot3d文件读取基本流导数
        procedure :: GetPartial_Xi    => Partial_Xi_3D                               !<读取基本流流向一阶导数
        procedure :: GetPartial_Eta   => Partial_Eta_3D                              !<读取基本流法向一阶导数
        procedure :: GetPartial_ZEta  => Partial_ZEta_3D                             !<读取基本流展向一阶导数
        procedure :: GetPartial_Eta2  => Partial_Eta2_3D                             !<读取基本流法向二阶导数
        procedure :: GetPartial_Zeta2 => Partial_Zeta2_3D                            !<读取基本流展向二阶导数
        procedure :: GetPartial_EtaZeta=> Partial_EtaZeta_3D                         !<读取基本流法向展向交叉二阶导数
        procedure :: finalize         => finalize_baseflow3d_diff                 !<析构类型

    end type Baseflow3D_Diff_type


contains

    !>  从plot3d文件创建并读取基本流导数.
    !!
    !!  @param[in] filename 基本流导数文件名
    subroutine CreateDiffFromFile(this, filename)

        use mod_difference
        use mod_vector
        use stringifor
        implicit none
        class(Baseflow2D_Diff_type), intent(inout) :: this
        type(string)                             :: filename
        integer                                  :: in, jn
        real(R_P), dimension(:, :), allocatable  :: rbx, ubx, vbx, wbx, tbx, &
            rby, uby, vby, wby, tby, rbyy, ubyy, vbyy, wbyy, tbyy
        type(vector_type)                        :: vel
        integer                                  :: i, j

        !call diff%PartialXI(1, jn, BF%Flux_org, this%DxFlux_org)
        !call diff%PartialEta(1, in, BF%Flux_org, this%DyFlux_org)
        !call diff%PartialEta2(1, in, BF%Flux_org, this%DyyFlux_org)
        open(1911, file=filename//'', form='unformatted')
        read(1911)
        read(1911)in, jn
        write(*, *) in, jn
        if(.not. allocated(this%DxFlux_org))then
          allocate(this%DxFlux_org(in, jn))
          allocate(this%DyFlux_org(in, jn))
          allocate(this%DyyFlux_org(in, jn))
          call this%DxFlux_org%Zeros()
          call this%DyFlux_org%Zeros()
          call this%DyyFlux_org%Zeros()
        endif
        if(.not. allocated(rbx)) &
        &  allocate(rbx(in, jn), ubx(in, jn), vbx(in, jn), wbx(in, jn), &
              tbx(in, jn), rby(in, jn), uby(in, jn), vby(in, jn), wby(in, jn), &
              tby(in, jn), rbyy(in, jn), ubyy(in, jn), vbyy(in, jn), &
              wbyy(in, jn), tbyy(in, jn))
        read(1911)rbx, ubx, vbx, wbx, tbx, &
            rby, uby, vby, wby, tby, rbyy, ubyy, vbyy, wbyy, tbyy
        close(1911)

        do j=1, jn
            do i=1, in
                call vel%Set(ubx(i, j), vbx(i, j), wbx(i, j))
                call this%DxFlux_org(i, j)%Set(rbx(i, j), vel, tbx(i, j))
                call vel%Set(uby(i, j), vby(i, j), wby(i, j))
                call this%DyFlux_org(i, j)%Set(rby(i, j), vel, tby(i, j))
                call vel%Set(ubyy(i, j), vbyy(i, j), wbyy(i, j))
                call this%DyyFlux_org(i, j)%Set(rbyy(i, j), vel, tbyy(i, j))
            end do
        end do

        deallocate(rbx, ubx, vbx, wbx, tbx, &
            rby, uby, vby, wby, tby, rbyy, ubyy, vbyy, wbyy, tbyy)

    end subroutine CreateDiffFromFile

    !>  从plot3d文件创建并读取基本流导数.
    !!
    !!  @param[in] filename 基本流导数文件名
    subroutine CreateDiffFromFile_3D(this, filename)

        use mod_difference
        use mod_vector
        use stringifor
        implicit none
        class(Baseflow3D_Diff_type), intent(inout) :: this
        type(string)                             :: filename
        integer                                  :: in, jn, kn
        real(R_P), dimension(:, :, :), allocatable  :: rbx, ubx, vbx, wbx, tbx, &
            rby, uby, vby, wby, tby, rbz, ubz, vbz, wbz, tbz, &
            rbyy, ubyy, vbyy, wbyy, tbyy, &
            rbyz, ubyz, vbyz, wbyz, tbyz, &
            rbzz, ubzz, vbzz, wbzz, tbzz
        type(vector_type)                        :: vel
        integer                                  :: i, j, k

        !call diff%PartialXI(1, jn, BF%Flux_org, this%DxFlux_org)
        !call diff%PartialEta(1, in, BF%Flux_org, this%DyFlux_org)
        !call diff%PartialEta2(1, in, BF%Flux_org, this%DyyFlux_org)
        open(1911, file=filename//'', form='unformatted')
        read(1911)
        read(1911)in, jn, kn
        write(*, *) in, jn, kn
        if(.not. allocated(this%DxFlux_org))then
          allocate(this%DxFlux_org(in, jn, kn))
          allocate(this%DyFlux_org(in, jn, kn))
          allocate(this%DzFlux_org(in, jn, kn))
          allocate(this%DyyFlux_org(in, jn, kn))
          allocate(this%DyzFlux_org(in, jn, kn))
          allocate(this%DzzFlux_org(in, jn, kn))
          call this%DxFlux_org%Zeros()
          call this%DyFlux_org%Zeros()
          call this%DzFlux_org%Zeros()
          call this%DyyFlux_org%Zeros()
          call this%DyzFlux_org%Zeros()
          call this%DzzFlux_org%Zeros()
        endif
        if(.not. allocated(rbx)) &
        &  allocate(rbx(in, jn, kn), ubx(in, jn, kn), vbx(in, jn, kn), wbx(in, jn, kn), &
              tbx(in, jn, kn), rby(in, jn, kn), uby(in, jn, kn), vby(in, jn, kn), wby(in, jn, kn), &
              tby(in, jn, kn), rbyy(in, jn, kn), ubyy(in, jn, kn), vbyy(in, jn, kn), &
              wbyy(in, jn, kn), tbyy(in, jn, kn), &
              rbz(in, jn, kn), ubz(in, jn, kn), vbz(in, jn, kn), wbz(in, jn, kn), tbz(in, jn, kn), &
              rbyz(in, jn, kn), ubyz(in, jn, kn), vbyz(in, jn, kn), wbyz(in, jn, kn), tbyz(in, jn, kn), &
              rbzz(in, jn, kn), ubzz(in, jn, kn), vbzz(in, jn, kn), wbzz(in, jn, kn), tbzz(in, jn, kn))
        read(1911)rbx, ubx, vbx, wbx, tbx, &
                  rby, uby, vby, wby, tby, &
                  rbz, ubz, vbz, wbz, tbz, &
                  rbyy, ubyy, vbyy, wbyy, tbyy, &
                  rbyz, ubyz, vbyz, wbyz, tbyz, &
                  rbzz, ubzz, vbzz, wbzz, tbzz
        close(1911)

        do k=1, kn
            do j=1, jn
                do i=1, in
                    call vel%Set(ubx(i, j, k), vbx(i, j, k), wbx(i, j, k))
                    call this%DxFlux_org(i, j, k)%Set(rbx(i, j, k), vel, tbx(i, j, k))
                    call vel%Set(uby(i, j, k), vby(i, j, k), wby(i, j, k))
                    call this%DyFlux_org(i, j, k)%Set(rby(i, j, k), vel, tby(i, j, k))
                    call vel%Set(ubz(i, j, k), vbz(i, j, k), wbz(i, j, k))
                    call this%DzFlux_org(i, j, k)%Set(rbz(i, j, k), vel, tbz(i, j, k))
                    call vel%Set(ubyy(i, j, k), vbyy(i, j, k), wbyy(i, j, k))
                    call this%DyyFlux_org(i, j, k)%Set(rbyy(i, j, k), vel, tbyy(i, j, k))
                    call vel%Set(ubyz(i, j, k), vbyz(i, j, k), wbyz(i, j, k))
                    call this%DyzFlux_org(i, j, k)%Set(rbyz(i, j, k), vel, tbyz(i, j, k))
                    call vel%Set(ubzz(i, j, k), vbzz(i, j, k), wbzz(i, j, k))
                    call this%DzzFlux_org(i, j, k)%Set(rbzz(i, j, k), vel, tbzz(i, j, k))
                end do
            end do
        enddo

        deallocate(rbx, ubx, vbx, wbx, tbx, &
            rby, uby, vby, wby, tby, &
            rbz, ubz, vbz, wbz, tbz, &
            rbyy, ubyy, vbyy, wbyy, tbyy, &
            rbyz, ubyz, vbyz, wbyz, tbyz, &
            rbzz, ubzz, vbzz, wbzz, tbzz)

    end subroutine CreateDiffFromFile_3D

    !>   利用in, jn信息初始化2D流场.
    !!
    !!   \param[in] in 流向网格点数
    !!   \param[in] jn 法向网格点数
    subroutine create_baseflow2d(this, in, jn)

        implicit none
        class(Baseflow2D_type), intent(inout) :: this
        integer, intent(in)                  :: in, jn

        if(.not. allocated(this%flux_org))then
            this%in=in
            this%jn=jn
            allocate(this%Flux_org(in, jn))
            call this%Flux_org%Zeros()
        endif

    end subroutine create_baseflow2d

    !>   利用in, jn, kn信息初始化3D流场.
    !!
    !!   \param[in] in 流向网格点数
    !!   \param[in] jn 法向网格点数
    !!   \param[in] kn 展向网格点数
    subroutine create_baseflow3d(this, in, jn, kn)

        implicit none
        class(Baseflow3D_type), intent(inout) :: this
        integer, intent(in)                  :: in, jn, kn

        if(.not. allocated(this%flux_org))then
            this%in=in
            this%jn=jn
            this%kn=kn
            allocate(this%Flux_org(in, jn, kn))
            call this%Flux_org%Zeros()
        endif

    end subroutine create_baseflow3d

    !>   利用2D网格信息初始化2D流场.
    !!
    !!   \param[in] grid 2D网格
    subroutine create_baseflow2d_from_grid(this, grid)

        implicit none
        class(Baseflow2D_type), intent(inout) :: this
        type(grid2d_type)     , intent(in)   :: grid
        integer                              :: in,  jn

        if(.not. allocated(this%flux_org) )then
            Call Grid%GetSize(in, jn)
            this%in=in
            this%jn=jn
            allocate(this%Flux_org(in, jn))
            call this%Flux_org%Zeros()
        endif

    end subroutine create_baseflow2d_from_grid

    !>   利用3D网格信息初始化3D流场.
    !!
    !!   \param[in] grid 3D网格
    subroutine create_baseflow3d_from_grid(this, grid)

        implicit none
        class(Baseflow3D_type), intent(inout) :: this
        type(grid3d_type)     , intent(in)   :: grid
        integer                              :: in,  jn, kn

        if(.not. allocated(this%flux_org) )then
            Call Grid%GetSize(in, jn, kn)
            this%in=in
            this%jn=jn
            this%kn=kn
            allocate(this%Flux_org(in, jn, kn))
            call this%Flux_org%Zeros()
        endif

    end subroutine create_baseflow3d_from_grid

    !>  从基本流创建并计算基本流导数.
    !!
    !!  @param[in] BF 2D基本流
    !!  @param[in] diff 基本流导数计算的基函数
    !!  @return 创建基本流导数
    subroutine create_diff(this, BF, diff)

        use mod_difference
        implicit none
        class(Baseflow2D_Diff_type), intent(inout) :: this
        class(Baseflow2D_type)     , intent(in)   :: BF
        type(difference_2D_type)   , intent(in)   :: diff
        integer                                   :: in, jn

        in=BF%in; jn=BF%jn
        if(.not. allocated(this%DxFlux_org)) then
            allocate(this%DxFlux_org (in, jn))
            allocate(this%DyFlux_org (in, jn))
            allocate(this%DyyFlux_org(in, jn))
            call this%DxFlux_org%Zeros()
            call this%DyFlux_org%Zeros()
            call this%DyyFlux_org%Zeros()
        endif
        call diff%PartialXI  (1, jn, BF%Flux_org, this%DxFlux_org)
        call diff%PartialEta (1, in, BF%Flux_org, this%DyFlux_org)
        call diff%PartialEta2(1, in, BF%Flux_org, this%DyyFlux_org)

    end subroutine create_diff

    !>  从基本流创建并计算基本流导数.
    !!
    !!  @param[in] BF 3D基本流
    !!  @param[in] diff 基本流导数计算的基函数
    !!  @return 创建基本流导数
    subroutine create_diff_3d(this, BF, diff)      !TODO the comment must be deleted

        use mod_difference
        implicit none
        class(Baseflow3D_Diff_type), intent(inout) :: this
        class(Baseflow3D_type)     , intent(in)   :: BF
        type(difference_3D_type)   , intent(in)   :: diff
        integer                                   :: in, jn, kn
        integer :: k, i, j

        in=BF%in; jn=BF%jn; kn=BF%kn
        if(.not. allocated(this%DxFlux_org)) then
            allocate(this%DxFlux_org (in, jn, kn))
            allocate(this%DyFlux_org (in, jn, kn))
            allocate(this%DzFlux_org (in, jn, kn))
            allocate(this%DyyFlux_org(in, jn, kn))
            allocate(this%DyzFlux_org(in, jn, kn))
            allocate(this%DzzFlux_org(in, jn, kn))
            call this%DxFlux_org%Zeros()
            call this%DyFlux_org%Zeros()
            call this%DzFlux_org%Zeros()
            call this%DyyFlux_org%Zeros()
            call this%DyzFlux_org%Zeros()
            call this%DzzFlux_org%Zeros()
        endif
        do k=1, kn
            call diff%PartialXI  (1, jn, BF%Flux_org(:,:,k), this%DxFlux_org(:,:,k))
            call diff%PartialEta (1, in, BF%Flux_org(:,:,k), this%DyFlux_org(:,:,k))
            call diff%PartialEta2(1, in, BF%Flux_org(:,:,k), this%DyyFlux_org(:,:,k))
        enddo
        do i=1, in
            do j=1, jn
                call diff%PartialZeta (BF%Flux_org(i, j, :), this%DzFlux_org(i, j, :))
                call diff%PartialZeta (this%DyFlux_org(i, j, :), this%DyzFlux_org(i, j, :))
                call diff%PartialZEta2(BF%Flux_org(i, j, :), this%DzzFlux_org(i, j, :))
            enddo
        enddo

    end subroutine create_diff_3d

    !>  从plot3d文件读取基本流流场信息.
    !!
    !!  @param[in] filename plot3d流场文件名
    subroutine create_from_PLOT3D(this, filename)

        use stringifor
        implicit none
        type(string),        intent(in)   :: filename
        class(Baseflow2D_type), intent(inout) :: this
        real(R_P), allocatable              :: Q(:, :, :, :)
        integer                             :: in, jn, kn, ln
        integer                             :: i, j, l

        open(11, file=filename//'', form='unformatted')
        read(11)
        read(11)in, jn, ln
        kn=1
        print*, in, jn, kn, ln
        if(kn .ne. 1) stop "kn must be 1 for 2DPSE"
        call this%Create(in, jn)
        if(.not. allocated(Q)) allocate(Q(in, jn, kn, ln))
        read(11)Q
        close(11)
        call this%Set(Q(:, :, 1, 1), Q(:, :, 1, 2), Q(:, :, 1, 3), Q(:, :, 1, 4), &
            &  Q(:, :, 1, 5))
        deallocate(Q)

    end subroutine create_from_PLOT3D

    !>  从plot3d文件读取基本流流场信息.
    !!
    !!  @param[in] filename plot3d流场文件名
    subroutine create_from_PLOT3D_3d(this, filename)

        use stringifor
        implicit none
        type(string),        intent(in)   :: filename
        class(Baseflow3D_type), intent(inout) :: this
        real(R_P), allocatable              :: Q(:, :, :, :)
        integer                             :: in, jn, kn, ln
        integer                             :: i, j, k, l

        open(11, file=filename//'', form='unformatted')
        read(11)
        read(11)in, jn, ln
        kn=1
        print*, in, jn, kn, ln
        if(kn .eq. 1) stop "kn must be greater than 1 for 3DPSE"
        call this%Create(in, jn, kn)
        if(.not. allocated(Q)) allocate(Q(in, jn, kn, ln))
        read(11)Q
        close(11)
        call this%Set(Q(:, :, :, 1), Q(:, :, :, 2), Q(:, :, :, 3), Q(:, :, :, 4) , Q(:, :, :, 5))
        deallocate(Q)

    end subroutine create_from_PLOT3D_3D

    !>  读取全部基本流流场通量.
    !!
    !!  @param[out] rho 密度
    !!  @param[out] Vel 速度矢量
    !!  @param[out] T   温度
    !!  @return 全部流场信息
    subroutine get_baseflow2d(this, rho, Vel, T)

        implicit none
        class(Baseflow2D_type) , intent(in)                  :: this
        real(R_P)              , dimension(:, :), intent(inout) :: rho,  T
        type(vector_type), dimension(:, :), intent(inout)     :: Vel
        integer :: i, j

        do j=1, this%jn
            do i=1, this%in
                call this%Flux_org(i, j)%Get(rho(i, j), Vel(i, j), T(i, j))
            enddo
        enddo

    end subroutine get_baseflow2d

    !>  读取全部基本流流场通量.
    !!
    !!  @param[out] rho 密度
    !!  @param[out] Vel 速度矢量
    !!  @param[out] T   温度
    !!  @return 全部流场信息
    subroutine get_baseflow3d(this, rho, Vel, T)

        implicit none
        class(Baseflow3D_type) , intent(in)                  :: this
        real(R_P)              , dimension(:, :, :), intent(inout) :: rho,  T
        type(vector_type), dimension(:, :, :), intent(inout)     :: Vel
        integer :: i, j, k

        do k=1, this%kn
            do j=1, this%jn
                do i=1, this%in
                    call this%Flux_org(i, j, k)%Get(rho(i, j, k), Vel(i, j, k), T(i, j, k))
                enddo
            enddo
        enddo

    end subroutine get_baseflow3d

    !>  读取部分基本流流场通量.
    !!
    !!  @param[in] is 开始流向位置序号
    !!  @param[in] ie 结束流向位置序号
    !!  @param[in] js 开始法向位置序号
    !!  @param[in] je 结束法向位置序号
    !!  @param[out] Flux 流场通量
    !!  @return 返回部分流场通量
    subroutine get_baseflow2dPartFlux(this, is, ie, js, je, Flux)

        implicit none
        class(Baseflow2D_type), intent(in) :: this
        integer, intent(in) :: is, ie, js, je
        type(bf_flux_org_ij_type), intent(inout) :: Flux(is:ie, js:je)

        Flux(is:ie, js:je)=this%Flux_org(is:ie, js:je)

    end subroutine get_baseflow2dPartFlux

    !>  读取部分基本流流场通量.
    !!
    !!  @param[in] is 开始流向位置序号
    !!  @param[in] ie 结束流向位置序号
    !!  @param[in] js 开始法向位置序号
    !!  @param[in] je 结束法向位置序号
    !!  @param[in] ks 开始展向位置序号
    !!  @param[in] ke 结束展向位置序号
    !!  @param[out] Flux 流场通量
    !!  @return 返回部分流场通量
    subroutine get_baseflow3dPartFlux(this, is, ie, js, je, ks, ke, Flux)

        implicit none
        class(Baseflow3D_type), intent(in) :: this
        integer, intent(in) :: is, ie, js, je, ks, ke
        type(bf_flux_org_ij_type), intent(inout) :: Flux(is:ie, js:je, ks:ke)

        Flux(is:ie, js:je, ks:ke)=this%Flux_org(is:ie, js:je, ks:ke)

    end subroutine get_baseflow3dPartFlux

    !!>  获得某流向站位的基本流流场通量.
    !!!
    !!!  @param[in] iloc 流向位置
    !!!  @param[out] Flux 流场通量
    !!!  @return 某站位一列的基本流流场法分布
    !subroutine get_iloc(this, iloc, Flux)
    !
    !    implicit none
    !    class(Baseflow2D_type), intent(in) :: this
    !    integer, intent(in) :: iloc
    !    type(bf_flux_org_ij_type), intent(inout) :: Flux(1:this%jn)
    !    integer :: j
    !    do j=1, this%jn
    !        Flux(j)=this%Flux_org(iloc, j)
    !    enddo
    !
    !end subroutine get_iloc

    !>  读取某点流场的基本流通量.
    !!
    !!  @param[in] iloc 流向位置
    !!  @param[in] jloc 法向位置
    !!  @return 某点的基本流通量
    !!  @retval Flux 该点基本流通量
    function get_PointFlux(this, iloc, jloc) result(Flux)

        implicit none
        class(Baseflow2D_type), intent(in) :: this
        integer, intent(in) :: iloc, jloc
        type(bf_flux_org_ij_type) :: Flux

        Flux=this%Flux_org(iloc, jloc)

    end function get_PointFlux

    !>  读取某点流场的基本流通量.
    !!
    !!  @param[in] iloc 流向位置
    !!  @param[in] jloc 法向位置
    !!  @param[in] kloc 展向位置
    !!  @return 某点的基本流通量
    !!  @retval Flux 该点基本流通量
    function get_PointFlux_3d(this, iloc, jloc, kloc) result(Flux)

        implicit none
        class(Baseflow3D_type), intent(in) :: this
        integer, intent(in) :: iloc, jloc, kloc
        type(bf_flux_org_ij_type) :: Flux

        Flux=this%Flux_org(iloc, jloc, kloc)

    end function get_PointFlux_3d

    !>  获得流场的流向和法向点数.
    !!
    !!  @param[out] in 流向点数
    !!  @param[out] jn 法向点数
    !!  @return 基本流流向和法向点数
    subroutine get_size(this, in, jn)

        implicit none
        class(Baseflow2D_type), intent(in) :: this
        integer , intent(inout) :: in, jn

        in=this%in; jn=this%jn

    end subroutine get_size

    !>  获得流场的流向和法向点数.
    !!
    !!  @param[out] in 流向点数
    !!  @param[out] jn 法向点数
    !!  @param[out] jn 展向点数
    !!  @return 基本流流向和法向点数
    subroutine get_size_3d(this, in, jn, kn)

        implicit none
        class(Baseflow3D_type), intent(in) :: this
        integer , intent(inout) :: in, jn, kn

        in=this%in; jn=this%jn; kn=this%kn

    end subroutine get_size_3d

    !>  检查是否创建基本流流场.
    !!
    !!  @return 流场是否创建
    !!  @retval true 已创建
    !!  @retval false 未创建
    function IsCreate_BF(this) result(Flag)

        implicit none
        class(Baseflow2D_type), intent(in) :: this
        logical :: Flag

        Flag=allocated(this%Flux_org)

    end function IsCreate_BF

    !>  检查是否创建基本流流场.
    !!
    !!  @return 流场是否创建
    !!  @retval true 已创建
    !!  @retval false 未创建
    function IsCreate_BF_3d(this) result(Flag)

        implicit none
        class(Baseflow3D_type), intent(in) :: this
        logical :: Flag

        Flag=allocated(this%Flux_org)

    end function IsCreate_BF_3d

    !>  检查是否创建基本流流场导数.
    !!
    !!  @return 流场导数是否创建
    !!  @retval true 已创建
    !!  @retval false 未创建
    function IsCreate_BFDiff(this) result(is)

        implicit none
        class(Baseflow2D_Diff_type), intent(in) :: this
        logical :: is

        is=allocated(this%DxFlux_org)

    end function IsCreate_BFDiff

    !>  检查是否创建基本流流场导数.
    !!
    !!  @return 流场导数是否创建
    !!  @retval true 已创建
    !!  @retval false 未创建
    function IsCreate_BFDiff_3d(this) result(is)

        implicit none
        class(Baseflow3D_Diff_type), intent(in) :: this
        logical :: is

        is=allocated(this%DxFlux_org)

    end function IsCreate_BFDiff_3d

    !>  读取某点的基本流法向一阶导数通量.
    !!
    !!  @param[in] iloc 流向位置
    !!  @param[in] jloc 法向位置
    !!  @return 基本流在该点的法向一阶导数
    !!  @retval Partial_Eta 此点法向一阶导数
    type(bf_flux_org_ij_type) function Partial_Eta_2d(this, iloc, jloc)

        use mod_difference
        implicit none
        class(Baseflow2D_Diff_type), intent(in) :: this
        integer, intent(in) :: iloc, jloc

        Partial_Eta_2d=this%DyFlux_org(iloc, jloc)

    end function Partial_Eta_2d

    !>  读取某点的基本流法向二阶导数通量.
    !!
    !!  @param[in] iloc 流向位置
    !!  @param[in] jloc 法向位置
    !!  @return 基本流在该点的法向二阶导数
    !!  @retval Partial_Eta 此点法向二阶导数
    type(bf_flux_org_ij_type) function Partial_Eta2_2d(this, iloc, jloc)

        use mod_difference
        implicit none
        class(Baseflow2D_Diff_type), intent(in) :: this
        integer, intent(in) :: iloc, jloc

        Partial_Eta2_2d=this%DyyFlux_org(iloc, jloc)

    end function Partial_Eta2_2d

    !>  读取某点的基本流流向一阶导数通量.
    !!
    !!  @param[in] iloc 流向位置
    !!  @param[in] jloc 法向位置
    !!  @return 基本流在该点的流向一阶导数
    !!  @retval Partial_Xi 此点流向一阶导数
    type(bf_flux_org_ij_type) function Partial_Xi_2d(this, iloc, jloc)

        use mod_difference
        implicit none
        class(Baseflow2D_Diff_type), intent(in) :: this
        integer, intent(in) :: iloc, jloc

        Partial_Xi_2d=this%DxFlux_org(iloc, jloc)

    end function Partial_Xi_2d

    !>  读取某点的基本流法向一阶导数通量.
    !!
    !!  @param[in] iloc 流向位置
    !!  @param[in] jloc 法向位置
    !!  @param[in] kloc 展向位置
    !!  @return 基本流在该点的法向一阶导数
    !!  @retval Partial_Eta 此点法向一阶导数
    type(bf_flux_org_ij_type) function Partial_Eta_3d(this, iloc, jloc, kloc)

        use mod_difference
        implicit none
        class(Baseflow3D_Diff_type), intent(in) :: this
        integer, intent(in) :: iloc, jloc, kloc

        Partial_Eta_3d=this%DyFlux_org(iloc, jloc, kloc)

    end function Partial_Eta_3d

    !>  读取某点的基本流法向二阶导数通量.
    !!
    !!  @param[in] iloc 流向位置
    !!  @param[in] jloc 法向位置
    !!  @param[in] kloc 展向位置
    !!  @return 基本流在该点的法向二阶导数
    !!  @retval Partial_Eta 此点法向二阶导数
    type(bf_flux_org_ij_type) function Partial_Eta2_3d(this, iloc, jloc, kloc)

        use mod_difference
        implicit none
        class(Baseflow3D_Diff_type), intent(in) :: this
        integer, intent(in) :: iloc, jloc, kloc

        Partial_Eta2_3d=this%DyyFlux_org(iloc, jloc, kloc)

    end function Partial_Eta2_3d

    !>  读取某点的基本流流向一阶导数通量.
    !!
    !!  @param[in] iloc 流向位置
    !!  @param[in] jloc 法向位置
    !!  @param[in] kloc 法向位置
    !!  @return 基本流在该点的流向一阶导数
    !!  @retval Partial_Xi 此点流向一阶导数
    type(bf_flux_org_ij_type) function Partial_Xi_3d(this, iloc, jloc, kloc)

        use mod_difference
        implicit none
        class(Baseflow3D_Diff_type), intent(in) :: this
        integer, intent(in) :: iloc, jloc, kloc

        Partial_Xi_3d=this%DxFlux_org(iloc, jloc, kloc)

    end function Partial_Xi_3d

    !>  读取某点的基本流展向一阶导数通量.
    !!
    !!  @param[in] iloc 流向位置
    !!  @param[in] jloc 法向位置
    !!  @param[in] kloc 法向位置
    !!  @return 基本流在该点的展向一阶导数
    !!  @retval Partial_Xi 此点展向一阶导数
    type(bf_flux_org_ij_type) function Partial_Zeta_3d(this, iloc, jloc, kloc)

        use mod_difference
        implicit none
        class(Baseflow3D_Diff_type), intent(in) :: this
        integer, intent(in) :: iloc, jloc, kloc

        Partial_Zeta_3d=this%DzFlux_org(iloc, jloc, kloc)

    end function Partial_Zeta_3d

    !>  读取某点的基本流展向二阶导数通量.
    !!
    !!  @param[in] iloc 流向位置
    !!  @param[in] jloc 法向位置
    !!  @param[in] kloc 展向位置
    !!  @return 基本流在该点的展向二阶导数
    !!  @retval Partial_Eta 此点展向二阶导数
    type(bf_flux_org_ij_type) function Partial_Zeta2_3d(this, iloc, jloc, kloc)

        use mod_difference
        implicit none
        class(Baseflow3D_Diff_type), intent(in) :: this
        integer, intent(in) :: iloc, jloc, kloc

        Partial_Zeta2_3d=this%DzzFlux_org(iloc, jloc, kloc)

    end function Partial_Zeta2_3d

    !>  读取某点的基本流法向展向二阶交叉导数通量.
    !!
    !!  @param[in] iloc 流向位置
    !!  @param[in] jloc 法向位置
    !!  @param[in] kloc 展向位置
    !!  @return 基本流在该点的法向展向二阶交叉导数
    !!  @retval Partial_Eta 此点法向展向二阶交叉导数
    type(bf_flux_org_ij_type) function Partial_EtaZeta_3d(this, iloc, jloc, kloc)

        use mod_difference
        implicit none
        class(Baseflow3D_Diff_type), intent(in) :: this
        integer, intent(in) :: iloc, jloc, kloc

        Partial_EtaZeta_3d=this%DyzFlux_org(iloc, jloc, kloc)

    end function Partial_EtaZeta_3d

    !>  设置基本流流场通量.
    !!
    !!  @param[in] rho 密度
    !!  @param[in] Vel 速度矢量
    !!  @param[in] T   温度
    !!  @return 全部2D流场信息
    subroutine set_baseflow2d(this, rho, Vel, T)

        implicit none
        class(Baseflow2D_type) , intent(inout)            :: this
        real(R_P)              , dimension(:, :), intent(in) :: rho,  T
        type(vector_type), dimension(:, :), intent(in) :: Vel
        integer :: i, j

        call this%Flux_org%Set(rho, Vel, T)

        do j=1, this%jn
            do i=1, this%in
                call this%Flux_org(i, j)%SetIJ(i, j)
            end do
        end do

    end subroutine set_baseflow2d

    !>  设置基本流流场通量.
    !!
    !!  @param[in] rho 密度
    !!  @param[in] Vel 速度矢量
    !!  @param[in] T   温度
    !!  @return 全部3D流场信息
    subroutine set_baseflow3d(this, rho, Vel, T)

        implicit none
        class(Baseflow3D_type) , intent(inout)            :: this
        real(R_P)              , dimension(:, :, :), intent(in) :: rho,  T
        type(vector_type), dimension(:, :, :), intent(in) :: Vel
        integer :: i, j, k

        call this%Flux_org%Set(rho, Vel, T)

        do k=1, this%kn
            do j=1, this%jn
                do i=1, this%in
                    call this%Flux_org(i, j, k)%SetIJ(i, j)
                end do
            end do
        enddo

    end subroutine set_baseflow3d

    !>  设置基本流流场通量.
    !!
    !!  @param[in] rho 密度
    !!  @param[in] U   流向速度
    !!  @param[in] V   法向速度
    !!  @param[in] W   展向速度
    !!  @param[in] T   温度
    !!  @return 全部2D流场信息
    subroutine set_baseflow2d_dble(this, rho, U, V, W, T)

        implicit none
        class(Baseflow2D_type) , intent(inout)            :: this
        real(R_P)              , dimension(:, :), intent(in) :: rho,  T
        real(R_P), dimension(:, :), intent(in) :: U, V, W
        integer :: i, j

        do j=1, this%jn
            do i=1, this%in
                call this%Flux_org(i, j)%Set(rho(i, j), Vector(U(i, j), V(i, j), W(i, j)), T(i, j))
                call this%Flux_org(i, j)%SetIJ(i, j)
            end do
        end do

    end subroutine set_baseflow2d_dble

    !>  设置基本流流场通量.
    !!
    !!  @param[in] rho 密度
    !!  @param[in] U   流向速度
    !!  @param[in] V   法向速度
    !!  @param[in] W   展向速度
    !!  @param[in] T   温度
    !!  @return 全部3D流场信息
    subroutine set_baseflow3d_dble(this, rho, U, V, W, T)

        implicit none
        class(Baseflow3D_type) , intent(inout)            :: this
        real(R_P), dimension(:, :, :), intent(in) :: rho,  T
        real(R_P), dimension(:, :, :), intent(in) :: U, V, W
        integer :: i, j, k

        do k=1, this%kn
            do j=1, this%jn
                do i=1, this%in
                    call this%Flux_org(i, j, k)%Set(rho(i, j, k), Vector(U(i, j, k), V(i, j, k), W(i, j, k)), T(i, j, k))
                    call this%Flux_org(i, j, k)%SetIJ(i, j)
                end do
            end do
        enddo

    end subroutine set_baseflow3d_dble

    !>  二维基本流类的析构函数.
    !!

    Subroutine finalize_baseflow2d(this)

        Implicit None
        class(Baseflow2D_type), intent(Inout) :: this

        this%in=0; this%jn=0
        deallocate(this%Flux_org)

    End Subroutine finalize_baseflow2d

    !>  三维基本流类的析构函数.
    !!
    Subroutine finalize_baseflow3d(this)

        Implicit None
        class(Baseflow3D_type), intent(Inout) :: this

        this%in=0; this%jn=0; this%kn=0
        deallocate(this%Flux_org)

    End Subroutine finalize_baseflow3d

    !>  二维基本流导数类的析构函数.
    !!
    Subroutine finalize_baseflow2d_diff(this)

        Implicit None
        class(Baseflow2D_Diff_type), intent(Inout) :: this

        deallocate(this%DxFlux_org)
        deallocate(this%DyFlux_org)
        deallocate(this%DyyFlux_org)

    End Subroutine finalize_baseflow2d_diff

    !>  二维基本流导数类的析构函数.
    !!
    Subroutine finalize_baseflow3d_diff(this)

        Implicit None
        class(Baseflow3D_Diff_type), intent(Inout) :: this

        deallocate(this%DxFlux_org)
        deallocate(this%DyFlux_org)
        deallocate(this%DzFlux_org)
        deallocate(this%DyyFlux_org)
        deallocate(this%DyzFlux_org)
        deallocate(this%DzzFlux_org)

    End Subroutine finalize_baseflow3d_diff


end module mod_baseflow
