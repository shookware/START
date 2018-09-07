!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_dis.f90
!> @file
!> @breif 扰动模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: dis
!> @breif 扰动模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-03 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-03
!------------------------------------------------------------------------------`
module mod_dis

    use mod_lpse_dis_normal
    use mod_dis_flux
    use mod_vector_cmplx
    use mod_grid
    use mod_dis_wavenum
    use penf, only: R_P

    implicit none

    private

    !> 扰动类型
    type, public :: dis_type ! linear PSE dis

        private
        integer, private :: in !< 流向点数
        integer, private :: jn !< 法向点数
        type(lpse_dis_normal_type), allocatable, private :: dis_normal(:) !<全场扰动,以流向排布的法向形函数分布体现

        Contains

          Procedure :: Create !< 创建类型,分配内存
          procedure :: Set !< 设置某站位处的形函数分布
          procedure :: finalize => finalize_dis !< 析构函数
          generic :: Print => print_dis_location, print_dis_whole !< 输出扰动场到文件
          generic :: Get => Get_disNorm, Get_wavenum !<获得扰动信息

          procedure, private :: print_dis_location !< 输出某站位处的扰动
          procedure, private :: print_dis_whole !< 输出全场扰动
          procedure, private :: Get_disNorm !< 获得某流向站位处扰动信息
          procedure, private :: Get_wavenum !< 获得某流向站位处的色散关系

    end type dis_type

    contains

    !> 输出全场扰动到文件
    !! @param[in] fn_surf 输出文件文件名前缀
    subroutine print_dis_whole(this, fn_surf)

        use mod_print
        use stringifor
        implicit none
        class(dis_type), intent(in) :: this
        type(string), intent(in) :: fn_surf
        type(lpse_dis_normal_type) :: DisNorm
        type(dis_flux_ij_type) :: FluxPoint
        complex(R_P), dimension(:, :, :), allocatable :: DisShape
        integer :: i, j, iloc
        type(string) :: fn

        allocate(DisShape(this%in, this%jn, 5))
        do iloc=1, this%in
          call this%Get(iloc, DisNorm)
          do j=1, this%jn
            call DisNorm%Get(j, FluxPoint)
            call FluxPoint%get(DisShape(iloc, j, 1), DisShape(iloc, j, 2), DisShape(iloc, j, 3), &
                               DisShape(iloc, j, 4), DisShape(iloc, j, 5))
          enddo
        enddo

        fn=fn_surf//'_DisNorm.dat'

        call PrintComplexArray2D(this%in, this%jn, 5, DisShape, fn)

        deallocate(DisShape)

    end subroutine print_dis_whole

    !>输出某流向站位处的扰动到文件
    !! @param[in] iloc 流向站位序号
    !! @param[in] grid 计算网格
    subroutine print_dis_location(this, iloc, grid)

        implicit none
        class(dis_type), intent(in) :: this
        type(grid2d_type), pointer, intent(in) :: Grid
        integer, intent(in) :: iloc
        type(lpse_dis_normal_type) :: DisNorm

        type(dis_flux_ij_type) :: FluxPoint
        type(vector_cmplx_type) :: Vel
        complex(R_P), dimension(this%jn) :: Rho, U, V, W, T
        integer :: j
        real(R_P), dimension(this%jn) :: yy
        type(dis_wavenum_lpse_type) :: wavenum
        complex(R_P) :: alf
        character(10) :: loc_char

        call this%Get(iloc, DisNorm)
        wavenum=DisNorm%GetWaveNum()
        alf=wavenum%getAlpha()
        do j=1, this%jn
          call DisNorm%Get(j, FluxPoint)
          call FluxPoint%get(rho(j), Vel, T(j))
          call Vel%Get(U(j), V(j), W(j))
        enddo
        call grid%get_iy(iloc, yy)
        write(loc_char, '(I4.4)')iloc
        open(999, file='dis_locationx'//trim(loc_char)//'_DisNorm.plt', form='formatted')
        write(999, *)"variables='y', 'rho', 'u', 'v', 'w', 'T'"
        do j=1, this%jn
          write(999, *)yy(j), abs(rho(j)), abs(u(j)), abs(v(j)), abs(w(j)), abs(T(j))
        end do
        close(999)

        open(999, file='dis_locationx'//trim(loc_char)//'_DisNorm_imag.plt', form='formatted')
        write(999, *)"variables='y', 'rho', 'u', 'v', 'w', 'T'"
        do j=1, this%jn
          write(999, *)yy(j), aimag(rho(j)), aimag(u(j)), aimag(v(j)), aimag(w(j)), aimag(T(j))
        end do
        close(999)

        open(999, file='dis_locationx'//trim(loc_char)//'_DisNorm_real.plt', form='formatted')
        write(999, *)"variables='y', 'rho', 'u', 'v', 'w', 'T'"
        do j=1, this%jn
          write(999, *)yy(j), real(rho(j)), real(u(j)), real(v(j)), real(w(j)), real(T(j))
        end do
        close(999)

        open(997, file='dis_locationx'//trim(loc_char)//'_alf.dat', form='formatted')
        write(997, *)"variables= 'alfr', 'alfi'"

          write(997, *)real(alf), aimag(alf)

        close(997)

    end subroutine print_dis_location

    !> 设置某流向站位处的扰动信息
    !! @param[in] iloc 流向站位序号
    !! @param[in] DisNorm 该流向站位处的扰动信息
    subroutine Set(this, iloc, DisNorm)

        implicit none
        type(lpse_dis_normal_type), intent(in) :: DisNorm
        integer, intent(in) :: iloc
        class(dis_type), intent(inout) :: this

!        call DisNorm%
        this%dis_normal(iloc)=DisNorm

    end subroutine Set

    !> 获得某流向站位处的扰动信息
    !! @param[in] iloc 流向站位序号
    !! @param[out] DisNorm 扰动信息
    subroutine Get_disNorm(this, iloc, DisNorm)

        implicit none
        type(lpse_dis_normal_type), intent(inout) :: DisNorm
        class(dis_type), intent(in) :: this
        integer, intent(in) :: iloc

        call DisNorm%Create(this%jn)
        DisNorm=this%dis_normal(iloc)

    end subroutine Get_disNorm

    !> 获得某流向站位处的扰动色散信息
    !! @param[in] iloc 流向站位序号
    !! @param[out] wavenum 扰动色散信息
    subroutine Get_wavenum(this, iloc, wavenum)

        implicit none
        type(dis_wavenum_lpse_type), intent(inout) :: wavenum
        class(dis_type), intent(in) :: this
        integer, intent(in) :: iloc

        wavenum=this%dis_normal(iloc)%GetWaveNum()

    end subroutine Get_wavenum

    !> 创建类型,分配内存
    !! @param[in]  grid 计算网格
    subroutine create(this, grid)

        use mod_grid
        implicit none
        type(grid2d_type), intent(in) :: grid
        class(dis_type), intent(inout) :: this
        integer :: in, jn
        integer :: i

        call grid%getsize(in, jn)
        this%in=in; this%jn=jn
        if(.not. allocated(this%dis_normal)) &
        &   allocate(this%dis_normal(in))

        do i=1, in
            call this%dis_normal(i)%Create(jn)
            !call this%dis_normal(i)%CreateDiff(jn)
            call this%dis_normal(i)%SetILoc(i)
            !call this%dis_normal(i)%SetDiffIloc(i)
        end do

    end subroutine create

    !> 析构函数
    Subroutine finalize_dis(this)

        Implicit None
        class(dis_type), intent(Inout) :: this
        integer :: i

        do i=1, this%in
          call this%dis_normal(i)%finalize()
        enddo

        deallocate(this%dis_normal)
        this%in=0; this%jn=0

    End Subroutine finalize_dis

end module mod_dis
