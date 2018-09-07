!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_dis_shape.f90
!> @file
!> @breif 扰动形函数模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: dis_shape
!> @breif 扰动形函数模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-03 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-03
!------------------------------------------------------------------------------`
module mod_dis_shape

    use mod_dis_flux
    use penf, only: R_P

    implicit none

    private

    !> 扰动形函数沿法向分布类型
    type, public :: dis_shape_type !!linear used

        private
        integer, public :: jn=0  !< 法向点数
        integer, private :: iloc  !< 流向站位
        type(dis_flux_ij_type), allocatable, private :: flux(:) !< 扰动形函数法向分布

        Contains

          Procedure :: Create   !< 根据法向点数创建法向形函数法分布,分配内存
          procedure :: SetILoc => set_iloc !< 设置流向站位
          procedure :: GetDisShape !< 获得扰动形函数分布
          procedure :: SetDisShape !< 设置扰动形函数分布
          procedure :: GetFlux !< 获得整个法向的扰动形函数通量
          procedure :: getflux_array !< 获得整个法向的扰动形函数通量(数组)
          procedure :: GetFlux_array_5 !<获得整个法向扰动形函数通量(数组，变量分开)
          procedure :: GetILoc => get_i_index !<获得流向站位
          procedure :: getJn => get_jn !<获得法向点数
          procedure :: SetFluxPoint !< 设置法向某点处的扰动形函数通量值
          procedure :: IsCreate !< 判断是否已经分配内存,完成初始化
          procedure :: ReadDisShape=>Read !< 读取扰动
          procedure :: WriteDisShape=>Write !< 写入扰动
          procedure :: TAmp !<温度扰动幅值
          procedure :: UAmp !<速度扰动幅值

          generic :: Get => GetFluxPoint_flux, &
                          & GetFluxPoint_cmplx, & !< 获得法向某点处的扰动形函数通量值
                          & GetFlux, &
                          & getflux_array, &
                          & GetFlux_array_5
          generic :: Set => set_array, set_flux  !< 设置形函数法向分布的值

          generic :: operator(*) => mul_scalar_disshape  !< 标量乘
          generic :: operator(+) => add_disshape  !< 形函数相加

          procedure, pass(this), private :: mul_scalar_disshape !< 标量乘
          procedure, pass(this), private :: add_disshape !< 形函数相加
          procedure, pass(This), private :: GetFluxPoint_flux !<获得法向某点处的形函数通量值
          procedure, pass(this), private :: GetFluxPoint_cmplx !< 获得法向某点处的形函数的通量值
          procedure, pass(This), private :: set_flux !< 设置形函数法向分布的值(通量)
          procedure, pass(This), private :: set_array !< 设置形函数法向分布的值(矩阵)

          procedure :: finalize => finalize_dis_shape !< 析构函数

    end type dis_shape_type

    !!> 沿法向分布的扰动形函数通量的导数值类型
    !type, public :: dis_shape_diff_type     !!linear used
    !
    !    private
    !    integer, private :: jn=0 !< 法向点数
    !    integer, private :: iloc !< 流向站位
    !    type(dis_flux_ij_type), allocatable, private :: Dyflux(:) !< 形函数沿法向的一阶导数法分布
    !    type(dis_flux_ij_type), allocatable, private :: Dyyflux(:) !< 形函数沿法向的二阶导数法分布
    !
    !    Contains
    !
    !      Procedure :: Create => Create_Diff !< 根据法向点数创建法向形函数导数法分布,分配内存
    !      procedure :: Set => Set_diff      !< 设置形函数导数的法分布
    !      procedure :: SetILoc => SetILoc_Diff !< 设置流向站位
    !      procedure :: GetPartialEta !< 获得形函数某法向位置的法向的一阶导数
    !      procedure :: GetPartialEta2 !< 获得形函数某法向位置的法向的二阶导数
    !      procedure :: finalize => finalize_dis_shape_diff !< 析构函数
    !
    !end type dis_shape_diff_type

    contains

    !> 温度扰动幅值
    real function TAmp(this)

      implicit none
      class(dis_shape_type), intent(in) :: this
      integer :: j
      complex(R_P) :: rho, u, v, w, t

      TAmp=0.0d0
      do j=1, this%jn
        call this%flux(j)%Get(rho, u, v, w, T)
        if(abs(T)>=TAmp) TAmp=abs(T)
      enddo

    end function TAmp

    !> 速度扰动幅值
    real function UAmp(this)

      implicit none
      class(dis_shape_type), intent(in) :: this
      integer :: j
      complex(R_P) :: rho, u, v, w, t
      complex(R_P) :: Usqrt

      UAmp=0.0d0
      do j=1, this%jn
        call this%flux(j)%Get(rho, u, v, w, T)
        Usqrt=U*conjg(U)+v*conjg(v)+w*conjg(w)
        if(sqrt(real(Usqrt))>=UAmp) &
        & UAmp=sqrt(real(Usqrt))
      enddo

    end function UAmp

    !> 读取文件
    subroutine Read(this, unit)

      implicit none
      class(dis_shape_type), intent(inout) :: this
      integer, intent(in) :: unit
      complex(R_P), allocatable, save :: flux(:)

      read(unit)this%jn, this%iloc
      if(.not. allocated(flux)) allocate(flux(5*this%jn))
      read(unit)flux
      if(.not. this%IsCreate()) call this%Create(this%jn)
      call this%set_array(flux, this%iloc)

    end subroutine Read

    !> 写入文件
    subroutine Write(this, unit)

      implicit none
      class(dis_shape_type), intent(in) :: this
      integer, intent(in) :: unit
      complex(R_P) :: flux(5*this%jn)

      write(unit)this%jn, this%iloc
      call this%getflux_array(flux)
      write(unit)flux

    end subroutine write
    !-----------------------------------------------------------------------
    !> 判断变量内存是否已经分配,完成初始化过程
    !-----------------------------------------------------------------------
    logical function  IsCreate(this)

        implicit none
        class(dis_shape_type), intent(in) :: this

        IsCreate=allocated(this%flux)

   end function IsCreate
    !> 类型相加
    function add_disshape(this, obj) result(add)

        implicit none
        class(dis_shape_type), intent(in) :: this
        type(dis_shape_type), intent(in) :: obj
        type(dis_shape_type) :: add
        integer :: i

        add%jn=this%jn
        call add%Create(this%jn)

        add%iloc=this%iloc
        do i=1, this%jn
          add%flux(i)=this%flux(i)+obj%flux(i)
        end do

    end function add_disshape

    !> 类型标量乘
    function mul_scalar_disshape(scalar, this) result(mul)

        implicit none
        real(R_P), intent(in) :: scalar
        class(dis_shape_type), intent(in), target :: this
        type(dis_shape_type), target :: mul
        integer :: i

        call mul%create(this%jn)
        mul%jn=this%jn
        mul%iloc=this%iloc

        do i=1, this%jn
          mul%flux(i)=scalar*this%flux(i)
        enddo

    end function mul_scalar_disshape


    !> 获得形函数法分布
    !! @return 形函数法分布类型
    function GetDisShape(this) result(DisShape)

        implicit none
        class(dis_shape_type), intent(in) :: this
        type(dis_shape_type) :: Disshape

        DisShape%jn   = this%jn
        call Disshape%Create(this%jn)
        DisShape%flux = this%flux
        DisShape%iloc = this%iloc

    end function GetDisShape

    !> 设置扰动形函数法分布
    !! @param[in] Disshape 扰动形函数的法分布
    subroutine SetDisShape(this, Disshape)

        implicit none
        class(dis_shape_type), intent(inout) :: this
        type(dis_shape_type), intent(in) :: Disshape

        this%jn   = DisShape%jn
        this%iloc = Disshape%iloc
        call this%Create(DisShape%jn)
        this%flux = DisShape%flux

    end subroutine SetDisShape

    !< 根据法向点数创建法向形函数法分布,分配内存
    !! @param[in] jn 法向点数
    elemental subroutine create(this, jn)

        use mod_grid
        implicit none
        integer , intent(in) :: jn
        class(dis_shape_type), intent(inout) :: this
        integer :: j

        this%jn=jn
        this%iloc=0
        if(.not. allocated(this%flux)) &
        &   allocate(this%flux(this%jn))
        do j=1, this%jn
            this%flux(j)=DIS_FLUXIJ_NULL
        enddo

    end subroutine create

    !!< 根据法向点数创建法向形函数导数法分布,分配内存
    !!! @param[in] jn 法向点数
    !elemental subroutine create_diff(this, jn)
    !
    !    use mod_grid
    !    implicit none
    !    integer , intent(in) :: jn
    !    class(dis_shape_diff_type), intent(inout) :: this
    !
    !    this%jn=jn
    !    this%iloc=0
    !    if(.not. allocated(this%Dyflux))then
    !      allocate(this%Dyflux(jn))
    !      allocate(this%Dyyflux(jn))
    !    endif
    !
    !end subroutine create_diff

    !> 设置流向站位
    !! @param[in] iloc 流向站位序号
    subroutine set_iloc(this, iloc)

        implicit none
        class(dis_shape_type), intent(inout) :: this
        integer, intent(in) :: iloc

        this%iloc=iloc

    end subroutine set_iloc

    !!> 设置流向站位
    !!! @param[in] iloc 流向站位序号
    !subroutine setiloc_diff(this, iloc)
    !
    !    implicit none
    !    class(dis_shape_diff_type), intent(inout) :: this
    !    integer, intent(in) :: iloc
    !
    !    this%iloc=iloc
    !
    !end subroutine setiloc_diff

    !> 设置形函数法向分布的值
    !! @param[in] flux 形函数法分布
    !! @param[in] iloc 流向站位序号
    subroutine set_flux(this, flux, iloc)

        implicit none
        class(dis_shape_type), intent(inout) :: this
        type(dis_flux_ij_type), intent(in) :: flux(:)
        integer, intent(in) :: iloc
        integer :: j

        if(.not. this%IsCreate()) call this%Create(size(flux))
        this%iloc=iloc
        do j=1, this%jn
            this%flux(j)=flux(j)
            call this%Flux(j)%setIJ(iloc, j)
        enddo

    end subroutine set_flux

    !> 设置形函数法向分布的值(矩阵）
    !! @param[in] flux 形函数法分布
    subroutine set_array(this, flux, iloc)

        implicit none
        class(dis_shape_type), intent(inout) :: this
        complex(R_P), intent(in) :: flux(:)
        integer, intent(in) :: iloc
        integer :: j

        if(.not. this%IsCreate()) call this%Create(size(flux)/5)
        this%iloc=iloc
        do j=1, this%jn
            call this%flux(j)%set(flux((j-1)*5+1:(j-1)*5+5))
            call this%Flux(j)%setIJ(iloc, j)
        enddo

    end subroutine set_array

    !< 获得整个法向的扰动形函数通量
    !! @param[in] flux 整个法向的扰动形函数通量
    subroutine getflux(this, flux)

        implicit none
        class(dis_shape_type), intent(in) :: this
        type(dis_flux_ij_type), intent(inout) :: flux(:)
        integer :: j

        do j=1, this%jn
            flux(j)=this%flux(j)
        enddo

    end subroutine getflux

    !< 获得整个法向的扰动形函数通量
    !! @param[inout] flux 整个法向的扰动形函数通量(数组形式)
    subroutine getflux_array(this, flux)

        implicit none
        class(dis_shape_type), intent(in) :: this
        complex(R_P), intent(inout) :: flux(this%jn*5)
        integer :: j

#IFDEF DEBUG
        ! print*, 'getflux_array', this%jn
#ENDIF

        do j=0, this%jn-1
            call this%flux(j+1)%get(flux(j*5+1), flux(j*5+2), flux(j*5+3), flux(j*5+4), flux(j*5+5))
        enddo

    end subroutine getflux_array

    !< 获得整个法向的扰动形函数通量
    !! @param[inout] flux 整个法向的扰动形函数通量(数组形式，但是数组分开)
    subroutine getflux_array_5(this, flux)

        implicit none
        class(dis_shape_type), intent(in) :: this
        complex(R_P), intent(inout) :: flux(5, this%jn)
        integer :: j

#IFDEF DEBUG
        ! print*, 'getflux_array5'
        ! print*, this%jn
#ENDIF
        do j=1, this%jn
            call this%flux(j)%get(flux(1, j), flux(2, j), flux(3, j), flux(4, j), flux(5, j))
        enddo


    end subroutine getflux_array_5

    !> 获得流向站位
    !! @return 流向站位
    integer function get_i_index(this)

        implicit none
        class(dis_shape_type), intent(in) :: this

        get_i_index=this%iloc

    end function get_i_index

    !!< 设置形函数导数的法分布
    !!! @param[in] Dis 扰动形函数在某流向站位上的分布
    !!! @param[in] Diff 扰动差分算子的基函数
    !subroutine Set_diff(this, Dis, diff)
    !
    !    use mod_difference
    !    implicit none
    !    class(dis_shape_diff_type), intent(inout) :: this
    !    class(dis_shape_type), intent(in) :: Dis
    !    type(difference_2D_type), intent(in) :: diff
    !    integer :: iloc, jn
    !
    !    iloc=Dis%iloc; jn=Dis%jn
    !    call diff%PartialEta(iloc, iloc, Dis%Flux, this%DyFlux)
    !    call diff%PartialEta2(iloc, iloc, Dis%Flux, this%DyyFlux)
    !
    !end subroutine Set_diff

    !!< 获得形函数沿法向的一阶导数
    !!! @param[in] jloc 法向位置序号
    !!! @return 在jloc点处的形函数法向一阶导数
    !function GetPartialEta(this, jloc) result(diff)
    !
    !    implicit none
    !    class(dis_shape_diff_type), intent(in) :: this
    !    integer, intent(in) :: jloc
    !    type(dis_flux_ij_type) :: Diff
    !
    !    Diff=this%Dyflux(jloc)
    !
    !end function GetPartialEta

    !!< 获得形函数沿法向的二阶导数
    !!! @param[in] jloc 法向位置序号
    !!! @return 在jloc点处的形函数法向二阶导数
    !function GetPartialEta2(this, jloc) result(diff)
    !
    !    implicit none
    !    class(dis_shape_diff_type), intent(in) :: this
    !    integer, intent(in) :: jloc
    !    type(dis_flux_ij_type) :: Diff
    !
    !    Diff=this%Dyyflux(jloc)
    !
    !end function GetPartialEta2

    !> 获得法向点数
    !! @return 法向点数
    function get_jn(this) result(jn)

        class(dis_shape_type), intent(in) :: this
        integer :: jn

        jn=this%jn

    end function get_jn

    !<获得法向某点处的形函数通量值
    !! @param[out] j 法向站位
    !! @param[out] FluxPoint 在j点处的扰动形函数通量值
    subroutine GetFluxPoint_flux(this, j, FluxPoint)

        implicit none
        class(dis_shape_type), intent(in) :: this
        integer, intent(in) :: j
        type(dis_flux_ij_type), intent(inout) :: fluxPoint

        fluxPoint=this%Flux(j)

    end subroutine GetFluxPoint_flux

    !<获得法向某点处的形函数通量值(复数数组形式)
    !! @param[out] j 法向站位
    !! @param[out] FluxPoint 在j点处的扰动形函数通量值
    subroutine GetFluxPoint_cmplx(this, j, Flux)

        implicit none
        class(dis_shape_type), intent(in) :: this
        integer, intent(in) :: j
        complex(R_P), intent(out) :: flux(5)
        type(dis_flux_ij_type) :: FluxPoint

        fluxpoint=this%Flux(j)
        flux=(0.0d0, 0.0d0)
        call Fluxpoint%get(flux(1), flux(2), flux(3), flux(4), flux(5))

    end subroutine GetFluxPoint_cmplx

    !< 设置法向某点处的扰动形函数通量值
    !! @param[in] j 法向站位
    !! @param[in] FluxPoint 在j点处的扰动形函数通量值
    subroutine SetFluxPoint(this, j, FluxPoint)

        implicit none
        class(dis_shape_type), intent(inout) :: this
        integer, intent(in) :: j
        type(dis_flux_ij_type), intent(in) :: fluxPoint

        this%Flux(j)=fluxPoint

    end subroutine SetFluxPoint

    !> 析构函数
    Subroutine finalize_dis_shape(this)

        Implicit None
        class(dis_shape_type), intent(Inout) :: this

        this%jn=0
        deallocate(this%flux)

    End Subroutine finalize_dis_shape

    !!> 析构函数
    !Subroutine finalize_dis_shape_diff(this)
    !
    !    Implicit None
    !    class(dis_shape_diff_type), intent(Inout) :: this
    !
    !    this%jn=0
    !    deallocate(this%Dyflux)
    !    deallocate(this%Dyyflux)
    !
    !End Subroutine finalize_dis_shape_diff

end module mod_dis_shape
