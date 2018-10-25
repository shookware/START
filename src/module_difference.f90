!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_difference.f90
!> @file
!> @breif 差分算子模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: difference
!> @breif 差分算子模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-02 - Initial Version
!  TODO_2017-08-02
!> @author
!> Liu Jianxin
!> \date 2017-20-02
!------------------------------------------------------------------------------
module mod_difference

   use mod_grid
   use mod_baseflow_org
   use mod_dis_flux
   use penf, only: R_P

   implicit none

   private

   type, public :: Norm_Diff_Coef_type

      integer :: jn
      real(R_P), allocatable :: Coef_di(:)
      real(R_P), allocatable :: Coef_dj(:, :)
      real(R_P), allocatable :: Coef_djj(:, :)

    contains

      generic :: dx=> dx_dis_flux
      generic :: dy=> dy_dis_flux, dy_complex
      generic :: dyy=> dyy_dis_flux, dyy_complex
      procedure, private :: dx_dis_flux
      procedure, private :: dy_dis_flux
      procedure, private :: dyy_dis_flux
      procedure, private :: dy_complex
      procedure, private :: dyy_complex

   end type Norm_Diff_Coef_type

   type, public :: difference_object_type

       private

      integer, private :: in !< 流向点数
      integer, private :: jn !< 法向点数
      integer, private :: XIorder !<流向精度
      integer, private :: ETAorder !<法向精度
      integer, private :: difftype !<差分类型:基本流还是扰动
      real(R_P), dimension(:, :, :, :), allocatable :: Coef_di !<流向一阶差分基函数
      real(R_P), dimension(:, :, :, :), allocatable :: Coef_dj !<法向一阶差分基函数
      real(R_P), dimension(:, :, :, :), allocatable :: Coef_djj !<法向儿阶差分基函数

       Contains
      procedure :: GetXIOrder => get_xiorder      !<获得流向差分精度
      procedure :: GetETAOrder => get_etaorder    !<获得法向差分精度
      procedure :: GetXICoef => get_Coef_XI       !<获得流向插值基函数
      procedure :: GetETACoef => Get_Coef_ETA     !<获得法向插值基函数
      procedure :: IsCreate                       !<判断是否已经分配空间

      generic, public :: PartialXI => Partial_XI_BFFluxIJ, &
                                    & Partial_XI_BFFluxIJ_iloc !<流向求一阶导数
      generic, public :: PartialETA => Partial_Eta_BFFluxIJ, &
                                    &  Partial_Eta_DisFluxIJ, &
                                    & Partial_Eta_BFFluxIJ_iloc !<法向求一阶导数
      generic, public :: PartialETA2 => Partial_Eta2_BFFluxIJ, &
                                    &   Partial_Eta2_DisFluxIJ, &
                                    &   Partial_Eta2_BFFluxIJ_iloc !<法向求二阶导数

      procedure, public :: Partial_Eta_DisFluxIJ
      procedure, public :: Partial_Eta2_DisFluxIJ
      procedure, public :: Partial_Eta_BFFluxIJ_iloc
      procedure, public :: Partial_Eta2_BFFluxIJ_iloc
      procedure, public :: Partial_Eta_BFFluxIJ
      procedure, public :: Partial_Eta2_BFFluxIJ
      procedure, public :: Partial_XI_BFFluxIJ
      procedure, public :: Partial_XI_BFFluxIJ_iloc

   end type difference_object_type


   !> 二维差分算子类
   type, extends(difference_object_type), public :: difference_2D_type

      private
      type(grid2d_type), Pointer, private :: grid=>null() !<待差分二维网格

   Contains

      procedure :: GetEta                         !<获得法向坐标
      procedure :: GetLocalNormDiffCoef           !<获得当地法向插值基函数
      Procedure :: InitDiff => initial_difference_2d !<初始化差分算子
      procedure :: Finalize => finalize_difference_2D !<析构函数

   end type difference_2D_type

   type, extends(difference_object_type), public :: difference_3d_type

      private
      integer, private :: kn !<展向点数
      type(grid3d_type), Pointer, private :: grid !<待差分三维网格
      real(R_P), private :: beta0 !<谱方法中的基本波数

      contains

      generic :: PartialZeta => Partial_ZETA_BFFluxIJ, Partial_ZETA_DISFluxIJ
      generic :: PartialZeta2 => Partial_ZETA2_BFFluxIJ, Partial_ZETA2_DISFluxIJ
!      generic :: PartialEtaZeta2 => Partial_ETAZETA2_BFFluxIJ, Partial_ETAZETA2_DISFluxIJ

      procedure, private :: Partial_ZETA_BFFluxIJ
      procedure, private :: Partial_ZETA_DISFluxIJ
      procedure, private :: Partial_ZETA2_BFFluxIJ
      procedure, private :: Partial_ZETA2_DISFluxIJ

   end type  difference_3d_type


   real(R_P), allocatable :: abx(:, :), bby(:, :), bbyy(:, :)

  contains

    function dx_dis_flux(this, disflux)

      class(Norm_Diff_Coef_type), intent(in) :: this
      type(dis_flux_ij_type), intent(in) :: disflux(:, :)
      type(dis_flux_ij_type), allocatable :: dx_dis_flux(:)
      integer :: dimsize
      integer :: i, j

      dimsize=size(this%Coef_di, dim=1)
      dx_dis_flux=DIS_FLUXIJ_NULL

      if(size(disflux, dim=2)/=this%jn) stop 'dx disflux j num is wrong!'

      allocate(dx_dis_flux(this%jn))

      do j=1, this%jn
        do i=1, dimsize
          dx_dis_flux(j)=dx_dis_flux(j)+this%Coef_di(i)*disflux(i, j)
        enddo
      enddo

    end function  dx_dis_flux

    function dy_dis_flux(this, disflux)
      class(Norm_Diff_Coef_type), intent(in) :: this
      type(dis_flux_ij_type), intent(in) :: disflux(:)
      type(dis_flux_ij_type), allocatable :: dy_dis_flux(:)
      integer :: dimsize

      integer :: j, l

      dimsize=size(this%Coef_dj, dim=1)

#ifdef DEBUG
      ! print*, size(disflux), this%jn
#endif

      if(size(disflux)/=this%jn) stop 'dy disflux j num is wrong!'
      allocate(dy_dis_flux(this%jn))

      dy_dis_flux=DIS_FLUXIJ_NULL
      do j=1, dimsize/2
        do l=1, dimsize
          dy_dis_flux(j)=dy_dis_flux(j)+this%Coef_dj(l, j)*disflux(l)
        enddo
      enddo
      do j=dimsize/2+1, this%jn-dimsize/2
        do l=1, dimsize
          dy_dis_flux(j)=dy_dis_flux(j)+this%Coef_dj(l, j)*disflux(j-dimsize/2-1+l)
        enddo
      enddo
      do j=this%jn-dimsize/2+1, this%jn
        do l=1, dimsize
          dy_dis_flux(j)=dy_dis_flux(j)+this%Coef_dj(l, j)*disflux(this%jn-dimsize+l)
        enddo
      enddo

    end function dy_dis_flux

    function dy_complex(this, f) result(dy_f)
      class(Norm_Diff_Coef_type), intent(in) :: this
      complex(R_P), intent(in) :: f(:)
      complex(R_P), allocatable :: dy_f(:)
      integer :: dimsize

      integer :: j, l

      dimsize=size(this%Coef_dj, dim=1)

#ifdef DEBUG
      ! print*, f
      ! print*, size(f), this%jn
      ! print*, dimsize/2

#endif
      if(size(f)/=this%jn) stop 'dy complex f j num is wrong!'
      allocate(dy_f(this%jn))

      dy_f=0.0d0
      do j=1, dimsize/2
        do l=1, dimsize
          dy_f(j)=dy_f(j)+this%Coef_dj(l, j)*f(l)
        enddo
      enddo
      do j=dimsize/2+1, this%jn-dimsize/2
        do l=1, dimsize
          dy_f(j)=dy_f(j)+this%Coef_dj(l, j)*f(j-dimsize/2-1+l)
        enddo
      enddo
      do j=this%jn-dimsize/2+1, this%jn
        do l=1, dimsize
          dy_f(j)=dy_f(j)+this%Coef_dj(l, j)*f(this%jn-dimsize+l)
        enddo
      enddo

    end function dy_complex

    function dyy_dis_flux(this, disflux)
      class(Norm_Diff_Coef_type), intent(in) :: this
      type(dis_flux_ij_type), intent(in) :: disflux(:)
      type(dis_flux_ij_type), allocatable :: dyy_dis_flux(:)
      integer :: dimsize

      integer :: j, l

      dimsize=size(this%Coef_dj, dim=1)

      if(size(disflux)/=this%jn) stop 'dyy disflux j num is wrong!'
      allocate(dyy_dis_flux(this%jn))

      dyy_dis_flux=DIS_FLUXIJ_NULL
      do j=1, dimsize/2
        do l=1, dimsize
          dyy_dis_flux(j)=dyy_dis_flux(j)+this%Coef_djj(l, j)*disflux(l)
        enddo
      enddo
      do j=dimsize/2+1, this%jn-dimsize/2
        do l=1, dimsize
          dyy_dis_flux(j)=dyy_dis_flux(j)+this%Coef_djj(l, j)*disflux(j-dimsize/2-1+l)
        enddo
      enddo
      do j=this%jn-dimsize/2+1, this%jn
        do l=1, dimsize
          dyy_dis_flux(j)=dyy_dis_flux(j)+this%Coef_djj(l, j)*disflux(this%jn-dimsize+l)
        enddo
      enddo


    end function dyy_dis_flux

    function dyy_complex(this, f) result(dyy_f)
      class(Norm_Diff_Coef_type), intent(in) :: this
      complex(R_P), intent(in) :: f(:)
      complex(R_P), allocatable :: dyy_f(:)
      integer :: dimsize

      integer :: j, l

      dimsize=size(this%Coef_dj, dim=1)

      if(size(f)/=this%jn) stop 'dyy complex f j num is wrong!'
      allocate(dyy_f(this%jn))

      dyy_f=0.0d0
      do j=1, dimsize/2
        do l=1, dimsize
          dyy_f(j)=dyy_f(j)+this%Coef_djj(l, j)*f(l)
        enddo
      enddo
      do j=dimsize/2+1, this%jn-dimsize/2
        do l=1, dimsize
          dyy_f(j)=dyy_f(j)+this%Coef_djj(l, j)*f(j-dimsize/2-1+l)
        enddo
      enddo
      do j=this%jn-dimsize/2+1, this%jn
        do l=1, dimsize
          dyy_f(j)=dyy_f(j)+this%Coef_djj(l, j)*f(this%jn-dimsize+l)
        enddo
      enddo

    end function dyy_complex

    subroutine GetLocalNormDiffCoef(this, iloc, NormDiffCoef)

      implicit none
      class(difference_2D_type), intent(in) :: this
      integer, intent(in) :: iloc
      type(Norm_Diff_Coef_type), intent(out) :: NormDiffCoef
      integer :: iSize, jSize
      integer :: iCount, jCount
      integer :: err

      NormDiffCoef%jn=this%jn
      iSize=size(this%Coef_di, dim=1)
      jSize=size(this%Coef_dj, dim=1)

      if(.not. allocated(NormDiffCoef%Coef_di)) then
        allocate(NormDiffCoef%Coef_di(iSize), stat=err)
        if ( err/= 0) print *, "NormDiffCoef%Coef_di: Allocation request denied"
        allocate(NormDiffCoef%Coef_dj(jSize, this%jn), NormDiffCoef%Coef_djj(jSize, this%jn), stat=err)
        if ( err/= 0) print *, "NormDiffCoef%Coef_di: Allocation request denied"
      endif

      NormDiffCoef%Coef_di=this%Coef_di(:, iloc, 1, 1)
      NormDiffCoef%Coef_dj=this%Coef_dj(:, iloc, :, 1)
      NormDiffCoef%Coef_djj=this%Coef_djj(:, iloc, :, 1)

    end subroutine GetLocalNormDiffCoef

    !> 判断是否已经分配空间.
    !! @retval TRUE 已经分配空间
    !! @retval FALSE 未分配空间
    function IsCreate(this) result(Flag)

        implicit none
        class(difference_object_type), intent(in) :: this
        logical :: Flag

        Flag=allocated(this%Coef_di)

    end function IsCreate

   !> 获得某站位法向坐标
   !! @param[in] iloc 流向站位
   !! @param[out] Eta 法向坐标
   subroutine GetEta(this, iloc, Eta)

       implicit none
       class(difference_2D_type), intent(in) :: this
       integer , intent(in) :: iloc
       real(R_P), intent(inout) :: Eta(this%jn)

       call this%Grid%Get_iy(iloc, Eta)

   end subroutine GetEta

   !> 初始化差分算子
   !! @param[in] grid 待差分二维网格
   !! @param[in] xi_order 流向差分精度
   !! @param[in] eta_order 法向差分精度
   !! @param[in] difftype 差分类型
   subroutine initial_difference_2d(this, grid, xi_order, eta_order, difftype)

      implicit none
      type(grid2d_type), Pointer, intent(in) :: grid
      integer, intent(in) :: xi_order, eta_order
      integer, intent(in) :: difftype
      class(difference_2D_type), intent(inout) :: this
      integer :: in, jn
      integer, parameter :: DIFFBF=1000, DIFFDIS=1002

      if (.not. Grid%HasBodyFixed()) stop "The Grid must be a body-fixed one! Use Trans Wrapper to get one!"
      call grid%GetSize(in, jn)
      this%in = in; this%jn = jn
      if((associated(this%grid)))  this%grid=>null()
      this%grid => grid
      if (eta_order .ne. 4) stop "Normal order must be 4 now!"
      !if (xi_order .ne. 1) stop 'streamwise order must be 1 now!'
      this%XIorder = xi_order
      this%ETAorder= eta_order
      if (difftype .ne. DIFFBF .and. difftype .ne. DIFFDIS) stop "Please input a correct difference type"
      this%difftype = difftype

      call initial_grid_coef(this)

      call compute_Coef_eta(this)
      call compute_Coef_xi (this)

      call finalize_grid_coef()

   end subroutine initial_difference_2d

   !>计算计算网格的导数 \f$\frac {\partial a} {\partial x},\frac {\partial b} {\partial y}\f$
   subroutine initial_grid_coef(this)

     implicit none
     type(difference_2D_type), intent(in) :: this
     integer :: in, jn
     integer :: i, j
     real(R_P) :: x(this%in), y(this%jn)

     in=this%In
     jn=this%Jn

     allocate(abx(in, jn)); allocate(bby(in, jn)); allocate(bbyy(in, jn))
     do j=1, jn
       call this%grid%Get_jx(j, x)
       call df4center(in, x, abx(:, j))
       abx(:, j)=1.0d0/abx(:, j)
     enddo

     do i=1, in
       call this%grid%Get_iy(i, y)
       call df4center(jn, y, bby(i, :))
       bby(i, :)=1.0d0/bby(i, :)
       call df4center(jn, bby(i, :), bbyy(i, :))
       bbyy(i, :)=bbyy(i, :)*bby(i, :)
     enddo

   end subroutine initial_grid_coef

   subroutine df4center(n, f, df)

     implicit none
     integer::N
    real(R_P) :: f(0:N-1), df(0:N-1)
   	real(R_P) :: dq(-2:2)
   	real(R_P) :: ry
   	integer :: n2, j

   !      hy=1.d0/dble(n2)
   	n2=N-1
   	ry=real(n2, R_P)
   !---------j=2,n-2
   	dq(-2)=ry/12.d0
   	dq(-1)=-ry*8.d0/12.d0
   	dq(0)=0.d0
   	dq(1)=ry*8.d0/12.d0
   	dq(2)=-ry/12.d0
   	do j=2,n2-2
         df(j)=f(j-2)*dq(-2)+f(j-1)*dq(-1)+f(j)*dq(0)+&
        &	   f(j+1)*dq(1)+f(j+2)*dq(2)
         end do
   !---------j=0
   	dq(0)=-ry*3.d0/2.d0
   	dq(1)=ry*2.d0
   	dq(2)=-ry/2.d0
         df(0)=f(0)*dq(0)+f(1)*dq(1)+f(2)*dq(2)
   !---------j=1
         dq(-1)=-ry/3.d0
   	dq(0)=-ry/2.d0
   	dq(1)=ry
   	dq(2)=-ry/6.d0
         df(1)=f(0)*dq(-1)+f(1)*dq(0)+f(2)*dq(1)+f(3)*dq(2)
   !---------j=n-1
         dq(-2)=ry/6.d0
   	dq(-1)=-ry
   	dq(0)=ry/2.d0
   	dq(1)=ry/3.d0
         df(n2-1)=f(n2-3)*dq(-2)+f(n2-2)*dq(-1)+f(n2-1)*dq(0)+&
        &          f(n2)*dq(1)
   !---------j=n
   	dq(-2)=ry/2.d0
   	dq(-1)=-ry*2.d0
   	dq(0)=ry*3.d0/2.d0
    df(n2)=f(n2-2)*dq(-2)+f(n2-1)*dq(-1)+f(n2)*dq(0)

   endsubroutine df4center

   !>释放计算网格导数空间
   subroutine finalize_grid_coef()

     implicit none

     deallocate(bby, bbyy, abx)


   end subroutine finalize_grid_coef


   !> 初始化差分算子
   !! @param[in] grid 待差分二维网格
   !! @param[in] xi_order 流向差分精度
   !! @param[in] eta_order 法向差分精度
   !! @param[in] zeta_order 展向差分精度  ! -1对应谱方法
   !! @param[in] difftype 差分类型
   subroutine initial_difference_3d(this, grid, xi_order, eta_order, zeta_order, difftype)

      use mod_parameter, only: PI
      implicit none
      class(difference_3D_type), intent(inout) :: this
      type(grid3d_type), Pointer, intent(in) :: grid
      integer, intent(in) :: xi_order, eta_order, zeta_order
      integer, intent(in) :: difftype
      integer :: in, jn, kn
      real(R_P), allocatable :: zz(:)
      integer, parameter :: DIFFBF=1000, DIFFDIS=1002

      if (.not. Grid%HasBodyFixed()) stop "The Grid must be a body-fixed one! Use Trans Wrapper to get one!"
      call grid%GetSize(in, jn, kn)
      this%in = in; this%jn = jn; this%kn=kn
      if(.not. (associated(this%grid, grid))) this%grid => grid
      if (eta_order .ne. 4) stop "Normal order must be 4 now!"
      !if (xi_order .ne. 1) stop 'streamwise order must be 1 now!'
      this%XIorder = xi_order
      this%ETAorder= eta_order
      if (difftype .ne. DIFFBF .and. difftype .ne. DIFFDIS) stop "Please input a correct difference type"
      this%difftype = difftype
      if( zeta_order == -1) then
        print*, 'Spectral method is used!'
        print*, 'The Grid corresponds to the whole period in the spanwise direction z(kn)-z(1)==L_z'
        if (zeta_order/2*2 - eta_order == 0) stop "With spectral method, kn must an odd number!"
        allocate(zz(kn))
        call grid%Get_ijz(1, 1, zz)
        this%beta0=2.0d0*PI/(zz(kn)-zz(1))
      end if

      !call initial_grid_coef(this) !TODO_fill in the initial programe

      call compute_Coef_eta(this)
      call compute_Coef_xi (this)

      !call finalize_grid_coef()


   end subroutine initial_difference_3d

   !> 计算基本流流向一阶导数
   !! @param[in] js 法向起始点序号
   !! @param[in] je 法向终止点序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_XI_BFFluxIJ(this, js, je, A, Diff)

       implicit none
       class(difference_object_type), intent(in) :: this
       integer, intent(in) :: js, je
       type(bf_flux_org_ij_type), intent(in) :: A(:, :)
       type(bf_flux_org_ij_type), intent(inout) :: Diff(:, :)
       integer :: i, j, l
       integer :: order
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002

       order=this%XIorder
       if(this%difftype==DIFFBF)then
         do j=js, je
            do i=1, order/2
              Diff(i,j)=BF_FLUXIJ_NULL
              do l=-order/2, order+2
                Diff(i, j)=Diff(i, j)+A(1+l+order/2, j)*this%Coef_di(l, i, j, 1)
              end do
            end do
            do i=order/2+1, this%in-order/2
              Diff(i,j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(i, j)=Diff(i, j)+A(i+l, j)*this%Coef_di(l, i, j, 1)
              enddo
            end do
            do i=this%in-order/2+1, this%in
              Diff(i,j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(i, j)=Diff(i, j)+A(this%in-order/2+l, j)*this%Coef_di(l, i, j, 1)
              end do
            end do
         end do
       elseif(this%difftype==DIFFDIS)then
         do j=js, je
            do i=1, order
              Diff(i,j)=BF_FLUXIJ_NULL
              do l=-order, 0
                Diff(i, j)=Diff(i, j)+A(i+l+order, j)*this%Coef_di(l, i, j, 1)
              end do
            end do
            do i=order+1, this%in
              Diff(i,j)=BF_FLUXIJ_NULL
              do l=-order, 0
                Diff(i, j)=Diff(i, j)+A(i+l, j)*this%Coef_di(l, i, j, 1)
              enddo
            end do
         end do
       endif

   end subroutine Partial_XI_BFFluxIJ

   !> 计算基本流在某一流向站位处的流向一阶导数
   !! @param[in] js 法向起始点序号
   !! @param[in] je 法向终止点序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_XI_BFFluxIJ_iloc(this, js, je, iloc, A, Diff)

       implicit none
       class(difference_object_type), intent(in) :: this
       integer, intent(in) :: js, je
       integer, intent(in) :: iloc
       type(bf_flux_org_ij_type), intent(in) :: A(:, :)
       type(bf_flux_org_ij_type), intent(inout) :: Diff(:)
       integer :: j, l
       integer :: order
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002

       order=this%XIorder
       if(this%difftype==DIFFBF)then
         do j=js, je
             Diff(j)=BF_FLUXIJ_NULL
             if( iloc>=1 .and. iloc<=order/2) then
                do l=-order/2, order/2
                    Diff(j)=Diff(j)+A(1+order/2+l, j)*this%Coef_di(l, iloc, j, 1)
                end do
             else if(iloc>=order/2+1 .and. iloc<=this%in-order/2)then
                do l=-order/2, order/2
                    Diff(j)=Diff(j)+A(iloc+l, j)*this%Coef_di(l, iloc, j, 1)
                enddo
             else if(iloc>=this%in-order/2+1 .and. iloc<=this%in)then
                do l=-order/2, order/2
                    Diff(j)=Diff(j)+A(this%in-order/2+l, j)*this%Coef_di(l, iloc, j, 1)
                end do
             endif
         enddo
       elseif(this%difftype==DIFFDIS)then
         do j=js, je
            Diff(j)=BF_FLUXIJ_NULL
            if(iloc>=1 .and. iloc<=order)then
              do l=-order, 0
                Diff(j)=Diff(j)+A(iloc+l+order, j)*this%Coef_di(l, iloc, j, 1)
              end do
            elseif(iloc>=order+1 .and. iloc<=this%in)then
              do l=-order, 0
                Diff(j)=Diff(j)+A(iloc+l, j)*this%Coef_di(l, iloc, j, 1)
              enddo
            endif
         end do
       endif

   end subroutine Partial_XI_BFFluxIJ_iloc

   !> 计算基本流法向一阶导数
   !! @param[in] is 流向起始点序号
   !! @param[in] ie 流向终止点序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_Eta_BFFluxIJ(this, is, ie, A, Diff)

       use mod_vector
       implicit none
       class(difference_object_type), intent(in) :: this
       integer, intent(in) :: is, ie
       type(bf_flux_org_ij_type), intent(in) :: A(:, :)
       type(bf_flux_org_ij_type), intent(inout) :: Diff(:, :)
       integer :: i, j, l
       integer :: order
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002
       real(R_P) :: rho, T, U, V, W
       type(vector_type) :: Vel

       order=this%ETAorder
         do i=is, ie
            do j=1, order/2
              Diff(i,j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(i, j)=Diff(i, j)+A(i, 1+order/2+l)*this%Coef_dj(l, i, j, 1)
              end do
            end do
            do j=order/2+1, this%jn-order/2
              Diff(i,j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(i, j)=Diff(i, j)+A(i, j+l)*this%Coef_dj(l, i, j, 1)
              enddo
            end do
            do j=this%jn-order/2+1, this%jn
              Diff(i,j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(i, j)=Diff(i, j)+A(i, this%jn-order/2+l)*this%Coef_dj(l, i, j, 1)
              end do
            end do
         end do

   end subroutine Partial_Eta_BFFluxIJ

   !> 计算基本流法向二阶导数
   !! @param[in] is 流向起始点序号
   !! @param[in] ie 流向终止点序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_Eta2_BFFluxIJ(this, is, ie, A, Diff)

       implicit none
       class(difference_object_type), intent(in) :: this
       integer, intent(in) :: is, ie
       type(bf_flux_org_ij_type), intent(in) :: A(:, :)
       type(bf_flux_org_ij_type), intent(inout) :: Diff(:, :)
       integer :: i, j, l
       integer :: order
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002

       order=this%ETAorder
         do i=is, ie
            do j=1, order/2
              Diff(i, j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(i, j)=Diff(i, j)+A(i, 1+order/2+l)*this%Coef_djj(l, i, j, 1)
              end do
            end do
            do j=order/2+1, this%jn-order/2
              Diff(i, j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(i, j)=Diff(i, j)+A(i, j+l)*this%Coef_djj(l, i, j, 1)
              enddo
            end do
            do j=this%jn-order/2+1, this%jn
              Diff(i, j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(i, j)=Diff(i, j)+A(i, this%jn-order/2+l)*this%Coef_djj(l, i, j, 1)
              end do
            end do
         end do

   end subroutine Partial_Eta2_BFFluxIJ

   !> 计算基本流法向一阶导数
   !! @param[in] iloc 流向站位序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_Eta_BFFluxIJ_iloc(this, iloc, A, Diff)

       use mod_vector
       implicit none
       class(difference_object_type), intent(in) :: this
       integer, intent(in) :: iloc
       type(bf_flux_org_ij_type), intent(in) :: A(:, :)
       type(bf_flux_org_ij_type), intent(inout) :: Diff(:)
       integer :: i, j, l
       integer :: order
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002
       real(R_P) :: rho, T, U, V, W
       type(vector_type) :: Vel

       order=this%ETAorder
         do i=iloc, iloc
            do j=1, order/2
              Diff(j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(j)=Diff(j)+A(i, 1+order/2+l)*this%Coef_dj(l, i, j, 1)
              end do
            end do
            do j=order/2+1, this%jn-order/2
              Diff(j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(j)=Diff(j)+A(i, j+l)*this%Coef_dj(l, i, j, 1)
              enddo
            end do
            do j=this%jn-order/2+1, this%jn
              Diff(j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(j)=Diff(j)+A(i, this%jn-order/2+l)*this%Coef_dj(l, i, j, 1)
              end do
            end do
         end do

   end subroutine Partial_Eta_BFFluxIJ_iloc

   !> 计算基本流法向二阶导数
   !! @param[in] iloc 流向序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_Eta2_BFFluxIJ_iloc(this, iloc, A, Diff)

       implicit none
       class(difference_object_type), intent(in) :: this
       integer, intent(in) :: iloc
       type(bf_flux_org_ij_type), intent(in) :: A(:, :)
       type(bf_flux_org_ij_type), intent(inout) :: Diff(:)
       integer :: i, j, l
       integer :: order
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002

       order=this%ETAorder
         do i=iloc, iloc
            do j=1, order/2
              Diff(j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(j)=Diff(j)+A(i, 1+order/2+l)*this%Coef_djj(l, i, j, 1)
              end do
            end do
            do j=order/2+1, this%jn-order/2
              Diff(j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(j)=Diff(j)+A(i, j+l)*this%Coef_djj(l, i, j, 1)
              enddo
            end do
            do j=this%jn-order/2+1, this%jn
              Diff(j)=BF_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(j)=Diff(j)+A(i, this%jn-order/2+l)*this%Coef_djj(l, i, j, 1)
              end do
            end do
         end do

   end subroutine Partial_Eta2_BFFluxIJ_iloc


   !> 计算扰动法向一阶导数
   !! @param[in] is 流向起始点序号
   !! @param[in] ie 流向终止点序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_Eta_DisFluxIJ(this, iloc, A, Diff)

       implicit none
       class(difference_object_type), intent(in) :: this
       integer, intent(in) :: iloc
       type(dis_flux_ij_type), intent(in) :: A(:)
       type(dis_flux_ij_type), intent(inout) :: Diff(:)
       integer :: i, j, l
       integer :: order
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002

       order=this%ETAorder
            do j=1, order/2
              Diff(j)=DIS_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(j)=Diff(j)+A(l+order/2+1)*this%Coef_dj(l, iloc, j, 1)
              end do
            end do
            do j=order/2+1, this%jn-order/2
              Diff(j)=DIS_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(j)=Diff(j)+A(j+l)*this%Coef_dj(l, iloc, j, 1)
              enddo
            end do
            do j=this%jn-order/2+1, this%jn
              Diff(j)=DIS_FLUXIJ_NULL
              do l=-order/2, order/2
                Diff(j)=Diff(j)+A(this%jn+l-order/2)*this%Coef_dj(l, iloc, j, 1)
              end do
            end do

   end subroutine Partial_Eta_DisFluxIJ

   !> 计算扰动法向二阶导数
   !! @param[in] is 流向起始点序号
   !! @param[in] ie 流向终止点序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_Eta2_DisFluxIJ(this, iloc, A, Diff)

       implicit none
       class(difference_object_type), intent(in) :: this
       integer, intent(in) :: iloc
       type(dis_flux_ij_type), intent(in) :: A(:)
       type(dis_flux_ij_type), intent(inout) :: Diff(:)
       integer :: i, j, l
       integer :: order
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002

       order=this%ETAorder
        do j=1, order/2
            Diff(j)=DIS_FLUXIJ_NULL
            do l=-order/2, order/2
            Diff(j)=Diff(j)+A(1+l+order/2)*this%Coef_djj(l, iloc, j, 1)
            end do
        end do
        do j=order/2+1, this%jn-order/2
            Diff(j)=DIS_FLUXIJ_NULL
            do l=-order/2, order/2
            Diff(j)=Diff(j)+A(j+l)*this%Coef_djj(l, iloc, j, 1)
            enddo
        end do
        do j=this%jn-order/2+1, this%jn
            Diff(j)=DIS_FLUXIJ_NULL
            do l=order/2, order/2
            Diff(j)=Diff(j)+A(this%jn+l-order/2)*this%Coef_djj(l, iloc, j, 1)
            end do
        end do

   end subroutine Partial_Eta2_DisFluxIJ

   !> 计算基本流展向一阶导数
   !! @param[in] is 流向起始点序号
   !! @param[in] ie 流向终止点序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_Zeta_BFFluxIJ(this, A, Diff)

       use mod_parameter, only: ZI=>CPLI
       implicit none
       class(difference_3D_type), intent(in) :: this
       type(bf_flux_org_ij_type), intent(in) :: A(this%kn)
       type(bf_flux_org_ij_type), intent(inout) :: Diff(this%kn)
       type(dis_flux_ij_type) :: zA(this%kn)
       type(dis_flux_ij_type) :: zDiff(this%kn)
       integer :: i, k, l
       integer :: m
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002

       m=(this%kn-1)/2

       call FFT(A, zA)
         do l=0, m
            zDiff(l+1)=ZI*this%beta0*real(l, 8)*zA(l+1)
         end do
        if(this%kn-2*m==2) zDiff(m+2)=DIS_FLUXIJ_NULL

         do l=-m, -1
            zDiff(this%kn+1+l)=ZI*this%beta0*real(l, 8)*zA(this%kn+1+l)
         end do
       call iFFT(zDiff, Diff)

   end subroutine Partial_Zeta_BFFluxIJ

   !> 计算扰动展向一阶导数
   !! @param[in] is 流向起始点序号
   !! @param[in] ie 流向终止点序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_Zeta_DisFluxIJ(this, A, Diff)

       use mod_parameter, only: ZI=>CPLI
       implicit none
       class(difference_3D_type), intent(in) :: this
       type(dis_flux_ij_type), intent(in) :: A(this%kn)
       type(dis_flux_ij_type), intent(inout) :: Diff(this%kn)
       type(dis_flux_ij_type) :: zA(this%kn), zDiff(this%kn)
       integer :: i, k, l
       integer :: m
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002

       m=(this%kn-1)/2

       call FFT(A, zA)
        do l=0, m
            zDiff(l+1)=ZI*this%beta0*real(l, 8)*zA(l+1)
        end do
        if(this%kn-2*m==2) zDiff(m+2)=DIS_FLUXIJ_NULL
        do l=-m, -1
            zDiff(this%kn+1+l)=ZI*this%beta0*real(l, 8)*zA(this%kn+1+l)
        end do
        call iFFT(zDiff, Diff)

   end subroutine Partial_Zeta_DisFluxIJ

   !> 计算基本流展向二阶导数
   !! @param[in] is 流向起始点序号
   !! @param[in] ie 流向终止点序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_Zeta2_BFFluxIJ(this, A, Diff)

       use mod_parameter, only: ZI=>CPLI
       implicit none
       class(difference_3D_type), intent(in) :: this
       type(bf_flux_org_ij_type), intent(in) :: A(this%kn)
       type(bf_flux_org_ij_type), intent(inout) :: Diff(this%kn)
       type(dis_flux_ij_type) :: zA(this%kn), zDiff(this%kn)
       integer :: i, k, l
       integer :: m
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002

       m=(this%kn-1)/2

       call FFT(A, zA)
         do l=0, m
            zDiff(l+1)=(-this%beta0**2)*real(l, 8)**2*zA(l+1)
         end do
        if(this%kn-2*m==2) zDiff(m+2)=DIS_FLUXIJ_NULL
         do l=-m, -1
            zDiff(this%kn+1+l)=(-this%beta0**2)*real(l, 8)**2*zA(this%kn+1+l)
         end do
       call iFFT(zDiff, Diff)

   end subroutine Partial_Zeta2_BFFluxIJ

   !> 计算扰动展向二阶导数
   !! @param[in] is 流向起始点序号
   !! @param[in] ie 流向终止点序号
   !! @param[in] A 待计算导数的原函数矩阵
   !! @param[out] Diff 差分结果
   subroutine Partial_Zeta2_DisFluxIJ(this, A, Diff)

       use mod_parameter, only: ZI=>CPLI
       implicit none
       class(difference_3D_type), intent(in) :: this
       type(dis_flux_ij_type), intent(in) :: A(this%kn)
       type(dis_flux_ij_type), intent(inout) :: Diff(this%kn)
       type(dis_flux_ij_type) :: zA(this%kn), zDiff(this%kn)
       integer :: i, k, l
       integer :: m
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002

       m=(this%kn-1)/2
       call FFT(A, zA)
        do l=0, m
            zDiff(l+1)=ZI*this%beta0*real(l, 8)*zA(l+1)
        end do
        if(this%kn-2*m==2) zDiff(m+2)=DIS_FLUXIJ_NULL

        do l=-m, -1
            zDiff(this%kn+1+l)=ZI*this%beta0*real(l, 8)*zA(this%kn+1+l)
        end do
       call iFFT(zDiff, Diff)

   end subroutine Partial_Zeta2_DisFluxIJ

   !> 计算流向导数基函数
   subroutine compute_Coef_xi(this)     !TODO  the difference 3d part

       implicit none
       class(difference_object_type), intent(inout) :: this
       integer :: order
       integer :: in, jn
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002

       order=this%XIorder;       in=this%in ; jn=this%jn
       if(this%difftype==DIFFBF) then
          allocate (this%Coef_di(-order/2:order/2, in, jn, 1))
          select type (self=>this)
          type is (difference_2d_type)
            call get_Coef_di_bf_2d(self%grid, order, self%Coef_di(:,:,:,1))
          type is (difference_3d_type)
            call get_Coef_di_bf_3d(self%grid, order, self%Coef_di(:,:,:,1))
          end select

       else if(this%difftype==DIFFDIS) then
          allocate (this%Coef_di(-order:0, in, jn, 1))
          select type (self=>this)
          type is (difference_2d_type)
            call get_Coef_di_dis_2d(self%grid, order, self%Coef_di(:,:,:,1))
          type is (difference_3d_type)
            call get_Coef_di_dis_3d(self%grid, order, self%Coef_di(:,:,:,1))
          end select
       else
          stop "Please input a correct difference type"
       end if

   end subroutine compute_Coef_xi


   !< 计算法向导数基函数
   subroutine compute_Coef_eta(this)

      implicit none
      class(difference_object_type), intent(inout) :: this
      integer :: order
      integer :: in, jn

      order= this%ETAorder
      in= this%in; jn= this%jn

!      allocate (this%Coef_di(-this%order/2:this%order/2, this%in, this%jn))
!      call get_Coef_di(this, this%difftype)
      allocate (this%Coef_dj (-order/2:order/2, in, jn, 1))
      allocate (this%Coef_djj(-order/2:order/2, in, jn, 1))
      call get_Coef_djdjj(this)

   end subroutine compute_Coef_eta

   !< 获得扰动流向插值基函数
   !! @param[in] order 离散精度
   !! @param[out] Coef_di 流向一阶导数插值基函数
   subroutine get_Coef_di_dis_2d(grid, order, Coef_di)

      implicit none
      integer, intent(in)   :: order
      real(R_P), intent(inout) :: Coef_di(-1:, :, :)
      type(grid2d_type), intent(in)   :: grid
      real(R_P), allocatable  :: xx(:)
      integer :: in, jn
      integer :: i, j
      integer :: order_bd
      integer, parameter :: one=1
      real(R_P) :: dx

      Coef_di=0.0d0
      call grid%GetSize(in, jn)
      allocate (xx(in))
      do j = 1, jn
         call grid%Get_jx(j, xx)
         do i = 1,  1
           dx=xx(2)-xx(1)
           Coef_di(-1, i, j)=-1.0d0/dx
           Coef_di(0, i, j)=1.0d0/dx
         enddo
         do i = 2, in
           dx=xx(i)-xx(i-1)
           Coef_di(-1, i, j)=-1.0d0/dx
           Coef_di(0, i, j)=1.0d0/dx
         end do
      end do
      deallocate (xx)

   end subroutine get_Coef_di_dis_2d

   !< 获得扰动流向插值基函数
   !! @param[in] order 离散精度
   !! @param[out] Coef_di 流向一阶导数插值基函数
   subroutine get_Coef_di_dis_3d(grid, order, Coef_di)

      implicit none
      integer, intent(in)   :: order
      real(R_P), intent(inout) :: Coef_di(:, :, :)
      type(grid3d_type), intent(in)   :: grid
      real(R_P), allocatable  :: xx(:)
      integer :: in, jn, kn
      integer :: i, j
      integer :: order_bd
      integer, parameter :: one=1
      real(R_P) :: dx

      !TODO_write subroutine

      ! Coef_di=0.0d0
      ! call grid%GetSize(in, jn, kn)
      ! allocate (xx(in))
      ! do j = 1, jn
      !    call grid%Get_jx(j, xx)
      !    do i = 1,  1
      !      dx=xx(2)-xx(1)
      !      Coef_di(-1, i, j)=1.0d0/dx
      !      Coef_di(0, i, j)=-1.0d0/dx
      !    enddo
      !    do i = 2, in
      !      dx=xx(i)-xx(i-1)
      !      Coef_di(-1, i, j)=1.0d0/dx
      !      Coef_di(0, i, j)=1.0d0/dx
      !    end do
      ! end do
      ! deallocate (xx)

   end subroutine get_Coef_di_dis_3d

   !> 获得基本流流向插值基函数
   !! @param[in] order 离散精度
   !! @param[out] Coef_di 流向一阶导数插值基函数
   subroutine get_Coef_di_bf_2d(grid, order, Coef_di)

      implicit none
      integer, intent(in)   :: order
      real(R_P), intent(inout) :: Coef_di(:, :, :)
      class(grid2d_type), intent(in)   :: grid
      integer :: in, jn
      integer :: i, j
      real(R_P) :: ry, dq(-2:2)

      Coef_di=0.0d0

      call grid%GetSize(in, jn)

      ry=real(in-1, R_P)
      do j = 1, jn
         do i = 1, 1
            dq=0.0d0
            dq(-2)=ry*(-3.d0/2.d0)
            dq(-1)=ry*2.0d0
            dq( 0)=ry*(-1.d0/2.d0)
            dq( 1)=0.0d0
            dq( 2)=0.0d0
            Coef_di(:, i, j)=abx(i, j)*dq
         enddo
         do i=2, 2
            dq=0.0d0
            dq(-2)=ry*(-2.d0/6.d0)
       	    dq(-1)=ry*(-3.d0/6.d0)
       	    dq( 0)=ry
            dq( 1)=ry*(-1.d0/6.d0)
            dq( 2)=0.0d0
            Coef_di(:, i, j)=abx(i, j)*dq
         enddo
         do i = 3, in - 2
            dq(-2:2)=0.d0
            dq(-2)=ry/12.d0
       	    dq(-1)=ry*(-8.d0/12.d0)
       	    dq( 0)=0.d0
       	    dq( 1)=ry*(8.d0/12.d0)
            dq( 2)=ry*(-1.d0/12.d0)
            Coef_di(:, i, j)=abx(i, j)*dq
         end do
         do i = in - 1, in-1
           dq=0.0d0
           dq(-2)=0.0d0
           dq(-1)=ry*(1.d0/6.d0)
           dq( 0)=-1.0d0*ry
           dq( 1)=ry*(3.d0/6.d0)
           dq( 2)=ry*(2.d0/6.d0)
           Coef_di(:, i, j)=abx(i, j)*dq
         end do
         do i= in, in
           dq=0.0d0
           dq(-2)=0.0d0
           dq(-1)=0.0d0
           dq( 0)=ry*(1.d0/2.d0)
           dq( 1)=ry*(-2.0d0)
           dq( 2)=ry*(3.d0/2.d0)
           Coef_di(:, i, j)=abx(i, j)*dq
         enddo
      end do

   end subroutine get_Coef_di_bf_2d

   !> 获得基本流流向插值基函数
   !! @param[in] order 离散精度
   !! @param[out] Coef_di 流向一阶导数插值基函数
   subroutine get_Coef_di_bf_3d(grid, order, Coef_di)

      implicit none
      integer, intent(in)   :: order
      real(R_P), intent(inout) :: Coef_di(:, :, :)
      class(grid3d_type), intent(in)   :: grid
      real(R_P), allocatable  :: xx(:)
      integer :: in, jn, kn
      integer :: i, j
      integer :: order_bd

      Coef_di=0.0d0

      !TODO_write the subroutine

      ! order_bd=order-2
      !
      ! call grid%GetSize(in, jn, kn)
      ! allocate (xx(in))
      ! do j = 1, jn
      !    call grid%Get_jkx(j, 1, xx)
      !    do i = 1, order_bd
      !       call dslad_order_d1(order_bd, xx(1:order_bd + 1), -order_bd + i, &
      !                           Coef_di(:, i, j))
      !    enddo
      !    do i = 1 + order_bd, in - order_bd
      !       call dslad_order_d1(order, xx(i - order/2:i + order/2), 0, Coef_di(:, i, j))
      !    end do
      !    do i = in - order_bd + 1, in
      !       call dslad_order_d1(order_bd, xx(in - order_bd:in), i - in + order_bd-1, &
      !                           Coef_di(order-order_bd+1:order+1, i, j))
      !    end do
      ! end do
      ! deallocate (xx)

   end subroutine get_Coef_di_bf_3d

   !> 获得法向插值基函数
   subroutine get_Coef_djdjj(this)

       implicit none
       class(difference_object_type), intent(inout) :: this
       integer, parameter :: DIFFBF=1000, DIFFDIS=1002
       integer :: order
       Integer :: difftype

       order=this%ETAorder
       Difftype=this%Difftype

       select case (difftype)
           case (DIFFBF)
              select type (self=>this)
              type is (difference_2D_type)
                 call get_Coef_djdjj_bf_2d(self%grid, order, self%Coef_dj(:,:,:,1), self%Coef_djj(:,:,:,1))
              type is (difference_3d_type)
                 call get_Coef_djdjj_bf_3d(self%grid, order, self%Coef_dj(:,:,:,1), self%Coef_djj(:,:,:,1))
              end select
           case (DIFFDIS)
              select type (self=>this)
              type is (difference_2D_type)
                call get_Coef_djdjj_dis_2d(self%grid, order, self%Coef_dj(:,:,:,1), self%Coef_djj(:,:,:,1))
              type is (difference_3d_type)
                call get_Coef_djdjj_dis_3d(self%grid, order, self%Coef_dj(:,:,:,1), self%Coef_djj(:,:,:,1))
              end select
           case default
             stop "Please input a correct difference type"
       end select

   end subroutine get_Coef_djdjj

   !> 获得基本流法向插值基函数
   !! @param[in] order 法向离散精度
   !! @param[out] Coef_dj 法向一阶导数插值基函数
   !! @param[out] Coef_djj 法向二阶导数插值基函数
   subroutine get_Coef_djdjj_bf_2d(grid, order, Coef_dj, Coef_djj)

      implicit none
      integer, intent(in) :: order
      real(R_P), intent(inout) :: Coef_dj(:, :, :), Coef_djj(:, :, :)
      type(grid2d_type), Pointer, intent(in) :: grid
      real(R_P) :: ry, ry2, dq(-2:2), dr(-2:2)
      integer :: in, jn
      integer :: i, j
      integer, parameter :: zero = 0

      Coef_dj=0.0d0; Coef_djj=0.0d0

      call grid%GetSize(in, jn)
      ry=real(jn-1, R_P)
      ry2=ry*ry
      do i = 1, in
         do j = 1, 1
            dq=0.0d0; dr=0.0d0
            dq(-2)=ry*(-3.d0/2.d0)
            dq(-1)=ry*2.0d0
            dq( 0)=ry*(-1.d0/2.d0)
            dq( 1)=0.0d0
            dq( 2)=0.0d0
            Coef_dj(:, i, j)=bby(i, j)*dq
            Coef_djj(:, i, j)=bby(i, j)**2*dq+bbyy(i, j)*dr
         enddo
         do j=2, 2
            dq=0.0d0; dr=0.0d0
            dq(-2)=ry*(-2.d0/6.d0)
       	    dq(-1)=ry*(-3.d0/6.d0)
       	    dq( 0)=ry
            dq( 1)=ry*(-1.d0/6.d0)
            dq( 2)=0.0d0
            dr(-2)=ry2
        	  dr(-1)=ry2*(-2.d0)
        	  dr( 0)=ry2
            dr( 1)=0.0d0
            dr( 2)=0.0d0
            Coef_dj(:, i, j)=bby(i, j)*dq
            Coef_djj(:, i, j)=bby(i, j)**2*dq+bbyy(i, j)*dr
         enddo
         do j = 3, jn - 2
            dq(-2:2)=0.d0; dr=0.0d0
            dq(-2)=ry/12.d0
       	    dq(-1)=ry*(-8.d0/12.d0)
       	    dq( 0)=0.d0
       	    dq( 1)=ry*(8.d0/12.d0)
            dq( 2)=ry*(-1.d0/12.d0)
            dr(-2)=ry2*(-1.d0/12.d0)
        	  dr(-1)=ry2*(16.d0/12.d0)
        	  dr( 0)=ry2*(-30.d0/12.d0)
        	  dr( 1)=ry2*(16.d0/12.d0)
        	  dr( 2)=ry2*(-1.d0/12.d0)
            Coef_dj(:, i, j)=bby(i, j)*dq
            Coef_djj(:, i, j)=bby(i, j)**2*dr+bbyy(i, j)*dq
         end do
         do j = jn - 1, jn-1
           dq=0.0d0; dr=0.0d0
           dq(-2)=0.0d0
           dq(-1)=ry*(1.d0/6.d0)
           dq( 0)=-1.0d0*ry
           dq( 1)=ry*(3.d0/6.d0)
           dq( 2)=ry*(2.d0/6.d0)
           dr(-2)=0.0d0
           dr(-1)=0.0d0
           dr( 0)=ry2
           dr( 1)=ry2*(-2.d0)
           dr( 2)=ry2
           Coef_dj(:, i, j)=bby(i, j)*dq
           Coef_djj(:, i, j)=bby(i, j)**2*dr+bbyy(i, j)*dq
         end do
         do j= jn, jn
           dq=0.0d0; dr=0.0d0
           dq(-2)=0.0d0
           dq(-1)=0.0d0
           dq( 0)=ry*(1.d0/2.d0)
           dq( 1)=ry*(-2.0d0)
           dq( 2)=ry*(3.d0/2.d0)
           Coef_dj(:, i, j)=bby(i, j)*dq
           Coef_djj(:, i, j)=bby(i, j)**2*dr+bbyy(i, j)*dq
         enddo
      end do

   end subroutine get_Coef_djdjj_bf_2d

   !> 获得基本流法向插值基函数
   !! @param[in] order 法向离散精度
   !! @param[out] Coef_dj 法向一阶导数插值基函数
   !! @param[out] Coef_djj 法向二阶导数插值基函数
   subroutine get_Coef_djdjj_bf_3d(grid, order, Coef_dj, Coef_djj)

      implicit none
      integer, intent(in) :: order
      real(R_P), intent(inout) :: Coef_dj(:, :, :), Coef_djj(:, :, :)
      type(grid3d_type), Pointer, intent(in) :: grid
      real(R_P), allocatable :: yy(:)
      integer :: in, jn, kn
      integer :: i, j
      integer :: order_bd
      integer, parameter :: zero = 0

      Coef_dj=0.0d0; Coef_djj=0.0d0
      !TODO_write the subroutine

   end subroutine get_Coef_djdjj_bf_3d

   !> 获得扰动法向插值基函数
   !! @param[in] order 法向离散精度
   !! @param[out] Coef_dj 法向一阶导数插值基函数
   !! @param[out] Coef_djj 法向二阶导数插值基函数
   subroutine get_Coef_djdjj_dis_2d(grid, order, Coef_dj, Coef_djj)

      implicit none
      integer, intent(in) :: order
      real(R_P), intent(inout) :: Coef_dj(:, :, :), Coef_djj(:, :, :)
      type(grid2d_type), Pointer, intent(in) :: grid
      real(R_P) :: ry, ry2, dq(-2:2), dr(-2:2)
      integer :: in, jn
      integer :: i, j
      integer, parameter :: zero = 0

      Coef_dj=0.0d0; Coef_djj=0.0d0

      call grid%GetSize(in, jn)
      ry=real(jn-1, R_P)
      ry2=ry*ry
      do i = 1, in
         do j = 1, 1
            dq=0.0d0
            dq(-2)=ry*(-3.d0/2.d0)
            dq(-1)=ry*2.0d0
            dq( 0)=ry*(-1.d0/2.d0)
            dq( 1)=0.0d0
            dq( 2)=0.0d0
            Coef_dj(:, i, j)=bby(i, j)*dq
            Coef_djj(:, i, j)=bby(i, j)**2*dr+bbyy(i, j)*dq
         enddo
         do j=2, 2
            dq=0.0d0; dr=0.0d0
            dq(-2)=ry*(-2.d0/6.d0)
       	    dq(-1)=ry*(-3.d0/6.d0)
       	    dq( 0)=ry
            dq( 1)=ry*(-1.d0/6.d0)
            dq( 2)=0.0d0
            dr(-2)=ry2
        	  dr(-1)=ry2*(-2.d0)
        	  dr( 0)=ry2
            dr( 1)=0.0d0
            dr( 2)=0.0d0
            Coef_dj(:, i, j)=bby(i, j)*dq
            Coef_djj(:, i, j)=bby(i, j)**2*dr+bbyy(i, j)*dq
         enddo
         do j = 3, jn - 2
            dq(-2:2)=0.0d0; dr=0.0d0
            dq(-2)=ry/12.d0
       	    dq(-1)=ry*(-8.d0/12.d0)
       	    dq( 0)=0.0d0
       	    dq( 1)=ry*(8.d0/12.d0)
            dq( 2)=ry*(-1.d0/12.d0)
            dr(-2)=ry2*(-1.d0/12.d0)
        	  dr(-1)=ry2*(16.d0/12.d0)
        	  dr( 0)=ry2*(-30.d0/12.d0)
        	  dr( 1)=ry2*(16.d0/12.d0)
        	  dr( 2)=ry2*(-1.d0/12.d0)
            Coef_dj(:, i, j)=bby(i, j)*dq
            Coef_djj(:, i, j)=bby(i, j)**2*dr+bbyy(i, j)*dq
         end do
         do j = jn - 1, jn-1
           dq=0.0d0; dr=0.0d0
           dq(-2)=0.0d0
           dq(-1)=ry*(1.d0/6.d0)
           dq( 0)=-1.0d0*ry
           dq( 1)=ry*(3.d0/6.d0)
           dq( 2)=ry*(2.d0/6.d0)
           dr(-2)=0.0d0
           dr(-1)=0.0d0
           dr( 0)=ry2
           dr( 1)=ry2*(-2.d0)
           dr( 2)=ry2
           Coef_dj(:, i, j)=bby(i, j)*dq
           Coef_djj(:, i, j)=bby(i, j)**2*dr+bbyy(i, j)*dq
         end do
         do j= jn, jn
           dq=0.0d0; dr=0.0d0
           dq(-2)=0.0d0
           dq(-1)=0.0d0
           dq( 0)=ry*(1.d0/2.d0)
           dq( 1)=ry*(-2.0d0)
           dq( 2)=ry*(3.d0/2.d0)
           Coef_dj(:, i, j)=bby(i, j)*dq
           Coef_djj(:, i, j)=bby(i, j)**2*dr+bbyy(i, j)*dq
         enddo
      end do

   end subroutine get_Coef_djdjj_dis_2d

   !> 获得扰动法向插值基函数
   !! @param[in] order 法向离散精度
   !! @param[out] Coef_dj 法向一阶导数插值基函数
   !! @param[out] Coef_djj 法向二阶导数插值基函数
   subroutine get_Coef_djdjj_dis_3d(grid, order, Coef_dj, Coef_djj)

      implicit none
      integer, intent(in) :: order
      real(R_P), intent(inout) :: Coef_dj(:, :, :), Coef_djj(:, :, :)
      type(grid3d_type), intent(in) :: grid
      real(R_P), allocatable :: yy(:)
      integer :: in, jn, kn
      integer :: i, j
      integer :: order_bd
      integer, parameter :: zero = 0

      Coef_dj=0.0d0; Coef_djj=0.0d0
      !TODO_write the subroutine

   end subroutine get_Coef_djdjj_dis_3d

   !> 获得流向导数精度
   !! @return 流向导数精度
   integer function get_xiorder(this)

      implicit none
      class(difference_object_type), intent(in) :: this

      get_xiorder = this%xiorder

   end function get_xiorder

   !> 获得法向导数精度
   !! @return 法向导数精度
   integer function get_etaorder(this)

      implicit none
      class(difference_object_type), intent(in) :: this

      get_etaorder = this%etaorder

   end function get_etaorder

   !> 获得某点处的流向导数插值模板
   !!
   !! 导数类型|对应常量
   !! :-----------:|:------:
   !! 法向一阶导数|ETA1
   !! 法向二阶导数|ETA2
   !! @param[in] i_index 流向位置序号
   !! @param[in] j_index 法向位置序号
   !! @param[in] Flag 导数类型
   !! @param[out] Coef 插值模板系数
   subroutine Get_Coef_ETA(this, i_index, j_index, Flag, Coef)

      implicit none
      class(difference_object_type), intent(in) :: this
      integer, intent(in) :: i_index, j_index
      integer :: Flag
      real(R_P), intent(inout) :: Coef(1:this%ETAorder + 1)
      integer, parameter :: XI1 = 1, ETA1 = 2, ETA2 = 3

      select case (Flag)
      case (ETA1)
         Coef = this%Coef_dj(:, i_index, j_index, 1)
      case (ETA2)
         Coef = this%Coef_djj(:, i_index, j_index, 1)
      case default
         write (*, *) "the Coef parameter must be the one of 'xi1', 'eta1', 'eta2'"
      end select

   end subroutine get_Coef_ETA

   !> 获得某点处的流向导数插值模板
   !!
   !! 导数类型|对应常量
   !! :-----------|:------:
   !! 流向一阶导数|XI1
   !! @param[in] i_index 流向位置序号
   !! @param[in] j_index 法向位置序号
   !! @param[in] Flag 导数类型
   !! @param[out] Coef 插值模板系数
   subroutine Get_Coef_XI(this, i_index, j_index, Coef)

      implicit none
      class(difference_object_type), intent(in) :: this
      integer, intent(in) :: i_index, j_index
      real(R_P), intent(inout) :: Coef(1:this%XIorder + 1)
      integer, parameter :: XI1 = 1, ETA1 = 2, ETA2 = 3

      Coef = this%Coef_di(:, i_index, j_index, 1)

   end subroutine get_Coef_XI

   !> 析构函数
   Subroutine finalize_difference_2D(this)

      Implicit None
      class(difference_2D_type), intent(Inout) :: this

      this%in = 0; this%jn = 0
      deallocate (this%Coef_di)
      deallocate (this%Coef_dj)
      deallocate (this%Coef_djj)

   End Subroutine finalize_difference_2D

end module mod_difference
