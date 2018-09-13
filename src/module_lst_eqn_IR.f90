!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_lst_IR.f90
!> @file
!> @breif 反幂法求LST模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: LST_IR
!> @breif 反幂法求LST模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-03 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2017-08-03
!------------------------------------------------------------------------------
module mod_lst_eqn_IR

    use mod_BF_normal
    use mod_local_normal_coordinate
    use mod_lns_op_normal
    use mod_lst_dis_OP_normal
    use mod_difference
    use mod_baseflow
    use mod_curvature_2d
    use mod_grid
    use mod_dis_wavenum
    use mod_inverse_rayleigh
    use mod_sparse_matrix
    use mod_lpse_dis_normal
    use penf, only: R_P

    implicit none

    private

    !> 反幂法求LST解算器类
    type, public :: lst_eqn_IR_type

        integer, private :: iloc !< 流向站位
        integer, private :: jn !< 法向网格点数
        type(dis_wavenum_type), private :: wavenum !<扰动色散信息
        type(BF_normal_type), pointer, private :: BFNorm !<基本流法向分布信息
        type(lpse_dis_normal_type), Pointer, private :: DisNorm !<扰动流向某站位法向分布
        type(local_normal_coordinate_type), pointer, private :: NormCoord !<当地法向坐标系
        type(Norm_Diff_Coef_type), pointer, private :: NormCoef
        type(lst_bf_op_normal_type), private :: BFOPNorm !<基本流法向算子
        type(lst_dis_op_normal_type), private :: DisOPNorm !<扰动法向算子
        !type(grid2d_type), pointer, private :: Grid !< 场网格
        !type(difference_2D_type), pointer, private :: Diff !<扰动差分基函数
        !type(Baseflow2D_type), pointer, private :: BaseFlow !<2D基本流
!        type(Baseflow2D_Diff_type), pointer, private :: DBaseFlow !< 2D基本流导数
        !type(curvature_2d_type), pointer, private :: Curvature !<2D曲率
        type(inverse_rayleigh_type), private :: solver !<反幂法解算器
        type(mat_coo_type), private :: A !>广义特征值问题左矩阵
        type(mat_coo_type), private :: B !>广义特征值问题右矩阵
        real(R_P), allocatable, private :: Eta(:) !< 当地位置法向网格分布

        Contains

          !generic :: Solve => solve_iloc, solve_domain !<求解过程

          Procedure :: Create !<初始化求解器,分配内存
          !procedure :: Set !< 设置解算器需要的基本信息
          procedure :: SolveLST !<求解LST

          procedure, private :: SetInitialWavenum !<设置初始猜测的扰动色散信息(特征值)
          procedure, private :: GetShapefun !< 获得扰动形函数
          procedure, private :: GetShapefunAdjoint !<获得扰动伴随形函数
          !procedure :: GetIntoDis !< 将扰动信息返回到2D全场扰动
          procedure, private :: GetWavenum !< 获得扰动色散信息
          procedure, private :: SetBFOPNorm !< 设置基本流算子法向分布
          procedure, private :: SetMatrix !< 设置特征值问题的左右矩阵
          procedure, private :: GetEigenFunc !<获得特征函数到形函数
          procedure, private :: GetEigenFuncAdjoint !<获得特征函数的伴随到伴随形函数
          procedure, private :: GetEigenVal !<获得特征值
          procedure, private :: Solve_iloc !< 求解某个站位
          !procedure, private :: Solve_domain !< 求解整个域

    end type lst_eqn_IR_type

    contains

    !>初始化求解器,分配内存
    !! @param[in] jn 法向点数
    subroutine Create(this, jn)

        implicit none
        class(lst_eqn_IR_type), intent(inout) :: this
        integer, intent(in) :: jn
        integer :: eqkind

        this%jn=jn
        eqkind=1
        call this%BFOPNorm%Create(this%jn, eqkind)
        call this%DisOPNorm%Create(this%jn)
        allocate(this%eta(this%jn))

    end subroutine Create

    !> LST求解
    subroutine SolveLST(this, DisNorm, iloc, &
    &               BFNorm, NormCoord, Eta, NormCoef, isAdjoint)

        use mod_lpse_dis_normal

        implicit none
        class(lst_eqn_ir_type), intent(inout) :: this
        type(lpse_dis_normal_type), target, intent(inout) :: DisNorm
        integer, intent(in) :: iloc
        type(BF_normal_type), target, intent(in) :: BFNorm
        type(local_normal_coordinate_type), target, intent(in) :: NormCoord
        type(Norm_Diff_Coef_type), target, intent(in) :: NormCoef
        real(R_P), intent(in) :: Eta(:)
        logical, optional, intent(in) :: isAdjoint
        type(dis_wavenum_lpse_type) :: wavenum

        this%Iloc=iloc

        this%BFNorm=>BFNorm
        this%NormCoord=>NormCoord
        this%NormCoef=>NormCoef
        this%DisNorm=>DisNorm

        this%Eta=Eta

        call this%SetBFOPNorm(this%iloc)

        Wavenum=this%DisNorm%GetWaveNum()
        call this%SetInitialWavenum(wavenum)
        call this%solve_iloc()

        if(present(isAdjoint) .and. isAdjoint) then
          call this%GetShapefunAdjoint(DisNorm)
        else
          call this%GetShapefun(DisNorm)
        endif
        call DisNorm%SetILoc(iloc)

        wavenum=this%GetWavenum()
        call DisNorm%SetWaveNum(wavenum)

    end subroutine SolveLST

!     !> 设置解算器需要的基本信息
!     !! @param[in] iloc 流向站位
!     !! @param[in] Grid 2D网格
!     !! @param[in] BaseFlow 2D基本流
!     !! @param[in] Curvature 2D区率
!     !! @param[in] Diff 扰动差分基函数
!     subroutine Set(this, iloc, Grid, BaseFlow, Curvature, Diff)
!
!         implicit none
!         class(lst_eqn_IR_type), intent(inout) :: this
!         integer, intent(in) :: iloc
!         type(grid2d_type), target, intent(in) :: Grid
!         type(Baseflow2D_type), target, intent(in) :: Baseflow
! !        type(Baseflow2D_Diff_type), target, intent(in) :: DBaseflow
!         type(difference_2D_type), target, intent(in) :: Diff
!         type(curvature_2d_type), target, intent(in) :: Curvature
!
!         this%iloc=iloc
!         !if(.not. (associated(this%Grid, Grid))) this%Grid => Grid
!         !if(.not. (associated(this%BaseFlow, BaseFlow))) &
!         !&   this%BaseFlow => Baseflow
!         !if(.not. (associated(this%DBaseFlow, DBaseFlow))) &
!         !&   this%DBaseFlow => DBaseflow
!         !if(.not. (associated(this%Diff, Diff))) &
!         !&   this%Diff => Diff
!         !if(.not. (associated(this%Curvature, Curvature))) &
!         !&   this%Curvature => Curvature
!
!         call this%SetBFOPNorm(iloc)
!
!     end subroutine Set

    !> 设置基本流算子法向分布
    subroutine SetBFOPNorm(this, iloc)

        implicit none
        class(lst_eqn_IR_type), intent(inout) :: this
        integer, intent(in) :: iloc

!        call this%BFNorm%Set(iloc, this%BaseFlow, This%Diff)
!        call this%NormCoord%SetLameFromGrid(this%grid, iloc)
!        call this%NormCoord%SetCurvature(This%Curvature%GetPoint(iloc))
        call this%BFOPNorm%Set(this%BFNorm, this%NormCoord)

    end subroutine SetBFOPNorm

    !> 设置初始猜测的扰动色散信息(特征值)
    !! @param[in] wavenum 初始扰动色散信息(特征值和参数)
    subroutine SetInitialWavenum(this, wavenum)

        implicit none
        class(lst_eqn_IR_type), intent(inout) :: this
        class(dis_wavenum_type) :: wavenum
        complex(R_P) :: alpha, beta, omega

        call wavenum%Get(alpha, beta, omega)

        call this%wavenum%Set(alpha, beta, omega)

    end subroutine SetInitialWavenum

    !> 求解过程(单独Iloc站位)
    subroutine solve_iloc(this)

        use mod_sparse_matrix
        implicit none
        class(lst_eqn_IR_type), intent(inout) :: this

        call this%DisOPNorm%Set(this%BFOPNorm, this%NormCoord, this%wavenum)
        call this%SetMatrix()
        call this%solver%SetMatrix(this%jn*5*2, this%A, this%B)
        call this%solver%SetEigenValueGuess(this%wavenum%getAlpha())
        call this%solver%solve()

    end subroutine solve_iloc

    ! !> 求解过程(整个计算域）
    ! subroutine solve_domain(this, Istart, Iend, Dis)
    !
    !     use mod_sparse_matrix
    !     use mod_dis
    !     implicit none
    !     class(lst_eqn_IR_type), intent(inout) :: this
    !     integer, intent(inout) :: Istart
    !     integer, intent(inout) :: Iend
    !     type(dis_type), intent(inout) :: Dis
    !     type(dis_wavenum_type) :: wavenumBak
    !     integer :: i
    !     integer :: in, jn
    !
    !     wavenumBak=this%wavenum
    !     call this%Grid%GetSize(in, jn)
    !     if(Istart==-1) then
    !         Istart=1
    !     endif
    !
    !     if(Iend==-1) then
    !         Iend=in
    !     endif
    !
    !     do i=this%iloc, Iend
    !         call this%SetBFOPNorm(i)
    !         call this%SetInitialWavenum(this%wavenum)
    !         call this%solve()
    !         call this%GetIntoDis(i, Dis)
    !     enddo
    !     this%wavenum=wavenumBak
    !     do i=this%iloc-1, Istart, -1
    !         call this%SetBFOPNorm(i)
    !         call this%SetInitialWavenum(this%wavenum)
    !         call this%solve()
    !         call this%GetIntoDis(i, Dis)
    !     enddo
    !
    ! end subroutine solve_domain


    !> 获得特征函数
    !! @param[out] DisShape 特征函数(按形函数形式排列)
    subroutine GetEigenFunc(this, DisShape)

        implicit none
        class(lst_eqn_IR_type), intent(inout) :: this
        complex(R_P), intent(out) :: DisShape(5, this%jn)
        complex(R_P) :: Eigenfunc(10*this%jn)
        integer :: j
        real(R_P) :: AbsU(this%jn)

        call this%solver%GetEigenFunc(Eigenfunc)

        Disshape=reshape(eigenfunc(1:5*this%jn), [5, this%jn])

        absu=abs(DisShape(2, :))
        DisShape=DisShape/DisShape(2, maxloc(absu, 1))

    end subroutine GetEigenFunc

    !> 获得特征函数的伴随
    !! @param[out] DisShape 特征函数(按形函数形式排列)
    subroutine GetEigenFuncAdjoint(this, DisShape)

        implicit none
        class(lst_eqn_IR_type), intent(inout) :: this
        complex(R_P), intent(out) :: DisShape(5, this%jn)
        complex(R_P) :: Eigenfunc(10*this%jn)
        integer :: j
        real(R_P) :: AbsU(this%jn)

        call this%solver%GetEigenFuncAdjoint(Eigenfunc)

        Disshape=reshape(eigenfunc(1:5*this%jn), [5, this%jn])

        absu=abs(DisShape(2, :))
        DisShape=DisShape/DisShape(2, maxloc(absu, 1))

    end subroutine GetEigenFuncAdjoint

    !> 获得特征值
    !! @param[out] 特征值
    subroutine GetEigenVal(this, EigenVal)

        implicit none
        class(lst_eqn_IR_type), intent(inout) :: this
        complex(R_P), intent(out) :: EigenVal

        call this%solver%GetEigenValue(EigenVal)

    end subroutine GetEigenVal

    !> 获得扰动形函数
    !! @param[in, out] shapefun 扰动形函数
    subroutine GetShapefun(this, shapefun)

        use mod_dis_shape
        use mod_dis_flux
        use mod_vector_cmplx
        implicit none
        class(lst_eqn_IR_type), intent(inout) :: this
        class(dis_shape_type), intent(inout) :: shapefun
        type(dis_flux_ij_type) :: Flux(this%jn)
        type(vector_cmplx_type) :: aVel
        integer :: j
        complex(R_P) :: disshape(5, this%jn)

        call this%GetEigenFunc(disshape)

        call shapefun%Create(this%jn)

        do j=1, this%jn
          call aVel%set(disshape(2, j), disshape(3, j), disshape(4, j))
          call Flux(j)%Set(disshape(1, j), aVel, disshape(5, j))
       end do
       call shapefun%set(flux, 1)

    end subroutine GetShapefun

    !> 获得扰动伴随形函数
    !! @param[in, out] shapefun 扰动形函数的伴随
    subroutine GetShapefunAdjoint(this, shapefun)

        use mod_dis_shape
        use mod_dis_flux
        use mod_vector_cmplx
        implicit none
        class(lst_eqn_IR_type), intent(inout) :: this
        class(dis_shape_type), intent(inout) :: shapefun
        type(dis_flux_ij_type) :: Flux(this%jn)
        type(vector_cmplx_type) :: aVel
        integer :: j
        complex(R_P) :: disshape(5, this%jn)

        call this%GetEigenFuncAdjoint(disshape)

        call shapefun%Create(this%jn)

        do j=1, this%jn
          call aVel%set(disshape(2, j), disshape(3, j), disshape(4, j))
          call Flux(j)%Set(disshape(1, j), aVel, disshape(5, j))
       end do
       call shapefun%set(flux, 1)

    end subroutine GetShapefunAdjoint

    !> 获得扰动色散信息
    !! @return 扰动的色散信息
    function GetWavenum(this) result(wavenum)

        implicit none
        class(lst_eqn_IR_type), intent(inout) :: this
        type(dis_wavenum_lpse_type) :: wavenum
        complex(R_P) :: alpha, beta, omega

        call this%wavenum%Get(alpha, beta, omega)
        call this%solver%GetEigenValue(alpha)
        call wavenum%Set(alpha, beta, omega)

    end function GetWavenum

    ! !> 将扰动信息返回到2D全场扰动
    ! !! @param[in, out] Dis 2D全场扰动
    ! subroutine GetIntoDis(this, iloc, Dis)
    !
    !     use mod_dis
    !     use mod_lpse_dis_normal
    !     implicit none
    !     class(lst_eqn_IR_type), intent(inout) :: this
    !     integer, intent(in) :: iloc
    !     type(dis_type), intent(inout) :: Dis
    !     type(lpse_dis_normal_type) :: DisNorm
    !     type(dis_wavenum_lpse_type) :: wavenum
    !
    !     call this%GetShapefun(DisNorm)
    !     call DisNorm%SetILoc(iloc)
    !
    !     wavenum=this%GetWavenum()
    !     call DisNorm%SetWaveNum(wavenum)
    !     print*, '============', 'iloc=', iloc, '============='
    !     call wavenum%print()
    !     call this%SetInitialWavenum(wavenum)
    !     call Dis%set(iloc, DisNorm)
    !     call Dis%Print(iloc, this%Grid)
    !
    ! end subroutine GetIntoDis

    !> 设置特征值问题的左右矩阵
    subroutine SetMatrix(this)

        use mod_lst_dis_OP_point
        use mod_sparse_matrix
        use mod_difference
        implicit none

        class(lst_eqn_IR_type), intent(inout), target :: this
        real(R_P) :: coef_dy(5), coef_dyy(5)
        type(lst_dis_op_point_type) :: DisOP
        complex(R_P), dimension(5, 5) :: BL, BR, DL, DR1, DR2, Vyy
        integer :: j, l
        complex(R_P) :: tmpij(-2:2, 5, 5), tmpMatrix(5, 5)
        type(mat_coo_type), pointer :: A, B

        integer, parameter :: XI1 = 1, ETA1 = 2, ETA2 = 3

        if(.not. (associated(A, this%A)))then
          A => this%A; B => this%B
        endif

        if( .not.A%ISCreate()) then
            call A%Create(this%jn*10, (this%jn*5*5*5-3*2*5*5)+this%jn*5)
            call B%Create(this%jn*10, (this%jn*5*5*5-3*2*5*5+this%jn*5+5*5*this%jn))
        else
            call A%zeros()
            call B%zeros()
        endif

        associate(iloc => this%iloc, DisOPNorm => this%DisOPNorm, jn => this%jn)

            j=1
            !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
            !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
            Coef_dy=this%NormCoef%Coef_dj(:, j)
            Coef_dyy=this%NormCoef%Coef_djj(:, j)
            DisOP=DisOPNorm%GetPoint(j)
            call DisOP%Get(BL, BR, DL, DR1, DR2, Vyy)
            tmpIJ=0.0d0
            tmpIJ(0, :, :)= DL+BL*Coef_dy(1)-Vyy*Coef_dyy(1)
            tmpIJ(1, :, :)=    BL*Coef_dy(2)-Vyy*Coef_dyy(2)
            tmpIJ(2, :, :)=    BL*Coef_dy(3)-Vyy*Coef_dyy(3)
            do l=0, 2
                print*, tmpIJ(l, :, :)
                call A%set(5, j, l+1, tmpIJ(l, :, :))
            end do
            read(*, *)l
!            do l=1, jn*5
!                call A%Set(l+jn*5, l+jn*5, (1.0d0, 0.0d0))
!            end do
            tmpIJ=0.0d0
            tmpIJ(0, :, :)= BR*Coef_dy(1)+DR1
            tmpIJ(1, :, :)= BR*Coef_dy(2)
            tmpIJ(2, :, :)= BR*Coef_dy(3)
            do l=0, 2
                call B%set(5, j, l+1, tmpIJ(l, :, :))
            end do
            tmpIJ=0.0d0
            tmpIJ(0, :, :)= DR2
            call B%Set(5, j, j+jn, tmpIJ(0, :, :))
!            do l=1, jn*5
!                call B%set(l+jn*5, l, (1.0d0, 0.0d0))
!            end do

            j=2
            !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
            !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
            Coef_dy=this%NormCoef%Coef_dj(:, j)
            Coef_dyy=this%NormCoef%Coef_djj(:, j)
            DisOP=DisOPNorm%GetPoint(j)
            call DisOP%Get(BL, BR, DL, DR1, DR2, Vyy)
            tmpIJ=0.0d0
            tmpIJ(-1, :, :)=   BL*Coef_dy(1)-Vyy*Coef_dyy(1)
            tmpIJ( 0, :, :)=DL+BL*Coef_dy(2)-Vyy*Coef_dyy(2)
            tmpIJ( 1, :, :)=   BL*Coef_dy(3)-Vyy*Coef_dyy(3)
            tmpIJ( 2, :, :)=   BL*Coef_dy(4)-Vyy*Coef_dyy(4)
            do l=-1, 2
                call A%set(5, j, l+j, tmpIJ(l, :, :))
            end do
!            do l=1, jn*5
!                call A%Set(l+jn*5, l+jn*5, (1.0d0, 0.0d0))
!            end do
            tmpIJ=0.0d0
            tmpIJ(-1, :, :)= BR*Coef_dy(1)
            tmpIJ( 0, :, :)= BR*Coef_dy(2)+DR1
            tmpIJ( 1, :, :)= BR*Coef_dy(3)
            tmpIJ( 2, :, :)= BR*Coef_dy(4)
            do l=-1, 2
                call B%set(5, j, l+j, tmpIJ(l, :, :))
            end do
            tmpIJ=0.0d0
            tmpIJ(1, :, :)= DR2
            call B%Set(5, j, j+jn, tmpIJ(1, :, :))
!            do l=1, jn*5
!                call B%set(l+jn*5, l, (1.0d0, 0.0d0))
!            end do

            do j=3, jn-2
              !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
              !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
              Coef_dy=this%NormCoef%Coef_dj(:, j)
              Coef_dyy=this%NormCoef%Coef_djj(:, j)
              DisOP=DisOPNorm%GetPoint(j)
                call DisOP%Get(BL, BR, DL, DR1, DR2, Vyy)
                tmpIJ=0.0d0
                tmpIJ(-2, :, :)=BL*Coef_dy(1)-Vyy*Coef_dyy(1)
                tmpIJ(-1, :, :)=BL*Coef_dy(2)-Vyy*Coef_dyy(2)
                tmpIJ( 0, :, :)=BL*Coef_dy(3)-Vyy*Coef_dyy(3)+DL
                tmpIJ( 1, :, :)=BL*Coef_dy(4)-Vyy*Coef_dyy(4)
                tmpIJ( 2, :, :)=BL*Coef_dy(5)-Vyy*Coef_dyy(5)
                do l=-2, 2
                    call A%set(5, j, j+l, tmpIJ(l, :, :))
                enddo
!                do l=1, jn*5
!                    call A%Set(l+jn*5, l+jn*5, (1.0d0, 0.0d0))
!                end do
                tmpIJ=0.0d0
                tmpIJ(-2, :, :)=BR*Coef_dy(1)
                tmpIJ(-1, :, :)=BR*Coef_dy(2)
                tmpIJ( 0, :, :)=BR*Coef_dy(3)+DR1
                tmpIJ( 1, :, :)=BR*Coef_dy(4)
                tmpIJ( 2, :, :)=BR*Coef_dy(5)
                do l=-2, 2
                    call B%set(5, j, j+l, tmpIJ(l, :, :))
                enddo
                tmpIJ=0.0d0
                tmpIJ(0, :, :)= DR2
                call B%Set(5, j, j+jn, tmpIJ(0, :, :))
!                do l=1, jn*5
!                    call B%set(l+jn*5, l, (1.0d0, 0.0d0))
!                end do
            enddo

            j=jn-1
            !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
            !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
            Coef_dy=this%NormCoef%Coef_dj(:, j)
            Coef_dyy=this%NormCoef%Coef_djj(:, j)
            DisOP=DisOPNorm%GetPoint(j)
            call DisOP%Get(BL, BR, DL, DR1, DR2, Vyy)
            tmpIJ=0.0d0
            tmpIJ(-2, :, :)=BL*Coef_dy(2)-Vyy*Coef_dyy(2)
            tmpIJ(-1, :, :)=BL*Coef_dy(3)-Vyy*Coef_dyy(3)
            tmpIJ( 0, :, :)=BL*Coef_dy(4)-Vyy*Coef_dyy(4)+DL
            tmpIJ( 1, :, :)=BL*Coef_dy(5)-Vyy*Coef_dyy(5)
            do l=-2, 1
                call A%set(5, j, j+l, tmpIJ(l, :, :))
            enddo
!            do l=1, jn*5
!                call A%Set(l+jn*5, l+jn*5, (1.0d0, 0.0d0))
!            end do
            tmpIJ=0.0d0
            tmpIJ(-2, :, :)=BR*Coef_dy(2)
            tmpIJ(-1, :, :)=BR*Coef_dy(3)
            tmpIJ( 0, :, :)=BR*Coef_dy(4)+DR1
            tmpIJ( 1, :, :)=BR*Coef_dy(5)
            do l=-2, 1
                call B%set(5, j, j+l, tmpIJ(l, :, :))
            enddo
            tmpIJ=0.0d0
            tmpIJ(0, :, :)= DR2
            call B%Set(5, j, j+jn, tmpIJ(0, :, :))
!            do l=1, jn*5
!                call B%set(l+jn*5, l, (1.0d0, 0.0d0))
!            end do

            j=jn
            !call DisDiff%GetETACoef(iloc, j, ETA1 , Coef_dy)
            !call DisDiff%GetETACoef(iloc, j, ETA2 , Coef_dyy)
            Coef_dy=this%NormCoef%Coef_dj(:, j)
            Coef_dyy=this%NormCoef%Coef_djj(:, j)
            DisOP=DisOPNorm%GetPoint(j)
            call DisOP%Get(BL, BR, DL, DR1, DR2, Vyy)
            tmpIJ(-2, :, :)=BL*Coef_dy(3)
            tmpIJ(-1, :, :)=BL*Coef_dy(4)
            tmpIJ( 0, :, :)=BL*Coef_dy(5)+DL
            do l=-2, 0
                call A%set(5, j, jn+l, tmpIJ(l, :, :))
            enddo

            tmpIJ(-2, :, :)=BR*Coef_dy(3)
            tmpIJ(-1, :, :)=BR*Coef_dy(4)
            tmpIJ( 0, :, :)=BR*Coef_dy(5)+DR1
            do l=-2, 0
                call B%set(5, j, jn+l, tmpIJ(l, :, :))
            enddo
            tmpIJ(0, :, :)= DR2
            call B%Set(5, j, j+jn, tmpIJ(0, :, :))
            do l=1, jn*5
                call A%Set(l+jn*5, l+jn*5, (1.0d0, 0.0d0))
            end do
            do l=1, jn*5
                call B%set(l+jn*5, l, (1.0d0, 0.0d0))
            end do

          end associate

    end subroutine SetMatrix

end module mod_lst_eqn_IR
