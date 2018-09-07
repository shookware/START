!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_lns_op_point.f90
!> @file
!> @breif LNS方程的基本流在法向一点处相关算子文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  Module: lns_OP_point
!
!  DESCRIPTION:
!> @breif LNS方程的基本流在法向一点处相关算子模块.
!>
!! 相当于LNS方程的算子
!  REVISION HISTORY:
!  2017-08-09 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-08-09
!------------------------------------------------------------------------------
    module mod_lns_OP_point

    use mod_op_mat
    use mod_baseflow_org
    use mod_BF_normal
    use mod_local_coordinate
    use mod_parameter
    use mod_op_mat_cmplx

    implicit none

    private

    !> LNS方程在法向一点处相关算子(PSE用)类.
    !!
    !! 算子定义见下：
    type, public :: lns_OP_point_type

        private

        integer, private :: jloc
        type(op_mat_type), private :: G !< 时间导数系数算子
        type(op_mat_type), private :: A !< 流向一阶导数系数算子
        type(op_mat_type), private :: B !< 法向一阶导数系数算子
        type(op_mat_type), private :: C !< 展向一阶导数系数算子
        type(op_mat_type), private :: D !< 扰动系数算子
        type(op_mat_type), private :: Vxx !< 流向二阶导数系数算子
        type(op_mat_type), private :: Vxy !< 流向法向交叉二阶导数系数算子
        type(op_mat_type), private :: Vxz !< 流向展向交叉二阶导数系数算子
        type(op_mat_type), private :: Vyy !< 法向二阶导数系数算子
        type(op_mat_type), private :: Vyz !< 法向展向二阶导数系数算子
        type(op_mat_type), private :: Vzz !< 展向二阶导数系数算子
        type(local_coordinate_type), private :: Coord !< 该点的局部坐标

    Contains

    generic :: Get => Get_op_mat, getCmplx, getMat !< 获得系数LNS方程系数

!    Procedure :: Set
    procedure :: SetFromFlow => set_from_flow_lns !< 直接从基本流信息设置算子的值
    procedure :: SetCoord => set_coord !<设置局部坐标系
    procedure :: GetCoord !< 获得局部坐标系信息
    procedure :: GetSpt !< 获得谱半径
    procedure, private :: Get_op_mat !< 获得系数算子(实数算子形式)
    procedure, private :: getCmplx !<获得系数算子(复数算子形式)
    procedure, private :: getMat !< 获得系数矩阵(矩阵形式)

    end type lns_OP_point_type

    !> LST用LNS算子
    type, extends(lns_OP_point_type), public :: lst_bf_OP_point_type

    end type lst_bf_OP_point_type

    !> LPSE用LNS算子
    type, extends(lns_OP_point_type), public :: lpse_bf_op_point_type

        type(op_mat_type), private :: DxA !< 流向一阶导数系数算子沿流向的一阶导数(只在伴随方程中用到)
        real(R_P) :: VigneronA(2) !< Vigneron技术处理下流向动量方程的密度和温度导数流向一阶导数系数A(2,1)和A(2,5)

    contains

        generic :: GetDxA => getDxA_opmat_cmplx, getDxA_realmat !< 获得流向一阶导数系数算子沿流向的一阶导数

        procedure, pass(this) :: SetFromFlow => set_from_flow_lpse !< 直接从基本流信息设置算子的值
        procedure :: GetVigneronA => GetVigneronA !<获得Vigneron技术处理下流向动量方程的密度和温度导数流向一阶导数系数
        procedure, private :: getDxA_opmat_cmplx !< 获得\f$\mathrm{d}A/\mathrm{d}x\f$
        procedure, private :: getDxA_realmat !< 获得\f$\mathrm{d}A/\mathrm{d}x\f$(矩阵形式)

    end type lpse_bf_op_point_type

    contains

!    subroutine Set(this, G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz, Coord)
!
!    implicit none
!    real(R_P), dimension(5, 5), intent(in) :: G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
!    type(local_coordinate_type), intent(in) :: Coord
!    class(lns_OP_point_type), intent(inout) :: this
!
!    this%G  = (G  )
!    this%A  = (A  )
!    this%B  = (B  )
!    this%C  = (C  )
!    this%D  = (D  )
!    this%Vxx= (Vxx)
!    this%Vxy= (Vxy)
!    this%Vxz= (Vxz)
!    this%Vyy= (Vyy)
!    this%Vyz= (Vyz)
!    this%Vzz= (Vzz)
!    this%Coord = Coord
!
!    end subroutine Set

    !> 获得局部坐标系
    !! @return 局部坐标系
    function GetCoord(this) result(Coord)

        implicit none
        class(lns_OP_point_type), intent(in) :: this
        type(local_coordinate_type) :: Coord

        Coord=this%Coord

    end function GetCoord

    !> 获得系数算子
    !!@param[out] G     G  算子
    !!@param[out] A     A  算子
    !!@param[out] B     B  算子
    !!@param[out] C     C  算子
    !!@param[out] D     D  算子
    !!@param[out] Vxx   \f$V_{xx}\f$算子
    !!@param[out] Vxy   \f$V_{xy}\f$算子
    !!@param[out] Vxz   \f$V_{xz}\f$算子
    !!@param[out] Vyy   \f$V_{yy}\f$算子
    !!@param[out] Vyz   \f$V_{yz}\f$算子
    !!@param[out] Vzz   \f$V_{zz}\f$算子
    !!@note 算子系数形式具体见@ref mod_lns_op_point模块的详细说明
    subroutine Get_op_mat(this, G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz)

    implicit none
    type(op_mat_type), intent(inout) :: G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
    class(lns_OP_point_type), intent(in) :: this

    G   = this%G
    A   = this%A
    B   = this%B
    C   = this%C
    D   = this%D
    Vxx = this%Vxx
    Vxy = this%Vxy
    Vxz = this%Vxz
    Vyy = this%Vyy
    Vyz = this%Vyz
    Vzz = this%Vzz

    end subroutine Get_op_mat

    !> 获得系数算子(复数形式)
    !!@param[out] G     G  算子
    !!@param[out] A     A  算子
    !!@param[out] B     B  算子
    !!@param[out] C     C  算子
    !!@param[out] D     D  算子
    !!@param[out] Vxx   \f$V_{xx}\f$算子
    !!@param[out] Vxy   \f$V_{xy}\f$算子
    !!@param[out] Vxz   \f$V_{xz}\f$算子
    !!@param[out] Vyy   \f$V_{yy}\f$算子
    !!@param[out] Vyz   \f$V_{yz}\f$算子
    !!@param[out] Vzz   \f$V_{zz}\f$算子
    !!@note 算子系数形式具体见@ref mod_lns_op_point模块的详细说明
    subroutine getCmplx(this, G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz)

    implicit none
    type(op_mat_cmplx_type), intent(inout) :: G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
    class(lns_OP_point_type), intent(in) :: this

    G   = this%G
    A   = this%A
    B   = this%B
    C   = this%C
    D   = this%D
    Vxx = this%Vxx
    Vxy = this%Vxy
    Vxz = this%Vxz
    Vyy = this%Vyy
    Vyz = this%Vyz
    Vzz = this%Vzz

    end subroutine getCmplx

    !> 获得系数算子(矩阵形式)
    !!@param[out] G     G  算子
    !!@param[out] A     A  算子
    !!@param[out] B     B  算子
    !!@param[out] C     C  算子
    !!@param[out] D     D  算子
    !!@param[out] Vxx   \f$V_{xx}\f$算子
    !!@param[out] Vxy   \f$V_{xy}\f$算子
    !!@param[out] Vxz   \f$V_{xz}\f$算子
    !!@param[out] Vyy   \f$V_{yy}\f$算子
    !!@param[out] Vyz   \f$V_{yz}\f$算子
    !!@param[out] Vzz   \f$V_{zz}\f$算子
    !!@note 算子系数形式具体见@ref mod_lns_op_point模块的详细说明
    subroutine getMat(this, G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz)

    implicit none
    real(R_P), dimension(5, 5), intent(inout) :: G, A, B, C, D, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
    class(lns_OP_point_type), intent(in) :: this

    G   = this%G
    A   = this%A
    B   = this%B
    C   = this%C
    D   = this%D
    Vxx = this%Vxx
    Vxy = this%Vxy
    Vxz = this%Vxz
    Vyy = this%Vyy
    Vyz = this%Vyz
    Vzz = this%Vzz

    end subroutine getMat

    !> 设置局部坐标系
    !! @param[in] 局部坐标系
    subroutine set_coord(this, Coord)

    implicit none
    type(local_coordinate_type), intent(in) :: Coord
    class(lns_OP_point_type), intent(inout) :: this

    this%coord=Coord

    end subroutine set_coord

    !> 直接从基本流信息设置算子的值(LST用)
    !!
    !! @param[in] Flowpoint 基本流信息(含导数)
    !! @param[in] isParrallel 是否流动平行
    !! @note 具体表达形式请参照详细注释
    !! @todo 需要写一些公式作为详细注释
    subroutine set_from_flow_lns(this, FlowPoint, isParallel)

    use mod_gas
    use mod_parameter
    use mod_vector

    implicit none
    class(lns_OP_point_type), intent(inout) :: this
    type(BF_point_type), intent(in) :: FlowPoint
    logical, intent(in), optional :: isParallel
    type(BF_flux_org_ij_type) :: Flux, DxFlux, DyFlux, DDyFlux !, DxDyFlux! rho, u, v, w, T
    type(lame_type) :: lame
    type(lame_grad_type) :: LameGrad
    type(basis_type) :: basis
    type(vector_type) :: Vel, Velx, Vely, Velyy, Velxy
    real(R_P) :: Pe, gf, g1, g2
    real(R_P), parameter :: d1d3=1.d0/3.d0, d2d3=2.d0/3.d0, d4d3=4.d0/3.d0, d7d3=7.d0/3.d0
    real(R_P) :: hx, hy, hz
    real(R_P) :: rho, U, V, W, T
    real(R_P) :: rhox, Ux, Vx, Wx, Tx
    real(R_P) :: Uxx, Vxx, Wxx, Txx
    real(R_P) :: rhoy, Uy, Vy, Wy, Ty
    real(R_P) :: rhoyy, Uyy, Vyy, Wyy, Tyy
    real(R_P) :: Uxy, Vxy, Wxy!, Txy, rhoxy,
    real(R_P) :: d12, d32, d31               !没有

    real(R_P) :: G(5, 5), A(5, 5), B(5, 5), C(5, 5), D(5, 5), DxA(5, 5)
    real(R_P) :: MVxx(5, 5), MVxy(5, 5), MVxz(5, 5), MVyy(5, 5), MVyz(5, 5), MVzz(5, 5)
    real(R_P) :: Vig, acoustic

    G=0.0d0; A=0.0d0; B=0.0d0; C=0.0d0; D=0.0d0; DxA=0.0d0
    MVxx=0.0d0; MVxy=0.0d0; MVxz=0.0d0; MVyy=0.0d0; MVyz=0.0d0; MVzz=0.0d0

    Pe=1.d0/(GAMMA*Ma*Ma)
    gf=1.d0/GAMMA
    g1=(1.d0-GAMMA)/GAMMA
    g2=(GAMMA-1.d0)*Ma*Ma


    call this%coord%Get(basis, lame)
    call lame%Get(hx, hy, hz)
    call this%coord%GetGrad(LameGrad)
    call LameGrad%Get(d12, d32, d31)

    call FlowPoint%Get(flux, dxflux, dyflux, ddyflux)
    call Flux%get(rho, Vel, T)
    call Vel%Get(U, V, W)

    call DxFlux%get(rhox, Velx, Tx)
    call Velx%Get(Ux, Vx, Wx)

    call DyFlux%get(rhoy, Vely, Ty)
    call Vely%Get(Uy, Vy, Wy)

    call DDyFlux%get(rhoyy, Velyy, Tyy)
    call Velyy%Get(Uyy, Vyy, Wyy)

    if(.not. present(isParallel) .or. isParallel)then
        Ux=0.0d0; Wx=0.0d0; Rhox=0.0d0; Tx=0.0d0
        V=0.0d0; Vx=0.0d0; Vy=0.0d0; Vyy=0.0d0
    end if

    Uxy=0.0d0; Vxy=0.0d0; Wxy=0.0d0
    Uxx=0.0d0; Vxx=0.0d0; Wxx=0.0d0; Txx=0.0d0

    G(1, 1)=1;G(2, 2)=rho;G(3, 3)=rho;G(4, 4)=rho
    G(5, 1)=(1-GAMMA)/(GAMMA)*T; G(5, 5)=Rho/(GAMMA)

    call this%G%set(G)

    A(1, 1)=U
    A(1, 2)=rho
    A(2, 1)=Pe*T
    A(2, 2)=rho*U-d4d3*(miuT(T)*Tx/hx+miu(T)*d31)/Re
    A(2, 3)=-1.d0/Re*(d7d3*miu(T)*d12+d1d3*miu(T)*d32+miuT(T)*Ty)
    A(2, 5)=rho*Pe-d4d3*miuT(T)*(Ux/hx+V*d12)/Re+d2d3*miuT(T)*(Vy+d31*U+V*d32)/Re
    A(3, 2)=(d7d3*miu(T)*d12+d2d3*miuT(T)*Ty)/Re
    A(3, 3)=rho*U-(miuT(T)*Tx/hx+miu(T)*d31)/Re
    A(3, 5)=-1.d0*miuT(T)/Re*(Vx/hx+Uy-U*d12)
    A(4, 4)=rho*U-(miuT(T)*Tx/hx+miu(T)*d31)/re
    A(4, 5)=-miuT(T)/re*(Wx/hx-W*d31)
    A(5, 1)=g1*T*U
    A(5, 2)=-4.d0*g2*miu(T)/re*(Ux/hx+V*d12-d1d3*(Ux/hx+Vy+U*d31+V*(d12+d32)))
    A(5, 3)=-2.d0*g2*miu(T)/re*(Vx/hx+Uy-U*d12)
    A(5, 4)=-2.d0*g2*miu(T)/re*(Wx/hx-W*d31)
    A(5, 5)=gf*rho*U-2.d0*miuT(T)*Tx/hx/re/pr-miu(T)*d31/re/pr

    A=A/hx

    call this%A%set(A)

    B(1, 1)=V
    B(1, 3)=rho
    B(2, 2)=rho*V-(miu(T)*(d12+d32)+miuT(T)*Ty)/Re
    B(2, 3)=d2d3/hx*miuT(T)*Tx/Re
    B(2, 5)=-miuT(T)/Re*(Vx/hx+Uy-U*d12)
    B(3, 1)=Pe*T
    B(3, 2)=-1.d0*(miuT(T)*Tx/hx+d1d3*miu(T)*d31)/Re
    B(3, 3)=rho*V-d4d3*(miuT(T)*Ty+miu(T)*(d12+d32))/Re
    B(3, 5)=rho*pe+d2d3*miuT(T)/Re*(Ux/hx-2.d0*Vy+U*d31+V*(d12+d32))
    B(4, 4)=rho*V-(miuT(T)*Ty+miu(T)*(d32+d12))/re
    B(4, 5)=-miuT(T)/re*(Wy-W*d32)
    B(5, 1)=g1*T*V
    B(5, 2)=-2.d0*g2*miu(T)/re*(Vx/hx+Uy-U*d12)
    B(5, 3)=-4.d0*g2*miu(T)/re*(Vy-d1d3*(Ux/hx+Vy+U*d31+V*(d12+d32)))
    B(5, 4)=-2.d0*g2*miu(T)/re*(Wy-W*d32)
    B(5, 5)=gf*rho*V-2.d0*miuT(T)*Ty/re/pr-(d12+d32)*miu(T)/re/pr
    call this%B%set(B)

    C(1, 1)=W
    C(1, 4)=rho
    C(2, 2)=rho*W
    C(2, 4)=(d2d3*miuT(T)*Tx/hx+d7d3*miu(T)*d31)/Re
    C(2, 5)=-miuT(T)*(Wx/hx-W*d31)/Re
    C(3, 3)=rho*W
    C(3, 4)=(d2d3*miuT(T)*Ty+d7d3*miu(T)*d32)/Re
    C(3, 5)=-1.d0*miuT(T)*(Wy-W*d32)/Re
    C(4, 1)=T*Pe
    C(4, 2)=-(miuT(T)*Tx/hx+d7d3*miu(T)*d31)/re
    C(4, 3)=-(miuT(T)*Ty+d7d3*miu(T)*d32+d1d3*miu(T)*d12)/re
    C(4, 4)=rho*W
    C(4, 5)=rho*Pe-miuT(T)/re*(-d2d3*Ux/hx-d2d3*Vy+d4d3*d31*U+d4d3*d32*V-d2d3*d12*V)
    C(5, 1)=g1*T*W
    C(5, 2)=-2.d0*g2*miu(T)/re*(Wx/hx-W*d31)
    C(5, 3)=-2.d0*g2*miu(T)/re*(Wy-W*d32)
    C(5, 4)=-4.d0*g2*miu(T)/re*(U*d31+V*d32-d1d3*(Ux/hx+Vy+U*d31+V*(d12+d32)))
    C(5, 5)=gf*rho*W
    C=C/hz
    call this%C%set(C)


    D(1, 1)=Ux/hx+Vy+V*(d12+d32)+U*d31
    D(1, 2)=rhox/hx+rho*d31
    D(1, 3)=rhoy+rho*(d12+d32)
    D(2, 1)=(U*Ux/hx+V*Uy+U*V*d12-W*W*d31)+Tx/hx*Pe
    D(2, 2)=rho*Ux/hx+rho*d12*V+miuT(T)*Ty*d12/Re+d2d3*miuT(T)*Tx*d31/hx/Re
    D(2, 3)=rho*Uy+rho*U*d12-2.d0*d12*miuT(T)*Tx/hx/Re+d2d3*(d12+d32)*miuT(T)*Tx/hx/Re
    D(2, 5)=rhox/hx*Pe+d2d3*miuTT(T)*Tx/hx/re*(Ux/hx+Vy+U*d31+V*(d12+d32))-&
        &miuTT(T)*Ty/Re*(Vx/hx+Uy-U*d12)-2.d0*miuTT(T)*Tx/hx/re*(Ux/hx+V*d12)-&
        &miuT(T)/Re*(d4d3*Uxx/hx/hx+d1d3*Vxy/Hx+Uyy+d4d3*d31*Ux/hx+(d7d3*d12+d1d3*d32)*Vx/hx+&
        &(d12+d32)*Uy)
    D(3, 1)=(U*Vx/hx+V*Vy-U*U*d12-W*W*d32)+Pe*Ty
    D(3, 2)=rho*Vx/hx-2.d0*rho*d12*U+miuT(T)/Re*(Tx*d12/hx+d2d3*Ty*d31)
    D(3, 3)=rho*Vy+d2d3*miuT(T)*Ty/Re*(d12+d32)
    D(3, 4)=-2.d0*rho*W*d32
    D(3, 5)=rhoy*Pe-2.d0*miuT(T)/re*(0.5d0*(Vxx/hx/hx+Uxy/hx-d12*Ux/hx)+Vyy+(d12+d32)*Vy+&
        &0.5d0*d31*(Vx/hx+Uy)-d12*Ux/hx)-miuTT(T)*Tx/hx/re*(Vx/hx+Uy-U*d12)-2.d0*miuTT(T)/re*Ty*Vy+&
        &d2d3/Re*(miuTT(T)*Ty*(Ux/hx+Vy+U*d31+V*(d12+d32))+miuT(T)*(Uxy/hx-d12*Ux/hx+Vyy+d31*Uy+&
        &(d12+d32)*Vy))
    D(4, 1)=U*Wx/hx+V*Wy+W*U*d31+W*V*d32
    D(4, 2)=rho*Wx/hx+rho*d31*W
    D(4, 3)=rho*Wy+rho*d32*W
    D(4, 4)=rho*d31*U+rho*d32*V+(miuT(T)*Tx*d31/hx+miuT(T)*Ty*d32)/re
    D(4, 5)=-miuT(T)/re*(Wxx/hx/hx-d31*Wx/hx+Wyy+(d12+d32)*Wy+2.d0*d31*Wx/hx)-miuTT(T)*Tx/hx/re*(Wx/hx-d31*W)-&
        &miuTT(T)*Ty/re*(Wy-W*d32)
    D(5, 1)=gf*(U*Tx/hx+V*Ty)
    D(5, 2)=gf/hx*rho*Tx+g1*T*rhox/hx+d4d3*g2*miu(T)*d31*(Ux/hx+Vy)/re+2.d0*g2*miu(T)*d12*(Vx/hx+Uy)/re
    D(5, 3)=gf*rho*Ty+g1*T*rhoy-4.d0*g2*d12*miu(T)*Ux/hx/re+d4d3*g2*miu(T)/re*(d12+d32)*(Ux/hx+Vy)
    D(5, 4)=-2.d0*g2*miu(T)/re*(Wx/hx*d31+Wy*d32)
    D(5, 5)=g1*(U*rhox/hx+V*rhoy)-miuT(T)/re/pr*(Txx/hx/hx+Tyy+d12*Ty+d31*Tx/hx+d32*Ty)-miuTT(T)/re/pr*&
        &(Tx*Tx/hx/hx+Ty*Ty)-2.d0*g2*miuT(T)/re*((Ux/hx+V*d12)**2.d0+Vy**2.d0+(U*d31+V*d32)**2.d0+&
        &0.5d0*(Vx/hx+Uy-U*d12)**2.d0+0.5d0*(Wx/hx-W*d31)**2.d0+0.5d0*(Wy-W*d32)**2.d0)+&
        &d2d3*g2*miuT(T)/re*(Ux/hx+Vy+U*d31+V*(d12+d32))**2.d0

    call this%D%set(D)

    MVxx(2, 2)=d4d3*miu(T)/re
    MVxx(3, 3)=1.d0*miu(T)/re
    MVxx(4, 4)=1.d0*miu(T)/re
    MVxx(5, 5)=miu(T)/re/pr
    MVxx=MVxx/hx/hx
    call this%Vxx%set(MVxx)

    MVxy(2, 3)=1.d0*d1d3*miu(T)/re
    MVxy(3, 2)=1.d0*d1d3*miu(T)/re
    MVxy=MVxy/hx
    call this%Vxy%set(MVxy)

    MVxz(2, 4)=1.d0*d1d3*miu(T)/re
    MVxz(4, 2)=1.d0*d1d3*miu(T)/re
    MVxz=MVxz/hx/hz
    call this%Vxz%set(MVxz)

    MVyy(2, 2)=1.d0*miu(T)/re
    MVyy(3, 3)=d4d3*miu(T)/re
    MVyy(4, 4)=1.d0*miu(T)/re
    MVyy(5, 5)=miu(T)/re/pr
    call this%Vyy%set(MVyy)

    MVyz(3, 4)=d1d3*miu(T)/re
    MVyz(4, 3)=d1d3*miu(T)/re
    MVyz=MVyz/hz
    call this%Vyz%set(MVyz)

    MVzz(2, 2)=1.d0*miu(T)/re
    MVzz(3, 3)=1.d0*miu(T)/re
    MVzz(4, 4)=d4d3*miu(T)/re
    MVzz(5, 5)=miu(T)/re/pr
    MVzz=MVzz/hz/hz
    call this%Vzz%set(MVzz)

    end subroutine set_from_flow_lns

    !> 直接从基本流信息设置算子的值(LPSE用)
    !!
    !! @param[in] Flowpoint 基本流信息(含导数)
    !! @param[in] isParrallel 是否流动平行
    !! @note 具体表达形式请参照详细注释
    !! @todo 需要写一些公式作为详细注释
    subroutine set_from_flow_lpse(this, FlowPoint, isParallel)

    use mod_gas
    use mod_parameter
    use mod_vector

    implicit none
    class(lpse_BF_OP_point_type), intent(inout) :: this
    type(BF_point_type), intent(in) :: FlowPoint
    logical, intent(in), optional :: isParallel
    type(BF_flux_org_ij_type) :: Flux, DxFlux, DyFlux, DDyFlux !, DxDyFlux! rho, u, v, w, T
    type(lame_type) :: lame
    type(lame_grad_type) :: LameGrad
    type(basis_type) :: basis
    type(vector_type) :: Vel, Velx, Vely, Velyy, Velxy
    real(R_P) :: Pe, gf, g1, g2
    real(R_P), parameter :: d1d3=1.d0/3.d0, d2d3=2.d0/3.d0, d4d3=4.d0/3.d0, d7d3=7.d0/3.d0
    real(R_P) :: hx, hy, hz
    real(R_P) :: rho, U, V, W, T
    real(R_P) :: rhox, Ux, Vx, Wx, Tx
    real(R_P) :: Uxx, Vxx, Wxx, Txx
    real(R_P) :: rhoy, Uy, Vy, Wy, Ty
    real(R_P) :: rhoyy, Uyy, Vyy, Wyy, Tyy
    real(R_P) :: Uxy, Vxy, Wxy!, Txy, rhoxy,
    real(R_P) :: d12, d32, d31               !没有

    real(R_P) :: G(5, 5), A(5, 5), B(5, 5), C(5, 5), D(5, 5), DxA(5, 5)
    real(R_P) :: MVxx(5, 5), MVxy(5, 5), MVxz(5, 5), MVyy(5, 5), MVyz(5, 5), MVzz(5, 5)
    real(R_P) :: Vig, acoustic

    G=0.0d0; A=0.0d0; B=0.0d0; C=0.0d0; D=0.0d0; DxA=0.0d0
    MVxx=0.0d0; MVxy=0.0d0; MVxz=0.0d0; MVyy=0.0d0; MVyz=0.0d0; MVzz=0.0d0

    Pe=1.d0/(GAMMA*Ma*Ma)
    gf=1.d0/GAMMA
    g1=(1.d0-GAMMA)/GAMMA
    g2=(GAMMA-1.d0)*Ma*Ma

    call FlowPoint%Get(flux, dxflux, dyflux, ddyflux)

    call this%coord%Get(basis, lame)
    call lame%Get(hx, hy, hz)
    call this%coord%GetGrad(LameGrad)
    call LameGrad%Get(d12, d32, d31)

    !!!!!

    call Flux%get(rho, Vel, T)
    call Vel%Get(U, V, W)

    call DxFlux%get(rhox, Velx, Tx)
    call Velx%Get(Ux, Vx, Wx)

    call DyFlux%get(rhoy, Vely, Ty)
    call Vely%Get(Uy, Vy, Wy)

    call DDyFlux%get(rhoyy, Velyy, Tyy)
    call Velyy%Get(Uyy, Vyy, Wyy)

    if(.not. present(isParallel) .or. isParallel)then
        Ux=0.0d0; Wx=0.0d0; rhox=0.0d0; Tx=0.0d0
        V=0.0d0; Vx=0.0d0; Vy=0.0d0; Vyy=0.0d0
    end if

    !call DxDyFlux%get(rhoxy, Velxy, Txy)
    !call Velxy%Get(Uxy, Vxy, Wxy)
    Uxy=0.0d0; Vxy=0.0d0; Wxy=0.0d0
    Uxx=0.0d0; Vxx=0.0d0; Wxx=0.0d0; Txx=0.0d0

    !! debug use

!    v=0.0d0; vy=0.0d0; vyy=0.0d0; vx=0.0d0

    G(1, 1)=1;G(2, 2)=rho;G(3, 3)=rho;G(4, 4)=rho
    G(5, 1)=(1-GAMMA)/(GAMMA)*T; G(5, 5)=Rho/(GAMMA)

    call this%G%set(G)

    A(1, 1)=U
    A(1, 2)=rho
    A(2, 1)=Pe*T
    A(2, 2)=rho*U-d4d3*(miuT(T)*Tx/hx+miu(T)*d31)/Re
    A(2, 3)=-1.d0/Re*(d7d3*miu(T)*d12+d1d3*miu(T)*d32+miuT(T)*Ty)
    A(2, 5)=rho*Pe-d4d3*miuT(T)*(Ux/hx+V*d12)/Re+d2d3*miuT(T)*(Vy+d31*U+V*d32)/Re
    A(3, 2)=(d7d3*miu(T)*d12+d2d3*miuT(T)*Ty)/Re
    A(3, 3)=rho*U-(miuT(T)*Tx/hx+miu(T)*d31)/Re
    A(3, 5)=-1.d0*miuT(T)/Re*(Vx/hx+Uy-U*d12)
    A(4, 4)=rho*U-(miuT(T)*Tx/hx+miu(T)*d31)/re
    A(4, 5)=-miuT(T)/re*(Wx/hx-W*d31)
    A(5, 1)=g1*T*U
    A(5, 2)=-4.d0*g2*miu(T)/re*(Ux/hx+V*d12-d1d3*(Ux/hx+Vy+U*d31+V*(d12+d32)))
    A(5, 3)=-2.d0*g2*miu(T)/re*(Vx/hx+Uy-U*d12)
    A(5, 4)=-2.d0*g2*miu(T)/re*(Wx/hx-W*d31)
    A(5, 5)=gf*rho*U-2.d0*miuT(T)*Tx/hx/re/pr-miu(T)*d31/re/pr

    acoustic=U/sqrt(T)*Ma

    A=A/hx

    Vig=min(1.0d0, acoustic **2)
    this%VigneronA(1)=pe*T/hx*(1.0d0-Vig)
    this%VigneronA(2)=pe*rho/hx*(1.0d0-Vig)

    call this%A%set(A)

    B(1, 1)=V
    B(1, 3)=rho
    B(2, 2)=rho*V-(miu(T)*(d12+d32)+miuT(T)*Ty)/Re
    B(2, 3)=d2d3/hx*miuT(T)*Tx/Re
    B(2, 5)=-miuT(T)/Re*(Vx/hx+Uy-U*d12)
    B(3, 1)=Pe*T
    B(3, 2)=-1.d0*(miuT(T)*Tx/hx+d1d3*miu(T)*d31)/Re
    B(3, 3)=rho*V-d4d3*(miuT(T)*Ty+miu(T)*(d12+d32))/Re
    B(3, 5)=rho*pe+d2d3*miuT(T)/Re*(Ux/hx-2.d0*Vy+U*d31+V*(d12+d32))
    B(4, 4)=rho*V-(miuT(T)*Ty+miu(T)*(d32+d12))/re
    B(4, 5)=-miuT(T)/re*(Wy-W*d32)
    B(5, 1)=g1*T*V
    B(5, 2)=-2.d0*g2*miu(T)/re*(Vx/hx+Uy-U*d12)
    B(5, 3)=-4.d0*g2*miu(T)/re*(Vy-d1d3*(Ux/hx+Vy+U*d31+V*(d12+d32)))
    B(5, 4)=-2.d0*g2*miu(T)/re*(Wy-W*d32)
    B(5, 5)=gf*rho*V-2.d0*miuT(T)*Ty/re/pr-(d12+d32)*miu(T)/re/pr
    call this%B%set(B)

    C(1, 1)=W
    C(1, 4)=rho
    C(2, 2)=rho*W
    C(2, 4)=(d2d3*miuT(T)*Tx/hx+d7d3*miu(T)*d31)/Re
    C(2, 5)=-miuT(T)*(Wx/hx-W*d31)/Re
    C(3, 3)=rho*W
    C(3, 4)=(d2d3*miuT(T)*Ty+d7d3*miu(T)*d32)/Re
    C(3, 5)=-1.d0*miuT(T)*(Wy-W*d32)/Re
    C(4, 1)=T*Pe
    C(4, 2)=-(miuT(T)*Tx/hx+d7d3*miu(T)*d31)/re
    C(4, 3)=-(miuT(T)*Ty+d7d3*miu(T)*d32+d1d3*miu(T)*d12)/re
    C(4, 4)=rho*W
    C(4, 5)=rho*Pe-miuT(T)/re*(-d2d3*Ux/hx-d2d3*Vy+d4d3*d31*U+d4d3*d32*V-d2d3*d12*V)
    C(5, 1)=g1*T*W
    C(5, 2)=-2.d0*g2*miu(T)/re*(Wx/hx-W*d31)
    C(5, 3)=-2.d0*g2*miu(T)/re*(Wy-W*d32)
    C(5, 4)=-4.d0*g2*miu(T)/re*(U*d31+V*d32-d1d3*(Ux/hx+Vy+U*d31+V*(d12+d32)))
    C(5, 5)=gf*rho*W
    C=C/hz
    call this%C%set(C)


    D(1, 1)=Ux/hx+Vy+V*(d12+d32)+U*d31
    D(1, 2)=rhox/hx+rho*d31
    D(1, 3)=rhoy+rho*(d12+d32)
    D(2, 1)=(U*Ux/hx+V*Uy+U*V*d12-W*W*d31)+Tx/hx*Pe
    D(2, 2)=rho*Ux/hx+rho*d12*V+miuT(T)*Ty*d12/Re+d2d3*miuT(T)*Tx*d31/hx/Re
    D(2, 3)=rho*Uy+rho*U*d12-2.d0*d12*miuT(T)*Tx/hx/Re+d2d3*(d12+d32)*miuT(T)*Tx/hx/Re
    D(2, 5)=rhox/hx*Pe+d2d3*miuTT(T)*Tx/hx/re*(Ux/hx+Vy+U*d31+V*(d12+d32))-&
        &miuTT(T)*Ty/Re*(Vx/hx+Uy-U*d12)-2.d0*miuTT(T)*Tx/hx/re*(Ux/hx+V*d12)-&
        &miuT(T)/Re*(d4d3*Uxx/hx/hx+d1d3*Vxy/Hx+Uyy+d4d3*d31*Ux/hx+(d7d3*d12+d1d3*d32)*Vx/hx+&
        &(d12+d32)*Uy)
    D(3, 1)=(U*Vx/hx+V*Vy-U*U*d12-W*W*d32)+Pe*Ty
    D(3, 2)=rho*Vx/hx-2.d0*rho*d12*U+miuT(T)/Re*(Tx*d12/hx+d2d3*Ty*d31)
    D(3, 3)=rho*Vy+d2d3*miuT(T)*Ty/Re*(d12+d32)
    D(3, 4)=-2.d0*rho*W*d32
    D(3, 5)=rhoy*Pe-2.d0*miuT(T)/re*(0.5d0*(Vxx/hx/hx+Uxy/hx-d12*Ux/hx)+Vyy+(d12+d32)*Vy+&
        &0.5d0*d31*(Vx/hx+Uy)-d12*Ux/hx)-miuTT(T)*Tx/hx/re*(Vx/hx+Uy-U*d12)-2.d0*miuTT(T)/re*Ty*Vy+&
        &d2d3/Re*(miuTT(T)*Ty*(Ux/hx+Vy+U*d31+V*(d12+d32))+miuT(T)*(Uxy/hx-d12*Ux/hx+Vyy+d31*Uy+&
        &(d12+d32)*Vy))
    D(4, 1)=U*Wx/hx+V*Wy+W*U*d31+W*V*d32
    D(4, 2)=rho*Wx/hx+rho*d31*W
    D(4, 3)=rho*Wy+rho*d32*W
    D(4, 4)=rho*d31*U+rho*d32*V+(miuT(T)*Tx*d31/hx+miuT(T)*Ty*d32)/re
    D(4, 5)=-miuT(T)/re*(Wxx/hx/hx-d31*Wx/hx+Wyy+(d12+d32)*Wy+2.d0*d31*Wx/hx)-miuTT(T)*Tx/hx/re*(Wx/hx-d31*W)-&
        &miuTT(T)*Ty/re*(Wy-W*d32)
    D(5, 1)=gf*(U*Tx/hx+V*Ty)
    D(5, 2)=gf/hx*rho*Tx+g1*T*rhox/hx+d4d3*g2*miu(T)*d31*(Ux/hx+Vy)/re+2.d0*g2*miu(T)*d12*(Vx/hx+Uy)/re
    D(5, 3)=gf*rho*Ty+g1*T*rhoy-4.d0*g2*d12*miu(T)*Ux/hx/re+d4d3*g2*miu(T)/re*(d12+d32)*(Ux/hx+Vy)
    D(5, 4)=-2.d0*g2*miu(T)/re*(Wx/hx*d31+Wy*d32)
    D(5, 5)=g1*(U*rhox/hx+V*rhoy)-miuT(T)/re/pr*(Txx/hx/hx+Tyy+d12*Ty+d31*Tx/hx+d32*Ty)-miuTT(T)/re/pr*&
        &(Tx*Tx/hx/hx+Ty*Ty)-2.d0*g2*miuT(T)/re*((Ux/hx+V*d12)**2.d0+Vy**2.d0+(U*d31+V*d32)**2.d0+&
        &0.5d0*(Vx/hx+Uy-U*d12)**2.d0+0.5d0*(Wx/hx-W*d31)**2.d0+0.5d0*(Wy-W*d32)**2.d0)+&
        &d2d3*g2*miuT(T)/re*(Ux/hx+Vy+U*d31+V*(d12+d32))**2.d0
    D(2, 1)=D(2, 1)-pe*Tx/hx*(1.0d0-Vig)
    D(2, 5)=D(2, 5)-pe*rhox/hx*(1.0d0-Vig)

    call this%D%set(D)

    MVxx(2, 2)=d4d3*miu(T)/re
    MVxx(3, 3)=1.d0*miu(T)/re
    MVxx(4, 4)=1.d0*miu(T)/re
    MVxx(5, 5)=miu(T)/re/pr
    MVxx=MVxx/hx/hx
    call this%Vxx%set(MVxx)

    MVxy(2, 3)=1.d0*d1d3*miu(T)/re
    MVxy(3, 2)=1.d0*d1d3*miu(T)/re
    MVxy=MVxy/hx
    call this%Vxy%set(MVxy)

    MVxz(2, 4)=1.d0*d1d3*miu(T)/re
    MVxz(4, 2)=1.d0*d1d3*miu(T)/re
    MVxz=MVxz/hx/hz
    call this%Vxz%set(MVxz)

    MVyy(2, 2)=1.d0*miu(T)/re
    MVyy(3, 3)=d4d3*miu(T)/re
    MVyy(4, 4)=1.d0*miu(T)/re
    MVyy(5, 5)=miu(T)/re/pr
    call this%Vyy%set(MVyy)

    MVyz(3, 4)=d1d3*miu(T)/re
    MVyz(4, 3)=d1d3*miu(T)/re
    MVyz=MVyz/hz
    call this%Vyz%set(MVyz)

    MVzz(2, 2)=1.d0*miu(T)/re
    MVzz(3, 3)=1.d0*miu(T)/re
    MVzz(4, 4)=d4d3*miu(T)/re
    MVzz(5, 5)=miu(T)/re/pr
    MVzz=MVzz/hz/hz
    call this%Vzz%set(MVzz)

    DxA(1, 1)=Ux
    DxA(1, 2)=rhox
    DxA(2, 1)=0.0d0
    DxA(2, 2)=rhox*U+rho*Ux
    DxA(2, 3)=0.0d0
    DxA(2, 5)=0.0d0
    DxA(3, 2)=0.0d0
    DxA(3, 3)=rhox*U+rho*Ux
    DxA(3, 5)=0.0d0
    DxA(4, 4)=rhox*U+rho*Ux
    DxA(4, 5)=0.0d0
    DxA(5, 1)=g1*Tx*U+g1*T*Ux
    DxA(5, 2)=0.0d0
    DxA(5, 3)=0.0d0
    DxA(5, 4)=0.0d0
    DxA(5, 5)=gf*rhox*U+gf*rho*Ux

    DxA=DxA/hx
    call this%DxA%set(DxA)


            !write(ENO_BF, *)'matA'
            !call this%A%Print()
            !write(ENO_BF, *)'matB'
            !call this%B%Print()
            !write(ENO_BF, *)'matC'
            !call this%C%Print()
            !write(ENO_BF, *)'matD'
            !call this%D%Print()
            !write(ENO_BF, *)'matG'
            !call this%G%Print()

    end subroutine set_from_flow_lpse

    !> 直接从基本流信息设置算子的值(LPSE用) !内能方程
    !!
    !! @param[in] Flowpoint 基本流信息(含导数)
    !! @param[in] isParrallel 是否流动平行
    !! @note 具体表达形式请参照详细注释
    !! @todo 需要写一些公式作为详细注释
    subroutine set_from_flow_lpse_ie(this, FlowPoint, isParallel)

    use mod_gas
    use mod_parameter
    use mod_vector

    implicit none
    class(lpse_BF_OP_point_type), intent(inout) :: this
    type(BF_point_type), intent(in) :: FlowPoint
    logical, intent(in), optional :: isParallel
    type(BF_flux_org_ij_type) :: Flux, DxFlux, DyFlux, DDyFlux !, DxDyFlux! rho, u, v, w, T
    type(BF_flux_org_ij_type) :: DzFlux, DzzFlux
    type(lame_type) :: lame
    type(lame_grad_type) :: LameGrad
    type(basis_type) :: basis
    type(vector_type) :: Vel, Velx, Vely, Velyy, Velxy, Velz, Velzz
    real(R_P) :: Pe, gf, g1, g2
    real(R_P), parameter :: d1d3=1.d0/3.d0, d2d3=2.d0/3.d0, d4d3=4.d0/3.d0, d7d3=7.d0/3.d0
    real(R_P) :: hx, hy, hz
    real(R_P) :: rho, U, V, W, T
    real(R_P) :: rhox, Ux, Vx, Wx, Tx
    real(R_P) :: Uxx, Vxx, Wxx, Txx
    real(R_P) :: rhoy, Uy, Vy, Wy, Ty
    real(R_P) :: rhoyy, Uyy, Vyy, Wyy, Tyy
    real(R_P) :: Uxy, Vxy, Wxy!, Txy, rhoxy,
    real(R_P) :: d12, d32, d31               !没有
    real(R_P) :: tmp

    real(R_P) :: G(5, 5), A(5, 5), B(5, 5), C(5, 5), D(5, 5), DxA(5, 5)
    real(R_P) :: MVxx(5, 5), MVxy(5, 5), MVxz(5, 5), MVyy(5, 5), MVyz(5, 5), MVzz(5, 5)
    real(R_P) :: Vig, acoustic

    real(R_P) :: e11, e22, e33, e12, e13, e21, e23, e31, e32
    real(R_P) :: p11, p12, p13, p22, p23, p33
    real(R_P) :: divu
    real(R_P) :: ugradux, ugraduy, ugraduz
    real(R_P) :: divpx, divpy, divpz
    real(R_P) :: divgradux, divgraduy, divgraduz
    real(R_P) :: gradxdivu, gradydivu, gradzdivu
    real(R_P) :: rhoz, Uz, Vz, Wz, Tz, Tzz, Uzz, Vzz, Wzz
    real(R_P) :: Uxz, Vyz, Wxz, Wyz
    real(R_P) :: cv0, cv0t, cp0
    real(R_P) :: divgradkai, kaixx, kaiyy, kaizz
    real(R_P) :: phi


    G=0.0d0; A=0.0d0; B=0.0d0; C=0.0d0; D=0.0d0; DxA=0.0d0
    MVxx=0.0d0; MVxy=0.0d0; MVxz=0.0d0; MVyy=0.0d0; MVyz=0.0d0; MVzz=0.0d0

    Pe=1.d0/(GAMMA*Ma*Ma)
    gf=1.d0/GAMMA
    g1=(1.d0-GAMMA)/GAMMA
    g2=(GAMMA-1.d0)*Ma*Ma

!    call FlowPoint%Get(flux, dxflux, dyflux, ddyflux, dzFlux, DzzFlux)
    call FlowPoint%Get(flux, dxflux, dyflux, ddyflux)

    call this%coord%Get(basis, lame)
    call lame%Get(hx, hy, hz)
    call this%coord%GetGrad(LameGrad)
    call LameGrad%Get(d12, d32, d31)

    !!!!!

    call Flux%get(rho, Vel, T)
    call Vel%Get(U, V, W)

    call DxFlux%get(rhox, Velx, Tx)
    call Velx%Get(Ux, Vx, Wx)

    call DyFlux%get(rhoy, Vely, Ty)
    call Vely%Get(Uy, Vy, Wy)

    call DDyFlux%get(rhoyy, Velyy, Tyy)
    call Velyy%Get(Uyy, Vyy, Wyy)

!    call DzFlux%get(rhoz, Velz, Tz)
!    call Velz%Get(Uz, Vz, Wz)

    rhoz=0.0d0; Tz=0.0d0
    Uz=0.0d0; Vz=0.0d0; Wz=0.0d0
!    call DzzFlux%get(tmp, Velzz, Tzz)
!    call Velzz%Get(Uzz, Vzz, Wzz)
    Tzz=0.0d0
    Uzz=0.0d0; Vzz=0.0d0; Wzz=0.0d0

    if(.not. present(isParallel) .or. isParallel)then
        Ux=0.0d0; Wx=0.0d0; rhox=0.0d0; Tx=0.0d0
        V=0.0d0; Vx=0.0d0; Vy=0.0d0; Vyy=0.0d0
    end if

    !call DxDyFlux%get(rhoxy, Velxy, Txy)
    !call Velxy%Get(Uxy, Vxy, Wxy)
    Uxy=0.0d0; Vxy=0.0d0; Wxy=0.0d0
    Uxx=0.0d0; Vxx=0.0d0; Wxx=0.0d0; Txx=0.0d0
    Uxz=0.0d0; Vyz=0.0d0; Wxz=0.0d0; Wyz=0.0d0

    !! local gas parameteter
    CV0=1.0d0/(Gamma*(GAMMA-1.d0)*Ma*Ma) !< 定容比热
    CV0T=0.0d0 !< 定容比热随温度导数
    Cp0=GAMMA*CV0 !< 定压比热


    !! other variables (stress etc.)
    e11=Ux/hx+V*d12
    e22=Vy
    e33=Wz/hz+d32*v
    e12=Vx/hx-U*d12
    e13=wx/hx
    e21=uy
    e23=wy
    e31=uz/hz
    e32=vz/hz-d32*w

    divu=e11+e22+e33

    p11=e11+e11-d2d3*e11
    p12=e12+e21
    p13=e13+e31
    p22=e22+e22-d2d3*e22
    p23=e23+e32
    p33=e33+e33-d2d3*e33

    ugradux=u*e11+V*e21+w*e31
    ugraduy=u*e12+V*e22+w*e32
    ugraduz=u*e13+V*e23+w*e33

    divgradux=uxx/hx**2+uyy+uzz/hz**2+(d32+d12)*uy+2.0d0*d12*vx/hx
    divgraduy=Vxx/hx**2+Vyy+Vzz/hz**2+(d32+d12)*Vy-2.0d0*d12*ux/hx-2.0d0*d32*wz/hz-v*(d12**2+d32**2)
    divgraduz=wxx/hx**2+wyy+wzz/hz**2+(d32+d12)*wy+2.0d0*d32*vz/hz-w*(d32**2)

    gradxdivu=uxx/hx**2+vxy/hx+wxz/hx/hz-d12*v+(d12+d32)*vx/hx
    gradydivu=Uxy/hx+Vyy+Wyz/hz+(d12+d32)*Vy
    gradzdivu=Uxz/hx/hz+Vyz/hz+Wzz/hz**2-d32*Vy+(d12+d32)*Vz/hz

    divpx=divgradux+d1d3*gradxdivu
    divpy=divgraduy+d1d3*gradydivu
    divpz=divgraduz+d1d3*gradzdivu

    kaixx=(miuT(T)*Txx+miuTT(T)*Tx**2)*cp0/RE/PR
    kaiyy=(miuT(T)*Tyy+miuTT(T)*Ty**2)*cp0/RE/PR
    kaizz=(miuT(T)*Tzz+miuTT(T)*Tz**2)*cp0/RE/PR
    divgradkai=kaixx/hx**2+kaiyy+kaizz/hz**2+(d12+d32)*miuT(T)*Ty*Cp0/RE/PR

    phi=2.0d0*(p11*e11+p22*e22+p33*e33+p12*(e12+e21)+p13*(e13+e31)+p23*(e23+e32))

    !! debug use

!    v=0.0d0; vy=0.0d0; vyy=0.0d0; vx=0.0d0

    G(1, 1)=1;G(2, 2)=rho;G(3, 3)=rho;G(4, 4)=rho
    G(5, 5)=Rho*Cv0

    call this%G%set(G)

    A(1, 1)=U
    A(1, 2)=rho
    A(2, 1)=Pe*T
    A(2, 2)=rho*U-d4d3*(miuT(T)*Tx/hx)/Re
    A(2, 3)=-1.d0/Re*(d7d3*miu(T)*d12+d1d3*miu(T)*d32+miuT(T)*Ty)
    A(2, 4)=-1.d0/Re*miuT(T)*Tz/hz
    A(2, 5)=rho*Pe-miuT(T)/Re*p11
    A(3, 2)=(d7d3*miu(T)*d12+d2d3*miuT(T)*Ty)/Re
    A(3, 3)=rho*U-(miuT(T)*Tx/hx)/Re
    A(3, 5)=-1.d0*miuT(T)/Re*(p12)
    A(4, 2)=d2d3*miuT(T)/Re*Tz/hz
    A(4, 4)=rho*U-(miuT(T)*Tx/hx)/re
    A(4, 5)=-miuT(T)/re*(p13)
    A(5, 2)=pe*rho*T-2.0d0*miu(T)/Re*p11
    A(5, 3)=-2.0d0*miu(T)/Re*p12
    A(5, 4)=-2.0d0*miu(T)/Re*p13
    A(5, 5)=rho*U*cv0-2.d0*miuT(T)*Tx/hx*cp0/Re/Pr

    acoustic=U/sqrt(T)

    Vig=min(1.0d0, acoustic **2)
    this%VigneronA(1)=pe*T/hx*(1.0d0-Vig)
    this%VigneronA(2)=pe*rho/hx*(1.0d0-Vig)

    A=A/hx

    call this%A%set(A)

    B(1, 1)=V
    B(1, 3)=rho
    B(2, 2)=rho*V-(miu(T)*(d12+d32)+miuT(T)*Ty)/Re
    B(2, 3)=d2d3/hx*miuT(T)*Tx/Re
    B(2, 5)=-miuT(T)/Re*(p12)
    B(3, 1)=Pe*T
    B(3, 2)=-1.d0*(miuT(T)*Tx/hx)/Re
    B(3, 3)=rho*V-d4d3*(miuT(T)*Ty+miu(T)*(d12+d32))/Re
    B(3, 4)=-miuT(T)/Re*Tz/hz
    B(3, 5)=rho*pe-miuT(T)/Re*p22
    B(4, 3)=d2d3*miuT(T)*Tz/hz/Re
    B(4, 4)=rho*V-(miuT(T)*Ty+miu(T)*(d32+d12))/re
    B(4, 5)=-miuT(T)/re*(p23)
    B(5, 1)=0.0d0
    B(5, 2)=-2.0d0*miu(T)/Re*p12
    B(5, 3)=pe*rho*T-2.0d0*miu(T)/Re*p22
    B(5, 4)=-2.0d0*miu(T)/Re*p23
    B(5, 5)=rho*V*cv0-(2.d0*miuT(T)*Ty+miu(T)*(d12+d32))*cp0/Re/Pr
    call this%B%set(B)

    C(1, 1)=W
    C(1, 4)=rho
    C(2, 2)=rho*W-miuT(T)/Re*Tz/hz
    C(2, 4)=(d2d3*miuT(T)*Tx/hx)/Re
    C(2, 5)=-miuT(T)*(p13)/Re
    C(3, 3)=rho*W-miuT(T)*Tz/hz/Re
    C(3, 4)=(d2d3*miuT(T)*Ty+d7d3*miu(T)*d32)/Re
    C(3, 5)=-1.d0*miuT(T)*(p23)/Re
    C(4, 1)=T*Pe
    C(4, 2)=-(miuT(T)*Tx/hx)/re
    C(4, 3)=-(miuT(T)*Ty+d7d3*miu(T)*d32+d1d3*miu(T)*d12)/re
    C(4, 4)=rho*W-d4d3*miuT(T)*Tz/hz/Re
    C(4, 5)=rho*Pe-miuT(T)/re*(p33)
    C(5, 1)=0.0d0
    C(5, 2)=-2.0d0*miu(T)/Re*p13
    C(5, 3)=-2.0d0*miu(T)/Re*p23
    C(5, 4)=-2.0d0*miu(T)/Re*p33+pe*rho*T
    C(5, 5)=rho*W*cv0-2.d0*miuT(T)*Tz/hz*cp0/Re/Pr
    C=C/hz
    call this%C%set(C)


    D(1, 1)=divu
    D(1, 2)=rhox/hx
    D(1, 3)=rhoy+rho*(d12+d32)
    D(1, 4)=rhoz/hz
    D(2, 1)=(ugradux)+Tx/hx*Pe
    D(2, 2)=rho*e11+miuT(T)*Ty*d12/Re+miu(T)/Re*d12**2
    D(2, 3)=rho*(e21+U*d12)-d2d3*(2.0d0*d12-d32)*miuT(T)*Tx/hx/Re
    D(2, 4)=rho*e31
    D(2, 5)=rhox/hx*Pe-(miuTT(T)*Tx/hx*p11+miuTT(T)*Ty*p12+miuTT(T)*Tz/hz*p13)/Re-miuT(T)/Re*divpx
    D(3, 1)=(ugraduy)+Pe*Ty
    D(3, 2)=rho*(e12-U*d12)+miuT(T)/Re*(Tx*d12/hx)
    D(3, 3)=rho*e22-d2d3*miuT(T)*Ty/Re*(d12+d32)+miu(T)/Re*(d12**2+d32**2)
    D(3, 4)=rho*(e32-W*d32)+(d32*miuT(T)*Tz/hz/Re)
    D(3, 5)=rhoy*Pe-(miuTT(T)*Tx/hx*p12+miuTT(T)*Ty*p22+miuTT(T)*Tz/hz*p23)/Re-miuT(T)/Re*divpy
    D(4, 1)=(ugraduz)+Pe*Tz/hz
    D(4, 2)=rho*(e13)

    D(4, 3)=rho*(e23+w*d32)-d2d3*miuT(T)*Tz/hz/Re*(2.0d0*d32-d12)
    D(4, 4)=rho*(e33)+(miuT(T)*Ty*d32+miu(T)*d32**2)/re
    D(4, 5)=rhoz/hz*Pe-(miuTT(T)*Tx/hx*p13+miuTT(T)*Ty*p23+miuTT(T)*Tz/hz*p33)/Re-miuT(T)/Re*divpz
    D(5, 1)=cv0*(u*Tx/hx+v*Ty+w*Tz/hz)+pe*T*divu
    D(5, 2)=rho*cv0*Tx/hx+2.0d0*miu(T)/Re*d12*p12
    D(5, 3)=rho*cv0*Ty+pe*rho*T*(d12+d32)-2.d0*miu(T)/Re*(p11*d12+p33*d32)
    D(5, 4)=rho*cv0*Tz/hz+2.0d0*miu(T)/Re*d32*p23
    D(5, 5)=rho*cv0T*(u*Tx/hx+v*Ty+w*Tz/hz)+pe*rho*divu-divgradkai-miuT(T)/Re*phi

    D(2, 1)=D(2, 1)-pe*Tx/hx*(1.0d0-Vig)
    D(2, 5)=D(2, 5)-pe*rhox/hx*(1.0d0-Vig)

    call this%D%set(D)

    MVxx(2, 2)=d4d3*miu(T)/re
    MVxx(3, 3)=1.d0*miu(T)/re
    MVxx(4, 4)=1.d0*miu(T)/re
    MVxx(5, 5)=miu(T)*cp0/re/pr
    MVxx=MVxx/hx/hx
    call this%Vxx%set(MVxx)

    MVxy(2, 3)=1.d0*d1d3*miu(T)/re
    MVxy(3, 2)=1.d0*d1d3*miu(T)/re
    MVxy=MVxy/hx
    call this%Vxy%set(MVxy)

    MVxz(2, 4)=1.d0*d1d3*miu(T)/re
    MVxz(4, 2)=1.d0*d1d3*miu(T)/re
    MVxz=MVxz/hx/hz
    call this%Vxz%set(MVxz)

    MVyy(2, 2)=1.d0*miu(T)/re
    MVyy(3, 3)=d4d3*miu(T)/re
    MVyy(4, 4)=1.d0*miu(T)/re
    MVyy(5, 5)=miu(T)*cp0/re/pr
    call this%Vyy%set(MVyy)

    MVyz(3, 4)=d1d3*miu(T)/re
    MVyz(4, 3)=d1d3*miu(T)/re
    MVyz=MVyz/hz
    call this%Vyz%set(MVyz)

    MVzz(2, 2)=1.d0*miu(T)/re
    MVzz(3, 3)=1.d0*miu(T)/re
    MVzz(4, 4)=d4d3*miu(T)/re
    MVzz(5, 5)=miu(T)*cp0/re/pr
    MVzz=MVzz/hz/hz
    call this%Vzz%set(MVzz)

    DxA(1, 1)=Ux
    DxA(1, 2)=rhox
    DxA(2, 1)=Pe*Tx
    DxA(2, 2)=rhox*U+rho*Ux
    DxA(2, 3)=0.0d0
    DxA(2, 4)=0.0d0
    DxA(2, 5)=rhox*Pe
    DxA(3, 2)=0.0d0
    DxA(3, 3)=rhox*U+rho*Ux
    DxA(3, 5)=0.0d0
    DxA(4, 2)=0.0d0
    DxA(4, 4)=rhox*U+rho*Ux
    DxA(4, 5)=0.0d0
    DxA(5, 2)=pe*rhox*T+pe*rho*Tx
    DxA(5, 3)=0.0d0
    DxA(5, 4)=0.0d0
    DxA(5, 5)=rho*Ux*cv0+rhox*U*cv0
    DxA=DxA/hx
    call this%DxA%set(DxA)


            !write(ENO_BF, *)'matA'
            !call this%A%Print()
            !write(ENO_BF, *)'matB'
            !call this%B%Print()
            !write(ENO_BF, *)'matC'
            !call this%C%Print()
            !write(ENO_BF, *)'matD'
            !call this%D%Print()
            !write(ENO_BF, *)'matG'
            !call this%G%Print()

    end subroutine set_from_flow_lpse_ie

    !> 获得\f$\mathrm{d}A/\mathrm{d}x\f$(实数矩阵形式)
    !! @param[out] DxA \f$\mathrm{d}A/\mathrm{d}x\f$
    subroutine getDxA_realmat(this, DxA)

    implicit none
    real(R_P), dimension(5, 5), intent(inout) :: DxA
    class(lpse_BF_OP_point_type), intent(in) :: this

    DxA   = this%DxA

    end subroutine getDxA_realmat

    !> 获得\f$\mathrm{d}A/\mathrm{d}x\f$(复数算子形式)
    !! @param[out] DxA \f$\mathrm{d}A/\mathrm{d}x\f$
    subroutine getDxA_opmat_cmplx(this, DxA)

    implicit none
    type(op_mat_cmplx_type), intent(inout) :: DxA
    class(lpse_BF_OP_point_type), intent(in) :: this

    DxA   = this%DxA

    end subroutine getDxA_opmat_cmplx

    !>获得Vigneron技术处理下流向动量方程的密度和温度导数流向一阶导数系数
    !!@param[out] VigneronA 流向动量方程的密度和温度导数流向一阶导数系数
    subroutine GetVigneronA(this, VigneronA)

    implicit none
    class(lpse_BF_OP_point_type), intent(in) :: this
    real(R_P) , intent(out) :: VigneronA(2)

    VigneronA=this%VigneronA

    end subroutine GetVigneronA

    !> 获得谱半径
    subroutine GetSpt(this, rhox, rhoy, rhoz)

        use mod_gas, only: GAMMA
        implicit none
        class(lns_OP_point_type),intent(inout) :: this
        real(R_P), intent(out) :: rhox !< x方向谱半径
        real(R_P), intent(out) :: rhoy !< y方向谱半径
        real(R_P), intent(out) :: rhoz !< z方向谱半径
        real(R_P) :: acoustic
        real(R_P), dimension(5,5) :: A, B, C
        real(R_P) :: U, V, W
        real(R_P), parameter :: KAPPA=1.5d0

        call this%A%Get(A)
        call this%B%Get(B)
        call this%C%Get(C)
        acoustic=sqrt(B(3,1)*GAMMA)
        U=A(1,1)
        V=B(1,1)
        W=C(1,1)
        rhox=max(abs(U), abs(U-acoustic), abs(U+acoustic))*KAPPA
        rhoy=max(abs(V), abs(V-acoustic), abs(V+acoustic))*KAPPA
        rhoz=max(abs(W), abs(W-acoustic), abs(W+acoustic))*KAPPA

    end subroutine GetSpt



    end module mod_lns_OP_point
