!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: /media/psf/Home/Fortran/START/src/module_cfgio_adapter.f90
!> @file
!> @breif cfgio代码适配器文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: cfgio_adapter
!> @breif cfgio库的适配器模块.
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
module mod_cfgio_adapter

    use cfgio_mod
    use mod_parameter
    use mod_dis_wavenum
    use penf, only: R_P
    use stringifor
    implicit none

    public :: cfg_loader, prefix
    public :: wavenum_init, gridfile, flowfile, Curvaturefile
    public :: bctype
    public :: Istart
    public :: Iend
    public :: Iloc
    public :: kindSolver
    public :: omega0, beta0, m_index, n_index, mdim, ndim, amp
    public :: isCrossflow
    public :: Coef_Nonlinear !< 非线性系数
    public :: Coef_Non_sub_loc
    private

    type(string) :: gridfile   !<网格文件名
    type(string) :: flowfile   !<流场文件名
    type(String) :: Curvaturefile     !<曲率文件名
    type(string) :: prefix            !<文件名前缀
    integer :: bctype(5, 2)             !<边界条件类型
    type(dis_wavenum_type), allocatable :: wavenum_init(:) !<初始波数
    real(R_P) :: omega0 !<时间基本频率
    real(R_P) :: beta0 !<展向基本波数
    integer, allocatable :: m_index(:), n_index(:) !<非线性扰动序号
    real(R_P), allocatable :: amp(:) !<扰动初始幅值
    integer :: mdim, ndim !<频率和展向谱空间维度
    type(cfg_t), save :: cfg                            !<cfg类
    integer, parameter :: kd = kind('default')
    integer :: Istart         !< 流向计算域开始位置点序号
    integer :: Iend  !< 流向计算域结束位置点序号
    integer :: Iloc !< 当前流向位置
    integer :: kindSolver !< 求解器类型
    logical :: isCrossflow
    real(R_P) :: Coef_Nonlinear=1.0d0
    integer :: Coef_Non_sub_loc=50

    contains

    !> 从配置文件读取信息
    !! @param[in] cfgfn 配置文件文件名
    subroutine cfg_loader(cfgfn)

        implicit none
        character(kind=kd, len=*), intent(in) :: cfgfn

        cfg=parse_cfg(trim(cfgfn))
        call cfg_general(cfg)
        call cfg_domain(cfg)
        call cfg_filename(cfg)
        call cfg_freestream(cfg)
        call cfg_solver(cfg)
        call cfg_disturbances(cfg)
        call cfg_generic_bc(cfg)
        call cfg%write()
      !  pause

    end subroutine cfg_loader

    !> 通用信息提取函数
    !> @param[in] cfg cfg类
    !> @return 返回配置文件general字段相关信息
    subroutine cfg_general(cfg)

        implicit none
        type(cfg_t) :: cfg
        character(kind=kd, len=256) :: StringT

        if(cfg%has_key("general", "Title")) then
          call cfg%get("general", "Title", StringT)
        else
            StringT='Problem'
        end if

    end subroutine cfg_general

    !> 计算域信息提取函数
    !> @param[in] cfg cfg类
    !> @return 返回配置文件domain字段相关信息
    subroutine cfg_domain(cfg)

        implicit none
        type(cfg_t) :: cfg

        if(cfg%has_key("Domain", "Istart")) then
          call cfg%get("Domain", "Istart", istart)
        else
            istart=-1
        end if
        if(cfg%has_key("Domain", "ILoc")) then
          call cfg%get("Domain", "ILoc", iloc)
        else
            iloc=1
        end if

        if(cfg%has_key("Domain", "Iend")) then
          call cfg%get("Domain", "Iend", iend)
        else
            iend=-1
        end if

    end subroutine cfg_domain

    !> 文件信息提取函数
    !> @param[in] cfg cfg类
    !> @return 返回配置文件filnames字段相关信息
    subroutine cfg_filename(cfg)

        implicit none
        type(cfg_t) :: cfg
        character(kind=kd, len=256) :: grid
        character(kind=kd, len=256) :: flow
        character(kind=kd, len=256) :: Curvature
        character(kind=kd, len=256) :: Dir
        character(kind=kd, len=256) :: stringPrefix

        if(cfg%has_key("filenames", "Dir")) then
          call cfg%get("filenames", "Dir", dir)
        else
          Dir="./"
        endif
        call cfg%get("filenames", "Grid", grid)
        call cfg%get("filenames", "Flow", flow)
        call cfg%get("filenames", "Curvature", Curvature)
        if(cfg%has_key("filenames", "Output prefix") )then
          call cfg%get("filenames", "Output prefix", stringPrefix)
      prefix=trim(stringPrefix)
        else
          stringPrefix='UnNameCase'
        end if

        gridfile=trim(dir)//trim(grid)
        flowfile=trim(dir)//trim(flow)
        if(trim(Curvature)/='null')then
            Curvaturefile=trim(dir)//trim(Curvature)
        else
            Curvaturefile='null'
        end if

    end subroutine cfg_filename

    !> 来流信息提取函数
    !> @param[in] cfg cfg类
    !> @return 返回配置文件freestream字段相关信息
    subroutine cfg_freestream(cfg)

        implicit none
        type(cfg_t) :: cfg

        call cfg%get("Freestream", "Ma", Ma)
        call cfg%get("Freestream", "Re", Re)
        call cfg%get("Freestream", "Te", Te)

    end subroutine cfg_freestream

    !> 求解器信息提取函数
    !> @param[in] cfg cfg类
    !> @return 返回配置文件Solver字段相关信息
    subroutine cfg_solver(cfg)

        implicit none
        type(cfg_t) :: cfg
        character(kind=kd, len=256) :: solver, isCross
        integer, parameter :: LPSE=0, ALPSE=1, LST=2, ALPSE_LPSE=3, LPSE_ALPSE=4
        integer, parameter :: NPSE=5

        if(cfg%has_key("Solver", "Type")) then
            call cfg%get("Solver", "Type", solver)
        elseif(cfg%has_key("Solver", "type")) then
            call cfg%get("Solver", "type", solver)
        else
            pause 'no solver set and LPSE is adopted!'
            solver='LPSE'
        endif

        if(cfg%has_key("Solver", "isCrossflow")) then
          call cfg%get("Solver", "isCrossflow", isCross)
        else
          isCross='no'
        endif

        if(cfg%has_key("Solver", "Non_sub_loc")) then
          call cfg%get("Solver", "Non_sub_loc", Coef_Non_sub_loc)
          if(cfg%has_key("Solver", "Non_sub_factor")) then
            call cfg%get("Solver", "Non_sub_factor", Coef_Nonlinear)
          else
            Coef_Nonlinear=1.0d0
          endif
        else
          Coef_Non_sub_loc=100000
        endif

        select case (trim(solver))
            case ('ALPSE_LPSE', 'ALPSE-LPSE')
                kindsolver=ALPSE_LPSE
            case ('LPSE_ALPSE', 'LPSE-ALPSE')
                kindsolver=LPSE_ALPSE
            case ('LPSE')
                kindsolver=LPSE
            case ('ALPSE')
                kindsolver=ALPSE
            case ('LST')
                kindsolver=LST
            case ('NPSE')
                kindsolver=NPSE
            case default
                pause "The solver type is wrong, and LPSE is used."
                kindsolver=LPSE
         end select

         select case (trim(isCross))
             case ('yes')
               isCrossflow=.True.
             case ('no')
               isCrossflow=.False.
             case default
               write(*, *) 'The crossflow stability :', trim(isCross), &
               &    'is set wrong'
               pause 'no crossflow is adopted'
               isCrossflow=.False.
          end select



    end subroutine cfg_solver

    !> 扰动信息提取函数
    !> @param[in] cfg cfg类
    !> @return 返回配置文件Disturbances字段相关信息
    subroutine cfg_disturbances(cfg)

        implicit none
        type(cfg_t) :: cfg
        real(R_P) :: bat_in, ome_in
        real(R_P), allocatable :: alf_in(:)
        complex(R_P) :: alf, bat, ome
        integer :: DisNum, i
        character(3) :: chai
        character(kind=kd, len=256) :: optn
        integer :: npar
        integer :: midx, nidx
        logical :: isNPSE
        integer, parameter :: NPSE=5

        call cfg%get("Disturbances", "Disturbances Number", DisNum)

        if(kindSolver==NPSE) then
          isNPSE=.True.
          call cfg%get("Disturbances", "mdim", mdim)
          call cfg%get("Disturbances", "ndim", ndim)
          if(cfg%has_key("Disturbances", "omega0")) then
            call cfg%get("Disturbances", "omega0", omega0)
          else
            stop 'NPSE must set omega0'
          endif
          if(cfg%has_key("Disturbances", "beta0")) then
            call cfg%get("Disturbances", "beta0", beta0)
          else
            pause 'no beta0 can be found, NPSE set beta0=0'
          endif
        else
          isNPSE=.False.
        endif


        if(Disnum>1 .and. .not. isNPSE) &
          & stop 'Multiple disturbances must use NPSE.'
        if(.not. allocated(wavenum_init)) allocate(wavenum_init(DisNum))
        if(.not. allocated(m_index)) allocate(m_index(DisNum))
        if(.not. allocated(n_index)) allocate(n_index(DisNum))
        if(.not. allocated(amp)) allocate(amp(DisNum))

        if(isNPSE)then
          write(*, *)'Disturbances at inlet are:'
          do i=1, DisNum
            write(chai, "(I3.3)")i
            optn="Disturbances_"//trim(chai)
            npar=2
            if(.not. cfg%has_key(trim(optn), 'm_index').and. i==1) then
              pause 'no m_index is set in dis001, and m_index is set to 1'
              midx=1
            else
              call cfg%get(trim(optn), 'm_index', midx)
            endif
            if(.not. cfg%has_key(trim(optn), 'n_index').and. i==1) then
              pause 'no n_index is set in dis001, and n_index is set to 1'
              nidx=1
            else
              call cfg%get(trim(optn), 'n_index', nidx)
            endif
            call cfg%get(trim(optn), "amplitude", amp(i))
            call cfg%get(trim(optn), "alpha", alf_in, npar)

            alf=alf_in(1)+CPLI*alf_in(2)
            ome=real(midx, R_P)*omega0
            bat=real(nidx, R_P)*beta0
            m_index(i)=midx
            n_index(i)=nidx

            call wavenum_init(i)%Set(alf, bat, ome)
            print*, '['//trim(optn)//']'
            print*, 'omega=', wavenum_init(i)%getOmega()
            print*, 'alpha=', wavenum_init(i)%getAlpha()
            print*, 'beta=', wavenum_init(i)%getBeta()
          enddo
        else
            i=1
            write(chai, "(I3.3)")i
            optn="Disturbances_"//trim(chai)
            npar=2
            call cfg%get(trim(optn), "alpha", alf_in, npar)
            call cfg%get(trim(optn), "omega", ome_in)
            if( cfg%has_key(trim(optn), "beta")) then
              call cfg%get(trim(optn), "beta", bat_in)
            else
              bat_in=0.0d0
            endif

            alf=alf_in(1)+CPLI*alf_in(2)
            ome=ome_in
            bat=bat_in

            call wavenum_init(i)%Set(alf, bat, ome)
            print*, '['//trim(optn)//']'
            print*, 'omega=', wavenum_init(i)%getOmega()
            print*, 'alpha=', wavenum_init(i)%getAlpha()
            print*, 'beta=', wavenum_init(i)%getBeta()
        endif

        print*, 'Disturbances are done.'

        !pause

    end subroutine cfg_disturbances

    !> 边界条件信息提取函数
    !> @param[in] cfg cfg类
    !> @return 返回配置文件Boundary Conditions字段相关信息
   subroutine cfg_generic_bc(cfg)

        implicit none
        type(cfg_t) :: cfg

        call cfg_bc_wall(cfg)
        call cfg_bc_freestream(cfg)

    end subroutine cfg_generic_bc

    !> 壁面边界条件提取函数
    !> @param[in] cfg cfg类
    !> @return 返回配置文件Boundary Conditions字段Wall相关信息
    subroutine cfg_bc_wall(cfg)

        implicit none
        type(cfg_t) :: cfg

        character(kind=kd, len=256) :: bc_type

        if(cfg%has_key("Boundary Conditions", "Wall_Vel")) then
          call cfg%get("Boundary Conditions", "Wall_Vel", bc_type)
          select case (trim(bc_type))
              case ('no_slip')
                bctype(2:4, 1)=(/2000, 3000, 4000/)
              case default
                stop "Please input correct Wall_Vel boundary condition."
          end select
        else
          bctype(2:4, 1)=(/2000, 3000, 4000/)
        endif

        if(cfg%has_key("Boundary Conditions", "Wall_T")) then
          call cfg%get("Boundary Conditions", "Wall_T", bc_type)
          select case (trim(bc_type))
              case ('adiabatic')
                bctype(5, 1)=5001
              case ('isothermal')
                bctype(5, 1)=5000
              case default
                  stop "Please input correct Wall_T boundary condition."
          end select
        else
          bctype(5, 1)=5000
        endif

        bctype(1, 1)=1004

    end subroutine cfg_bc_wall

    !> 远场边界条件信息提取函数
    !> @param[in] cfg cfg类
    !> @return 返回配置文件Boundary Conditions字段Freestream相关信息
    subroutine cfg_bc_freestream(cfg)

        implicit none
        type(cfg_t) :: cfg

        character(kind=kd, len=256) :: bc_type

        if(cfg%has_key("Boundary Conditions", "Freestream"))then
          call cfg%get("Boundary Conditions", "Freestream", bc_type)
          select case (trim(bc_type))
              case ('Dirichlet')
                bctype(1:5, 2)=(/1004, 2000, 3000, 4000, 5000/)
              case default
                stop "Please input correct Freestream boundary condition."
          end select
        else
          bctype(1:5, 2)=(/1004, 2000, 3000, 4000, 5000/)
        endif

    end subroutine cfg_bc_freestream

end module mod_cfgio_adapter
