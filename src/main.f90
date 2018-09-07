!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_context.f90
!> @file
!> @breif 主程序文件.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2015-11-30 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-07-25
!------------------------------------------------------------------------------

Program Main

    use mod_cfgio_adapter
    !    use mod_cli_adapter
    use mod_context

    Implicit None

    integer, parameter :: LPSE=0, ALPSE=1, LST=2, ALPSE_LPSE=3, LPSE_ALPSE=4
    character(kind=1, len=256) :: cfg_file
    type(context_type) :: ctx
    !    call cli_parser(cfg_file)
    call get_command_argument(2, cfg_file)
    call cfg_loader(trim(cfg_file))
    call ctx%Initialize()

    !block
    !  use mod_debug
    !  call dbg_initial()
    !endblock
    call ctx%Solve(kindsolver)


End Program
