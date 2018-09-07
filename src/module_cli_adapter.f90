!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: /media/psf/Home/Fortran/foopse/src/module_cli_adapter.f90
!> @file
!> @breif 行命令提取器适配器文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: cli_adapter
!> @breif 行命令提取器适配器模块.
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
module mod_cli_adapter

    use flap
    implicit none

    public :: cli_parser
    private

    contains

    !> 行命令提取器.
    !!
    !! 从行命令提取配置文件文件名.
    !! @param[out] string  配置文件文件名
    subroutine cli_parser(string)

        implicit none

        type(command_line_interface) :: cli    ! Command Line Interface (CLI).
        character(kind=1, len=*), intent(out) :: string ! String value.
        integer                     :: error  ! Error trapping flag.

        call cli%init(progname = 'The General PSE Solver', &
                description = 'This a general PSE Solver for solving the &
                Parabolic Stability Equations (PSE) and its adjoint Equations', &
                authors = 'Liu Jianxin@TJU', &
                version = '0.7b', &
                examples = ['PSE_PetsC              ', &
                            'PSE_PetsC -h           ', &
                            'PSE_PetsC -f config.ini'])
        call cli%add(switch='--file', &
                     switch_ab='-f',    &
                     help='The configure file for PSE.',   &
                     required=.true.,   &
                     act='store',       &
                     error=error)
        if (error/=0) stop
        !call cli%parse(error=error)
        call cli%get(switch='-f', val=string, error=error)
        if (error/=0) stop
    !    print '(A)', cli%progname//' has been called with the following argument:'
    !    print '(A)', 'String = '//trim(adjustl(string))

    end subroutine cli_parser

end module mod_cli_adapter
