!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: /media/psf/Home/Fortran/START/src/module_parameter.f90
!> @file
!> @breif 基本参数文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: parameter
!> @breif 基本参数和常数模块.
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
module mod_parameter

    use penf, only: R_P
    implicit none

    real(R_P), parameter :: EPS=1.0d-8                !< 参差
    real(R_P) :: EPS_REL=1.0d-6
    complex(R_P), parameter :: CPLI=(0.0d0, 1.0d0)    !< 复数i
    real(R_P) :: Re                                 !< 雷诺数
    real(R_P) :: Te                                  !< 无穷远来流温度
    real(R_P) :: Ma                                  !< 来流马赫数

    real(R_P) :: PI=3.1415926535897932384626d0
    real(R_P) :: CriticalT=5.0d-3
    real(R_P) :: CriticalU=1.0d-2

    !! file no

    integer, parameter :: ENO_BF=99     !< 基本流文件号
    integer, parameter :: ENO_DIS=199   !< 扰动文件号

end module mod_parameter
