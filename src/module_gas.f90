!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: /media/psf/Home/Fortran/START/src/module_gas.f90
!> @file
!> @breif 气体属性相关模块文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: gas
!> @breif 气体属性模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-01 - Initial Version
!  TODO_2017-08-01 - 平衡气体模块
!> @author
!> Liu Jianxin
!> @date 2017-08-01
!------------------------------------------------------------------------------
module mod_gas

    use penf, only: R_P
    use mod_parameter

    implicit none
    real(R_P), parameter :: GAMMA=1.4d0  !< 比热比
    real(R_P), parameter :: PR=0.72       !< Pr数
    real(R_P), parameter :: CM=110.4d0    !< 粘性Surthland公式参考温度 \f$ C_m \f$
    real(R_P), parameter :: CK=110.4d0    !< 热传导系数Surthland公式参考温度 \f$ C_k \f$

    contains

    !> 气体粘性系数(Surthland公式).
    !!
    !! \f[ \frac \mu {\mu_0}= \frac {T_0+C_m} {T+C_m} \left(\frac T {T_0}\right)^{3/2} \f]
    !! @param[in] T 温度
    !! @return 气体粘性系数 \f$ \mu \f$
    real(R_P) function miu(T)

      real(R_P) :: T
      real(R_P) :: Cm
      real(R_P), parameter :: c1=110.4

      Cm=c1/Te

      miu=T**1.5*(1.0d0+cm)/(T+cm)

    end function miu

    !> 气体粘性系数随温度一阶导数(Surthland公式).
    !!
    !! @param[in] T 温度
    !! @return 气体粘性系数对温度一阶导数  \f$ \frac {\partial \mu} { \partial T} \f$
    real(R_P) function miuT(T)

      real(R_P) :: T
      real(R_P) :: Cm
      real(R_P), parameter :: c1=110.4

      Cm=c1/Te

      miuT=1.5d0*T**0.5d0*(1.0d0+cm)/(T+cm)
      miuT=miuT-T**1.5d0*(1.0d0+cm)/(T+cm)**2

    end function miuT

    !> 气体粘性系数随温度二阶导数(Surthland公式).
    !!
    !! @param[in] T 温度
    !! @return 气体粘性系数对温度二阶导数  \f$ \frac {\partial^2 \mu} { \partial T^2} \f$
    real(R_P) function miuTT(T)

      implicit none

      real(R_P) :: T
      real(R_P) :: Cm
      real(R_P), parameter :: c1=110.4

      Cm=c1/Te

      miuTT=0.25d0*(1.0d0+cm)*(-1.0d0*T**2.0d0-6.0d0*T*Cm+3.0d0*Cm**2.0d0) &
            /sqrt(T)/(T+Cm)**3.0d0

     end function miuTT

     ! !>无穷远压力
     ! real(R_P) function Pe()
     !
     !   use mod_parameter, only: Ma
     !   implicit none
     !
     !   Pe=1.d0/(GAMMA*Ma**2)
     !
     ! end function Pe

end module mod_gas
