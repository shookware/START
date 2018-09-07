module mod_debug

  implicit none

  integer :: jloc

  integer :: lloc
  integer :: mloc
  integer :: nloc

  integer :: mm=100
  integer :: nn=100

  integer, parameter :: file_log=2019



Contains

  subroutine dbg_initial()
    implicit none

    open(file_log, file='debug.log', action="write")

  end subroutine




end module mod_debug
