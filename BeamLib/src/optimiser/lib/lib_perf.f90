!->Module lib_time. Salvatore Maraniello. 10/Jul/2014
!
!->Description.-
!
!  This module includes few performance tools
!
!->Subroutines:
!  tic: starts a timer (cputime)
!  toc: stops a timer (cputime)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module lib_perf

    implicit none
    real(8) :: tstart, tstop

contains

    subroutine tic()
        call cpu_time(tstart)
    end subroutine tic

    subroutine toc()
      call cpu_time(tstop)
      print* ,"Done in... ", real(tstop-tstart), 'sec.!'
    end subroutine toc

end module lib_perf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Timer:
!real(8) :: e_time
!integer :: clock_rate, clock_start, clock_stop
!
!!! Start Timer:
!call system_clock(count_rate=clock_rate) !Find the time rate
!call system_clock(count=clock_start)     !Start Timer

!call system_clock(count=clock_stop)
!e_time = real(clock_stop-clock_start)/real(clock_rate)
!print *, 'job done in [sec.]: ', e_time
