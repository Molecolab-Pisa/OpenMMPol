#include "version.h"
#define OMMP_TIMING

module mod_profiling
    !! Unified Input/Output handling across the code.
    
    use mod_constants, only: OMMP_VERBOSE_DEBUG, &
                             OMMP_VERBOSE_HIGH, &
                             OMMP_VERBOSE_LOW, &
                             OMMP_VERBOSE_NONE, &
                             OMMP_VERBOSE_DEFAULT, &
                             OMMP_STR_CHAR_MAX, &
                             ip, rp
    use mod_io, only: fatal_error, ommp_message
    use mod_memory, only: mem_stat

    implicit none
    private

#ifdef OMMP_TIMING
    integer(ip), parameter :: ntimes = 128
    integer(ip) :: tcnt = 1
    real(rp) :: times(ntimes)
    integer(ip) :: maxmem(ntimes)
#endif

    public :: time_pull, time_push

    contains
    
    subroutine time_push()
        implicit none
#ifdef OMMP_TIMING
        real(rp) :: omp_get_wtime
        
        if(tcnt <= ntimes) then
            times(tcnt) = omp_get_wtime()
            ! Reset the memory counter, and save current value.
            maxmem(tcnt) = mem_stat()
            tcnt = tcnt + 1
        else
            call fatal_error('time_push Cannot push another time in the buffer.')
        end if
#endif
    end subroutine

    subroutine time_pull(s)
        implicit none

        character(len=*), intent(in) :: s
#ifdef OMMP_TIMING
        real(rp) :: elap, mm
        character(len=OMMP_STR_CHAR_MAX) :: msg

        real(rp) :: omp_get_wtime

        if(tcnt > 1) then
            elap = omp_get_wtime() - times(tcnt-1)
            !! Get maximum memory usage since last time push in
            !! GB, also make it ready for the next push/pull
            mm = mem_stat(maxmem(tcnt-1)) / 1e9
            tcnt = tcnt - 1
            write(msg, "(3a, ': ', e14.6E2, ' s')") repeat('-', tcnt), '> ', s, elap
            call ommp_message(msg, OMMP_VERBOSE_HIGH, 'time')
            write(msg, "(3a, ': ', e14.6E2, ' GB')") repeat('-', tcnt), '> ', s, mm
            call ommp_message(msg, OMMP_VERBOSE_HIGH, 'memory')
        else
            call fatal_error('time_pull Cannot pull any value.')
        end if
#endif
    end subroutine

end module mod_profiling
