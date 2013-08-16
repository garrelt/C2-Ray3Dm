!>
!! \short This module contains variables and routines for reporting
!!   the memory usage of a code
!!
!! This module uses the file /proc/self/status to establish the memory
!! use of the code using it. This works only on linux systems but Fortran
!! does not supply any internal routines for memory checking, so this
!! seems to be the best way.
!!
!! The file /proc/self/status contains a lot of information about a
!! running process. Here we only extract the memory information.
!! The relevant fields are (see
!! https://www.kernel.org/doc/Documentation/filesystems/proc.txt):
!! VmPeak                      peak virtual memory size
!! VmSize                      total program size
!! VmLck                       locked memory size
!! VmHWM                       peak resident set size ("high water mark")
!! VmRSS                       size of memory portions
!! VmData                      size of data, stack, and text segments
!! VmStk                       size of data, stack, and text segments
!! VmExe                       size of text segment
!! VmLib                       size of shared library code
!! VmPTE                       size of page table entries
!! VmSwap                      size of swap usage (the number of referred swapents)
!! Of these this module reports VmPeak, VmSize, VmRSS, VmHWM
!!
!! General module
!!
!! Dependencies: ISO_FORTRAN_ENV module (supplied by system).
!!
!! Author: Garrelt Mellema
!!
!! Date: 2013-08-07
!!
!! Version: 1.0
!!
!<

module report_memory_module

  ! Check the memory usage
  ! Uses /proc/self/status file
  
  use ISO_FORTRAN_ENV, only: IOSTAT_END

  implicit none

contains

  !> The reporting subroutine. Takes unit of report file as
  !! argument
  subroutine report_memory (report_file_unit)

    integer,intent(in) :: report_file_unit

    character(len=72) :: line
    integer :: read_status
    integer :: proc_unit !< unit to use for opening /proc/self/status
    integer :: open_status
    logical :: open_check
    integer :: lines_reported

    ! Initialize lines reported to zero
    lines_reported=0
    ! Set proc_unit to 99
    proc_unit=99
    ! Check if 99 is available, if not keep increasing proc_unit until
    ! an available unit is found.
    do
       INQUIRE(unit=proc_unit,OPENED=open_check)
       if (.not.open_check) exit
       proc_unit=proc_unit+1
    enddo

    write(report_file_unit,"(A)") "---------- Memory report ----------"

    ! Open /proc/self/status
    open(unit=proc_unit,file="/proc/self/status",status="old", &
         iostat=open_status)

    if (open_status == 0) then
       
       read(proc_unit,"(a)",iostat=read_status) line
       if (read_status == 0) &
            write(report_file_unit,"(2A)") "Memory report on Process ",line
       
       do 
          read(proc_unit,"(a)",iostat=read_status) line
          
          if (read_status == IOSTAT_END) exit ! if IOSTAT_END is not set, use -1
          
          if (read_status == 0) then ! read went well
             select case (line(1:6))
                ! List of lines to be reported, see preamble
             case ("VmPeak","VmSize","VmRSS:","VmHWM:")
                write(report_file_unit,"(a)") trim(line)
                lines_reported=lines_reported+1 ! count number of reported lines
             end select
          endif

       enddo
       
       close(proc_unit)

    else

       write(report_file_unit,"(A)") "Accessing /proc/self/status failed"

    endif
    
    ! Inform about problems during read
    if (lines_reported < 4) &
          write(report_file_unit,"(A)") &
          "Incomplete report due to read errors on /proc/self/status"
    write(report_file_unit,"(A)") "------- End of Memory report ------",line
    
  end subroutine report_memory

  !> Tests the module by allocating and deallocating several large arrays
  !! and reporting the memory usage before and after each event.
  subroutine test_module_report_memory

    !> The array to be allocated/deallocated
    real,dimension(:,:,:),allocatable :: array_buffer

    call report_memory(5)

    allocate(array_buffer(100,100,100))

    call report_memory(5)

    deallocate(array_buffer)

    call report_memory(5)

    allocate(array_buffer(200,200,200))

    call report_memory(5)

    deallocate(array_buffer)

    call report_memory(5)

  end subroutine test_module_report_memory

end module report_memory_module


  
