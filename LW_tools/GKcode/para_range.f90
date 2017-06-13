subroutine para_range (n1, n2, nprocs, irank, ista, iend)
  ! Assigns to each processor(irank) the starting index(ista) and the
  ! ending index(iend), if the full loop starts with index n1 and ends with
  ! n2, when nprocs (number of) processors are used.
  ! 
  ! Author: Kyungjin Ahn (27-May-2009)
  
  implicit none
  integer, intent(in)  :: n1, n2, nprocs, irank
  integer, intent(out) :: ista, iend
  integer :: iwork1, iwork2

  iwork1 =    (n2 -n1 +1)/nprocs
  iwork2 = MOD(n2 -n1 +1, nprocs)

  ista   = irank * iwork1 + n1 + MIN(irank, iwork2)
  iend   = ista + iwork1 - 1

  if (iwork2 > irank) iend = iend + 1
  
end subroutine para_range
