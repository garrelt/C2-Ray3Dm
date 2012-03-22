      program ifrit_format

!      include '../code/sizes.h'

      integer,parameter :: mesh=256  

      integer dummy
!     converts 2D map to a (quasi) 3D, for visualization with ifrit
      real,dimension(mesh,mesh) :: DeltaTtot

      open(unit=3,file='result_w_random_vel.bin',form='unformatted',status='unknown')
      read(3) dummy, dummy

      read(3) DeltaTtot

      open(unit=4,file='result_ifrit_random_vel.bin',form='unformatted',status='unknown')
      write(4) mesh,mesh,mesh
      write(4)(((real(DeltaTtot(i,j)*1e6),i=1,mesh),j=1,mesh),k=1,mesh)    

      stop
      end 
