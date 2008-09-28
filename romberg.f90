!>
!! \brief This module contains Romberg integration routines.
!!
!! \b Author: Frank Robijn / Garrelt Mellema
!!
!! \b Date: August 17, 1989	
!!

module romberg

  use precision
  implicit none

  integer,parameter :: maxpow=14 !< maximum number of integration points is 2^maxpow
  real(kind=dp),dimension(0:2**maxpow, -1:maxpow) :: romw !< Romberg weights
  
contains
  
  !========================================================================

  !> Romberg initialisation procedure (r and index are calculated)
  SUBROUTINE romberg_initialisation (nmax)
    
    !	In the initialisation procedure r and index are calculated	
    !____________________________________________________________________
    !	The argument is :						
    !	nmax	The maximum number of grid points in any direction	
    !_______________________________________________________________________

    !> The maximum number of grid points in any direction
    integer,intent(in) :: nmax

    integer :: pmax
    real(kind=dp),dimension(0:maxpow) :: a 
    real(kind=dp),dimension(0:maxpow) :: b 
    real(kind=dp),dimension(0:maxpow, 0:maxpow) :: s 
        
    integer :: i,j,k
    !	The main variables are :
    
    !	romw (i,j) Romberg weight for point # i if the integration grid has
    !		2^j+1 points. For a grid with one point, j = -1.
    !	a, b	Coefficients in Richardson's extrapolation


    pmax = nint(log (real(nmax,dp))/log(2.0))
    if (pmax > maxpow) stop  &
         'ROMBERG-Error : Too many grid points'
    if (nmax /= 2 ** pmax) stop &
         'ROMBERG-Error : # grid points has to be a power of 2'
    !					Set up extrapolation constants
    do k = 1, pmax
       b (k) = -1.0 / (4.0 ** k - 1.0)
       a (k) = - b (k) * 4.0 ** k
    enddo
    do i = 1, pmax
       s (i,0) = 0.0
    enddo
    do i = 0, pmax
       do j = 0, 2 ** i
          romw (j, i) = 0.
       enddo
    enddo
    !   Calculate weight of Integral (#grid=2^k)
    do k = 0, pmax
       !					Select Integ (#2^k)
       s (k,0) = 1.0
       !					Apply Romberg procedure
       do j = 1, pmax
          do i = pmax, j, -1
             s (i,j) = a(j)*s(i,j-1) + b(j)*s(i-1,j-1)
          enddo
       enddo
       ! s (i,i) is the weight of Integ(#2^k) if #grid = 2^1
       ! 2 ** (i-k) is weight of Sum(#2^k) relative to Sum(#2^i)
       do i = k, pmax
          do j = 0, 2 ** k
             romw(2**(i-k)*j,i)=s(i,i)*2**(i-k)+romw(2**(i-k)*j,i)
          enddo
       enddo
       s (k,0) = 0.0
    enddo
    romw (0, -1) = 1.0
    ! Weights at the edge have to be halved
    do i = 0, pmax
       romw (0, i) = 0.5 * romw (0, i)
       romw (2**i, i) = 0.5 * romw (2**i, i)
    enddo
    
  end SUBROUTINE romberg_initialisation
  
  !===========================================================================

  !>	The following procedure performs a two-dimensional integral     
  !>	using first a trapezoidal rule integral and then applying the	
  !>	Romberg integration technique.					
  !>	nx and ny have to be powers of two. For a one dimensional	
  !>	integral, ny = 0. First call Romberg_Initialisation before	
  !>	using the Romberg function					
  FUNCTION Scalar_Romberg (f, w, nc, nx, ny) result(scalar_romberg_result)
    
    !	The following procedure performs a two-dimensional integral     
    !	using first a trapezoidal rule integral and then applying the	
    !	Romberg integration technique.					
    !	nx and ny have to be powers of two. For a one dimensional	
    !	integral, ny = 0. First call Romberg_Initialisation before	
    !	using the Romberg function					
    
    !       Author  Frank Robijn    					
    !	Version August 17, 1989						
    !_______________________________________________________________________
    !	The arguments are :						
    !									
    !	f	Integrand values					
    !	w	Relative weights of grid points				
    !	nc	Leading dimension of f and w arrays			
    !	nx	Number of points - 1 in x direction			
    !	ny	Number of points - 1 in y direction			
    !_______________________________________________________________________
    
    real(kind=dp) :: scalar_romberg_result !< result of integral
    integer,intent(in) :: nc !< Leading dimension of f and w arrays
    integer,intent(in) :: nx !< Number of points - 1 in x direction
    integer,intent(in) :: ny !< Number of points - 1 in y direction
    real(kind=dp),dimension(0:nc,0:ny),intent(in) :: f !< integrand values
    real(kind=dp),dimension(0:nc,0:ny),intent(in) :: w !< relative weights of grid points

    integer :: px, py, x, y
    real(kind=dp) :: integral
    
    ! nx = 2^px, ny = 2^py
    px = nint(log (real(nx,dp)) / log (2.0))
    if (ny > 0) then
       py = nint(log (real(ny,dp)) / log (2.0))
    else
       ! One dimensional integral
       py = -1
    endif
    integral = 0.0
    ! Calculate integral
    do y = 0, ny
       do x = 0, nx
          integral = integral + f (x, y) * w (x, y) * &
               romw (x, px) * romw (y, py)
       enddo
    enddo
    scalar_romberg_result = integral
    
  end FUNCTION Scalar_Romberg
  
  !=============================================================================
  
    
  !>	This subroutin performs an one-dimensional integration in	
  !>	the x-direction for every y and stores the results in itg.	
  !>	In principle it works the same as the function Romberg.		
  !>	Extra argument : itg	the results of the integrations		
  SUBROUTINE Vector_Romberg (f, w, nc, nx, ny, itg)		
    
    
    !	This subroutin performs an one-dimensional integration in	
    !	the x-direction for every y and stores the results in itg.	
    !	In principle it works the same as the function Romberg.		
    !_______________________________________________________________________
    !	Extra argument : itg	the results of the integrations		
    !_______________________________________________________________________
    
    integer,intent(in) ::  nc !< Leading dimension of f and w arrays
    integer,intent(in) ::  nx !< number of x points -1 
    integer,intent(in) ::  ny !< number of y points -1
    real(kind=dp),dimension(0:nc, 0:ny),intent(in) :: f !< integrand values
    real(kind=dp),dimension(0:nc, 0:ny),intent(in) :: w !< relative weights of grid points
    real(kind=dp),dimension(0:ny),intent(out) :: itg !< results of the integrations

    integer :: px, x, y
    
    ! nx = 2^px
    px = nint(log (real(nx,dp)) / log (2.0))
    
    do y = 0, ny
       ! One dimensional integration
       itg (y) = 0.0
       do x = 0, nx
          itg (y) = itg(y)+f(x, y)*w(x, y)*romw(x, px)
       enddo
    enddo
  end SUBROUTINE Vector_Romberg
  
end module romberg

