!>
!! \brief This module contains routines for calculating the column density
!! for a point on a 3D grid
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2013-09-05
!!
!! \b Version: 3D, H only

module column_density

  use precision, only: dp
  use sizes, only: Ndim, mesh
  use evolve_data, only: coldensh_out

  implicit none

contains

  ! ===========================================================================

  !> Finds the column density at pos as seen from the source point srcpos
  !! through interpolation. The interpolation
  !! depends on the orientation of the ray. The ray crosses either
  !! a z-plane, a y-plane or an x-plane.
    subroutine cinterp (pos,srcpos,cdensi,path)
  
    ! Author: Garrelt Mellema
    
    ! Date: 21-Mar-2006 (06-Aug-2004)
    
    ! History:
    ! Original routine written by Alex Raga, Garrelt Mellema, Jane Arthur
    ! and Wolfgang Steffen in 1999.
    ! This version: Modified for use with a grid based approach.
    ! Better handling of the diagonals.
    ! Fortran90
    
    ! does the interpolation to find the column density at pos
    ! as seen from the source point srcpos. the interpolation
    ! depends on the orientation of the ray. The ray crosses either
    ! a z-plane, a y-plane or an x-plane.
    
    integer,dimension(Ndim),intent(in) :: pos !< cell position (mesh)
    integer,dimension(Ndim),intent(in) :: srcpos !< source position (mesh)
    real(kind=dp),intent(out) :: cdensi !< column density to cell
    real(kind=dp),intent(out) :: path !< path length over cell

    real(kind=dp),parameter :: sqrt3=sqrt(3.0)
    real(kind=dp),parameter :: sqrt2=sqrt(2.0)

    integer :: i,j,k,i0,j0,k0

    integer :: idel,jdel,kdel
    integer :: idela,jdela,kdela
    integer :: im,jm,km
    integer :: ip,imp,jp,jmp,kp,kmp
    integer :: sgni,sgnj,sgnk
    real(kind=dp) :: alam,xc,yc,zc,dx,dy,dz,s1,s2,s3,s4
    real(kind=dp) :: c1,c2,c3,c4
    real(kind=dp) :: dxp,dyp,dzp
    real(kind=dp) :: w1,w2,w3,w4
    real(kind=dp) :: di,dj,dk


    !DEC$ ATTRIBUTES FORCEINLINE :: weightf
    ! map to local variables (should be pointers ;)
    i=pos(1)
    j=pos(2)
    k=pos(3)
    i0=srcpos(1)
    j0=srcpos(2)
    k0=srcpos(3)
 
    ! calculate the distance between the source point (i0,j0,k0) and 
    ! the destination point (i,j,k)
    idel=i-i0
    jdel=j-j0
    kdel=k-k0
    idela=abs(idel)
    jdela=abs(jdel)
    kdela=abs(kdel)
     
    ! Find coordinates of points closer to source
    sgni=sign(1,idel)
!      if (idel == 0) sgni=0
    sgnj=sign(1,jdel)
!      if (jdel == 0) sgnj=0
    sgnk=sign(1,kdel)
!      if (kdel == 0) sgnk=0
    im=i-sgni
    jm=j-sgnj
    km=k-sgnk
    di=real(idel)
    dj=real(jdel)
    dk=real(kdel)
 
    ! Z plane (bottom and top face) crossing
    ! we find the central (c) point (xc,xy) where the ray crosses 
    ! the z-plane below or above the destination (d) point, find the 
    ! column density there through interpolation, and add the contribution
    ! of the neutral material between the c-point and the destination
    ! point.
    
    if (kdela >= jdela.and.kdela >= idela) then
       
       ! alam is the parameter which expresses distance along the line s to d
       ! add 0.5 to get to the interface of the d cell.
       alam=(real(km-k0)+sgnk*0.5)/dk
              
       xc=alam*di+real(i0) ! x of crossing point on z-plane 
       yc=alam*dj+real(j0) ! y of crossing point on z-plane
       
       dx=2.0*abs(xc-(real(im)+0.5*sgni)) ! distances from c-point to
       dy=2.0*abs(yc-(real(jm)+0.5*sgnj)) ! the corners.
       
       s1=(1.-dx)*(1.-dy)    ! interpolation weights of
       s2=(1.-dy)*dx         ! corner points to c-point
       s3=(1.-dx)*dy
       s4=dx*dy
        
       ip =modulo(i-1, mesh(1))+1
       imp=modulo(im-1,mesh(1))+1
       jp =modulo(j-1, mesh(2))+1
       jmp=modulo(jm-1,mesh(2))+1
       kmp=modulo(km-1,mesh(3))+1
       c1=     coldensh_out(imp,jmp,kmp)    !# column densities at the
       c2=     coldensh_out(ip,jmp,kmp)     !# four corners
       c3=     coldensh_out(imp,jp,kmp)
       c4=     coldensh_out(ip,jp,kmp)

       ! extra weights for better fit to analytical solution
       w1=   s1*weightf(c1,0)
       w2=   s2*weightf(c2,0)
       w3=   s3*weightf(c3,0)
       w4=   s4*weightf(c4,0)

       ! column density at the crossing point
       cdensi   =(c1   *w1   +c2   *w2   +c3   *w3   +c4   *w4   )/(w1+w2+w3+w4) 
       ! Take care of diagonals
       ! if (kdela == idela.or.kdela == jdela) then
       ! if (kdela == idela.and.kdela == jdela) then
       ! cdensi=sqrt3*cdensi
       !else
       !cdensi=sqrt2*cdensi
       !endif
       !endif

       if (kdela == 1.and.(idela == 1.or.jdela == 1)) then
          if (idela == 1.and.jdela == 1) then
             cdensi=   sqrt3*cdensi
          else
             cdensi=   sqrt2*cdensi
          endif
       endif
       ! if (kdela == 1) then
       ! if ((w3 == 1.0).or.(w2 == 1.0)) cdensi=sqrt(2.0)*cdensi
       ! if (w1 == 1.0) cdensi=sqrt(3.0)*cdensi
       ! write(logf,*) idela,jdela,kdela
       !endif

       ! Path length from c through d to other side cell.
       !dxp=di/dk
       !dyp=dj/dk
       path=sqrt((di*di+dj*dj)/(dk*dk)+1.0) ! pathlength from c to d point  
 

       ! y plane (left and right face) crossing
       ! (similar approach as for the z plane, see comments there)
    elseif (jdela >= idela.and.jdela >= kdela) then
           
       alam=(real(jm-j0)+sgnj*0.5)/dj
       zc=alam*dk+real(k0)
       xc=alam*di+real(i0)
       dz=2.0*abs(zc-(real(km)+0.5*sgnk))
       dx=2.0*abs(xc-(real(im)+0.5*sgni))
       s1=(1.-dx)*(1.-dz)
       s2=(1.-dz)*dx
       s3=(1.-dx)*dz
       s4=dx*dz
       ip=modulo(i-1,mesh(1))+1
       imp=modulo(im-1,mesh(1))+1
       jmp=modulo(jm-1,mesh(2))+1
       kp=modulo(k-1,mesh(3))+1
       kmp=modulo(km-1,mesh(3))+1

       c1=  coldensh_out(imp,jmp,kmp)

       c2=  coldensh_out(ip,jmp,kmp)

       c3=  coldensh_out(imp,jmp,kp)

       c4=  coldensh_out(ip,jmp,kp)

       ! extra weights for better fit to analytical solution
       w1=s1*weightf(c1,0)
       w2=s2*weightf(c2,0)
       w3=s3*weightf(c3,0)
       w4=s4*weightf(c4,0)

       cdensi=   (c1   *w1   +c2   *w2   +c3   *w3   +c4   *w4   )/(w1+w2+w3+w4)

       ! Take care of diagonals
       if (jdela == 1.and.(idela == 1.or.kdela == 1)) then
          if (idela == 1.and.kdela == 1) then
             !write(logf,*) 'error',i,j,k
             cdensi=   sqrt3*cdensi
          else
             !write(logf,*) 'diagonal',i,j,k
             cdensi=   sqrt2*cdensi
          endif
       endif

       !dxp=di/dj
       !dzp=dk/dj
       !path=sqrt(dxp*dxp+1.0+dzp*dzp)
       path=sqrt((di*di+dk*dk)/(dj*dj)+1.0)
       

       ! x plane (front and back face) crossing
       ! (similar approach as with z plane, see comments there)

    elseif(idela >= jdela.and.idela >= kdela) then

       alam=(real(im-i0)+sgni*0.5)/di
       zc=alam*dk+real(k0)
       yc=alam*dj+real(j0)
       dz=2.0*abs(zc-(real(km)+0.5*sgnk))
       dy=2.0*abs(yc-(real(jm)+0.5*sgnj))
       s1=(1.-dz)*(1.-dy)
       s2=(1.-dz)*dy
       s3=(1.-dy)*dz
       s4=dy*dz
  
       imp=modulo(im-1,mesh(1))+1
       jp= modulo(j-1,mesh(2))+1
       jmp=modulo(jm-1,mesh(2))+1
       kp= modulo(k-1,mesh(3))+1
       kmp=modulo(km-1,mesh(3))+1
       c1=  coldensh_out(imp,jmp,kmp)
       c2=  coldensh_out(imp,jp,kmp)
       c3=  coldensh_out(imp,jmp,kp)
       c4=  coldensh_out(imp,jp,kp)

       ! extra weights for better fit to analytical solution
       w1   =s1*weightf(c1,0)
       w2   =s2*weightf(c2,0)
       w3   =s3*weightf(c3,0)
       w4   =s4*weightf(c4,0)

       cdensi   =(c1   *w1   +c2   *w2   +c3   *w3   +c4   *w4   )/(w1+w2+w3+w4)

       if ( idela == 1 .and. ( jdela == 1 .or. kdela == 1 ) ) then
          if ( jdela == 1 .and. kdela == 1 ) then
             cdensi=   sqrt3*cdensi
          else
             cdensi   =sqrt2*cdensi
          endif
       endif
        
       !dyp=dj/di
       !dzp=dk/di
       !path=sqrt(1.0+dyp*dyp+dzp*dzp)
       path=sqrt(1.0+(dj*dj+dk*dk)/(di*di))
       
    end if
  
  end subroutine cinterp

  ! =========================================================================

  !> Weight function for interpolation in cinterp
  real(kind=dp) function weightf (cd,id)

    use cgsphotoconstants, only: sigma_HI_at_ion_freq
    real(kind=dp):: sig
    real(kind=dp),intent(in) :: cd
    integer,intent(in) :: id
    real(kind=dp),parameter :: minweight=1.0_dp/0.6_dp

    !weightf=1.0
    ! weightf=1.0/max(1.0d0,cd**0.54)
    ! weightf=exp(-min(700.0,cd*0.15*6.3d-18))
    sig=sigma_HI_at_ion_freq

    weightf=1.0/max(0.6_dp,cd*sig)

    ! weightf=1.0/log(max(e_ln,cd))

  end function weightf

end module column_density
