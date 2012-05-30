!>
!! \brief This module contains data and routines for permuting arrays randomly, but leaving elements close to their initial locations
!!
!! \b Author:  Michel Olagnon
!!
!! \b Date: May 2000
!!
!! \b Version: 
!!
Module m_ctrper
Use m_mrgrnk
Private
Integer, Parameter :: kdp = selected_real_kind(15)
public :: ctrper
!private :: kdp
private :: R_ctrper, I_ctrper, D_ctrper
interface ctrper
  module procedure d_ctrper, r_ctrper, i_ctrper
end interface ctrper
contains

Subroutine D_ctrper (XDONT, PCLS)
!   Permute array XVALT randomly, but leaving elements close
!   to their initial locations (nearbyness is controled by PCLS).
! _________________________________________________________________
!   The routine takes the 1...size(XVALT) index array as real
!   values, takes a combination of these values and of random
!   values as a perturbation of the index array, and sorts the
!   initial set according to the ranks of these perturbated indices.
!   The relative proportion of initial order and random order
!   is 1-PCLS / PCLS, thus when PCLS = 0, there is no change in
!   the order whereas the new order is fully random when PCLS = 1.
!   Michel Olagnon - May 2000.
! _________________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (InOut) :: XDONT
      Real, Intent (In) :: PCLS
! __________________________________________________________
!
      Real, Dimension (Size(XDONT)) :: XINDT
      Integer, Dimension (Size(XDONT)) :: JWRKT
      Real :: PWRK
      Integer :: I
!
      Call Random_number (XINDT(:))
      PWRK = Min (Max (0.0, PCLS), 1.0)
      XINDT = Real(Size(XDONT)) * XINDT
      XINDT = PWRK*XINDT + (1.0-PWRK)*(/ (Real(I), I=1,size(XDONT)) /)
      Call MRGRNK (XINDT, JWRKT)
      XDONT = XDONT (JWRKT)
!
End Subroutine D_ctrper

Subroutine R_ctrper (XDONT, PCLS)
!   Permute array XVALT randomly, but leaving elements close
!   to their initial locations (nearbyness is controled by PCLS).
! _________________________________________________________________
!   The routine takes the 1...size(XVALT) index array as real
!   values, takes a combination of these values and of random
!   values as a perturbation of the index array, and sorts the
!   initial set according to the ranks of these perturbated indices.
!   The relative proportion of initial order and random order
!   is 1-PCLS / PCLS, thus when PCLS = 0, there is no change in
!   the order whereas the new order is fully random when PCLS = 1.
!   Michel Olagnon - May 2000.
! _________________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (InOut) :: XDONT
      Real, Intent (In) :: PCLS
! __________________________________________________________
!
      Real, Dimension (Size(XDONT)) :: XINDT
      Integer, Dimension (Size(XDONT)) :: JWRKT
      Real :: PWRK
      Integer :: I
!
      Call Random_number (XINDT(:))
      PWRK = Min (Max (0.0, PCLS), 1.0)
      XINDT = Real(Size(XDONT)) * XINDT
      XINDT = PWRK*XINDT + (1.0-PWRK)*(/ (Real(I), I=1,size(XDONT)) /)
      Call MRGRNK (XINDT, JWRKT)
      XDONT = XDONT (JWRKT)
!
End Subroutine R_ctrper
Subroutine I_ctrper (XDONT, PCLS)
!   Permute array XVALT randomly, but leaving elements close
!   to their initial locations (nearbyness is controled by PCLS).
! _________________________________________________________________
!   The routine takes the 1...size(XVALT) index array as real
!   values, takes a combination of these values and of random
!   values as a perturbation of the index array, and sorts the
!   initial set according to the ranks of these perturbated indices.
!   The relative proportion of initial order and random order
!   is 1-PCLS / PCLS, thus when PCLS = 0, there is no change in
!   the order whereas the new order is fully random when PCLS = 1.
!   Michel Olagnon - May 2000.
! _________________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (InOut)  :: XDONT
      Real, Intent (In) :: PCLS
! __________________________________________________________
!
      Real, Dimension (Size(XDONT)) :: XINDT
      Integer, Dimension (Size(XDONT)) :: JWRKT
      Real :: PWRK
      Integer :: I
!
      Call Random_number (XINDT(:))
      PWRK = Min (Max (0.0, PCLS), 1.0)
      XINDT = Real(Size(XDONT)) * XINDT
      XINDT = PWRK*XINDT + (1.0-PWRK)*(/ (Real(I), I=1,size(XDONT)) /)
      Call MRGRNK (XINDT, JWRKT)
      XDONT = XDONT (JWRKT)
!
End Subroutine I_ctrper
end module m_ctrper
