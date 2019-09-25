!>
!! \brief This module contains routines for calculating the ionization and 
!! temperature evolution of a single point on the grid (3D).
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2013-09-05
!!
!! \b Version: 3D, MPI & OpenMP

module evolve_point

  ! This version has been adapted for efficiency in order to be able
  ! to calculate large meshes.
    
  ! - evolve0D : calculate the rates for one grid point, possibly
  !              apply them
  ! - evolve0D_global: apply the rates to a single point
  ! - do_chemistry: apply the rates

  ! Needs:
  ! doric : ionization calculation for one point + photo-ionization rates
  ! tped : temperature,pressure,electron density calculation

  use precision, only: dp
  use my_mpi ! supplies all the MPI and OpenMP definitions and variables
  use file_admin, only: logf,timefile,iterdump, results_dir, dump_dir
  use mathconstants, only: pi
  use cgsconstants, only: ini_rec_colion_factors
  use abundances, only: abu_he
  use c2ray_parameters, only: minimum_fractional_change
  use c2ray_parameters, only: minimum_fraction_of_atoms
  use c2ray_parameters, only: epsilon
  use c2ray_parameters, only: type_of_clumping, use_LLS,type_of_LLS
  use c2ray_parameters, only: add_photon_losses
  use sizes, only: Ndim, mesh
  use grid, only: vol,dr
  use density_module, only: ndens
  use ionfractions_module, only: xh
  use temperature_module, only: temper, temperature_grid
  use temperature_module, only: temperature_states_dbl
  use temperature_module, only: get_temperature_point, set_temperature_point
  use temperature_module, only: set_final_temperature_point, isothermal
  use ionfractions_module, only: ionstates
  use clumping_module, only: clumping_point
  use lls_module, only: coldensh_LLS, LLS_point
  use sourceprops, only: srcpos
  use radiation_photoionrates, only: photrates, photoion_rates
  !use thermalevolution, only: thermal
  use photonstatistics, only: photon_loss, total_LLS_loss
  use tped, only: electrondens
  use doric_module, only: doric, coldens

  use evolve_data, only: phih_grid, phiheat
  use evolve_data, only: xh_av, xh_intermed
  use evolve_data, only: coldensh_out
  use evolve_data, only: photon_loss_src_thread
  use evolve_data, only: last_l,last_r
  use evolve_data, only: tn

  use column_density, only: cinterp

  implicit none

  save

  private

  ! Flag to know whether the do_chemistry routine was called with local option
  logical,public ::local_chemistry=.false.

  public evolve0d, evolve0d_global

contains

  !=======================================================================

  !> Calculates the photo-ionization rate for one cell due to one source
  !! and adds this contribution to the collective rate.
  subroutine evolve0D(dt,rtpos,ns,niter)
    
    ! Note for multiple sources
    ! We call this routine for every grid point and for every source (ns).
    ! The photo-ionization rates for each grid point are found and added
    ! to phih_grid, but the ionization fractions are not updated.
    ! For the first pass (niter = 1) it makes sense to DO update the
    ! ionization fractions since this will increase convergence speed
    ! in the case of isolated sources.

    ! column density for stopping chemistry !***how should this criterion be for including he and more than one freq bands?
    ! for the moment, leave it as it is, it's probably ok. 
    real(kind=dp),parameter :: max_coldensh=2e29!2.0e22_dp!2e19_dp 
    
    logical :: falsedummy ! always false, for tests
    parameter(falsedummy=.false.)

    ! subroutine arguments
    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: rtpos ! cell position (for RT)
    integer,intent(in)      :: ns ! source number 
    integer,intent(in)      :: niter ! global iteration number
    
    integer :: nx,nd,idim ! loop counters
    integer,dimension(Ndim) :: pos
    real(kind=dp) :: dist2,path,vol_ph
    real(kind=dp) :: xs,ys,zs
    real(kind=dp) :: coldensh_in
    real(kind=dp),dimension(0:1) :: coldenshe_in,coldenshe_out_temp
    real(kind=dp) :: ndens_p
    
    type(photrates) :: phi, dummiphi
    type(ionstates) :: ion

    ! Map pos to mesh pos, assuming a periodic mesh
    do idim=1,Ndim
       pos(idim)=modulo(rtpos(idim)-1,mesh(idim))+1
    enddo

    ! If coldensh_out is zero, we have not done this point
    ! yet, so do it. Otherwise do nothing. (grid is set to 0 for every source)
    if (coldensh_out(pos(1),pos(2),pos(3)) == 0.0) then
       ! Initialize local ionization states to the global ones
#ifdef ALLFRAC       
       do nx=0,1
          ion%h_av(nx)=max(xh_av(pos(1),pos(2),pos(3),nx),epsilon)
          ion%h(nx)=max(xh_intermed(pos(1),pos(2),pos(3),nx),epsilon)
          ion%h_old(nx)=max(xh(pos(1),pos(2),pos(3),nx),epsilon)         
       enddo
#else
       ion%h_av(1)=max(xh_av(pos(1),pos(2),pos(3)),epsilon)
       ion%h(1)=max(xh_intermed(pos(1),pos(2),pos(3)),epsilon)
       ion%h_old(1)=max(xh(pos(1),pos(2),pos(3)),epsilon)
       ion%h_av(0)=max(1.0-ion%h_av(1),epsilon)
       ion%h(0)=max(1.0-ion%h(1),epsilon)
       ion%h_old(0)=max(1.0-ion%h_old(1),epsilon)
#endif

       ! Initialize local density and temperature
       ndens_p=ndens(pos(1),pos(2),pos(3))
!       call get_temperature_point(pos(1),pos(2),pos(3),temper)
       ! Find the column density at the entrance point of the cell (short
       ! characteristics)
       
       if ( all( rtpos(:) == srcpos(:,ns) ) ) then
          ! Do not call cinterp for the source point.
          ! Set coldensh and path by hand
          coldensh_in=0.0
          path=0.5*dr(1)

          ! Find the distance to the source (average?)
          !dist2=0.5*dr(1) !NOT NEEDED         ! this makes vol=dx*dy*dz
          !vol_ph=4.0/3.0*pi*dist2**3
          vol_ph=dr(1)*dr(2)*dr(3)
          
       else
          
          ! For all other points call cinterp to find the column density
          call cinterp(rtpos,srcpos(:,ns),coldensh_in,path)

          path=path*dr(1)
          
          ! Find the distance to the source
          xs=dr(1)*real(rtpos(1)-srcpos(1,ns))
          ys=dr(2)*real(rtpos(2)-srcpos(2,ns))
          zs=dr(3)*real(rtpos(3)-srcpos(3,ns))
          dist2=xs*xs+ys*ys+zs*zs
          
          ! Find the volume of the shell this cell is part of 
          ! (dilution factor).
          vol_ph=4.0*pi*dist2*path

          ! Add LLS opacity
          ! GM/110224: previously we added this to coldensh_out,
          ! but this gives funny results for the photo-ionization
          ! rate which is based on the difference between the
          ! in and out column density. To mimick a LLS fog, it
          ! may be better to add it here.
          ! Initialize local LLS (if type of LLS is appropriate)
          if (use_LLS) then
             if (type_of_LLS == 2) call LLS_point (pos(1),pos(2),pos(3))
             coldensh_in = coldensh_in + coldensh_LLS * path/dr(1)
          endif          
       endif

       ! Only ray trace and exit. Do not touch the ionization
       ! fractions. They are updated using phih_grid in evolve0d_global.
       ! Only do chemistry if this is the first pass over the sources,
       ! and if column density is below the maximum.
       ! On the first global iteration pass it may be beneficial to assume 
       ! isolated sources, but on later passes the effects of multiple sources 
       ! has to be taken into account. 
       ! Therefore no changes to xh, xh_av, etc. should happen on later passes!
       ! This option is temporarily disabled by testing niter == -1 which
       ! is always false.
       if (niter == 1 .and. coldensh_in < max_coldensh) then
          local_chemistry=.true.
          dummiphi%photo_cell_HI=0.0_dp

          call do_chemistry (dt, ndens_p, ion, dummiphi, &
               coldensh_in, path, vol_ph, pos, ns, local=.true.)

          ! Copy ion fractions to global arrays.
          ! This will speed up convergence if
          ! the sources are isolated and only ionizing up.
          ! In other cases it does not make a difference.
#ifdef ALLFRAC
          xh_intermed(pos(1),pos(2),pos(3),1)=max(ion%h(1), &
               xh_intermed(pos(1),pos(2),pos(3),1))
          xh_intermed(pos(1),pos(2),pos(3),0)=max(epsilon,1.0- &
               xh_intermed(pos(1),pos(2),pos(3),1))
          xh_av(pos(1),pos(2),pos(3),1)=max(ion%h_av(1), &
               xh_av(pos(1),pos(2),pos(3),1))
          xh_av(pos(1),pos(2),pos(3),0)=max(epsilon,1.0_dp- &
               xh_av(pos(1),pos(2),pos(3),1))
#else
          xh_intermed(pos(1),pos(2),pos(3))=max(ion%h(1), &
               xh_intermed(pos(1),pos(2),pos(3)))
          xh_av(pos(1),pos(2),pos(3))=max(ion%h_av(1), &
               xh_av(pos(1),pos(2),pos(3)))
#endif
       endif  !if niter==!

       ! Add the (time averaged) column density of this cell
       ! to the total column density (for this source)
       ! and add the LLS column density to this.
       ! GM/110224: No! This messes up phi since phi is based
       !  upon the difference between the in and out column density.
       !  Instead add the LLS to coldensh_in, see above
       coldensh_out(pos(1),pos(2),pos(3))=coldensh_in + &
            coldens(path,ion%h_av(0),ndens_p)!,(1.0_dp-abu_he))

       ! Calculate (photon-conserving) photo-ionization rate from the
       ! column densities.

       ! Limit the calculation to a certain maximum column density (hydrogen)
       if (coldensh_in < max_coldensh) then 

          ! photoion_rates the structure of rates (photo and heating)
          phi=photoion_rates(coldensh_in,coldensh_out(pos(1),pos(2),pos(3)), &
               vol_ph,ns,ion%h_av(1))

          ! Divide the photo-ionization rates by the appropriate neutral density
          ! (part of the photon-conserving rate prescription)
          phi%photo_cell_HI=phi%photo_cell_HI/(ion%h_av(0)*ndens_p)
          
          ! Calculate the losses due to LLSs.
          ! GM/110224: Add the factor vol/vol_ph to phi, just as we do for
          ! the photon losses.
          ! GM/110302: Use phi%h_in (not out) here since this is where we draw
          ! the photons from, above.
          if (use_LLS) call total_LLS_loss(phi%photo_in_HI*vol/vol_ph, &
               coldensh_LLS * path/dr(1))          
       else
          ! If the H0 column density is above the maximum, set rates to zero
          phi%photo_cell_HI = 0.0_dp
          phi%photo_out_HI = 0.0_dp
          phi%heat = 0.0_dp
          phi%photo_in = 0.0_dp
          phi%photo_out = 0.0_dp
       endif
       
       ! Add photo-ionization rate to the global array 
       ! (this array is applied in evolve0D_global)
       phih_grid(pos(1),pos(2),pos(3))= &
            phih_grid(pos(1),pos(2),pos(3))+phi%photo_cell_HI
       if (.not. isothermal) &
            phiheat(pos(1),pos(2),pos(3))=phiheat(pos(1),pos(2),pos(3))+phi%heat

       ! Photon statistics: register number of photons leaving the grid
       ! Note: This is only the H0 photo-ionization rate
       if ( (any(rtpos(:) == last_l(:))) .or. &
            (any(rtpos(:) == last_r(:))) ) then
          photon_loss_src_thread(tn)=photon_loss_src_thread(tn) + &
               phi%photo_out*vol/vol_ph
          !photon_loss_src(1,tn)=photon_loss_src(1,tn) + phi%h_out*vol/vol_ph
       endif

    endif ! end of coldens test
    
  end subroutine evolve0D

  ! =======================================================================

  !> Calculates the evolution of the ionization state for
  !! one cell (mesh position pos) and multiple sources.
  subroutine evolve0D_global(dt,pos,conv_flag)

    ! Calculates the evolution of the hydrogen + helium ionization state for
    ! one cell (pos) and multiple sources.

    ! Author: Garrelt Mellema

    ! Date: 11-Feb-2008 (20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: Multiple sources (global update, no ray tracing)

    ! Multiple sources
    ! Global update: the collected rates are applied and the new ionization 
    ! fractions and temperatures are calculated.
    ! We check for convergence.

    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: pos ! position on mesh
    integer,intent(inout) :: conv_flag ! convergence counter

    integer :: nx,nit ! loop counters
    real(kind=dp) :: de ! electron density
!    real(kind=dp),dimension(0:1) :: yh,yh_av,yh0 ! ionization fractions
    real(kind=dp) :: yh0_av_old, yh1_av_old 
    real(kind=dp) :: ndens_p ! local number density
    type(temperature_states_dbl) :: temperature_start, temperature_end

    real(kind=dp) :: phih ! local H photo-ionization rate (only non-zero when local=.false.!)
    real(kind=dp) :: phih_total ! local total photo-ionization rate (including
                                ! photon loss term)
    real(kind=dp),dimension(0:1) :: dummy=(/0.0_dp,0.0_dp/) ! dummy array
    real(kind=dp) :: convergence
    type(ionstates) :: ion    
    type(photrates) :: phi 

    ! Initialize local ionization states to global ones
#ifdef ALLFRAC
    do nx=0,1
       ion%h(nx)=max(epsilon,xh_intermed(pos(1),pos(2),pos(3),nx))
       ion%h_old(nx)=max(epsilon,xh(pos(1),pos(2),pos(3),nx))
       ion%h_av(nx)=max(epsilon,xh_av(pos(1),pos(2),pos(3),nx)) ! use calculated xh_av
    enddo
#else
       ion%h(1)=max(epsilon,xh_intermed(pos(1),pos(2),pos(3)))
       ion%h_old(1)=max(epsilon,xh(pos(1),pos(2),pos(3)))
       ion%h_av(1)=max(epsilon,xh_av(pos(1),pos(2),pos(3))) ! use calculated xh_av
       ion%h(0)=1.0-ion%h(1)
       ion%h_old(0)=1.0-ion%h_old(1)
       ion%h_av(0)=1.0-ion%h_av(1)
#endif

    ! Initialize local scalars for density and temperature
    ndens_p=ndens(pos(1),pos(2),pos(3))
    call get_temperature_point (pos(1),pos(2),pos(3),temperature_start)
    !temper_inter,temp_av_old,temper_old)
    !avg_temper=temper

    ! Use the collected photo-ionization rates
    phi%photo_cell_HI=phih_grid(pos(1),pos(2),pos(3))
    if(.not.isothermal) phi%heat=phiheat(pos(1),pos(2),pos(3))

    ! I think instead of calling here twice get_temp, it is perhaps better to pass t_new
    ! and t_old as arguments from/to do_chemistry. (?)
    call do_chemistry (dt, ndens_p, ion, &
         phi, 0.0_dp, 1.0_dp, 0.0_dp, pos, 0, local=.false.)
    
    ! Test for global convergence using the time-averaged neutral fraction.
    ! For low values of this number assume convergence
    ! use previously calculated xh_av
#ifdef ALLFRAC
    yh0_av_old=xh_av(pos(1),pos(2),pos(3),0)
    yh1_av_old=xh_av(pos(1),pos(2),pos(3),1)
#else
    yh1_av_old=max(epsilon,xh_av(pos(1),pos(2),pos(3)))
    yh0_av_old=1.0-yh1_av_old
#endif
    call get_temperature_point (pos(1),pos(2),pos(3),temperature_end)
    !call get_temperature_point (pos(1),pos(2),pos(3),temper_inter,temp_av_new,temper_old)

    if ( (abs((ion%h_av(0)-yh0_av_old)) > minimum_fractional_change                .and. &
          abs((ion%h_av(0)-yh0_av_old)/ion%h_av(0)) > minimum_fractional_change   .and. &
              (ion%h_av(0) > minimum_fraction_of_atoms)  ).or.                       &
         (abs((temperature_start%average-temperature_end%average)/temperature_end%average) > 1.0e-1_dp).and.              &
         (abs(temperature_start%average-temperature_end%average) >     100.0_dp)                          &                  
                                                                      ) then
       conv_flag=conv_flag+1
    endif

    ! Copy ion fractions to the global arrays.
#ifdef ALLFRAC
    do nx=0,1
       xh_intermed(pos(1),pos(2),pos(3),nx)=ion%h(nx)
       xh_av(pos(1),pos(2),pos(3),nx)=ion%h_av(nx)
    enddo
#else
    xh_intermed(pos(1),pos(2),pos(3))=ion%h(1)
    xh_av(pos(1),pos(2),pos(3))=ion%h_av(1)
#endif
    ! This was already done in do_chemistry
    !if (.not.isothermal) call set_temperature_point (pos(1),pos(2),pos(3),temper1,av)
    
  end subroutine evolve0D_global

  ! ===========================================================================

  subroutine do_chemistry (dt, ndens_p, ion, &
       phi, coldensh_in, path, vol_ph, pos, ns, local)

    real(kind=dp),intent(in) :: dt !< time step
    real(kind=dp),intent(in) :: ndens_p
    real(kind=dp),intent(in) :: coldensh_in
    real(kind=dp),intent(in) :: path
    real(kind=dp),intent(in) :: vol_ph
    integer,dimension(Ndim),intent(in) :: pos !< position on mesh
    integer,intent(in)      :: ns !< source number 
    logical,intent(in) :: local !< true if doing a non-global calculation.

    real(kind=dp) :: avg_temper, temper0, temper1,temper2,temper_inter
    type(temperature_states_dbl) :: temperature_start, temperature_end
    real(kind=dp) :: yh0_av_old,oldhe1av,oldhe0av,oldhav
    real(kind=dp) :: yh1_av_old
    real(kind=dp) :: de
    real(kind=dp) :: coldensh_cell

    real(kind=dp) :: ionh0old,ionh1old
    integer :: nx
    integer :: nit

    type(photrates) :: phi
    type(ionstates) :: ion

    ! Initialize local temperature
    call get_temperature_point (pos(1),pos(2),pos(3),temperature_start)
    !avg_temper=temper1
    temper0 = temperature_start%current
    !temper0   =temper1
  
    ! Initialize local clumping (if type of clumping is appropriate)
    if (type_of_clumping == 5 .or. type_of_clumping == 6) call clumping_point (pos(1),pos(2),pos(3))
    
    nit=0
    do 
       nit=nit+1
       temper2   =temper1  ! This is the temperature from last iteration

       ! Save the values of yh_av found in the previous iteration
       yh0_av_old=ion%h_av(0)
       yh1_av_old=ion%h_av(1)

       ! Copy ionic abundances back to initial values (doric assumes
       ! that it contains this) !*** this is not longer needed because I have the
       ! *** old values in ion%
       !yh(:)=yh0(:)
              
       ! Calculate (mean) electron density
       de=electrondens(ndens_p,ion%h_av)

       ! Find total photo-ionization rate
       if (local) then
          
          ! Calculate (time averaged) column density of cell
          coldensh_cell=coldens(path,ion%h_av(0),ndens_p)!,1.0_dp-abu_he)

          ! Calculate (photon-conserving) photo-ionization rate
          phi=photoion_rates(coldensh_in,coldensh_in+coldensh_cell, &
               vol_ph,ns,ion%h_av(1))
          
          phi%photo_cell_HI=phi%photo_cell_HI/(ion%h_av(0)*ndens_p)

       else

          ! (direct plus photon losses)
          ! DO THIS HERE, yh_av is changing
          ! (if the cell is ionized, add a fraction of the lost photons)
          !if (xh_intermed(pos(1),pos(2),pos(3),1) > 0.5)
          !phi%photo_cell_HI=phih_cell
          !phi%heat=phihv_cell

          ! initialize the collisional ionization and recombinations rates 
          ! (temperature dependent)
          if (.not.isothermal) call ini_rec_colion_factors(temperature_start%average) 
          
          ! Add photon losses to the photo-ionization rates
          if (add_photon_losses) then
             call distribute_photon_losses(ion,phi,ndens_p,vol)
          endif

       endif     ! local if/else

       !      Calculate the new and mean ionization states
       !*** DO THIS ONLY IF PHI IS NOT ==0
       ! if (phi%h.ne.0.0_dp) then
       
       coldensh_cell    =coldens(path,ion%h(0),ndens_p)!,(1.0_dp-abu_he))
       
       call doric(dt, temperature_start%average, de, ndens_p, &
            ion%h, ion%h_av, phi%photo_cell_HI)!,local)! 
       de=electrondens(ndens_p,ion%h_av)
       
       temper1=temper0 
       if (.not.isothermal) &
            !GM/141021 Change thermal so that it takes old values and outputs new
            !values, but not overwrites...
            call thermal(dt,temperature_start%current, &
            temperature_start%average,de,ndens_p, &
            ion,phi)    
       
       ! Test for convergence on time-averaged neutral fraction
       ! For low values of this number assume convergence
       if ((abs((ion%h_av(0)-yh0_av_old)/ion%h_av(0)) < &
            minimum_fractional_change .or. &
            (ion%h_av(0) < minimum_fraction_of_atoms)).and. &
            (abs((temper1-temper2)/temper1) < minimum_fractional_change) & 
            ) then  
          exit
       endif
       
       ! Warn about non-convergence and terminate iteration
       if (nit > 400) then
          if (rank == 0) then   
             write(logf,*) 'Convergence failing (global) nit=', nit
             write(logf,*) 'x',ion%h_av(0)
             write(logf,*) 'h',yh0_av_old
             write(logf,*) abs(ion%h_av(0)-yh0_av_old)
          endif
          exit
       endif
    enddo

    ! Update temperature
    ! GM/130815: Why is this done here?
    if (.not. isothermal) call set_temperature_point (pos(1),pos(2),pos(3),temperature_start)
    
  end subroutine do_chemistry
  
end module evolve_point
