program main
	!Module declaration
	!------------------
	use function_module
	use statistic_module
	
	!Variable declaration
	!--------------------
	implicit none
	
	include "../input/parameters0.h"
	!Simulation parameters
	!---------------------
	double precision,allocatable::positions(:,:),velocities(:,:)
	double precision,allocatable::force(:,:),force_post(:,:)
	double precision::rc,L,sigma_therm
	
	!Loops
	!-----
	integer::i,j,k

	!------------------------------------------------------------------------------------------------------
	!First configuration
	!------------------------------------------------------------------------------------------------------
	allocate(positions(num,3), velocities(num,3),force(num,3),force_post(num,3))
	
	!------------------------------------------------------------------------------------------------------
	!Initializing
	!------------------------------------------------------------------------------------------------------
	L= (num/density)**(1.d0/3.d0)
	rc = L/4.d0

	!Initialize positions in simple cubic
	!------------------------------------
	call sc(positions,density, 3)
		
	!Initialise velocities
	!---------------------
	velocities=0.0d0
		
	write(*,*) "Steps:",steps, " Time step:", dt
	
	!Run velocity verlet simulation
	!------------------------------
	!Previous Force calculation
	!--------------------------
	call lennard_jones_force_pb(force,positions,rc,L)
		
	do i=1, steps, 1
			
		!New positions
		!-------------
		call v_verlet_step_t(positions,positions,velocities,force,dt)
			
		!Boundary conditions
		!-------------------
		call pbc2_total(positions,L,L/2.d0)
			
		!New force
		!---------
		call lennard_jones_force_pb (force_post,positions,rc,L)
			
		!New velocity
		!------------
		call v_verlet_step_v(velocities,velocities,force,force_post,dt)
		
		!Thermostat
		!----------
		call thermo_andersen(velocities,nu,sigma_therm)
			
		!--------
		force(:,:) = force_post(:,:)
			
	enddo

	deallocate(positions,velocities,force,force_post)
!--------------
endprogram main
