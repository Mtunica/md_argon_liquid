
program main
	
	!Module declaration
	!------------------
	use function_module
	use statistic_module
	
	!Variable declaration
	!--------------------
	implicit none
	
	include "../input/parameters.h"
	!Simulation parameters
	!---------------------
	double precision,allocatable::positions(:,:),velocities(:,:)
	double precision,allocatable::force(:,:),force_post(:,:),vel_mod(:)
	integer,allocatable::histogram(:)
	double precision::rc,L,sigma_therm
	
	!Histogram
	!---------
	integer::bin,maximum,minimum,thermo
	double precision::h
	
	!Loops
	!-----
	integer::i,j,k
	
	!Thermodynamics properties
	!-------------------------
	double precision::kinetic,potential,total,momentum
	double precision::temper_avg

	write(*,*) "To go directly to Euler press 1 and enter."
	write(*,*) "Els, press any other integer number and enter."
	read(*,*) i
	
	if (i.eq.1) then
		go to 56
	endif
	!------------------------------------------------------------------------------------------------------
	!First configuration
	!------------------------------------------------------------------------------------------------------
	allocate(positions(num,3), velocities(num,3),force(num,3),force_post(num,3))
	allocate(vel_mod(num))

	!------------------------------------------------------------------------------------------------------
	!Initializing
	!------------------------------------------------------------------------------------------------------
	L= (num/density)**(1.d0/3.d0)
	rc = L/4.d0

	open(14,file="data_results/temperature.dat",status="unknown")
	open(13,file="data_results/histogram.dat",status="unknown")
	open(12,file="data_results/estimation_simulation_1.dat",status="unknown")
	
	do j=5,3,-1
		
		!Initialize positions in simple cubic
		!------------------------------------
		call sc(positions,density, 3)
		
		!Initialise velocities
		!---------------------
		call binit_velocities(velocities,T,1.0d0,1.0d0)
		
		temper_avg=0.0d0

		allocate(histogram(500))
		maximum=100; minimum=0; bin = 500
		
		histogram(:)=0
		h= (maximum-minimum)/dble(bin)
		do k=1, num
			vel_mod(k) = dsqrt(velocities(k,1)**2+velocities(k,2)**2+velocities(k,3)**2)
		enddo
		histogram=freq_vector(vel_mod,maximum,minimum,bin)
		do k=1, 200
			write(13,*) (2.0d0*minimum+k*h+(k-1)*h)/2.0d0, histogram(k)/dble(num)
		enddo 
		write(13,*)
		write(13,*)
		
		!Steps of the simulation
		steps=10**j
		dt=1.0d0/(10.0d0**j)
		thermo=steps/10
		write(*,*) "Steps:",steps, " Time step:", dt
		!Run velocity verlet simulation
		!------------------------------
		do k=1,bin,1
			histogram(k)=0
		enddo
		
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
			
			!Potential calculation
			!---------------------
			call lennard_jones_potential_pb (potential,num,positions, sigma, rc,L)

			!Kinetic energy calculation
			!--------------------------
			kinetic=kinetic_energy(velocities)

			!total energy
			!------------
			Total = kinetic+potential
			
		
			!Momentum
			!--------
			momentum=dsqrt(sum(velocities(:,1))**2+sum(velocities(:,2))**2+sum(velocities(:,3))**2)
			!--------
			force(:,:) = force_post(:,:)
			
			!Write
			!-----
			write(12,'(e25.17,1x,e25.17,1x,e25.17,1x,e25.17,1x,e25.17)') i*dt,potential/dble(num),kinetic/dble(num),total/dble(num),momentum
			write(14,'(e25.17,1x,e25.17)') i*dt, 2.0*kinetic/(3.d0*dble(num)-3.d0)
			
			if (mod(i,thermo)==0) then
				do k=1, num
					vel_mod(k) = dsqrt(velocities(k,1)**2+velocities(k,2)**2+velocities(k,3)**2)
				enddo
				maximum=100; minimum=0; bin=500
				histogram=histogram+freq_vector(vel_mod,maximum,minimum,bin)
				
				temper_avg = temper_avg+ 2.0*kinetic/(3.d0*dble(num)-3.d0)
			endif

		enddo
		write(12,*)
		write(12,*)
		
		write(*,*) "Temperature in reduced units:", temper_avg/(int(steps/thermo))
		
		!Histogram end	
		h= (maximum-minimum)/dble(bin)
		do k=1, 200
			write(13,*) (2.0d0*minimum+k*h+(k-1)*h)/2.0, histogram(k)/dble(num*(steps/thermo))
		enddo 
		write(13,*)
		write(13,*)
		write(14,*)
		write(14,*)
		deallocate(histogram)
		
	enddo
	close(12)
	close(13)
	close(14)
	deallocate(positions,velocities,force,force_post,vel_mod)
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!--Euler-----------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------------------------
	!First configuration
	!------------------------------------------------------------------------------------------------------
56	allocate(positions(num,3), velocities(num,3),force(num,3))
	allocate(vel_mod(num))

	!------------------------------------------------------------------------------------------------------
	!Initializing
	!------------------------------------------------------------------------------------------------------
	L= (num/0.7d0)**(1.d0/3.d0)
	rc = L/4.d0
	dt=0.0001d0
	epsilon_=1.0d0
	sigma=1.0d0
	T=100.0d0
	density = 0.7d0
	
	open(15,file="data_results/estimation_simulation_2.dat",status="unknown")
	
	do j=5,4,-1
		
		!Initialize positions in simple cubic
		!------------------------------------
		call sc(positions,density, 3)
		
		!Initialise velocities
		!---------------------
		call binit_velocities(velocities,T,1.0d0,1.0d0)
		
		!Steps of the simulation
		steps=10**j
		dt=1.0d0/(10.0d0**j)
		write(*,*) steps,dt
		!Run velocity verlet simulation
		!------------------------------
		
		
		do i=1, steps, 1

			!Force calculation
			!-----------------
			call lennard_jones_force_pb(force,positions,rc,L)
			
			!New positions
			!-------------
			call euler_step_t(positions,positions,velocities,force,dt)
			
			!Boundary conditions
			!-------------------
			call pbc2_total(positions,L,L/2.d0)
			
			!New velocity
			!------------
			call euler_step_v(velocities,velocities,force,dt)
			
			!Potential calculation
			!---------------------
			call lennard_jones_potential_pb (potential,num,positions, sigma, rc,L)

			!Kinetic energy calculation
			!--------------------------
			kinetic=kinetic_energy(velocities)

			!total energy
			!------------
			Total = kinetic+potential
			
			!Momentum
			!--------
			momentum=dsqrt(sum(velocities(:,1))**2+sum(velocities(:,2))**2+sum(velocities(:,3))**2)
			
			!Write
			!-----
			write(15,'(e25.17,1x,e25.17,1x,e25.17,1x,e25.17,1x,e25.17)') i*dt,potential,kinetic,total,momentum

		enddo
		write(15,*)
		write(15,*)
		
	enddo	
	close(15)
	
	deallocate(positions,velocities,force,vel_mod)
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
endprogram main
