
program main
	
	!Module declaration
	!------------------
	use function_module
	use statistic_module
	
	!Variable declaration
	!--------------------
	implicit none
	include "../input/parameters2.h"
	
	!Simulation parameters
	double precision,allocatable::positions(:,:),velocities(:,:),positions0(:,:),positionsnb(:,:)
	double precision,allocatable::force(:,:),force_post(:,:),displacement(:),histogram(:)
	double precision::rc,L,sigma_therm
	
	!Number of particles
	double precision:: k_boltz
	double precision,dimension(:)::densities(4)
	double precision::mean_displacement,mean,h,ek
	!Loops
	integer::i,j,k,index_,maximum,minimum,bin,count_,li
	integer::ii
	!Energies
	!--------
	double precision::kinetic,potential,total,momentum, potential_avg
	
	!Presion
	!--------
	double precision::pres11,pres22,pres_avg,pres2_avg,pres1_avg
	
	!Difussion
	!---------
	double precision,allocatable::mean_vector(:)
	integer::count_D,delta,time,time0,index_d
	
	!Radius
	!------
	double precision::dr,radius
	integer::num_rad
	integer,allocatable::rad_hist(:)
	double precision::d,dx,dy,dz
	!------------------------------------------------------------------------------------------------------
	!First configuration
	!------------------------------------------------------------------------------------------------------
	allocate(positions(num,3),positionsnb(num,3),positions0(num,3),velocities(num,3),force(num,3),force_post(num,3))
	allocate(displacement(num))
	

	!------------------------------------------------------------------------------------------------------
	!Initializing
	!------------------------------------------------------------------------------------------------------
	densities=[0.20d0,0.40d0,0.60d0,0.80d0]
	
	
	!Files
	!-----
	open(13,file="data_results/enegies.dat",status="unknown")
	open(14,file="data_results/mean_square_displacement.dat",status="unknown")
	open(15,file="data_results/histogram_b.dat",status="unknown")
	

	count_d=0
	
	!Difussivity
	!-----------
	allocate(mean_vector(thermo2))
	mean_vector=0.0d0
	!rad
	!---
	num_rad=5000
	allocate(rad_hist(num_rad))
	rad_hist=0
	
	
	do index_=1,4
	

	density =densities(index_)
	L= (num/density)**(1.d0/3.d0)
	rc = L/4.d0
	
	T=2.0d0
	
	if (index_ .eq. 4) then
		dr=L/dble(num_rad)
	endif
	
	!Initialize positions in simple cubic
	!------------------------------------
	call sc(positions,density, 3)

	!Initialise velocities
	!---------------------
	call binit_velocities(velocities,T,1.0d0,1.0d0)
	sigma_therm= dsqrt(T)

!Simulation
!------------------------------------------------------------------------------------------------------	
	!Equilibration
	!-------------
	
	
	!Previous force calculation
	!-----------------
	call lennard_jones_force_pb(force,positions,rc,L)
	
	do i=1, steps1, 1
		
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
		call  thermo_andersen(velocities,nu,sigma_therm)
		
		!New force
		!--------
		force(:,:) = force_post(:,:)
		
	enddo
	
	!Obtaining data after equilibrium
	!--------------------------------

	positions0(:,:)=positions(:,:)


	kinetic=0.0d0
	total=0.0d0
	potential_avg=0.0d0
	pres_avg=0.0d0
	pres2_avg=0.0d0
	pres1_avg=0.0d0
	count_=0
	
	!Difussion coeficient
	!--------------------
	index_d=0
	
	do i=1, steps2, 1
			
		!Force calculation
		!-----------------
		call lennard_jones_force_pb_pres(force,positions,rc,L,pres11,potential)
		
		!New positions
		!-------------
		call v_verlet_step_t(positions,positions,velocities,force,dt)
		
		positionsnb=positions
		
		!Boundary conditions
		!-------------------
		call pbc2_total(positions,L,L/2.d0)
		
		!New force
		!---------
		call lennard_jones_force_pb(force_post,positions,rc,L)
		
		!New velocity
		!------------
		call v_verlet_step_v(velocities,velocities,force,force_post,dt)
		
		!Thermostat
		!----------
		call  thermo_andersen(velocities,nu,sigma_therm)
		
		!New force
		!--------
		force(:,:) = force_post(:,:)
		
		
		
		if(index_.eq. 3) then
			index_d=index_d+1
			if(mod(i,thermo2) .eq. 0) then
				time0=time
				count_D=count_D+1
				positions0=positionsnb
				index_d=0
			else 
				do k=1, num
					displacement(k)=(positionsnb(k,1)-positions0(k,1))**2
					displacement(k)=displacement(k)+(positionsnb(k,2)-positions0(k,2))**2
					displacement(k)=displacement(k)+(positionsnb(k,3)-positions0(k,3))**2
				enddo
				
				mean_displacement =sum(displacement)/dble(num)
				mean_vector(index_d)=mean_displacement+mean_vector(index_d)
				
			endif
		endif
		
		!Compute radial function
		!-----------------------
		if (index_ .eq. 4) then
		
			do ii=1,num
				do j=ii+1,num
			
					dx = positions(ii,1)-positions(j,1)
					dx=pbc1(dx,L)
					dy = positions(ii,2)-positions(j,2)
					dy=pbc1(dy,L)
					dz = positions(ii,3)-positions(j,3)
					dz=pbc1(dz,L)
			
					d = dsqrt(dx**2+dy**2+dz**2)
				
									
					k=int(d/dr)+1
					rad_hist(k) = rad_hist(k)+2				
					
				enddo
			enddo
		
		endif
			
		!Compute thermodynamic average every thermo steps
		!------------------------------------------------
		if (mod(i,thermo)==0) then
			!Compute potential energy
			count_=count_+1
			potential_avg= potential_avg+potential		
			!Compute kinetic energy
			ek=kinetic_energy(velocities)
			kinetic = kinetic + ek
			pres22=2.0d0*ek/(3.d0*dble(num)-3.d0)*density
			
			pres1_avg=pres1_avg+pres11
			pres2_avg=pres2_avg+pres22
			pres_avg = Pres_avg + pres11 + pres22
			
		endif 	
		
	enddo
	
	kinetic=epsilon_*kinetic/dble(count_)/dble(num)
	potential_avg=epsilon_*potential_avg/dble(count_)/dble(num)
	pres_avg=pres_avg*(epsilon_/NA)/((sigma/1e10)**3.)/dble(count_)/dble(num)
	pres1_avg=pres1_avg*(epsilon_/NA)/((sigma/1e10)**3.)/dble(count_)/dble(num)
	pres2_avg=pres2_avg*(epsilon_/NA)/((sigma/1e10)**3.)/dble(count_)/dble(num)
	total= kinetic+potential_avg
	
	write(13,*) mass*density*(10.0d0**24)/(Na*(sigma**3)),kinetic,potential_avg,total,pres1_avg,pres2_avg,pres_avg
	
	
	
	enddo
	
	!Print r^2 (t)
	!-------------
	do k=1, thermo2-1
		write(14,*) 	sigma*dsqrt(mass/epsilon_)*k*dt, sigma*sigma*mean_vector(k)/count_d
	enddo
	
	do k=1, num_rad/2
		radius=(k*dr+(k-1)*dr)/2.0d0
		write(15,*)  sigma*radius, rad_hist(k)/(4.0d0*3.141592*radius*radius*dr)/density/dble(num)/dble(steps2)
	enddo
	
	!Close Files and deallocate
	!--------------------------
	close(13)
	close(14)
	close(15)
	deallocate(positions,velocities,force,force_post,positions0,displacement,positionsnb)
	deallocate(rad_hist)

endprogram main
