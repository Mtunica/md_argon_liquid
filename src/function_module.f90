module function_module

	implicit none

	contains
!-------------------------------------------------------------------------------
!===============================================================================
!Construction===================================================================
!===============================================================================
!*******************************************************************************
!SC-----------------------------------------------------------------------------
!Subroutine that given a matrix,  the density and the dimension of the 
!coordinates returns the matrix with the positions of a sc lattice. 
!The number of particles correspond to the number of rows of the matrix.
!*******************************************************************************
!Input: positions (matrix of dimension number of atoms, coordinates).
!		density of N_atoms/Volume
!*******************************************************************************
subroutine sc(positions,density, dim_)

	implicit none
	double precision,intent(out)::positions(:,:)
	integer::i,j,k,nx,ny,nz, M
	integer::num
	double precision::L,a
	double precision,intent(in)::density
	integer,intent(in)::dim_
	
	num=size(positions(:,1))
	!dim_=size(positions(1,:))
	
	M = (num)**(1.0/real(dim_))
	
	L = (num/density)**(1.0/real(dim_))
	a= dble(L/M)
	i=1
	do nx=0, M-1, 1
		do ny=0,M-1,1
			do nz=0,M-1,1
				positions(i,1) = 0.0d0 + a*nx
				positions(i,2) = 0.d0 + a*ny
				positions(i,3) = 0.d0 + a*nz
				i=i+1
			enddo
		enddo
	enddo

	return
endsubroutine sc

!*******************************************************************************
!Binit velocities---------------------------------------------------------------
!Subroutine that initializes the velocities by a bimodal given a temperature
!*******************************************************************************
!Input: T (temperature).
!		m (mass)
!		k (boltzmann constant)
!Output: vel (velocities)
!*******************************************************************************
subroutine binit_velocities(vel,T,k,m)

	implicit none
	double precision,intent(inout)::vel(:,:)
	double precision,intent(in)::T,k,m
	double precision::x1,x2,rnum
	integer::i,j,ind_pos,ind_neg,num
	
	x1= dsqrt(T*k/m)
	x2=-dsqrt(T*k/m)
	
	num= size(vel,1)
	ind_neg=0
	ind_pos=0
	do i = 1, num
		call random_number(rnum)
		if ((rnum.lt.0.5).and.(ind_pos.lt.int(num/2))) then
			
			ind_pos=ind_pos+1
			do j=1, size(vel,2)
				vel(i,j)=x1
			enddo
		elseif (ind_neg.lt.int(num/2)) then
			ind_neg=ind_neg+1
			do j=1, size(vel,2)
				vel(i,j)=x2
			enddo
		else
			ind_pos=ind_pos+1
			do j=1, size(vel,2)
				vel(i,j)=x1
			enddo
			
		endif
	enddo
	return
	
endsubroutine binit_velocities

!===============================================================================
!Potential and forces===========================================================
!===============================================================================

!*******************************************************************************
!Kinetic Energy-----------------------------------------------------------------
!Function that given a velocity matrix returns the kineticc energy
!*******************************************************************************
!Input: velocity (matrix of dimension number of atoms, coordinates).
!Output: Kinetic energy		
!*******************************************************************************
double precision function kinetic_energy(velocity)

	implicit none
	double precision,intent(in)::velocity(:,:)
	double precision::vel
	integer::i,j
	
	kinetic_energy = 0.d0
	
	do i=1, size(velocity(:,1)),1
		vel = 0.d0
		do j=1, size(velocity(1,:)), 1
			vel = vel + velocity(i,j)*velocity(i,j)
		enddo
		
		kinetic_energy= kinetic_energy + 0.5*vel
	enddo

endfunction kinetic_energy

!*******************************************************************************
!Lennard-Jones potential--------------------------------------------------------
!Soubroutine that computes the potential energy under a lennard-jones-potential.
!*******************************************************************************
!Input: positions(matrix of dimension number of atoms, coordinates).
!		num (number of particles)
!		sigma (typical sigma of lennard-jones potential).
!		rc  (cut-off radio) 
!Output: Upot	(potential energy)
!*******************************************************************************
subroutine lennard_jones_potential (Upot,num,positions, sigma, rc)
	implicit none
	double precision,intent(in)::positions(:,:)
	double precision::rc,dx,dy,dz,d2,d12,d6
	double precision,intent(out)::Upot
	integer::num,j,i
	double precision,intent(in)::sigma
	
	Upot = 0.0d0
	
	do i=1,num
		do j=i+1,num
		
			dx = positions(i,1)-positions(j,1)
			dy = positions(i,2)-positions(j,2)
			dz = positions(i,3)-positions(j,3)
		
			d2 = dx**2+dy**2+dz**2
			
			if (d2 .lt. ((rc)**2)) then
			
				d6=d2*d2*d2
				d12=d6*d6
				Upot = Upot + 4.d0*((sigma**12)/d12  - (sigma**6)/d6)
	
				if (rc .ne. 0.d0) then
					Upot=Upot + 4.d0*((sigma/rc)**12 - (sigma/rc)**6)
			
				endif
			
			endif
		enddo
	enddo
	
endsubroutine lennard_jones_potential

!*******************************************************************************
!Lennard-Jones potential--------------------------------------------------------
!Soubroutine that computes the potential energy under a lennard-jones-potential
!with boundaries conditions.
!*******************************************************************************
!Input: positions(matrix of dimension number of atoms, coordinates).
!		num (number of particles)
!		sigma (typical sigma of lennard-jones potential).
!		rc  (cut-off radio) 
!Output: Upot	(potential energy)
!*******************************************************************************
subroutine lennard_jones_potential_pb (Upot,num,positions, sigma, rc,L)
	implicit none
	double precision,intent(in)::positions(:,:)
	double precision::rc,dx,dy,dz,d2,d12,d6
	double precision,intent(out)::Upot
	integer::num,j,i
	double precision,intent(in)::sigma,L
	
	Upot = 0.0d0
	
	do i=1,num
		do j=i+1,num
		
			dx = positions(i,1)-positions(j,1)
			dx=pbc1(dx,L)
			dy = positions(i,2)-positions(j,2)
			dy=pbc1(dy,L)
			dz = positions(i,3)-positions(j,3)
			dz=pbc1(dz,L)
		
			d2 = dx*dx+dy*dy+dz*dz
			
			if (d2 .lt. ((rc)**2)) then
			
				d6=d2*d2*d2
				d12=d6*d6
				Upot = Upot + 4.d0*((sigma**12)/d12  - (sigma**6)/d6)
	
				if (rc .ne. 0.d0) then
					Upot=Upot - 4.d0*((sigma/rc)**12 - (sigma/rc)**6)
			
				endif
			
			endif
		enddo
	enddo
	
endsubroutine lennard_jones_potential_pb

!*******************************************************************************
!Lennard-Jones force------------------------------------------------------------
!Soubroutine that computes the force under lennard-jones-potential
!with boundaries conditions. computes also the pressure.
!*******************************************************************************
!Input: positions(matrix of dimension number of atoms, coordinates).
!		num (number of particles)
!		sigma (typical sigma of lennard-jones potential).
!		rc  (cut-off radio) 
!Output: pres (presion)
!		 Upot (potential energy)		 		
!*******************************************************************************
subroutine lennard_jones_force_pb_pres(Fpot,positions, rc,L,pres,pot)

	implicit none
	double precision,intent(in)::positions(:,:)
	double precision::rc,dx,dy,dz,d
	integer::num,j,i
	double precision,intent(in)::L
	double precision,intent(out)::Fpot(:,:)
	
	double precision::pres,volum,density,factp,facte
	double precision :: cutoff_pot, cutoff_pres, Fij_rij,pot
	double precision,parameter::pi= 3.14159d0

	Fpot(:,:) = 0.0d0
	
	num=size(Fpot(:,1))
		
	volum = L**3.0d0
	density = dble(num)/volum
	
	facte = (8.d0/3.d0)*pi*dble(num)*density
	
    factp = (16.d0/3.d0)*pi*(density*density)
    
    pres = 0.d0
	pot=0.0d0
	do i=1,num
		do j=i+1,num
			
			dx = positions(i,1)-positions(j,1)
			dx=pbc1(dx,L)
			dy = positions(i,2)-positions(j,2)
			dy=pbc1(dy,L)
			dz = positions(i,3)-positions(j,3)
			dz=pbc1(dz,L)
			
			d = dsqrt(dx**2+dy**2+dz**2)
				
			if (d .lt. rc) then
					  
					Fpot(i,1) = Fpot(i,1)+ (48/d**14-24/d**8)*dx
					Fpot(i,2) = Fpot(i,2) + (48/d**14-24/d**8)*dy
					Fpot(i,3) = Fpot(i,3) + (48/d**14-24/d**8)*dz
					
					Fpot(j,1) =Fpot(j,1) - (48/d**14-24/d**8)*dx
					Fpot(j,2) = Fpot(j,2) - (48/d**14-24/d**8)*dy
					Fpot(j,3) = Fpot(j,3) - (48/d**14-24/d**8)*dz
					
					Fij_rij = (Fpot(i,1)-Fpot(j,1))*dx
					Fij_rij =Fij_rij +(Fpot(2,i)-Fpot(2,j))*dy
					Fij_rij =Fij_rij +(Fpot(2,i)-Fpot(2,j))*dz
					 
					pot = pot + 4.d0*(1.d0/(d**12.) - 1.d0/(d**6.))
           			pres = pres + d*(48/d**14-24/d**8)
			endif
		enddo
	enddo
	
	 pot = pot + facte*((1.d0/3.d0)/(rc**9.) - 1.d0/(rc**3.))
     pres = (1.d0/(3.d0*volum))*pres
     pres = pres + factp*((2.d0/3.d0)/(rc**9.) - 1.d0/(rc**3.))
endsubroutine lennard_jones_force_pb_pres


!*******************************************************************************
!Lennard-Jones force------------------------------------------------------------
!Soubroutine that computes the force matrix under a lennard-jones-potential
!*******************************************************************************
!Input: positions(matrix of dimension number of atoms, coordinates).
!		num (number of particles)
!		rc  (cut-off radio) 
!Output: Fpot	(potential energy)
!*******************************************************************************
subroutine lennard_jones_force(Fpot,num,positions, rc)

	implicit none
	double precision,intent(in)::positions(:,:)
	double precision::rc,dx,dy,dz,d
	integer::num,j,i
	double precision,intent(inout)::Fpot(:,:)
	
	Fpot(:,:) = 0.0d0
	
	do i=1,num
		do j=i+1,num
		
				dx = positions(i,1) - positions(j,1) 
				dy = positions(i,2)-positions(j,2)
				dz = positions(i,3)-positions(j,3)
		
				d = (dx**2+dy**2+dz**2)**(1.d0/2.d0)
			
			if (d .lt. rc) then
			
				if (rc .ne. 0.d0) then
					Fpot(i,1) = Fpot(i,1) + (48/d**14-24/d**8)*dx
					Fpot(i,2) = Fpot(i,2) + (48/d**14-24/d**8)*dy
					Fpot(i,3) = Fpot(i,3) + (48/d**14-24/d**8)*dz
					
					Fpot(j,1) =Fpot(j,1) - (48/d**14-24/d**8)*dx
					Fpot(j,2) = Fpot(j,2) - (48/d**14-24/d**8)*dy
					Fpot(j,3) = Fpot(j,3) - (48/d**14-24/d**8)*dz
			
				endif
			
			endif
		enddo
	enddo

endsubroutine lennard_jones_force

!*******************************************************************************
!Lennard-Jones force boundary-----------------------------------------------
!Soubroutine that computes the force matrix under a lennard-jones-potential
!considering boundary conditions
!*******************************************************************************
!Input: positions(matrix of dimension number of atoms, coordinates).
!		L (longitude)
!		rc  (cut-off radio) 
!Output: Fpot	(potential energy)
!*******************************************************************************
subroutine lennard_jones_force_pb (Fpot,positions, rc,L)

	implicit none
	double precision,intent(in)::positions(:,:)
	double precision::rc,dx,dy,dz,d
	integer::num,j,i
	double precision,intent(in)::L
	double precision,intent(out)::Fpot(:,:)
	
	Fpot(:,:) = 0.0d0
	
	num=size(Fpot(:,1))
	do i=1,num
		do j=i+1,num
			
			dx = positions(i,1)-positions(j,1)
			dx=pbc1(dx,L)
			dy = positions(i,2)-positions(j,2)
			dy=pbc1(dy,L)
			dz = positions(i,3)-positions(j,3)
			dz=pbc1(dz,L)
			
			d = dsqrt(dx**2+dy**2+dz**2)
				
			if (d .lt. rc) then
			
					Fpot(i,1) = Fpot(i,1) + (48/d**14-24/d**8)*dx
					Fpot(i,2) = Fpot(i,2) + (48/d**14-24/d**8)*dy
					Fpot(i,3) = Fpot(i,3) + (48/d**14-24/d**8)*dz
					
					Fpot(j,1) =Fpot(j,1) - (48/d**14-24/d**8)*dx
					Fpot(j,2) = Fpot(j,2) - (48/d**14-24/d**8)*dy
					Fpot(j,3) = Fpot(j,3) - (48/d**14-24/d**8)*dz
			
			endif
		enddo
	enddo
	
endsubroutine lennard_jones_force_pb

!===============================================================================
!Periodic conditions============================================================
!===============================================================================

!*******************************************************************************
!Periodic boundary conditions----------------------------------------------------
!Function that correct the distance if we consider periodic boundary conditions
!*******************************************************************************
!Input:	L (longitude box)
!		x  (distance) 
!Output: Corrected dinstance
!*******************************************************************************
double precision function  pbc1(x,L)

	implicit none
	double precision,intent(in)::x,L
	
	if (x .gt. L/2.d0) then
	
		pbc1=x-L
	elseif (x .lt. -L/2.d0) then
		pbc1=x+L
	else
		pbc1=x
	endif
	
endfunction pbc1

!*******************************************************************************
!Periodic boundary conditions---------------------------------------------------
!Function that correct the position if we consider periodic boundary conditions
!*******************************************************************************
!Input:	L (longitude box)
!		x  (position) 
!		origin (center of the box)
!Output: Corrected position
!*******************************************************************************
double precision function  pbc2(x,L,origin)

	implicit none
	double precision,intent(in)::x,L,origin
	
	if (x .gt. origin+(L/2.d0)) then
	
		pbc2=x-L
	elseif (x .lt. origin-(L/2.d0)) then
		pbc2=x+L
	else
		pbc2=x
	endif
	
endfunction pbc2

!*******************************************************************************
!Periodic boundary conditions---------------------------------------------------
!Function that correct the positions if we consider periodic boundary conditions
!*******************************************************************************
!Input:	L (longitude box)
!		x  (position) 
!		origin (center of the box)
!Output: Corrected positions
!*******************************************************************************
subroutine pbc2_total(positions,L,origin)

	implicit none
	double precision,intent(inout)::positions(:,:)
	double precision,intent(in)::L,origin
	integer::i,j
	
	do i=1,size(positions(:,1)),1
		do j=1,size(positions(1,:)),1
			positions(i,j)=pbc2(positions(i,j),L,origin)
		enddo
	enddo
	
	return
endsubroutine


!===============================================================================
!Integration====================================================================
!===============================================================================

subroutine euler_step_v(velnew,velocity,force,deltat)

	implicit none
	double precision,intent(in)::force(:,:),deltat
	double precision,intent(out)::velnew(:,:),velocity(:,:)
	integer::dim_,i,j
	
	velnew(:,:) = velocity(:,:) + deltat*force(:,:)
	
endsubroutine
!-------------------------------------------------------------------------------
subroutine euler_step_t(trajectory,positions,velocity,force,deltat)

	implicit none
	double precision,intent(in)::positions(:,:),force(:,:),deltat
	double precision,intent(inout)::trajectory(:,:),velocity(:,:)

	trajectory(:,:) = positions(:,:) + deltat*velocity(:,:) + force(:,:)*0.5*deltat*deltat
	
endsubroutine
!-------------------------------------------------------------------------------
subroutine verlet_step_t(trajectory,positions,positions_prev,force,deltat)

	implicit none
	double precision,intent(in)::positions(:,:),positions_prev(:,:),force(:,:),deltat
	double precision,intent(inout)::trajectory(:,:)
	
	trajectory(:,:) = 2.d0*positions(:,:)-positions_prev(:,:) + force(:,:)*deltat*deltat
	
	return
endsubroutine
!-------------------------------------------------------------------------------
subroutine v_verlet_step_t(trajectory,positions,velocity,force,deltat)

	implicit none
	double precision,intent(in)::positions(:,:),velocity(:,:),force(:,:),deltat
	double precision,intent(out)::trajectory(:,:)
	
	trajectory(:,:) = positions(:,:)+velocity(:,:)*deltat + force(:,:)*0.5d0*deltat*deltat
	
	return
endsubroutine
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine v_verlet_step_v(velnew,velocity,force,force_post,deltat)

	implicit none
	double precision,intent(in)::velocity(:,:),force(:,:),deltat,force_post(:,:)
	double precision,intent(out)::velnew(:,:)
	
	velnew(:,:) = velocity(:,:) + (force(:,:)+force_post(:,:))*0.5d0*deltat
	
	return
endsubroutine

!===============================================================================
!Thermalization-----------------------------------------------------------------
!===============================================================================

!*******************************************************************************
!Andersen-----------------------------------------------------------------------
!*******************************************************************************
!Input:	vel (velocities)
!		sigma  (sigma) 
!		nu (nu)
!Output: Corrected velocities
!*******************************************************************************
subroutine thermo_andersen(vel,nu,sigma)

	implicit none
	double precision,intent(inout)::vel(:,:)
	double precision,intent(in)::nu,sigma
	double precision::x1,x2,rnum
	integer::i
	
	do i = 1, size(vel(:,1))
		
		call random_number(rnum)
		if (rnum .lt. nu) then
		
			call box_muller_anderson(x1,x2,sigma)
			vel(i,1)= x1
			vel(i,2) =x2
			
			call box_muller_anderson(x1,x2,sigma)
			vel(i,3) = x1
			
		endif
	enddo
	return
	
endsubroutine thermo_andersen

subroutine box_muller_anderson(x1,x2,sigma)
		
	implicit none
	double precision pi, sigma, x10, x20, x1, x2
	pi = 4d0*datan(1d0)

	call random_number(x10)
	call random_number(x20)
		       
	x1=sigma*dsqrt(-2d0*(dlog(1d0-x10)))*dcos(2d0*pi*x20)
	x2=sigma*dsqrt(-2d0*(dlog(1d0-x10)))*dsin(2d0*pi*x20)
	
	return
endsubroutine

!===============================================================================
!Statistics---------------------------------------------------------------------
!===============================================================================

function freq_vector(vector, maximum, minimum, bin) result(freq)
	
	implicit none	
	double precision, intent(in) :: vector(:)
	integer, allocatable:: freq(:)
	integer:: maximum, minimum
	integer bin,i,j
	double precision::h
		
		
	h=dble(maximum-minimum)/dble(bin)
	allocate(freq(bin))
		
	do i=1,bin
		freq(i)=0
	enddo
		
	do i= 1, size(vector), 1
		do j=1, bin, 1
			if ((vector(i).lt.(minimum+j*h))) then
				freq(j) = freq(j)+1
				exit
			endif
		enddo
			
		if ((vector(i).ge.maximum) .or. (vector(i).lt.minimum)) then
			write(*,*) "out of range histogram"
			stop
		endif
	enddo
		
endfunction freq_vector
!-----------------------------------------------------------
end module function_module
	
