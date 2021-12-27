! This module contains statistical analysis functions.
! average_func: Given a vector, return the average of its elements.	
!stand_desv_func: Given a vector, computes the standart desviation of their elements.
!random_parameter_func: Take a parameter x and returns a random number included in (-x,x).

module statistic_module
	! variables declarated
	implicit none
	integer::seed_const=25
	double precision::delta_metro,xlim
	double precision::histogram_metro(10000)
	
	
	double precision::p_value, sigma_module,error
	common/p_value/p_value, sigma_module, error
	common/histogram_metro/histogram_metro,delta_metro,xlim
	
	contains

!--------------------------------------------------------------------------------
!--Montecarlo--------------------------------------------------------------------
!--------------------------------------------------------------------------------
	subroutine metropolis2 (coordinates,X,new_coordinates,func)
		
		implicit none
		double precision,intent(in)::coordinates(:)
		double precision,intent(in)::X
		double precision, intent(out)::new_coordinates(:)
		double precision,external::func
		double precision::x1,x0
		integer::i
		
		if(size(coordinates).gt.size(new_coordinates)) then
			write(*,*) "Error in metropolis algorithm. With initial coordinate of ", coordinates(1)
			write(*,*) "Dimension of input bigger than output"
			write(*,*) "Forced to stop-----------------------"
			stop
		endif
		
		do i=1, size(coordinates),1
			
			x0 = coordinates(i)
			x1 = x0 + random_parameter_func(X)
			
			
			
			if (int(func(x1)/func(x0) + rand()).ne. 0) then
			
				new_coordinates(i)= x1			
			
			else	
			
				new_coordinates(i)= x0
			
			endif
		enddo

		return
	endsubroutine metropolis2
!--------------------------------------------------------------------------------------------
	double precision function metropolis_direct(x_ini,potential,max_iter,X,func) result(Epot)
		
		implicit none
		double precision,intent(in)::x_ini
		double precision,intent(in)::X
		double precision::potential,Epot_sqrt
		double precision,external::func
		double precision::x1,x0,f1,f0
		double precision::p_value,sigma,error
		integer::i
		integer::max_iter
		integer::count_=0
		
		double precision::T,k
		
		common/p_value/p_value,sigma,error
		common/temperature/T,k
		
		x0=x_ini	
		Epot = 0.d0
		Epot_sqrt=0.d0
		do i=1, max_iter,1
			
			x1 = x0 + random_parameter_func(X)
			f0 = func(x0)
			f1= func(x1)
			
			if (int(dexp((f0-f1)/(k*T)) + rand()).ne. 0) then
				count_=count_+1
				x0=x1
				Epot = Epot+f1
				Epot_sqrt = Epot_sqrt + f1*F1
			
			else
				Epot=Epot+f0
				Epot_sqrt = Epot_sqrt + f0*f0
			
			endif
			
			
		enddo
		
		Epot=Epot/dble(max_iter)
		Epot_sqrt=Epot_sqrt/dble(max_iter)
		p_value = dble(count_)/dble(max_iter)
		sigma = Epot_sqrt - Epot*Epot
		error= dsqrt(sigma)/dsqrt(dble(max_iter))
		return
	endfunction metropolis_direct
!--------------------------------------------------------------------------------------------
	double precision function metropolis_fermi(x_ini,potential,max_iter,x) result(Epot)
		
		implicit none
		double precision,intent(in)::x_ini(:)
		double precision,intent(in)::x
		double precision::potential,Epot_sqrt
		double precision,allocatable::x0(:)
		double precision::prod,x1
		double precision::p_value,sigma,error
		integer::i,j,k,dim_
		integer::max_iter
		integer::count_
		
		double precision::histogram_metro(10000)
		double precision::delta_metro
		double precision::xlim
		integer::index_
		
		common/histogram_metro/histogram_metro,delta_metro,xlim
		common/p_value/p_value,sigma,error
		
		dim_=size(x_ini)
		allocate(x0(dim_))
		
		histogram_metro= histogram_metro*0
		x0=x_ini	
		Epot = 0.d0
		Epot_sqrt=0.d0
		count_=0

		do i=1, max_iter,1
			
			do j=1,dim_,1
				x1 = x0(j) + random_parameter_func(x)
				
				prod = 1.d0
				
				do k=1,j-1,1
					prod = prod*(x1-x0(k))*(x1-x0(k))/((x0(j)-x0(k))*(x0(j)-x0(k)))
				enddo
				do k=j+1,dim_,1
					prod = prod*(x1-x0(k))*(x1-x0(k))/((x0(j)-x0(k))*(x0(j)-x0(k)))
				enddo
				
				if (int(dexp(-(x1*x1)+x0(j)*x0(j))* prod + rand()) .ne. 0) then
				
					count_=count_+1
					x0(j)=x1
					Epot = Epot+ x1*x1*0.5
					Epot_sqrt = Epot_sqrt + x1*x1*0.5*x1*x1*0.5 
					index_=int((x0(j)-xlim)/delta_metro)
					histogram_metro(index_)= histogram_metro(index_)+1
				
				else
					Epot=Epot+ x0(j)*x0(j)*0.5
					
					index_=int((x0(j)-xlim)/delta_metro)
					histogram_metro(index_)= histogram_metro(index_)+1
				
				endif
			enddo
			Epot_sqrt = Epot_sqrt + Epot*Epot
			
		enddo
		
		
		do i=1,size(histogram_metro),1
			histogram_metro(i) = histogram_metro(i)/(max_iter*dim_)
		enddo
		
		Epot=Epot/dble(max_iter)
		Epot_sqrt=Epot_sqrt/dble(max_iter)
		p_value = dble(count_)/dble(dim_*max_iter)
		sigma = Epot_sqrt - Epot*Epot
		error= dsqrt(sigma)/dsqrt(dble(max_iter))
		
		deallocate(x0)
		return
	endfunction metropolis_fermi

!=======================================================================================
!=Random Generator======================================================================
!=======================================================================================
	double precision function random_parameter_func(x)
	!This functions take a parameter x and returns a random number included in (-x,x)
	
		!Variable declaration
		implicit none
		double precision::x
		
		random_parameter_func = x*(2.d0*rand()-1.d0)
		
		return
	endfunction random_parameter_func
!---------------------------------------------------------------------------------------
	subroutine box_muller(x1,x2,sigma)
		
		implicit none
		double precision pi, sigma, x10, x20, x1, x2
       	pi = 4d0*datan(1d0)

		call random_number(x10)
		call random_number(x20)
		       
		x1=sigma*dsqrt(-2d0*(dlog(1d0-x10)))*dcos(2d0*pi*x20)
		x2=sigma*dsqrt(-2d0*(dlog(1d0-x10)))*dsin(2d0*pi*x20)
	
		return
	endsubroutine
!-----------------------------------------------------------------------------------------
!=======================================================================================
!=Statistics variables==================================================================
!=======================================================================================
	double precision function average(vec)
	
		implicit none
		double precision,intent(in)::vec(:)
		integer::i
		double precision::sum_
		
		sum_=0.0d0
		do i=1, size(vec), 1
			sum_=sum_+vec(i)
		enddo
		
		average= sum_/size(vec)
		
	endfunction average
!-----------------------------------------------------------------------------------------

	double precision function desv(vec)
	
		implicit none
		double precision,intent(in)::vec(:)
		integer::i
		double precision::aver,sum_
		
		sum_=0.0d0
		aver = average(vec)
		do i=1, size(vec),1
			sum_=sum_+ (vec(i)-aver)*(vec(i)-aver)
		enddo
		desv= dsqrt(sum_/size(vec))
	
	endfunction desv

end module statistic_module
