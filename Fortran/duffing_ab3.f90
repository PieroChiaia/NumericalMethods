!==================================================================
!==============      DUFFING DAMPED OSCILLATOR      ===============   
!==================================================================


PROGRAM duffing_oscillator

	IMPLICIT NONE
	
	REAL :: ti, tf, h, temp
	INTEGER, PARAMETER :: N=600000
	REAL, DIMENSION(N,2) :: sol
	REAL, DIMENSION(2) :: x0
	REAL :: omega=1

	ti=0
	tf=40000

	h=(tf-ti)/N
	
	x0(1) = 0
	x0(2) = 0
	
	CALL SYSTEM('clear')
	WRITE(*,*) ""
	WRITE(*,*) "======================================================="
	WRITE(*,*) "===========    DUFFING DAMPED OSCILLATOR    ==========="
	WRITE(*,*) "===========       AB3 EXPLICIT SCHEME       ==========="
	WRITE(*,*) "======================================================="
	WRITE(*,*) " ..        . "
	WRITE(*,*) " X + delta x + alpha x + beta x^3 = gamma cos(omega t) "
	WRITE(*,*) ""
	WRITE(*,*) "======================================================="
	
	
	WRITE(*,*) "Time interval:"
	WRITE(*,'(F12.2,F12.2)') ti,tf
	WRITE(*,*) "Number of time steps:"
	WRITE(*,*) N
	WRITE(*,*) "Initial conditions:"
	WRITE(*,'(F10.3,f10.3)') x0(1), x0(2)
	WRITE(*,*) ""
	
	CALL adams_bashforth(sol,h,x0,N)
	
	CALL poincare_map(sol,N,h,omega)
	
	WRITE(*,*) "-------------------------------------------------------"
	WRITE(*,*) "Solution WRITTEN on 'duffing_ab3.dat'"
	WRITE(*,*) "Poincarè Map WRITTEN on 'duffing_ab3_pm.dat'"
	
	


CONTAINS


!==================================================================
!	ADAMS-BASHFORTH 3 EXPLICIT SCHEME IMPLEMENTATION
	SUBROUTINE adams_bashforth(sol,h,x0,N)
	
		REAL, DIMENSION(N,2) :: sol
		REAL, DIMENSION(2) :: x0
		INTEGER :: N, i
		REAL :: h
		
		! INITIALIZE SOLUTION WITH
		
		OPEN(1, file='duffing_ab3.dat')
		
		sol(1,:)=x0
		sol(2,:)=sol(1,:)+h*duffing(h, sol(1,:))
		sol(3,:)=sol(2,:)+h*duffing(2*h, sol(2,:))
		
		WRITE(1,*) 0, sol(1,:)	
		WRITE(1,*) 1*h, sol(2,:)	
		WRITE(1,*) 2*h, sol(3,:)	
		
		DO i=3,(N-1)
			
			sol(i+1,:)=sol(i,:)+h/12*(23*duffing(i*h,sol(i,:)) &
			& -16*duffing((i-1)*h,sol(i-1,:))+5*duffing((i-2)*h,sol(i-2,:)))
			
			WRITE(1,*) (i*h), sol(i+1,:)	

		END DO
		CLOSE(1)
	
	
	END SUBROUTINE adams_bashforth


!==================================================================
!	CONSTRUCTION OF THE POINCARE MAP BY COMPUTED SOLUTION
	SUBROUTINE poincare_map(sol,N,h,omega)
	
		IMPLICIT NONE 
		
		REAL, DIMENSION(N,2) :: sol
		INTEGER :: N, i, ind
		REAL :: h, omega
		REAL :: pi=4*atan(1.d0)
		
		! INITIALIZE SOLUTION WITH
		
		OPEN(2, file='duffing_ab3_pm.dat')
		
		DO i=1,N
			ind=CEILING(1/h*(2*pi*i/omega))
			IF (ind.LE.N) THEN
				WRITE(2,*) i*h, sol(ind,:)
			END IF
		END DO
		
		CLOSE(2)
	
	
	END SUBROUTINE poincare_map


!==================================================================
!	IMPLEMENTING DUFFING OSCILLATOR I ORDER EQUIVALENT SYSTEM	
	FUNCTION duffing(t,x)
		
		REAL, DIMENSION(2) :: duffing
		REAL, DIMENSION(2) :: x
		REAL :: delta=0.15, alpha=-1, beta=1, gam=0.3, omega=1
		REAL :: t


		duffing(1)=x(2)
		duffing(2)=-delta*x(2)-alpha*x(1)-beta*x(1)**3+gam*cos(omega*t)

	END FUNCTION duffing





END PROGRAM duffing_oscillator

