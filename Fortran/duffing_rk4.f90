!==================================================================
!==============      DUFFING DAMPED OSCILLATOR      ===============   
!==================================================================


PROGRAM duffing_oscillator

	IMPLICIT NONE
	
	REAL :: ti, tf, h, temp
	INTEGER, PARAMETER :: N=4000
	
	REAL :: delta, alpha, beta, gam, omega
	
	
	REAL, DIMENSION(2) :: x0, x1
	REAL, DIMENSION(4,2) :: rk


	ti=0
	tf=50
	temp=0
	h=(tf-ti)/N
	
	x0(1) = 1.5
	x0(2)=1
	
	x1=duffing(ti,x0)
	
	CALL SYSTEM('clear')
	WRITE(*,*) ""
	WRITE(*,*) "======================================================="
	WRITE(*,*) "===========    DUFFING DAMPED OSCILLATOR    ==========="
	WRITE(*,*) "===========       RK4 EXPLICIT SCHEME       ==========="
	WRITE(*,*) "======================================================="
	WRITE(*,*) ""
	
	WRITE(*,*) "Time interval:"
	WRITE(*,'(F10.3,f10.3)') ti,tf
	WRITE(*,*) "Number of time steps:"
	WRITE(*,*) N
	WRITE(*,*) "Initial conditions:"
	WRITE(*,'(F10.3,f10.3)') x0(1), x0(2)

	WRITE(*,*) "Solution WRITTEN on 'duffing_rk4.dat'"

	

	!! WRITE SOLUTION ON FILE
	OPEN(1, file='duffing_rk4.dat')

	DO WHILE (temp <= tf )
		rk=rk4(temp,x0,h)
		x1(1)=x0(1)+h/6*(rk(1,1)+2*rk(2,1)+2*rk(3,1)+rk(4,1))
		x1(2)=x0(2)+h/6*(rk(1,2)+2*rk(2,2)+2*rk(3,2)+rk(4,2))
		WRITE(1,*) temp, x1
		x0 = x1
		temp=temp+h
	END DO
		
	CLOSE(1)

CONTAINS

	!! IMPLEMENTING VAN DER POL EQUIVALENT I ORDER SYSTEM	
	FUNCTION duffing(t,x)
		! Van Der Pol Oscillator
		REAL, DIMENSION(2) :: duffing
		REAL, DIMENSION(2) :: x
		REAL :: delta=0.02, alpha=1, beta=5, gam=8, omega=0.5
		REAL :: t


		duffing(1)=x(2)
		duffing(2)=-delta*x(2)-alpha*x(1)-beta*x(1)**3+gam*cos(omega*t)

	END FUNCTION duffing


	
	!! IMPLEMENTING CONSTANTS OF RK4
	FUNCTION rk4(t,x,h)
		REAL, DIMENSION(4,2) :: rk4
		REAL, DIMENSION(2) :: x
		REAL :: h,t

		rk4(1,:)=duffing(t,x)
		rk4(2,:)=duffing(t+h/2, x+h/2*rk4(1,:))
		rk4(3,:)=duffing(t+h/2, x+h/2*rk4(2,:))
		rk4(4,:)=duffing(t+h, x+h*rk4(3,:))
		

	END FUNCTION rk4



END PROGRAM duffing_oscillator

