!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! VAN DER POL OSCILLATOR: SOLUTION COMPUTATION !!!!
!!!!!!!  APPLICATION OF RUNGE-KUTTA 4 METHOD   !!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM vanderpol
	IMPLICIT NONE
	INTEGER :: i,j = 0
	REAL :: ti, tf, h, temp
	INTEGER, PARAMETER :: N=4000
	REAL, DIMENSION(2) :: x0, x1
	REAL, DIMENSION(4,2) :: coeff
	REAL, DIMENSION(N,2) :: vdp

	ti=0
	tf=50
	temp=0
	h=(tf-ti)/N
	x0(1) = 1.5
	x0(2)=1
	x1=vdp_eval(x0)
	
	WRITE(*,*) "Time interval:"
	WRITE(*,*) ti, tf
	WRITE(*,*) "Number of steps:"
	WRITE(*,*) N
	WRITE(*,*) "Initial conditions:"
	WRITE(*,*) x0

	WRITE(*,*) "Solution WRITTEN on 'sol.dat'"
	

	!! WRITE SOLUTION ON FILE
	OPEN(1, file='sol.dat')

	DO WHILE (temp <= tf )
		coeff=rk4(x0,h)
		x1(1)=x0(1)+h/6*(coeff(1,1)+2*coeff(2,1)+2*coeff(3,1)+coeff(4,1))
		x1(2)=x0(2)+h/6*(coeff(1,2)+2*coeff(2,2)+2*coeff(3,2)+coeff(4,2))
		WRITE(1,*) temp, x1
		x0 = x1
		temp=temp+h
	END DO
		
	CLOSE(1)
	CALL SYSTEM('gnuplot -p')

CONTAINS

	!! IMPLEMENTING VAN DER POL EQUIVALENT I ORDER SYSTEM	
	FUNCTION vdp_eval(x)
		! Van Der Pol Oscillator
		REAL, DIMENSION(2) :: vdp_eval
		REAL, DIMENSION(2) :: x
		REAL :: mu

		mu=3

		vdp_eval(1)=x(2)
		vdp_eval(2)=mu*(1-x(1)**2)*x(2)-x(1)

	END FUNCTION vdp_eval


	
	!! IMPLEMENTING CONSTANTS OF RK4
	FUNCTION rk4(x, h)
		REAL, DIMENSION(4,2) :: rk4
		REAL, DIMENSION(2) :: x
		REAL :: h

		rk4(1,:)=vdp_eval(x)
		rk4(2,:)=vdp_eval(x+h/2*rk4(1,:))
		rk4(3,:)=vdp_eval(x+h/2*rk4(2,:))
		rk4(4,:)=vdp_eval(x+h*rk4(3,:))

	END FUNCTION rk4



END PROGRAM vanderpol
