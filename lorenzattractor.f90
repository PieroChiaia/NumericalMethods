!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! LORENZ ATTRACTOR: SOLUTION COMPUTATION !!!!
!!!!  APPLICATION OF RUNGE-KUTTA 4 METHOD   !!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM vanderpol
	IMPLICIT NONE
	INTEGER :: i,j = 0
	REAL :: ti, tf, h, temp
	INTEGER, PARAMETER :: N=4000
	REAL, DIMENSION(3) :: x0, x1
	REAL, DIMENSION(4,3) :: coeff


	ti=0
	tf=50
	temp=0
	h=(tf-ti)/N
	x0(1) = 1.5
	x0(2)=1
	x0(3)=2
	x1=dynsys(x0)
	
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
		x1(3)=x0(3)+h/6*(coeff(1,3)+2*coeff(2,3)+2*coeff(3,3)+coeff(4,3))
		WRITE(1,*) temp, x1
		x0 = x1
		temp=temp+h
	END DO
		
	CLOSE(1)
	CALL SYSTEM('gnuplot -p')

CONTAINS

	!! IMPLEMENTING VAN DER POL EQUIVALENT I ORDER SYSTEM	
	FUNCTION dynsys(x)
		
		REAL, DIMENSION(3) :: dynsys
		REAL, DIMENSION(3) :: x
		REAL :: sigma, rho, beta

		
		!! PARAMETRI SISTEMA
		sigma=10
		rho=28
		beta=8/3

		dynsys(1)=sigma*(x(2)-x(1))
		dynsys(2)=rho*x(1)-x(2)-x(1)*x(3)
		dynsys(3)=-beta*x(3)+x(1)*x(2)

	END FUNCTION dynsys


	
	!! IMPLEMENTING CONSTANTS OF RK4
	FUNCTION rk4(x, h)
		REAL, DIMENSION(4,3) :: rk4
		REAL, DIMENSION(3) :: x
		REAL :: h

		rk4(1,:)=dynsys(x)
		rk4(2,:)=dynsys(x+h/2*rk4(1,:))
		rk4(3,:)=dynsys(x+h/2*rk4(2,:))
		rk4(4,:)=dynsys(x+h*rk4(3,:))

	END FUNCTION rk4



END PROGRAM vanderpol
