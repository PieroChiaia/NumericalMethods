!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!          ADVECTION EQUATION WITH          !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!            LAX-WENDROFF SCHEME            !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM advection_equation
	IMPLICIT NONE
	REAL, PARAMETER :: ti=0
	REAL, PARAMETER :: tf=10
	REAL, PARAMETER :: xs=-1
	REAL, PARAMETER :: xd=3

	REAL, PARAMETER :: CFL=0.8
	REAL, PARAMETER :: a=1

	REAL, PARAMETER :: l=0.5, h=l/20

	REAL, PARAMETER :: dt=1/a*CFL*h

	INTEGER, PARAMETER :: Nt=CEILING((tf-ti)/dt+1)
	INTEGER, PARAMETER :: Nx=CEILING((xd-xs)/h+1)

	REAL, DIMENSION(Nt,Nx) :: u
	REAL :: t_ask
	INTEGER :: ind

	CALL inputpar_print(ti,tf,xs,xd,CFL,dt,Nt,Nx)


	CALL ivp(u,Nt,Nx,l,xs,xd,h)
	CALL bvp(u,Nt,Nx,l,ti,tf,dt)
	
	CALL print_conditions(u,Nx,Nt)

	WRITE(*,*) "********************************************"
	WRITE(*,*) "   Computing solution..."

	CALL lax_wendroff(u,Nt,Nx,CFL)

	t_ask = 0.4
	ind = CEILING((t_ask-ti)/dt)
	
	CALL print_solution(u,Nt,Nx,ind,ti,dt)
			
	WRITE(*,*) "********************************************"
	WRITE(*,*) "   Solution computed at time", t_ask

CONTAINS



! ******************************************************************************
	
	! PRINT INPUT DATA PROBLEM
	SUBROUTINE inputpar_print(ti,tf,xs,xd,CFL,dt,Nt,Nx)
		REAL :: ti,tf,xs,xd,CFL,a,l,h,dt
		INTEGER :: Nt,Nx
		WRITE(*,*) ""
		WRITE(*,*) "********************************************"
		WRITE(*,*) "***********  ADVECTION EQUATION  ***********"
		WRITE(*,*) "********************************************"
		WRITE(*,'(A,F6.2,A,F6.2,A)') " Time Interval : [", ti, ", ", tf, " ] s"
		WRITE(*,'(A,F6.2,A,F6.2,A)') " Space Interval: [", xs, ", ", xd, " ] m"	
		WRITE(*,*) "CFL           :", CFL
		WRITE(*,*) "Time step     :", dt
		WRITE(*,*) "Grid points- t:", Nt
		WRITE(*,*) "Grid points- x:", Nx

	END SUBROUTINE inputpar_print



! ******************************************************************************
	
	! INITIAL VALUE PROBLEM CONDITIONS
	SUBROUTINE ivp(u,Nt,Nx,l,xs,xd,dx)
		INTEGER :: Nt, Nx
		REAL :: l,xs,xd,dx,step
		REAL, DIMENSION(Nt,Nx) :: u
		INTEGER :: k=1

		DO WHILE (k.LE.Nx)
			step=xs+(k-1)*dx
       			IF (step.GE.-l .AND. step.LE.l) THEN
				u(1,k) = f_ivp(step,l)
			ELSE
				u(1,k) = 0
			END IF
			
			k=k+1
		END DO

		WRITE(*,*) "********************************************"
		WRITE(*,*) "   Initial condition imposed"
		
	END SUBROUTINE ivp



! ******************************************************************************
	
	! BOUNDARY VALUE PROBLEM CONDITIONS
	SUBROUTINE bvp(u,Nt,Nx,l,ts,td,dt)
		INTEGER :: Nt, Nx
		REAL :: l,ts,td,dt,step
		REAL, DIMENSION(Nt,Nx) :: u
		INTEGER :: k=1

		DO WHILE (k.LE.Nt)
			step=ts+(k-1)*dt
       			IF (step.GE.1) THEN
				u(k,1) = f_bvp(step,l)
			ELSE
				u(k,1) = 0
			END IF

			k=k+1
		END DO

		WRITE(*,*) "********************************************"
		WRITE(*,*) "   Boundary condition imposed"

	END SUBROUTINE bvp



! ******************************************************************************
	
	! INITIAL CONDITION FUNCTION
	FUNCTION f_ivp(x,l)
		REAL :: x,l, f_ivp
		REAL :: pi = 3.1415926535

		f_ivp = sin(2*pi*x/l)

	END FUNCTION f_ivp



! ******************************************************************************
	
	! BOUNDARY CONDITION FUNCTION
	FUNCTION f_bvp(t,l)
		REAL :: t,l, f_bvp
		REAL :: pi = 3.1415926535

		f_bvp = -sin(2*pi*t/l)

	END FUNCTION f_bvp



! ******************************************************************************
	
	! CHECK IF IVP AND BVP ARE APPLIED
	SUBROUTINE print_conditions(u,Nx,Nt)
		INTEGER :: Nt, Nx
		REAL, DIMENSION(Nt,Nx) :: u
		INTEGER :: k=1
		
		OPEN(1, file='ivp.dat')

			DO WHILE (k.LE.Nx)
				WRITE(1,*) u(1,k)
				k=k+1
			END DO
			
		CLOSE(1)
		k=1
		OPEN(2, file='bvp.dat')

			DO WHILE (k.LE.Nt)
				WRITE(2,*) u(k,1)
				k=k+1
			END DO
			
		CLOSE(2)

	END SUBROUTINE print_conditions


! ******************************************************************************
	
	! LAX WENDROFF - SCHEME
	SUBROUTINE lax_wendroff(u,Nt,Nx,CFL)
		INTEGER :: Nt, Nx
		REAL :: CFL
		REAL, DIMENSION(Nt,Nx) :: u
		INTEGER :: n=1, j=2, k=1

		DO WHILE (n<(Nt-1))
			j=2
			DO WHILE (j.LE.(Nx-1))
				u(n+1,j)=u(n,j)-CFL/2*(u(n,j+1)-u(n,j-1))+CFL**2/2*(u(n,j+1)-2*u(n,j)+u(n,j-1))
				j=j+1
			END DO
			j=Nx
			u(n+1,j)=(1-CFL**2)*u(n,j)+CFL/2*((CFL-1)*(2*u(n,j)-u(n,j-1))+(CFL+1)*u(n,j-1))
			n=n+1
		END DO

	END SUBROUTINE lax_wendroff



! ******************************************************************************
	
	! PRINT THE SOLUTION
	SUBROUTINE print_solution(u,Nt,Nx,ind,ts,dt)
		INTEGER :: Nt, Nx
		REAL, DIMENSION(Nt,Nx) :: u
		REAL :: ts, dt, temp
		INTEGER :: ind, k=1
		
		OPEN(1, file='sol_adv.dat')
			
			DO WHILE (k.LE.Nx)
				temp=ts+(k-1)*dt
				WRITE(1,*) temp, u(ind,k)
				k=k+1
			END DO
			
		CLOSE(1)

	END SUBROUTINE print_solution



END PROGRAM advection_equation
